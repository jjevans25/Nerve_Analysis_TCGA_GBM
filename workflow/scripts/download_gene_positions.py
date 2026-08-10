"""Fetch Ensembl gene coordinates so CNV inference can order genes along the genome.

Biological question: which genes are physically adjacent on a chromosome? CNV
callers detect malignancy by finding *runs* of co-elevated or co-depressed
expression across neighbouring loci — a whole-arm gain or loss. That signal only
exists if genes are in genomic order.

The previous implementation ordered genes by the numeric suffix of their Ensembl
ID (`scrna_malignancy.py`, ENSG00000141510 -> 141510), which is an accession
counter, not a coordinate: it interleaves chromosomes arbitrarily, so the
sliding window averaged unrelated loci and the resulting "CNV score" was mostly
noise (18% recall / 48% precision against the Census annotation).

Nothing in this project carried real coordinates. The Census `.var` exposes only
feature_id / feature_name / feature_length, and `dataset_gene_symbol_map.py`
writes the literal string "unknown" into its `chromosome` column. Hence this
rule. Mirrors `download_msigdb_gmt.py`: retry + backoff, optional pinned SHA-256,
manifest sidecar, gitignored payload.
"""

# NOTE: no `from __future__ import annotations` here. Snakemake's `script:`
# directive prepends its own preamble to this file, so a __future__ import is no
# longer the first statement and raises SyntaxError at runtime. Python 3.12
# evaluates `str | None` natively, so it buys nothing.

import gzip
import json
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, "workflow/scripts")
from fair_utils import (
    file_sha256,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]

release: int = int(snakemake.params.release)
url: str = snakemake.params.url
timeout_s: int = int(snakemake.params.request_timeout)
max_retries: int = int(snakemake.params.max_retries)
expected_sha: str | None = snakemake.params.sha256 or None

gtf_path = Path(snakemake.output.gtf)
positions_path = Path(snakemake.output.positions)
manifest_path = Path(snakemake.output.manifest)

# Standard human assembly units only. Scaffolds and patches carry genes that are
# either duplicates of a primary-assembly locus or unplaced, and either way they
# cannot participate in a positional window.
PRIMARY_CONTIGS: list[str] = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def _download_with_retry(dest: Path, retries: int) -> None:
    """Fetch the GTF with exponential backoff and byte-range resume.

    Ensembl's HTTPS mirror stalls mid-stream on a 64 MB transfer often enough
    that a plain retry restarts from zero and stalls again. Resuming from the
    partial file makes progress monotonic. `timeout` is a (connect, read) pair so
    a stalled socket raises instead of hanging until the job is killed.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            have = tmp.stat().st_size if tmp.exists() else 0
            headers = {"Range": f"bytes={have}-"} if have else {}
            with requests.get(url, stream=True, timeout=(30, timeout_s), headers=headers) as resp:
                # 416 = the range starts at or past EOF, i.e. the partial file is
                # already the whole file. Retrying cannot improve on that.
                if resp.status_code == 416 and have:
                    log_transformation(log, "download_gene_positions",
                        f"Partial file is already complete ({have} bytes); accepting it")
                    tmp.replace(dest)
                    return
                # 206 = server honored the range; 200 = full body, so start over.
                if have and resp.status_code == 200:
                    have = 0
                resp.raise_for_status()
                with open(tmp, "ab" if have else "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=1 << 20):
                        fh.write(chunk)
            tmp.replace(dest)
            return
        except (requests.RequestException, OSError) as exc:
            last_exc = exc
            wait = 2 ** attempt
            log_transformation(log, "download_gene_positions",
                f"Attempt {attempt + 1}/{retries} failed ({exc}); retrying in {wait}s",
                status="WARNING")
            time.sleep(wait)
    raise RuntimeError(f"[FAIR-ALERT] could not download {url} after {retries} attempts: {last_exc}")


def parse_gene_records(path: Path) -> pd.DataFrame:
    """Extract one row per gene: Ensembl ID, contig, start, end, strand, biotype."""
    rows: list[tuple[str, str, int, int, str, str, str]] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            contig, start, end, strand, attrs = parts[0], int(parts[3]), int(parts[4]), parts[6], parts[8]
            gene_id = gene_name = biotype = ""
            for field in attrs.split(";"):
                field = field.strip()
                if field.startswith("gene_id "):
                    gene_id = field.split('"')[1]
                elif field.startswith("gene_name "):
                    gene_name = field.split('"')[1]
                elif field.startswith("gene_biotype "):
                    biotype = field.split('"')[1]
            if gene_id:
                # Ensembl IDs are versioned in the GTF (ENSG...\.3); the cohort
                # carries unversioned IDs, so strip the suffix to join on.
                rows.append((gene_id.split(".")[0], contig, start, end, strand, gene_name, biotype))
    return pd.DataFrame(
        rows, columns=["ensembl_id", "chromosome", "start", "end", "strand", "gene_symbol", "biotype"]
    )


# ---------------------------------------------------------------------------
# Download + verify
# ---------------------------------------------------------------------------
log_transformation(log, "download_gene_positions", f"Fetching Ensembl {release} GTF from {url}")
_download_with_retry(gtf_path, max_retries)
verify_artifact(gtf_path, min_size_bytes=1 << 20)

observed_sha = file_sha256(gtf_path)
if expected_sha and observed_sha != expected_sha:
    raise RuntimeError(
        f"[FAIR-ALERT] SHA-256 mismatch for {gtf_path.name}: expected {expected_sha}, "
        f"got {observed_sha}. The remote resource changed under a pinned release — "
        "do not proceed with an unverified reference."
    )
if not expected_sha:
    log_transformation(log, "download_gene_positions",
        f"No pinned sha256 in config; observed {observed_sha}. Copy it into "
        "config.gene_positions.sha256 and commit so future runs verify.",
        status="WARNING")

# ---------------------------------------------------------------------------
# Parse to a flat coordinate table
# ---------------------------------------------------------------------------
genes = parse_gene_records(gtf_path)
n_all = len(genes)
genes = genes[genes["chromosome"].isin(PRIMARY_CONTIGS)].copy()
genes = genes.drop_duplicates(subset="ensembl_id", keep="first")

# Sort into genome order once, here, so every consumer inherits the same
# ordering: chromosome 1..22, X, Y, MT, then ascending start within a contig.
genes["_contig_rank"] = genes["chromosome"].map({c: i for i, c in enumerate(PRIMARY_CONTIGS)})
genes = genes.sort_values(["_contig_rank", "start"]).drop(columns="_contig_rank").reset_index(drop=True)
genes["genome_order"] = range(len(genes))

positions_path.parent.mkdir(parents=True, exist_ok=True)
genes.to_csv(positions_path, sep="\t", index=False)
verify_artifact(positions_path, min_size_bytes=1024)

log_transformation(log, "download_gene_positions",
    f"Parsed {n_all} gene records → {len(genes)} on primary contigs "
    f"({genes['chromosome'].nunique()} contigs); written in genome order")

manifest = {
    "manifest_id":   str(uuid.uuid4()),
    "retrieved_at":  datetime.now(timezone.utc).isoformat(),
    "ensembl_release": release,
    "source_url":    url,
    "gtf_sha256":    observed_sha,
    "positions_sha256": file_sha256(positions_path),
    "n_genes_total": int(n_all),
    "n_genes_primary_contigs": int(len(genes)),
    "primary_contigs": PRIMARY_CONTIGS,
    "license":       "Ensembl data is released under the Apache 2.0 licence / no restriction",
}
manifest_path.parent.mkdir(parents=True, exist_ok=True)
with open(manifest_path, "w") as fh:
    json.dump(manifest, fh, indent=2)

prov = stamp_artifact(
    output_path=positions_path,
    rule_name="download_gene_positions",
    input_paths=[],
    tool_versions={"requests": requests.__version__, "pandas": pd.__version__},
    parameters={
        "ensembl_release": release,
        "source_url":      url,
        "gtf_sha256":      observed_sha,
        "n_genes":         int(len(genes)),
    },
    description="Ensembl gene coordinates (chromosome, start, end) in genome order for CNV inference",
    ontology_operation="operation:2422",  # EDAM: Data retrieval
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "download_gene_positions", "Complete", status="SUCCESS",
                   artifact_paths=[str(positions_path), str(manifest_path)])
