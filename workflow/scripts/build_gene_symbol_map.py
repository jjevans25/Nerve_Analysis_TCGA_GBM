"""Build a cached Ensembl-ID → HGNC-symbol map via MyGene.info for the GDC loom gene universe.

The GDC TCGA loom files use versioned Ensembl IDs (e.g. ENSG00000136492.9) as gene
identifiers. Downstream marker-based scoring operates on HGNC symbols, so every
analysis rule needs a persistent lookup. This rule fetches that lookup once per
Ensembl release and caches it as a TSV for reuse across samples and sessions.
"""

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import loompy
import pandas as pd
import requests

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
loom_paths: list[str] = list(snakemake.input.looms)
cache_tsv = Path(snakemake.output.cache)
ensembl_release: int = snakemake.params.ensembl_release
chunk_size: int = snakemake.params.chunk_size
request_timeout: int = snakemake.params.request_timeout


def _collect_ensembl_ids(loom_files: list[str]) -> tuple[list[str], list[str]]:
    """Return the sorted union of (versioned, unversioned) Ensembl IDs across loom files."""
    versioned: set[str] = set()
    for loom_file in loom_files:
        with loompy.connect(str(loom_file), mode="r") as ds:
            versioned.update(str(g) for g in ds.ra["Gene"])
    versioned_sorted = sorted(versioned)
    unversioned_sorted = sorted({vid.split(".", 1)[0] for vid in versioned_sorted})
    return versioned_sorted, unversioned_sorted


def _query_mygene(
    ensembl_ids: list[str],
    batch_size: int,
    timeout_s: int,
) -> pd.DataFrame:
    """Batch-POST to MyGene.info /v3/query and return hits as a DataFrame."""
    endpoint = "https://mygene.info/v3/query"
    hits: list[dict] = []
    n_batches = (len(ensembl_ids) + batch_size - 1) // batch_size
    for i in range(0, len(ensembl_ids), batch_size):
        chunk = ensembl_ids[i : i + batch_size]
        payload = {
            "q":       ",".join(chunk),
            "scopes":  "ensembl.gene",
            "fields":  "symbol,genomic_pos.chr",
            "species": "human",
        }
        response = requests.post(endpoint, data=payload, timeout=timeout_s)
        response.raise_for_status()
        hits.extend(response.json())
        log_transformation(log, "build_gene_symbol_map",
            f"MyGene batch {i // batch_size + 1}/{n_batches}: "
            f"sent {len(chunk)} IDs, received {len(response.json())} hits")
        time.sleep(0.1)  # courtesy rate-limit spacing
    return pd.DataFrame(hits)


def _extract_chromosome(row: pd.Series) -> str | None:
    """Normalize MyGene's genomic_pos field (list | dict | missing) to a single chr string."""
    pos = row.get("genomic_pos")
    if isinstance(pos, list) and pos:
        return pos[0].get("chr")
    if isinstance(pos, dict):
        return pos.get("chr")
    return None


log_transformation(log, "build_gene_symbol_map",
    f"Reading gene universe from {len(loom_paths)} loom files")
versioned_ids, unversioned_ids = _collect_ensembl_ids(loom_paths)
log_transformation(log, "build_gene_symbol_map",
    f"Gene universe: {len(versioned_ids)} versioned / {len(unversioned_ids)} unversioned Ensembl IDs")

log_transformation(log, "build_gene_symbol_map",
    f"Querying MyGene.info for {len(unversioned_ids)} IDs in batches of {chunk_size}")
try:
    hits_df = _query_mygene(unversioned_ids, chunk_size, request_timeout)
except requests.RequestException as exc:
    log_transformation(log, "build_gene_symbol_map",
        f"MyGene.info request failed: {exc}", status="ERROR")
    raise

if hits_df.empty:
    raise RuntimeError("MyGene.info returned zero hits — cannot build symbol map.")

# Handle multi-hits: MyGene can return >1 record per query if the Ensembl ID maps
# to multiple entries (rare for canonical genes). Keep the first per query.
hits_df = hits_df.drop_duplicates(subset="query", keep="first")
hits_df["chromosome"] = hits_df.apply(_extract_chromosome, axis=1)

# Build the cache DataFrame keyed on the unversioned ID with all versioned IDs retained.
versioned_map = pd.DataFrame({
    "ensembl_id":             versioned_ids,
    "ensembl_id_no_version":  [vid.split(".", 1)[0] for vid in versioned_ids],
})
cache_df = versioned_map.merge(
    hits_df.rename(columns={"query": "ensembl_id_no_version", "symbol": "gene_symbol"})[
        ["ensembl_id_no_version", "gene_symbol", "chromosome"]
    ],
    on="ensembl_id_no_version",
    how="left",
)

mapped = cache_df["gene_symbol"].notna().sum()
log_transformation(log, "build_gene_symbol_map",
    f"Mapped {mapped}/{len(cache_df)} Ensembl IDs to HGNC symbols "
    f"({100 * mapped / len(cache_df):.1f}%)")

# --- Write TSV cache ---------------------------------------------------------
cache_tsv.parent.mkdir(parents=True, exist_ok=True)
cache_df.to_csv(cache_tsv, sep="\t", index=False)
verify_artifact(cache_tsv, min_size_bytes=1024)

# --- Sidecar metadata (Ensembl release pin for FAIR reproducibility) --------
sidecar = cache_tsv.with_suffix(".meta.json")
sidecar.write_text(json.dumps({
    "ensembl_release":       ensembl_release,
    "mygene_endpoint":       "https://mygene.info/v3/query",
    "scopes":                "ensembl.gene",
    "species":               "human",
    "fields":                ["symbol", "genomic_pos.chr"],
    "chunk_size":            chunk_size,
    "n_versioned_ids":       len(versioned_ids),
    "n_unversioned_ids":     len(unversioned_ids),
    "n_mapped_to_symbol":    int(mapped),
    "retrieved_at_utc":      datetime.now(timezone.utc).isoformat(),
    "source_loom_files":     [str(Path(p).name) for p in loom_paths],
}, indent=2))

# --- Provenance --------------------------------------------------------------
prov = stamp_artifact(
    output_path=cache_tsv,
    rule_name="build_gene_symbol_map",
    input_paths=loom_paths,
    tool_versions={
        "requests": requests.__version__,
        "pandas":   pd.__version__,
        "loompy":   loompy.__version__,
    },
    parameters={
        "ensembl_release":   ensembl_release,
        "chunk_size":        chunk_size,
        "n_versioned_ids":   len(versioned_ids),
        "n_unversioned_ids": len(unversioned_ids),
        "n_mapped":          int(mapped),
    },
    description="Ensembl-ID → HGNC-symbol + chromosome cache built from MyGene.info",
    ontology_operation="operation:2497",  # EDAM: Gene ID conversion
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "build_gene_symbol_map", "Complete", status="SUCCESS",
                   artifact_paths=[str(cache_tsv), str(sidecar)])
