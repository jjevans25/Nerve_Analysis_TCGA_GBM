"""Download MSigDB .gmt collections for offline GSEA (gseapy.prerank).

Replaces the live Enrichr API call inside `nerve_cell_heterogeneity` with a
deterministic, FAIR-compliant local resource. Each .gmt is fetched with retry
+ exponential backoff, validated against an optional pinned SHA-256, and
recorded in a manifest sidecar (release, source URL, hash, retrieval time,
license).

The .gmt files themselves are gitignored; the manifest is committed so
collaborators can verify the integrity of any local copy they obtain.
"""

from __future__ import annotations

import json
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

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

release: str = snakemake.params.release
timeout_s: int = int(snakemake.params.request_timeout)
max_retries: int = int(snakemake.params.max_retries)
collections: dict = dict(snakemake.params.collections)

gmt_bp_path = Path(snakemake.output.gmt_bp)
gmt_mf_path = Path(snakemake.output.gmt_mf)
manifest_path = Path(snakemake.output.manifest)
provenance_path = Path(snakemake.output.provenance)

_OUTPUT_BY_KEY: dict[str, Path] = {
    "gmt_bp": gmt_bp_path,
    "gmt_mf": gmt_mf_path,
}


def _download_with_retry(
    url: str,
    dest: Path,
    timeout: int,
    retries: int,
) -> None:
    """Fetch `url` to `dest` with exponential backoff. Atomic via temp file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_exc: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=timeout) as resp:
                resp.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=65536):
                        if chunk:
                            fh.write(chunk)
            tmp.replace(dest)
            log_transformation(
                log,
                "download_msigdb_gmt",
                f"Downloaded {dest.name} on attempt {attempt}",
            )
            return
        except (requests.RequestException, OSError) as exc:
            last_exc = exc
            wait = 2 ** (attempt - 1)
            log_transformation(
                log,
                "download_msigdb_gmt",
                f"Attempt {attempt}/{retries} for {url} failed: {exc}; "
                f"sleeping {wait}s",
                status="WARNING",
            )
            if tmp.exists():
                tmp.unlink()
            if attempt < retries:
                time.sleep(wait)
    raise RuntimeError(
        f"[FAIR-ALERT] Could not download {url} after {retries} attempts: {last_exc}"
    )


def _gmt_term_count(path: Path) -> int:
    """Count gene-set lines (each line is one term)."""
    with open(path, "rt", encoding="utf-8") as fh:
        return sum(1 for line in fh if line.strip())


def _validate_gmt(path: Path) -> None:
    """Check the file looks like a GMT (≥3 tab-separated fields per line)."""
    with open(path, "rt", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"[FAIR-ALERT] {path.name}: line {line_no} has only "
                    f"{len(fields)} tab-separated fields (need ≥3)"
                )
            if line_no >= 5:
                break  # spot-check: format is consistent across MSigDB releases


manifest_entries: list[dict] = []
gmt_outputs: list[Path] = []

for key, meta in collections.items():
    if key not in _OUTPUT_BY_KEY:
        raise KeyError(
            f"[FAIR-ALERT] Unknown collection key '{key}' in config.msigdb.collections; "
            f"expected one of {sorted(_OUTPUT_BY_KEY)}"
        )
    dest = _OUTPUT_BY_KEY[key]
    expected_filename = meta["filename"]
    if dest.name != expected_filename:
        raise ValueError(
            f"[FAIR-ALERT] Output path {dest} does not match config filename {expected_filename}"
        )

    log_transformation(
        log,
        "download_msigdb_gmt",
        f"Fetching {key} ({meta['library_label']}) from {meta['url']}",
    )
    _download_with_retry(meta["url"], dest, timeout_s, max_retries)
    verify_artifact(dest, min_size_bytes=10_000)  # smallest MSigDB collection is ~40 KB
    _validate_gmt(dest)

    digest = file_sha256(dest)
    expected_sha = (meta.get("sha256") or "").strip().lower()
    if expected_sha:
        if digest != expected_sha:
            raise ValueError(
                f"[FAIR-ALERT] SHA-256 mismatch for {dest.name}: "
                f"expected {expected_sha}, got {digest}. "
                f"Pinned digest in config is stale or download is corrupt."
            )
        log_transformation(
            log,
            "download_msigdb_gmt",
            f"SHA-256 verified for {dest.name}",
        )
    else:
        log_transformation(
            log,
            "download_msigdb_gmt",
            f"No pinned SHA-256 in config for {key}; recording {digest} in manifest. "
            f"Copy this digest into config.msigdb.collections.{key}.sha256 to lock it.",
            status="WARNING",
        )

    n_terms = _gmt_term_count(dest)
    manifest_entries.append(
        {
            "key": key,
            "library_label": meta["library_label"],
            "filename": dest.name,
            "path": str(dest.resolve()),
            "url": meta["url"],
            "sha256": digest,
            "n_terms": n_terms,
            "size_bytes": dest.stat().st_size,
            "artifact_uuid": str(
                uuid.uuid5(uuid.NAMESPACE_URL, f"msigdb:{release}:{dest.name}")
            ),
        }
    )
    gmt_outputs.append(dest)
    log_transformation(
        log,
        "download_msigdb_gmt",
        f"{dest.name}: {n_terms} gene sets, sha256={digest[:12]}…",
    )

# --- Manifest -----------------------------------------------------------------
manifest = {
    "msigdb_release": release,
    "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
    "license": "CC-BY-4.0",
    "license_note": (
        "MSigDB v2024+ collections are released under CC-BY 4.0. "
        "Earlier releases require attribution per the Broad Institute terms."
    ),
    "ontology_data": "data_3753",  # EDAM: Gene set
    "source_root": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/",
    "collections": manifest_entries,
}
manifest_path.parent.mkdir(parents=True, exist_ok=True)
manifest_path.write_text(json.dumps(manifest, indent=2))
verify_artifact(manifest_path, min_size_bytes=128)

# --- Provenance --------------------------------------------------------------
prov = stamp_artifact(
    output_path=manifest_path,
    rule_name="download_msigdb_gmt",
    input_paths=[],
    tool_versions={"requests": requests.__version__},
    parameters={
        "release": release,
        "n_collections": len(manifest_entries),
        "request_timeout_s": timeout_s,
        "max_retries": max_retries,
        "sha256_per_collection": {
            entry["key"]: entry["sha256"] for entry in manifest_entries
        },
    },
    description=(
        "MSigDB gene-set collections (.gmt) for offline GSEA via gseapy.prerank"
    ),
    ontology_operation="operation:2422",  # EDAM: Data retrieval
)
write_provenance(prov, provenance_path)

log_transformation(
    log,
    "download_msigdb_gmt",
    "Complete",
    status="SUCCESS",
    artifact_paths=[str(p) for p in gmt_outputs] + [str(manifest_path)],
)
