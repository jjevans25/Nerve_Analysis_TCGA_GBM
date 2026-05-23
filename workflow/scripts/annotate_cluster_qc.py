"""Annotate downstream per-cluster tables with batch-QC verdict from purity table.

Question answered: when reading per-cluster GSEA, marker, or LIANA results,
which clusters failed the batch-correction QC and should not be interpreted
as cohort-level signatures?

Reads the v1.0.0 outputs of `nerve_batch_qc`, `nerve_cell_heterogeneity`, and
`nerve_tumor_interaction`; left-joins purity columns onto each downstream
table and writes annotated `_with_qc.csv` variants. Original v1.0.0 tables
stay byte-identical; this rule never modifies them in place.

Appended columns: `batch_qc_pass` (bool, from `pass_overall`),
`dominant_sample_fraction` (float), `n_contributing_samples` (int),
`normalised_entropy` (float).
"""

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import (
    artifact_id,
    file_sha256,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = "0"

PURITY_COLUMNS = [
    "cluster",
    "pass_overall",
    "dominant_sample_fraction",
    "n_contributing_samples",
    "normalised_entropy",
]
NERVE_CLUSTER_PREFIX = "nerve_c"


def _load_purity(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=PURITY_COLUMNS)
    df = df.rename(columns={"pass_overall": "batch_qc_pass"})
    df["cluster_key"] = df["cluster"].astype(str)
    return df.drop(columns=["cluster"])


def _annotate(
    source_path: str | Path,
    out_path: str | Path,
    purity: pd.DataFrame,
    join_column: str,
    *,
    strip_prefix: bool,
    exclude_clusters: set[str],
) -> tuple[int, int]:
    """Left-join purity onto a source table; drop excluded clusters; write CSV.

    Returns (rows_written, rows_dropped) for sanity logging. Source per-cluster
    tables are not modified; only the downstream `_with_qc.csv` variant has the
    excluded clusters filtered out.
    """
    src = pd.read_csv(source_path)
    if join_column not in src.columns:
        raise RuntimeError(
            f"[FAIR-ALERT] Expected join column '{join_column}' not in {source_path}"
        )
    raw_key = src[join_column].astype(str)
    if strip_prefix:
        key = raw_key.str.removeprefix(NERVE_CLUSTER_PREFIX)
        if key.equals(raw_key):
            raise RuntimeError(
                f"[FAIR-ALERT] None of the values in '{join_column}' carry "
                f"the expected '{NERVE_CLUSTER_PREFIX}' prefix in {source_path}"
            )
    else:
        key = raw_key
    src = src.assign(_join_key=key)

    merged = src.merge(purity, how="left", left_on="_join_key", right_on="cluster_key")
    if len(merged) != len(src):
        raise RuntimeError(
            f"[FAIR-ALERT] Row-count parity violated for {source_path}: "
            f"{len(src)} -> {len(merged)} after join"
        )
    n_unmatched = int(merged["batch_qc_pass"].isna().sum())
    if n_unmatched:
        unmatched_ids = sorted(set(merged.loc[merged["batch_qc_pass"].isna(), "_join_key"]))
        raise RuntimeError(
            f"[FAIR-ALERT] {n_unmatched} rows in {source_path} have a cluster "
            f"id not found in the purity table: {unmatched_ids}"
        )

    dropped_mask = merged["_join_key"].isin(exclude_clusters)
    n_dropped = int(dropped_mask.sum())
    if n_dropped:
        merged = merged.loc[~dropped_mask].reset_index(drop=True)

    merged = merged.drop(columns=["_join_key", "cluster_key"])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)
    verify_artifact(out_path, min_size_bytes=64)
    return len(merged), n_dropped


purity = _load_purity(snakemake.input.purity)
n_pass = int(purity["batch_qc_pass"].sum())
n_fail = int((~purity["batch_qc_pass"]).sum())
log_transformation(
    log,
    "annotate_cluster_qc",
    f"Loaded purity table: {n_pass} PASS / {n_fail} FAIL across {len(purity)} clusters",
)

exclude_clusters = {str(c) for c in getattr(snakemake.params, "exclude_clusters", []) or []}
if exclude_clusters:
    log_transformation(
        log,
        "annotate_cluster_qc",
        f"Excluding {len(exclude_clusters)} cluster(s) from _with_qc.csv outputs: "
        f"{sorted(exclude_clusters)}",
    )

jobs = [
    (snakemake.input.enrichment,      snakemake.output.enrichment,      "cluster",       False),
    (snakemake.input.markers,         snakemake.output.markers,         "cluster",       False),
    (snakemake.input.interactions,    snakemake.output.interactions,    "nerve_cluster", True),
    (snakemake.input.top_pairs,       snakemake.output.top_pairs,       "nerve_cluster", True),
]

row_counts: dict[str, int] = {}
rows_dropped: dict[str, int] = {}
for src, out, join_col, strip in jobs:
    n_rows, n_dropped = _annotate(
        src, out, purity, join_col,
        strip_prefix=strip, exclude_clusters=exclude_clusters,
    )
    row_counts[str(out)] = n_rows
    rows_dropped[str(out)] = n_dropped
    log_transformation(
        log,
        "annotate_cluster_qc",
        f"Annotated {src} -> {out} ({n_rows} rows; {n_dropped} dropped via exclude_clusters)",
    )

prov = stamp_artifact(
    output_path=snakemake.output.provenance,
    rule_name="annotate_cluster_qc",
    input_paths=[
        snakemake.input.purity,
        snakemake.input.enrichment,
        snakemake.input.markers,
        snakemake.input.interactions,
        snakemake.input.top_pairs,
    ],
    tool_versions={
        "pandas": pd.__version__,
    },
    parameters={
        "n_clusters": int(len(purity)),
        "n_clusters_pass": n_pass,
        "n_clusters_fail": n_fail,
        "row_counts": row_counts,
        "rows_dropped": rows_dropped,
        "exclude_clusters": sorted(exclude_clusters),
        "appended_columns": [
            "batch_qc_pass",
            "dominant_sample_fraction",
            "n_contributing_samples",
            "normalised_entropy",
        ],
    },
    description=(
        "Downstream annotation pass: left-joins batch-QC verdict from "
        "nerve_cluster_sample_purity.csv onto enrichment, marker, and LIANA "
        "tumor-nerve tables. Source tables are not modified."
    ),
    ontology_operation="operation:3436",  # EDAM: Aggregation
)
prov["outputs"] = [
    {
        "path": str(Path(p).resolve()),
        "artifact_id": artifact_id(p),
        "sha256": file_sha256(p),
    }
    for p in (
        snakemake.output.enrichment,
        snakemake.output.markers,
        snakemake.output.interactions,
        snakemake.output.top_pairs,
    )
]
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log,
    "annotate_cluster_qc",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.enrichment,
        snakemake.output.markers,
        snakemake.output.interactions,
        snakemake.output.top_pairs,
        snakemake.output.provenance,
    ],
)
