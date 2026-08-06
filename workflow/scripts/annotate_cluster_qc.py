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
from fair_utils import (  # noqa: E402
    NERVE_GROUP_PREFIX,
    nerve_group_key,
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
# Kept for reference; parsing now goes through fair_utils.nerve_group_key.
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
        # Via the shared helper, not a bare removeprefix: nerve groups are
        # `nerve_c{N}` for glial clusters AND `nerve_neuron` for the pooled
        # neuron group, which has a purity row of its own keyed "neuron".
        key = raw_key.map(nerve_group_key)
        if key.equals(raw_key):
            raise RuntimeError(
                f"[FAIR-ALERT] None of the values in '{join_column}' carry "
                f"the expected '{NERVE_GROUP_PREFIX}' prefix in {source_path}"
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


def _annotate_threeway(
    source_path: str | Path,
    out_path: str | Path,
    nerve_purity: pd.DataFrame,
    immune_purity: pd.DataFrame,
    *,
    nerve_exclude: set[str],
    immune_exclude: set[str],
) -> tuple[int, int]:
    """Join nerve + immune batch-QC verdicts onto the three-way interaction table.

    Each row involves at most one nerve cluster and at most one immune subtype
    (the third party, tumor, has no batch-QC verdict). Purity is joined on
    whichever compartment(s) the row spans; ``batch_qc_pass`` is the AND of the
    verdicts that apply (a not-involved compartment counts as pass). Rows whose
    nerve cluster or immune subtype is excluded are dropped, mirroring the
    per-cluster ``_annotate`` convention. Returns (rows_written, rows_dropped).
    """
    src = pd.read_csv(source_path)
    for col in ("nerve_cluster", "immune_subtype", "compartment_pair"):
        if col not in src.columns:
            raise RuntimeError(
                f"[FAIR-ALERT] Expected column '{col}' not in {source_path}"
            )
    n_src = len(src)
    src["nerve_cluster"] = src["nerve_cluster"].fillna("").astype(str)
    src["immune_subtype"] = src["immune_subtype"].fillna("").astype(str)
    src["_nerve_key"] = src["nerve_cluster"].str.removeprefix(NERVE_CLUSTER_PREFIX)
    src["_immune_key"] = src["immune_subtype"]

    nerve_cols = {
        "batch_qc_pass": "nerve_batch_qc_pass",
        "dominant_sample_fraction": "nerve_dominant_sample_fraction",
        "n_contributing_samples": "nerve_n_contributing_samples",
        "normalised_entropy": "nerve_normalised_entropy",
        "cluster_key": "_nerve_key",
    }
    immune_cols = {
        "batch_qc_pass": "immune_batch_qc_pass",
        "dominant_sample_fraction": "immune_dominant_sample_fraction",
        "n_contributing_samples": "immune_n_contributing_samples",
        "normalised_entropy": "immune_normalised_entropy",
        "cluster_key": "_immune_key",
    }
    merged = src.merge(nerve_purity.rename(columns=nerve_cols), on="_nerve_key", how="left")
    merged = merged.merge(immune_purity.rename(columns=immune_cols), on="_immune_key", how="left")
    if len(merged) != n_src:
        raise RuntimeError(
            f"[FAIR-ALERT] Row-count parity violated for {source_path}: "
            f"{n_src} -> {len(merged)} after join"
        )

    # Every *involved* cluster/subtype must resolve to a purity verdict.
    nerve_involved = merged["_nerve_key"] != ""
    immune_involved = merged["_immune_key"] != ""
    if (bad := merged.loc[nerve_involved & merged["nerve_batch_qc_pass"].isna(), "_nerve_key"]).any():
        raise RuntimeError(
            f"[FAIR-ALERT] nerve clusters in {source_path} missing from purity: "
            f"{sorted(set(bad))}"
        )
    if (bad := merged.loc[immune_involved & merged["immune_batch_qc_pass"].isna(), "_immune_key"]).any():
        raise RuntimeError(
            f"[FAIR-ALERT] immune subtypes in {source_path} missing from purity: "
            f"{sorted(set(bad))}"
        )

    # Combined verdict: AND of the verdicts that apply (not-involved = pass).
    nerve_ok = merged["nerve_batch_qc_pass"].fillna(True).astype(bool)
    immune_ok = merged["immune_batch_qc_pass"].fillna(True).astype(bool)
    merged["batch_qc_pass"] = nerve_ok & immune_ok

    dropped_mask = merged["_nerve_key"].isin(nerve_exclude) | merged["_immune_key"].isin(
        immune_exclude
    )
    n_dropped = int(dropped_mask.sum())
    merged = merged.loc[~dropped_mask].reset_index(drop=True)

    merged = merged.drop(columns=["_nerve_key", "_immune_key"])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)
    verify_artifact(out_path, min_size_bytes=64)
    return len(merged), n_dropped


purity = _load_purity(snakemake.input.purity)
immune_purity = _load_purity(snakemake.input.immune_purity)
n_pass = int(purity["batch_qc_pass"].sum())
n_fail = int((~purity["batch_qc_pass"]).sum())
log_transformation(
    log,
    "annotate_cluster_qc",
    f"Loaded purity table: {n_pass} PASS / {n_fail} FAIL across {len(purity)} clusters",
)

exclude_clusters = {str(c) for c in getattr(snakemake.params, "exclude_clusters", []) or []}
immune_exclude = {
    str(c) for c in getattr(snakemake.params, "immune_exclude_subtypes", []) or []
}
if exclude_clusters:
    log_transformation(
        log,
        "annotate_cluster_qc",
        f"Excluding {len(exclude_clusters)} nerve cluster(s) from _with_qc.csv outputs: "
        f"{sorted(exclude_clusters)}",
    )
if immune_exclude:
    log_transformation(
        log,
        "annotate_cluster_qc",
        f"Excluding {len(immune_exclude)} immune subtype(s) from three-way "
        f"_with_qc.csv outputs: {sorted(immune_exclude)}",
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

# --- Three-way nerve-tumor-immune interaction tables -------------------------
threeway_jobs = [
    (snakemake.input.tw_interactions, snakemake.output.tw_interactions),
    (snakemake.input.tw_top_pairs, snakemake.output.tw_top_pairs),
]
for src, out in threeway_jobs:
    n_rows, n_dropped = _annotate_threeway(
        src, out, purity, immune_purity,
        nerve_exclude=exclude_clusters, immune_exclude=immune_exclude,
    )
    row_counts[str(out)] = n_rows
    rows_dropped[str(out)] = n_dropped
    log_transformation(
        log,
        "annotate_cluster_qc",
        f"Annotated (3-way) {src} -> {out} ({n_rows} rows; {n_dropped} dropped)",
    )

prov = stamp_artifact(
    output_path=snakemake.output.provenance,
    rule_name="annotate_cluster_qc",
    input_paths=[
        snakemake.input.purity,
        snakemake.input.immune_purity,
        snakemake.input.enrichment,
        snakemake.input.markers,
        snakemake.input.interactions,
        snakemake.input.top_pairs,
        snakemake.input.tw_interactions,
        snakemake.input.tw_top_pairs,
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
        "immune_exclude_subtypes": sorted(immune_exclude),
        "appended_columns": [
            "batch_qc_pass",
            "dominant_sample_fraction",
            "n_contributing_samples",
            "normalised_entropy",
        ],
        "threeway_appended_columns": [
            "batch_qc_pass",
            "nerve_batch_qc_pass",
            "nerve_dominant_sample_fraction",
            "nerve_n_contributing_samples",
            "nerve_normalised_entropy",
            "immune_batch_qc_pass",
            "immune_dominant_sample_fraction",
            "immune_n_contributing_samples",
            "immune_normalised_entropy",
        ],
    },
    description=(
        "Downstream annotation pass: left-joins batch-QC verdict from "
        "nerve_cluster_sample_purity.csv and immune_subtype_sample_purity.csv onto "
        "enrichment, marker, and LIANA tumor-nerve and nerve-tumor-immune tables. "
        "Source tables are not modified."
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
        snakemake.output.tw_interactions,
        snakemake.output.tw_top_pairs,
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
        snakemake.output.tw_interactions,
        snakemake.output.tw_top_pairs,
        snakemake.output.provenance,
    ],
)
