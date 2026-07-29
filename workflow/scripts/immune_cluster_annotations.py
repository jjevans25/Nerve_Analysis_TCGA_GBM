"""Assign immune subtype labels to each immune Leiden cluster by marker scoring.

Biological question: which canonical immune subtype (homeostatic microglia,
tumor-associated macrophage, T cell, NK cell, dendritic cell) does each
re-clustered immune Leiden cluster correspond to?

Scores the canonical subtype marker panels (config ``immune_cells.subtype_markers``)
per cell with ``sc.tl.score_genes``, then labels each cluster by its highest
mean score (argmax) — the same cluster-level argmax scheme used by
``scrna_annotate``. Emits a labeled AnnData carrying ``immune_subtype`` and a
per-cluster annotation table.
"""

import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    H5AD_COMPRESSION,
    compute_cluster_purity,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]  # type: ignore[name-defined]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)  # type: ignore[name-defined]

log_transformation(
    log, "immune_cluster_annotations", f"Loading {snakemake.input.h5ad}"  # type: ignore[name-defined]
)
adata = ad.read_h5ad(snakemake.input.h5ad)  # type: ignore[name-defined]

# --- Subtype marker panels (HGNC symbols) from config ------------------------
marker_sets: dict[str, list[str]] = {
    k: list(v) for k, v in snakemake.params.subtype_markers.items()  # type: ignore[name-defined]
}

# --- Build HGNC-symbol → var-index lookup ------------------------------------
# adata.var_names are versioned Ensembl IDs; the panels are HGNC symbols.
# gene_symbol is carried in var (inherited from malignancy_labeled.h5ad).
if "gene_symbol" not in adata.var.columns:
    raise KeyError(
        "[FAIR-ALERT] gene_symbol missing from immune subset var — cannot map "
        "marker symbols to Ensembl IDs."
    )
symbol_to_index: dict[str, str] = (
    adata.var.dropna(subset=["gene_symbol"])
    .reset_index()
    .drop_duplicates(subset="gene_symbol", keep="first")
    .set_index("gene_symbol")["index"]
    .to_dict()
)
log_transformation(
    log,
    "immune_cluster_annotations",
    f"Built symbol→var-index lookup: {len(symbol_to_index)} unique symbols "
    f"({adata.var['gene_symbol'].notna().sum()}/{adata.n_vars} genes mapped)",
)

# --- Score each subtype panel ------------------------------------------------
scored_sets: list[str] = []
for label, genes in marker_sets.items():
    available = [symbol_to_index[s] for s in genes if s in symbol_to_index]
    available = [e for e in available if e in adata.var_names]
    unmapped = [s for s in genes if s not in symbol_to_index]
    if available:
        sc.tl.score_genes(
            adata,
            gene_list=available,
            score_name=f"score_{label}",
            random_state=snakemake.params.random_seed,  # type: ignore[name-defined]
        )
        scored_sets.append(label)
        log_transformation(
            log,
            "immune_cluster_annotations",
            f"Scored '{label}': {len(available)}/{len(genes)} markers mapped"
            + (f"; unmapped: {unmapped}" if unmapped else ""),
        )
    else:
        log_transformation(
            log,
            "immune_cluster_annotations",
            f"WARNING: no markers for '{label}' mapped (symbols tried: {genes})",
            status="WARNING",
        )

if not scored_sets:
    raise RuntimeError(
        "[FAIR-ALERT] No immune subtype marker set could be scored — check that "
        "gene_symbol mapping and marker panels are valid for this dataset."
    )

# --- Assign subtype per cluster by max mean score ----------------------------
score_cols = [f"score_{s}" for s in scored_sets]
cluster_means = (
    adata.obs[["immune_leiden", *score_cols]]
    .groupby("immune_leiden", observed=True)[score_cols]
    .mean()
)
label_map: dict[str, str] = {
    cluster: cluster_means.loc[cluster].idxmax().replace("score_", "")
    for cluster in cluster_means.index
}
adata.obs["immune_subtype"] = (
    adata.obs["immune_leiden"].map(label_map).astype("category")
)

# Confidence: gap between best and second-best cluster-level score, per cell.
if len(score_cols) >= 2:
    sorted_scores = np.sort(adata.obs[score_cols].values, axis=1)
    adata.obs["immune_subtype_confidence"] = (
        sorted_scores[:, -1] - sorted_scores[:, -2]
    ).astype(np.float32)
else:
    adata.obs["immune_subtype_confidence"] = np.float32(1.0)

adata.uns["immune_annotation_schema"] = {
    "method": "leiden_marker_scoring_argmax",
    "marker_sets": {k: marker_sets[k] for k in scored_sets},
    "edam_operation": "operation:3432",
}
log_transformation(
    log,
    "immune_cluster_annotations",
    f"Subtype distribution: {adata.obs['immune_subtype'].value_counts().to_dict()}",
)

# --- Per-cluster annotation table --------------------------------------------
sizes = adata.obs["immune_leiden"].value_counts().rename("n_cells")
annotations = (
    cluster_means.assign(
        immune_subtype=[label_map[c] for c in cluster_means.index],
        n_cells=[int(sizes.get(c, 0)) for c in cluster_means.index],
    )
    .reset_index()
    .rename(columns={"immune_leiden": "cluster"})
)
# Order columns: cluster, subtype, n_cells, then the score columns.
annotations = annotations[["cluster", "immune_subtype", "n_cells", *score_cols]]
annotations = annotations.sort_values("n_cells", ascending=False).reset_index(drop=True)
annotations.to_csv(snakemake.output.annotations, index=False)  # type: ignore[name-defined]

# --- Per-subtype batch purity (keyed on immune_subtype for downstream QC) -----
# The three-way interaction groups immune cells by subtype, so QC flags for the
# interaction table are joined on immune_subtype (not the finer immune_leiden).
subtype_purity = compute_cluster_purity(
    adata,
    cluster_col="immune_subtype",
    batch_col="sample_id",
    dominant_max=snakemake.params.dominant_fraction_max,  # type: ignore[name-defined]
    min_contributing_fraction=snakemake.params.min_contributing_fraction,  # type: ignore[name-defined]
    min_contributing_samples=snakemake.params.min_contributing_samples,  # type: ignore[name-defined]
)
subtype_purity.df.to_csv(snakemake.output.subtype_purity, index=False)  # type: ignore[name-defined]
log_transformation(
    log,
    "immune_cluster_annotations",
    f"Per-subtype batch purity: "
    f"{int(subtype_purity.df['pass_overall'].sum())}/{len(subtype_purity.df)} subtypes PASS",
)

# --- Write labeled AnnData + provenance --------------------------------------
adata.write_h5ad(snakemake.output.h5ad, compression=H5AD_COMPRESSION)  # type: ignore[name-defined]
verify_artifact(snakemake.output.h5ad, min_size_bytes=512)  # type: ignore[name-defined]

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,  # type: ignore[name-defined]
    rule_name="immune_cluster_annotations",
    input_paths=[snakemake.input.h5ad],  # type: ignore[name-defined]
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "scored_sets": scored_sets,
        "cluster_labels": label_map,
        "n_cells": int(adata.n_obs),
        "n_clusters": int(len(label_map)),
    },
    description="Immune subset AnnData with per-cluster subtype labels (marker argmax)",
    ontology_operation="operation:3432",
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]

log_transformation(
    log,
    "immune_cluster_annotations",
    "Complete",
    status="SUCCESS",
    artifact_paths=[snakemake.output.h5ad, snakemake.output.annotations],  # type: ignore[name-defined]
)
