"""Differential expression, GSEA, marker dot plot, and per-sample cluster abundance."""

import os
import sys
import warnings

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import (
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

log_transformation(log, "nerve_cell_heterogeneity", f"Loading {snakemake.input.h5ad}")
adata = ad.read_h5ad(snakemake.input.h5ad)

# --- Build HGNC-symbol ↔ Ensembl-ID lookup (mirrors scrna_annotate.py:51-77) -
# Markers in config.nerve_cells.markers are HGNC symbols; adata.var_names are
# versioned Ensembl IDs. We need both directions:
#   symbol_to_ensembl  → for §3.3 dot-plot (translate marker list)
#   ensembl_to_symbol  → for §3.1 markers CSV (annotate DE results)
if "gene_symbol" not in adata.var.columns:
    log_transformation(
        log,
        "nerve_cell_heterogeneity",
        "gene_symbol missing from var — re-joining from MyGene cache",
        status="WARNING",
    )
    _symbol_map = pd.read_csv(snakemake.input.symbol_map, sep="\t", dtype=str)
    if "ensembl_id" not in adata.var.columns:
        adata.var["ensembl_id"] = adata.var_names
    _var_joined = adata.var.merge(
        _symbol_map[["ensembl_id", "gene_symbol", "chromosome"]],
        on="ensembl_id",
        how="left",
    )
    _var_joined.index = adata.var.index
    adata.var = _var_joined

symbol_to_ensembl: dict[str, str] = (
    adata.var.dropna(subset=["gene_symbol"])
    .reset_index()
    .drop_duplicates(subset="gene_symbol", keep="first")
    .set_index("gene_symbol")["index"]
    .to_dict()
)
ensembl_to_symbol: dict[str, str] = adata.var.dropna(subset=["gene_symbol"])[
    "gene_symbol"
].to_dict()
log_transformation(
    log,
    "nerve_cell_heterogeneity",
    f"Built symbol↔Ensembl lookup: {len(symbol_to_ensembl)} unique symbols "
    f"({adata.var['gene_symbol'].notna().sum()}/{adata.n_vars} genes mapped)",
)

n_clusters = adata.obs["nerve_leiden"].nunique() if "nerve_leiden" in adata.obs else 0
log_transformation(
    log,
    "nerve_cell_heterogeneity",
    f"{adata.n_obs} cells across {n_clusters} nerve-cell clusters",
)

if adata.n_obs == 0 or n_clusters == 0:
    log_transformation(
        log,
        "nerve_cell_heterogeneity",
        "[FAIR-ALERT] Empty nerve-cell AnnData — writing placeholder outputs.",
        status="WARNING",
    )
    pd.DataFrame(
        columns=["cluster", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    ).to_csv(snakemake.output.markers, index=False)
    pd.DataFrame().to_csv(snakemake.output.enrichment, index=False)
    for fig_path in [snakemake.output.dotplot, snakemake.output.abundance]:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(
            0.5,
            0.5,
            "[FAIR-ALERT] No nerve cells found in subsampled data",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=11,
        )
        ax.set_axis_off()
        fig.savefig(fig_path, dpi=100, bbox_inches="tight")
        plt.close(fig)
    prov = stamp_artifact(
        output_path=snakemake.output.markers,
        rule_name="nerve_cell_heterogeneity",
        input_paths=[snakemake.input.h5ad],
        tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
        parameters={"n_clusters": 0, "n_cells": 0, "note": "no nerve cells found"},
        description="Placeholder: no nerve cells available for heterogeneity analysis",
        ontology_operation="operation:3223",
    )
    write_provenance(prov, snakemake.output.provenance)
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# Step 1: Differential expression (Wilcoxon) per nerve cluster
# ---------------------------------------------------------------------------
# Restore log-normalized layer from raw if available
if adata.raw is not None:
    adata_de = adata.raw.to_adata()
    adata_de.obs = adata.obs.copy()
else:
    adata_de = adata.copy()
    if not ("log1p" in adata_de.uns or adata_de.X.max() < 20):
        sc.pp.normalize_total(adata_de, target_sum=1e4)
        sc.pp.log1p(adata_de)

sc.tl.rank_genes_groups(
    adata_de,
    groupby="nerve_leiden",
    method="wilcoxon",
    n_genes=50,
    key_added="rank_genes_nerve",
)

# Flatten results to a DataFrame
marker_rows = []
for group in adata_de.obs["nerve_leiden"].unique():
    result = sc.get.rank_genes_groups_df(
        adata_de, group=group, key="rank_genes_nerve", pval_cutoff=0.05
    )
    result.insert(0, "cluster", group)
    marker_rows.append(result)

if marker_rows:
    markers_df = pd.concat(marker_rows, ignore_index=True)
else:
    markers_df = pd.DataFrame(
        columns=["cluster", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    )

# §3.1: attach human-readable HGNC symbol next to the Ensembl ID in `names`.
markers_df["gene_symbol"] = markers_df["names"].map(ensembl_to_symbol)
n_unmapped = int(markers_df["gene_symbol"].isna().sum())
_cols = list(markers_df.columns)
_cols.insert(_cols.index("names") + 1, _cols.pop(_cols.index("gene_symbol")))
markers_df = markers_df[_cols]
markers_df.to_csv(snakemake.output.markers, index=False)
log_transformation(
    log,
    "nerve_cell_heterogeneity",
    f"§3.1 fix: annotated markers CSV — {len(markers_df) - n_unmapped}/{len(markers_df)} "
    f"rows have a gene_symbol ({n_unmapped} unmapped)",
)
log_transformation(
    log,
    "nerve_cell_heterogeneity",
    f"DE complete: {len(markers_df)} significant marker genes across {n_clusters} clusters",
)

# ---------------------------------------------------------------------------
# Step 2: Gene set enrichment (gseapy) on top cluster markers
# ---------------------------------------------------------------------------
enrichment_rows = []
try:
    import gseapy as gp

    gene_sets = ["GO_Biological_Process_2023", "GO_Molecular_Function_2023"]

    for cluster in markers_df["cluster"].unique():
        cluster_markers = (
            markers_df[markers_df["cluster"] == cluster]
            .sort_values("scores", ascending=False)
            .head(100)["names"]
            .tolist()
        )
        if len(cluster_markers) < 5:
            continue
        for gs in gene_sets:
            try:
                enr = gp.enrichr(
                    gene_list=cluster_markers,
                    gene_sets=gs,
                    organism="human",
                    outdir=None,
                    verbose=False,
                )
                if enr.results is not None and len(enr.results) > 0:
                    top = enr.results.head(10).copy()
                    top.insert(0, "cluster", cluster)
                    top.insert(1, "gene_set_library", gs)
                    enrichment_rows.append(top)
            except Exception as exc:
                log_transformation(
                    log,
                    "nerve_cell_heterogeneity",
                    f"WARNING: enrichr failed for cluster {cluster}, {gs}: {exc}",
                    status="WARNING",
                )

    log_transformation(
        log,
        "nerve_cell_heterogeneity",
        f"GSEA complete: {len(enrichment_rows)} result blocks",
    )

except ImportError:
    warnings.warn("[FAIR-ALERT] gseapy not available; skipping GSEA")
    log_transformation(
        log,
        "nerve_cell_heterogeneity",
        "WARNING: gseapy not installed — GSEA skipped",
        status="WARNING",
    )

enrichment_df = (
    pd.concat(enrichment_rows, ignore_index=True) if enrichment_rows else pd.DataFrame()
)
enrichment_df.to_csv(snakemake.output.enrichment, index=False)

# ---------------------------------------------------------------------------
# Step 3: Marker dot plot for canonical nerve-cell markers
# ---------------------------------------------------------------------------
# §3.3: config.nerve_cells.markers are HGNC symbols; adata.var_names are
# Ensembl IDs. Translate symbols → Ensembl, then keep only those actually
# present in the dataset. Preserve the symbol for axis labelling.
all_markers_flat: list[str] = [
    g for genes in snakemake.params.markers.values() for g in genes
]
present_pairs: list[tuple[str, str]] = [
    (sym, symbol_to_ensembl[sym])
    for sym in all_markers_flat
    if sym in symbol_to_ensembl and symbol_to_ensembl[sym] in adata.var_names
]
unmapped_markers = [s for s in all_markers_flat if s not in symbol_to_ensembl]
present_markers = [eid for _sym, eid in present_pairs]
present_symbols = [sym for sym, _eid in present_pairs]
log_transformation(
    log,
    "nerve_cell_heterogeneity",
    f"§3.3 fix: {len(present_markers)}/{len(all_markers_flat)} canonical markers "
    f"present in dataset"
    + (f"; unmapped: {unmapped_markers}" if unmapped_markers else ""),
)

if len(present_markers) >= 1:
    fig, ax = plt.subplots(
        figsize=(max(8, len(present_markers) * 0.45), max(4, n_clusters * 0.6))
    )
    sc.pl.dotplot(
        adata_de,
        var_names=present_markers,
        groupby="nerve_leiden",
        ax=ax,
        show=False,
        title="Canonical Nerve-Cell Marker Expression per Cluster",
    )
    for sub_ax in fig.axes:
        if sub_ax.get_xticklabels():
            current = [t.get_text() for t in sub_ax.get_xticklabels()]
            relabeled = [ensembl_to_symbol.get(c, c) for c in current]
            sub_ax.set_xticklabels(relabeled, rotation=90)
    fig.tight_layout()
    fig.savefig(snakemake.output.dotplot, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log_transformation(
        log,
        "nerve_cell_heterogeneity",
        f"Dot plot: {len(present_markers)} markers across {n_clusters} clusters",
    )
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(
        0.5,
        0.5,
        "No canonical markers present in dataset",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=12,
    )
    fig.savefig(snakemake.output.dotplot, dpi=100, bbox_inches="tight")
    plt.close(fig)
    log_transformation(
        log,
        "nerve_cell_heterogeneity",
        "WARNING: no canonical markers present — empty dot plot written",
        status="WARNING",
    )

# ---------------------------------------------------------------------------
# Step 4: Per-sample cluster abundance heatmap
# ---------------------------------------------------------------------------
abundance = (
    adata.obs.groupby(["sample_id", "nerve_leiden"])
    .size()
    .unstack(fill_value=0)
    .apply(lambda row: row / row.sum(), axis=1)  # proportions per sample
)

fig, ax = plt.subplots(
    figsize=(max(6, n_clusters * 0.8), max(5, abundance.shape[0] * 0.4))
)
sns.heatmap(
    abundance,
    ax=ax,
    cmap="YlOrRd",
    linewidths=0.3,
    linecolor="white",
    vmin=0,
    vmax=abundance.values.max(),
    annot=abundance.shape[0] <= 20,
    fmt=".2f",
)
ax.set_xlabel("Nerve-cell cluster (leiden)", fontsize=11)
ax.set_ylabel("Sample UUID", fontsize=11)
ax.set_title("Nerve Cell Cluster Abundance per Sample (proportions)", fontsize=12)

# Annotate with primary_diagnosis if available
if "primary_diagnosis" in adata.obs.columns:
    diag_map = adata.obs.groupby("sample_id")["primary_diagnosis"].first().to_dict()
    new_labels = [f"{sid} | {diag_map.get(sid, 'N/A')[:25]}" for sid in abundance.index]
    ax.set_yticklabels(new_labels, rotation=0, fontsize=7)

fig.tight_layout()
fig.savefig(snakemake.output.abundance, dpi=150, bbox_inches="tight")
plt.close(fig)

log_transformation(
    log,
    "nerve_cell_heterogeneity",
    f"Abundance heatmap: {abundance.shape[0]} samples × {abundance.shape[1]} clusters",
)

# --- FAIR provenance ---------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.markers,
    rule_name="nerve_cell_heterogeneity",
    input_paths=[snakemake.input.h5ad],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "n_clusters": n_clusters,
        "n_cells": adata.n_obs,
        "n_marker_genes_total": len(markers_df),
        "n_markers_resolved_to_symbol": int(len(markers_df) - n_unmapped),
        "n_markers_unmapped": int(n_unmapped),
        "n_dotplot_markers_present": len(present_markers),
        "n_dotplot_markers_unmapped": len(unmapped_markers),
        "gsea_ran": bool(enrichment_rows),
    },
    description="Nerve cell DE markers, GSEA enrichment, marker dot plot, and sample abundance heatmap",
    ontology_operation="operation:3223",  # EDAM: Differential gene expression profiling
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log,
    "nerve_cell_heterogeneity",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.markers,
        snakemake.output.enrichment,
        snakemake.output.dotplot,
        snakemake.output.abundance,
    ],
)
