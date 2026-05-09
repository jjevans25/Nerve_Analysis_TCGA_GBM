"""Differential expression, GSEA, marker dot plot, and per-sample cluster abundance."""

import os
import sys

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

_TOP_N_MARKERS_PER_CLUSTER = 50

# Rank ALL genes per cluster: prerank GSEA needs the full ranking, while the
# exported markers CSV is truncated to the top-N significant genes (preserves
# prior CSV semantics).
sc.tl.rank_genes_groups(
    adata_de,
    groupby="nerve_leiden",
    method="wilcoxon",
    n_genes=None,
    key_added="rank_genes_nerve",
)

# Flatten results into two DataFrames:
#   markers_full_df → all genes per cluster (drives prerank rnk construction)
#   markers_df      → top-50 significant genes per cluster (written to CSV)
full_rows: list[pd.DataFrame] = []
marker_rows: list[pd.DataFrame] = []
for group in adata_de.obs["nerve_leiden"].unique():
    full = sc.get.rank_genes_groups_df(
        adata_de, group=group, key="rank_genes_nerve"
    )
    full.insert(0, "cluster", group)
    full_rows.append(full)

    sig = full[full["pvals_adj"] < 0.05].head(_TOP_N_MARKERS_PER_CLUSTER)
    marker_rows.append(sig)

_full_cols = ["cluster", "names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
markers_full_df = (
    pd.concat(full_rows, ignore_index=True) if full_rows else pd.DataFrame(columns=_full_cols)
)
markers_df = (
    pd.concat(marker_rows, ignore_index=True) if marker_rows else pd.DataFrame(columns=_full_cols)
)

# Annotate symbol on the full table too — used by prerank.
markers_full_df["gene_symbol"] = markers_full_df["names"].map(ensembl_to_symbol)

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
# Step 2: Offline GSEA via gseapy.prerank against local MSigDB GMTs
# (replaces Enrichr API — gap §3.2 in markdowns/next_steps_interpretation.md).
# Ranking metric: Wilcoxon z-score (`scores`) from rank_genes_groups, oriented
# so up-regulation = high. Genes are mapped Ensembl→HGNC symbol because MSigDB
# .gmt collections are keyed on symbols.
# ---------------------------------------------------------------------------
import gseapy as gp


def _read_gmt_set_sizes(gmt_path: str) -> dict[str, int]:
    """Map term name → number of genes in the set (used to format Overlap)."""
    sizes: dict[str, int] = {}
    with open(gmt_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            sizes[fields[0]] = len(fields[2:])
    return sizes


def _build_rnk(cluster_df: pd.DataFrame) -> pd.DataFrame:
    """Convert per-cluster markers into a 2-column [gene_symbol, score] rnk.

    Drops rows with no gene_symbol mapping. On duplicate symbols (rare, from
    Ensembl version collisions) keeps the entry with largest |score|.
    """
    sub = cluster_df.dropna(subset=["gene_symbol"]).copy()
    sub["abs_score"] = sub["scores"].abs()
    sub = (
        sub.sort_values("abs_score", ascending=False)
        .drop_duplicates(subset="gene_symbol", keep="first")
        .sort_values("scores", ascending=False)
    )
    return sub[["gene_symbol", "scores"]].rename(columns={"scores": "score"})


prerank_cfg = dict(snakemake.params.prerank)
top_n = int(prerank_cfg["top_n_per_cluster"])
min_size = int(prerank_cfg["min_size"])
max_size = int(prerank_cfg["max_size"])
permutation_num = int(prerank_cfg["permutation_num"])
min_ranked = int(prerank_cfg["min_ranked_genes"])

gmt_targets = [
    (snakemake.input.gmt_bp, snakemake.params.bp_label),
    (snakemake.input.gmt_mf, snakemake.params.mf_label),
]
set_size_lookup: dict[str, dict[str, int]] = {
    label: _read_gmt_set_sizes(path) for path, label in gmt_targets
}

enrichment_rows: list[pd.DataFrame] = []
n_clusters_skipped = 0
n_prerank_failures = 0

for cluster in markers_full_df["cluster"].unique():
    rnk = _build_rnk(markers_full_df[markers_full_df["cluster"] == cluster])
    if len(rnk) < min_ranked:
        n_clusters_skipped += 1
        log_transformation(
            log,
            "nerve_cell_heterogeneity",
            f"Skipping cluster {cluster}: only {len(rnk)} ranked symbols "
            f"(<{min_ranked} required)",
            status="WARNING",
        )
        continue

    for gmt_path, library_label in gmt_targets:
        try:
            pre = gp.prerank(
                rnk=rnk,
                gene_sets=str(gmt_path),
                threads=int(snakemake.threads),
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutation_num,
                seed=int(snakemake.params.random_seed),
                outdir=None,
                verbose=False,
            )
        except (ValueError, KeyError, RuntimeError) as exc:
            n_prerank_failures += 1
            log_transformation(
                log,
                "nerve_cell_heterogeneity",
                f"WARNING: prerank failed for cluster {cluster}, {library_label}: {exc}",
                status="WARNING",
            )
            continue

        res = getattr(pre, "res2d", None)
        if res is None or len(res) == 0:
            continue

        sizes = set_size_lookup[library_label]
        res = res.copy()
        res["__lead_size"] = res["Lead_genes"].fillna("").map(
            lambda s: 0 if not s else len([g for g in s.split(";") if g])
        )
        res["__set_size"] = res["Term"].map(sizes).fillna(0).astype(int)
        # Sort by FDR q-val ascending, then keep top-N for parity with prior
        # `enr.results.head(10)` slice.
        res = res.sort_values("FDR q-val", ascending=True).head(top_n)

        out = pd.DataFrame(
            {
                "cluster": cluster,
                "gene_set_library": library_label,
                "Term": res["Term"].values,
                "Adjusted P-value": res["FDR q-val"].values,
                "Overlap": [
                    f"{int(ls)}/{int(ss)}" if ss > 0 else f"{int(ls)}/NA"
                    for ls, ss in zip(res["__lead_size"], res["__set_size"])
                ],
            }
        )
        enrichment_rows.append(out)

if n_prerank_failures:
    log_transformation(
        log,
        "nerve_cell_heterogeneity",
        f"prerank had {n_prerank_failures} per-(cluster, library) failures "
        f"(see WARNING entries above)",
        status="WARNING",
    )
log_transformation(
    log,
    "nerve_cell_heterogeneity",
    f"GSEA complete: {len(enrichment_rows)} result blocks "
    f"({n_clusters_skipped} clusters skipped for insufficient ranked symbols)",
)

enrichment_df = (
    pd.concat(enrichment_rows, ignore_index=True) if enrichment_rows else pd.DataFrame(
        columns=["cluster", "gene_set_library", "Term", "Adjusted P-value", "Overlap"]
    )
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
    tool_versions={
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "gseapy": gp.__version__,
    },
    parameters={
        "n_clusters": n_clusters,
        "n_cells": adata.n_obs,
        "n_marker_genes_total": len(markers_df),
        "n_markers_resolved_to_symbol": int(len(markers_df) - n_unmapped),
        "n_markers_unmapped": int(n_unmapped),
        "n_dotplot_markers_present": len(present_markers),
        "n_dotplot_markers_unmapped": len(unmapped_markers),
        "gsea_ran": bool(enrichment_rows),
        "gsea_engine": "gseapy.prerank",
        "msigdb_release": snakemake.params.msigdb_release,
        "gsea_libraries": [
            snakemake.params.bp_label,
            snakemake.params.mf_label,
        ],
        "gsea_gmt_files": [
            str(snakemake.input.gmt_bp),
            str(snakemake.input.gmt_mf),
        ],
        "prerank_params": {
            "min_size": min_size,
            "max_size": max_size,
            "permutation_num": permutation_num,
            "top_n_per_cluster": top_n,
            "min_ranked_genes": min_ranked,
            "seed": int(snakemake.params.random_seed),
        },
        "n_prerank_failures": n_prerank_failures,
        "n_clusters_skipped_low_symbols": n_clusters_skipped,
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
