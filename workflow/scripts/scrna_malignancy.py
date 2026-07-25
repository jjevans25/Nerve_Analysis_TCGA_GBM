"""CNV-based tumor/normal separation using infercnvpy; labels each cell with is_malignant."""

import os
import sys
import warnings

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

log_transformation(log, "scrna_malignancy", f"Loading {snakemake.input.h5ad}")
adata = ad.read_h5ad(snakemake.input.h5ad)

# ---------------------------------------------------------------------------
# Step 1: Fetch gene chromosomal positions from CELLxGENE Census
# Needed by infercnvpy; tolerates genes absent from Census (marked NA).
# ---------------------------------------------------------------------------
log_transformation(log, "scrna_malignancy", "Fetching gene positions from CELLxGENE Census")

try:
    import cellxgene_census
    with cellxgene_census.open_soma(census_version=snakemake.params.census_version) as census:
        gene_df = (
            census["census_data"]["homo_sapiens"]["ms"]["RNA"]["var"]
            .read(column_names=["feature_name", "feature_id"])
            .concat()
            .to_pandas()
        )
    # Census doesn't expose chromosome positions directly; use feature_id (Ensembl) for ordering
    # Build a rough chromosomal order via Ensembl ID numeric suffix as a proxy
    gene_df["ensembl_numeric"] = (
        gene_df["feature_id"].str.extract(r"ENSG(\d+)")[0]
        .fillna("0")
        .astype(int)
    )
    gene_df = gene_df.sort_values("ensembl_numeric").reset_index(drop=True)
    gene_df["gene_order"] = range(len(gene_df))

    # Map to our dataset genes by Ensembl ID. var_names here are Ensembl IDs
    # (ENSG…), so join on feature_id, not feature_name (gene symbols) — matching
    # symbols against Ensembl IDs orders almost nothing.
    id_to_order = dict(zip(gene_df["feature_id"], gene_df["gene_order"]))
    gene_ids = adata.var["feature_id"] if "feature_id" in adata.var.columns \
        else pd.Series(adata.var_names, index=adata.var_names)
    adata.var["gene_order"] = gene_ids.map(id_to_order).fillna(-1).astype(int).values
    n_ordered = (adata.var["gene_order"] >= 0).sum()
    log_transformation(log, "scrna_malignancy",
        f"Mapped {n_ordered}/{adata.n_vars} genes to Ensembl genomic order")
except Exception as exc:
    warnings.warn(f"[FAIR-ALERT] Census gene order fetch failed: {exc}. Using alphabetical ordering.")
    log_transformation(log, "scrna_malignancy",
        f"WARNING: Census failed ({exc}); falling back to alphabetical gene order", status="WARNING")
    adata.var["gene_order"] = pd.Series(
        range(adata.n_vars), index=adata.var_names
    )

# ---------------------------------------------------------------------------
# Step 2: Identify reference cells (non-glial: T cells + endothelial)
# ---------------------------------------------------------------------------
reference_types = {"t_cell", "endothelial"}
adata.obs["is_reference"] = adata.obs["cell_type_predicted"].isin(reference_types)
n_ref = adata.obs["is_reference"].sum()
log_transformation(log, "scrna_malignancy",
    f"Reference cells (T cells + endothelial): {n_ref}")

if n_ref < 10:
    warnings.warn("[FAIR-ALERT] Fewer than 10 reference cells — CNV scoring may be unreliable.")
    log_transformation(log, "scrna_malignancy",
        "WARNING: very few reference cells for CNV baseline", status="WARNING")

# ---------------------------------------------------------------------------
# Step 3: CNV scoring via sliding-window expression smoothing.
# Streamed over cells in row-blocks to bound peak memory — every operation
# here is row-wise and independent per cell; the only cross-cell quantity is
# the per-gene reference mean, computed once. Chunking is numerically identical
# to a whole-matrix pass but never materializes the full dense matrix (which is
# ~59 GB for a 615k×24k float32 cohort, and OOM-kills at whole-matrix scale).
# Cells with high CNV variance relative to reference → malignant.
# ---------------------------------------------------------------------------
gene_order_idx = np.argsort(adata.var["gene_order"].values)
n_obs = adata.n_obs
n_genes_ord = adata.n_vars
ref_mask = adata.obs["is_reference"].values
chunk_size = int(snakemake.params.cnv_chunk_size)

# Sliding-window params (window centered per gene)
window = min(100, n_genes_ord)
pad_l, pad_r = window // 2, window - window // 2 - 1


def _load_ordered_block(lo, hi):
    """Densify cells [lo:hi), reorder genes to genomic order, log-normalize."""
    blk = adata.X[lo:hi][:, gene_order_idx].toarray().astype(np.float32, copy=False)
    np.log1p(blk, out=blk)
    return blk


# Per-gene reference baseline (mean of log-normalized reference cells)
if n_ref >= 1:
    ref_block = adata.X[ref_mask][:, gene_order_idx].toarray().astype(np.float32, copy=False)
    np.log1p(ref_block, out=ref_block)
    ref_mean = ref_block.mean(axis=0, keepdims=True)
    del ref_block
else:
    # No reference cells: baseline = mean over all cells (streamed).
    gene_sum = np.zeros((1, n_genes_ord), dtype=np.float64)
    for lo in range(0, n_obs, chunk_size):
        blk = _load_ordered_block(lo, min(lo + chunk_size, n_obs))
        gene_sum += blk.sum(axis=0, keepdims=True)
        del blk
    ref_mean = (gene_sum / n_obs).astype(np.float32)
    del gene_sum

# Pre-select heatmap cells before the loop so only 500 smoothed rows are kept.
n_plot = min(500, n_obs)
rng = np.random.default_rng(snakemake.params.random_seed)
plot_idx = np.sort(rng.choice(n_obs, size=n_plot, replace=False))
plot_buffer = np.empty((n_plot, n_genes_ord), dtype=np.float32)

# Stream CNV scores in row-blocks.
cnv_scores = np.empty(n_obs, dtype=np.float32)
for lo in range(0, n_obs, chunk_size):
    hi = min(lo + chunk_size, n_obs)
    blk = _load_ordered_block(lo, hi)
    blk -= ref_mean

    # Memory-efficient sliding window via cumsum: sum[i] = cs[i+w] - cs[i].
    padded = np.pad(blk, ((0, 0), (pad_l, pad_r)), mode="edge")
    del blk
    cs = np.empty((padded.shape[0], padded.shape[1] + 1), dtype=np.float32)
    cs[:, 0] = 0.0
    np.cumsum(padded, axis=1, out=cs[:, 1:])
    del padded
    smoothed = (cs[:, window : window + n_genes_ord] - cs[:, :n_genes_ord]) / window
    del cs

    cnv_scores[lo:hi] = smoothed.var(axis=1)

    # Retain smoothed rows for any pre-selected heatmap cells in this block.
    sel = plot_idx[(plot_idx >= lo) & (plot_idx < hi)]
    if sel.size:
        plot_buffer[np.searchsorted(plot_idx, sel)] = smoothed[sel - lo]
    del smoothed

adata.obs["cnv_score"] = cnv_scores

# Threshold: cells scoring > mean + 2 SD of reference CNV score → malignant
if n_ref >= 1:
    ref_cnv = cnv_scores[ref_mask]
    threshold = ref_cnv.mean() + 2 * ref_cnv.std()
else:
    threshold = np.percentile(cnv_scores, 75)

adata.obs["is_malignant"] = (cnv_scores > threshold).astype(bool)
n_malignant = adata.obs["is_malignant"].sum()

log_transformation(log, "scrna_malignancy",
    f"CNV threshold={threshold:.4f}; "
    f"{n_malignant}/{adata.n_obs} cells flagged as malignant "
    f"({100 * n_malignant / adata.n_obs:.1f}%)")

# ---------------------------------------------------------------------------
# Step 4: CNV heatmap (sampled cells × sorted genes), sorted by malignant flag
# ---------------------------------------------------------------------------
malignant_plot = adata.obs["is_malignant"].values[plot_idx]
order = np.argsort(malignant_plot.astype(int))
plot_buffer = plot_buffer[order]
n_normal_plot = int((~malignant_plot).sum())

fig, ax = plt.subplots(figsize=(14, 6))
im = ax.imshow(
    plot_buffer,
    aspect="auto",
    cmap="RdBu_r",
    vmin=-0.5,
    vmax=0.5,
    interpolation="nearest",
)
ax.axhline(n_normal_plot - 0.5, color="black", linewidth=1.5, linestyle="--")
ax.set_xlabel("Genes (genomic order)")
ax.set_ylabel(f"Cells (n={n_plot}; dashed = malignant boundary)")
ax.set_title("CNV Inference Heatmap — Smoothed Expression Signal")
plt.colorbar(im, ax=ax, label="Centered log-normalized expression")
fig.tight_layout()
fig.savefig(snakemake.output.cnv_plot, dpi=150, bbox_inches="tight")
plt.close(fig)

# --- Write output ------------------------------------------------------------
adata.write_h5ad(snakemake.output.h5ad)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="scrna_malignancy",
    input_paths=[snakemake.input.h5ad],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__, "numpy": np.__version__},
    parameters={
        "cnv_window_size":    window,
        "cnv_threshold":      float(threshold),
        "n_reference_cells":  int(n_ref),
        "n_malignant":        int(n_malignant),
        "n_cells":            adata.n_obs,
        "pct_malignant":      round(100 * n_malignant / adata.n_obs, 2),
    },
    description="CNV-scored AnnData with is_malignant label; sliding-window expression smoothing",
    ontology_operation="operation:3225",  # EDAM: Copy number variation detection
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "scrna_malignancy", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.cnv_plot])
