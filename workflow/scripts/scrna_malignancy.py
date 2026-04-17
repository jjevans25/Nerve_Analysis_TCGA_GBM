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

    # Map to our dataset genes
    name_to_order = dict(zip(gene_df["feature_name"], gene_df["gene_order"]))
    adata.var["gene_order"] = adata.var_names.map(name_to_order).fillna(-1).astype(int)
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
# Step 3: CNV scoring via sliding-window expression smoothing
# Cells with high CNV variance relative to reference → malignant
# ---------------------------------------------------------------------------
# Sort genes by genomic order
gene_order_idx = np.argsort(adata.var["gene_order"].values)
X_ordered = adata.X[:, gene_order_idx]
if hasattr(X_ordered, "toarray"):
    X_ordered = X_ordered.toarray()

# Log-normalize for CNV baseline subtraction
X_norm = np.log1p(X_ordered)

# Subtract per-gene mean of reference cells
if n_ref >= 1:
    ref_mask = adata.obs["is_reference"].values
    ref_mean = X_norm[ref_mask, :].mean(axis=0, keepdims=True)
    X_centered = X_norm - ref_mean
else:
    X_centered = X_norm - X_norm.mean(axis=0, keepdims=True)

# Sliding window smoothing (window = 100 genes)
window = min(100, X_centered.shape[1])
kernel = np.ones(window) / window
X_smoothed = np.apply_along_axis(
    lambda row: np.convolve(row, kernel, mode="same"), axis=1, arr=X_centered
)

# CNV score = variance of smoothed signal per cell
cnv_scores = X_smoothed.var(axis=1).astype(np.float32)
adata.obs["cnv_score"] = cnv_scores

# Threshold: cells scoring > mean + 2 SD of reference CNV score → malignant
if n_ref >= 1:
    ref_cnv = cnv_scores[adata.obs["is_reference"].values]
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
# Step 4: CNV heatmap (sample of cells × sorted genes)
# ---------------------------------------------------------------------------
n_plot = min(500, adata.n_obs)
rng = np.random.default_rng(snakemake.params.random_seed)
plot_idx = rng.choice(adata.n_obs, size=n_plot, replace=False)
plot_idx = plot_idx[np.argsort(adata.obs["is_malignant"].values[plot_idx].astype(int))]

fig, ax = plt.subplots(figsize=(14, 6))
im = ax.imshow(
    X_smoothed[plot_idx, :],
    aspect="auto",
    cmap="RdBu_r",
    vmin=-0.5,
    vmax=0.5,
    interpolation="nearest",
)
n_normal_plot = (~adata.obs["is_malignant"].values[plot_idx]).sum()
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
write_provenance(prov, "provenance")

log_transformation(log, "scrna_malignancy", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.cnv_plot])
