"""Batch-correction sanity check for the nerve-cell subspace.

Question answered: did scVI integration produce biologically-driven clusters,
or do clusters reflect patient identity? Per-cluster purity is summarised in
`nerve_cluster_sample_purity.csv` (one row per Leiden cluster) and visualised
as a UMAP coloured by sample plus a small-multiples panel highlighting each
sample's contribution.

Pass criterion (default): dominant single-sample fraction < 0.5 AND ≥ 3
samples contribute ≥ 1% of cluster cells. Both thresholds are tunable in
config.yaml under `nerve_cells.batch_qc`.
"""

import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import (
    compute_cluster_purity,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

DOMINANT_FRACTION_MAX: float = float(snakemake.params.dominant_fraction_max)
MIN_CONTRIBUTING_FRACTION: float = float(snakemake.params.min_contributing_fraction)
MIN_CONTRIBUTING_SAMPLES: int = int(snakemake.params.min_contributing_samples)

log_transformation(
    log,
    "nerve_batch_qc",
    f"Loading {snakemake.input.h5ad}",
)
# Backed read: this rule only reads obs and obsm["X_umap"]. Backed mode loads
# obs/var/obsm/obsp/uns but leaves .X (and its .raw twin) on disk, which is the
# entire multi-GB payload here.
adata = ad.read_h5ad(snakemake.input.h5ad, backed="r")
if adata.n_obs == 0 or "nerve_leiden" not in adata.obs:
    raise RuntimeError("[FAIR-ALERT] nerve_cells.h5ad is empty or missing nerve_leiden")
if "X_umap" not in adata.obsm:
    raise RuntimeError("[FAIR-ALERT] nerve_cells.h5ad has no obsm['X_umap']")

n_cells = adata.n_obs
samples = sorted(adata.obs["sample_id"].astype(str).unique().tolist())
n_samples = len(samples)
clusters = sorted(
    adata.obs["nerve_leiden"].astype(str).unique().tolist(),
    key=lambda c: int(c),
)
n_clusters = len(clusters)
log_transformation(
    log,
    "nerve_batch_qc",
    f"{n_cells} cells across {n_clusters} clusters and {n_samples} samples",
)

# ---------------------------------------------------------------------------
# Per-cluster purity table (shared helper — identical logic also used by the
# nerve_leiden_resolution_sweep rule to keep the comparison apples-to-apples)
# ---------------------------------------------------------------------------
purity = compute_cluster_purity(
    adata,
    cluster_col="nerve_leiden",
    batch_col="sample_id",
    dominant_max=DOMINANT_FRACTION_MAX,
    min_contributing_fraction=MIN_CONTRIBUTING_FRACTION,
    min_contributing_samples=MIN_CONTRIBUTING_SAMPLES,
)
purity_df = purity.df
expected_uniform = purity.expected_uniform_entropy
purity_df.to_csv(snakemake.output.purity, index=False)
verify_artifact(snakemake.output.purity, min_size_bytes=64)

pass_overall = purity_df["pass_overall"].to_numpy()
pass_dominant = purity_df["pass_dominant"].to_numpy()
pass_diversity = purity_df["pass_diversity"].to_numpy()
n_pass = int(pass_overall.sum())
n_fail_dominant = int((~pass_dominant).sum())
n_fail_diversity = int((~pass_diversity).sum())
median_dom = float(np.median(purity_df["dominant_sample_fraction"].to_numpy()))
median_norm_entropy = float(np.median(purity_df["normalised_entropy"].to_numpy()))

log_transformation(
    log,
    "nerve_batch_qc",
    f"Purity verdict: {n_pass}/{n_clusters} clusters PASS "
    f"(dominant<{DOMINANT_FRACTION_MAX:g}, ≥{MIN_CONTRIBUTING_SAMPLES} samples contributing). "
    f"Median dominant fraction = {median_dom:.3f}; "
    f"median normalised entropy = {median_norm_entropy:.3f} of {expected_uniform:.3f} bits.",
    status="SUCCESS" if n_pass == n_clusters else "WARNING",
)
if n_fail_dominant or n_fail_diversity:
    failed_clusters = purity_df.loc[~pass_overall, ["cluster", "dominant_sample", "dominant_sample_fraction", "n_contributing_samples"]]
    log_transformation(
        log,
        "nerve_batch_qc",
        f"Failing clusters:\n{failed_clusters.to_string(index=False)}",
        status="WARNING",
    )

# ---------------------------------------------------------------------------
# UMAP coloured by sample_id (single panel) + small-multiples per sample
# ---------------------------------------------------------------------------
umap = adata.obsm["X_umap"]
sample_arr = adata.obs["sample_id"].astype(str).to_numpy()
cluster_arr = adata.obs["nerve_leiden"].astype(str).to_numpy()

# Stable, colour-blind-friendly palette per sample
cmap_samples = plt.get_cmap("tab20", n_samples)
sample_to_color = {s: cmap_samples(i) for i, s in enumerate(samples)}
sample_short = {s: s[:8] for s in samples}  # first 8 chars of the GDC UUID

fig_main, axes_main = plt.subplots(1, 2, figsize=(18, 8))
ax_left, ax_right = axes_main

# Left panel: cells colored by sample
for sid in samples:
    mask = sample_arr == sid
    ax_left.scatter(
        umap[mask, 0],
        umap[mask, 1],
        s=1.5,
        c=[sample_to_color[sid]],
        alpha=0.55,
        linewidths=0,
        label=sample_short[sid],
    )
ax_left.set_xlabel("UMAP 1")
ax_left.set_ylabel("UMAP 2")
ax_left.set_title("Nerve-cell UMAP — coloured by sample")
ax_left.legend(
    markerscale=4,
    fontsize=7,
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    title="sample (first 8 chars)",
    title_fontsize=8,
)

# Right panel: cells colored by cluster (sanity reference)
cmap_clusters = plt.get_cmap("tab20", n_clusters)
for ci, cl in enumerate(clusters):
    mask = cluster_arr == cl
    ax_right.scatter(
        umap[mask, 0],
        umap[mask, 1],
        s=1.5,
        c=[cmap_clusters(ci % 20)],
        alpha=0.55,
        linewidths=0,
    )
ax_right.set_xlabel("UMAP 1")
ax_right.set_ylabel("UMAP 2")
ax_right.set_title("Nerve-cell UMAP — coloured by Leiden cluster")
fig_main.tight_layout()
fig_main.savefig(snakemake.output.umap_main, dpi=140, bbox_inches="tight")
plt.close(fig_main)
verify_artifact(snakemake.output.umap_main, min_size_bytes=10_000)

# Small multiples — one panel per sample, that sample's cells coloured, others greyed
ncols = 5
nrows = int(np.ceil(n_samples / ncols))
fig_sm, axes_sm = plt.subplots(
    nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows), sharex=True, sharey=True
)
axes_flat = axes_sm.ravel()
xlim = (umap[:, 0].min() - 0.5, umap[:, 0].max() + 0.5)
ylim = (umap[:, 1].min() - 0.5, umap[:, 1].max() + 0.5)
for i, sid in enumerate(samples):
    ax = axes_flat[i]
    other = sample_arr != sid
    this_sample = ~other
    ax.scatter(
        umap[other, 0], umap[other, 1], s=0.6, c="#dadada", alpha=0.4, linewidths=0
    )
    ax.scatter(
        umap[this_sample, 0],
        umap[this_sample, 1],
        s=1.6,
        c=[sample_to_color[sid]],
        alpha=0.85,
        linewidths=0,
    )
    ax.set_title(
        f"{sample_short[sid]}  (n={int(this_sample.sum())})", fontsize=8
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks([])
    ax.set_yticks([])
for j in range(n_samples, len(axes_flat)):
    axes_flat[j].set_axis_off()
fig_sm.suptitle(
    "Per-sample contribution to the nerve-cell UMAP (others greyed). "
    "Patient-pure clusters indicate residual batch effects.",
    fontsize=11,
    y=1.005,
)
fig_sm.tight_layout()
fig_sm.savefig(snakemake.output.umap_per_sample, dpi=130, bbox_inches="tight")
plt.close(fig_sm)
verify_artifact(snakemake.output.umap_per_sample, min_size_bytes=10_000)

# ---------------------------------------------------------------------------
# FAIR provenance
# ---------------------------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.purity,
    rule_name="nerve_batch_qc",
    input_paths=[snakemake.input.h5ad],
    tool_versions={
        "anndata": ad.__version__,
        "matplotlib": matplotlib.__version__,
        "pandas": pd.__version__,
    },
    parameters={
        "n_cells": int(n_cells),
        "n_clusters": int(n_clusters),
        "n_samples": int(n_samples),
        "dominant_fraction_max": DOMINANT_FRACTION_MAX,
        "min_contributing_fraction": MIN_CONTRIBUTING_FRACTION,
        "min_contributing_samples": MIN_CONTRIBUTING_SAMPLES,
        "n_clusters_pass": n_pass,
        "n_clusters_fail_dominant": n_fail_dominant,
        "n_clusters_fail_diversity": n_fail_diversity,
        "median_dominant_fraction": round(median_dom, 4),
        "median_normalised_entropy": round(median_norm_entropy, 4),
    },
    description=(
        "Per-cluster sample-purity sanity check on the nerve-cell scVI integration "
        "(higher dominance / lower entropy = residual batch effect)."
    ),
    ontology_operation="operation:3186",  # EDAM: Statistical calculation
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log,
    "nerve_batch_qc",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.purity,
        snakemake.output.umap_main,
        snakemake.output.umap_per_sample,
    ],
)
