"""Batch-correction sanity check on the v2 scANVI latent space.

Mirrors ``nerve_batch_qc.py`` exactly but operates on ``X_scANVI`` and
clusters via Leiden at the v1.0.0 resolution (1.0) so the v1 vs v2
comparison is apples-to-apples. Uses the shared
``fair_utils.compute_cluster_purity`` helper — identical purity criteria.
"""

import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import (
    H5AD_COMPRESSION,
    compute_cluster_purity,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]
SEED = int(snakemake.params.random_seed)
os.environ["PYTHONHASHSEED"] = str(SEED)

DOMINANT_FRACTION_MAX = float(snakemake.params.dominant_fraction_max)
MIN_CONTRIBUTING_FRACTION = float(snakemake.params.min_contributing_fraction)
MIN_CONTRIBUTING_SAMPLES = int(snakemake.params.min_contributing_samples)
LEIDEN_RESOLUTION = float(snakemake.params.leiden_resolution)

log_transformation(log, "nerve_batch_qc_v2", f"Loading {snakemake.input.h5ad}")
adata = ad.read_h5ad(snakemake.input.h5ad)
if "X_scANVI" not in adata.obsm:
    raise RuntimeError("[FAIR-ALERT] input is missing obsm['X_scANVI']")
if "neighbors" not in adata.uns:
    log_transformation(log, "nerve_batch_qc_v2", "Recomputing neighbors on X_scANVI")
    sc.pp.neighbors(adata, use_rep="X_scANVI", random_state=SEED)
if "X_umap" not in adata.obsm:
    log_transformation(log, "nerve_batch_qc_v2", "Recomputing UMAP")
    sc.tl.umap(adata, random_state=SEED)

log_transformation(log, "nerve_batch_qc_v2",
                   f"Running Leiden at resolution={LEIDEN_RESOLUTION}")
sc.tl.leiden(
    adata, resolution=LEIDEN_RESOLUTION, random_state=SEED,
    key_added="nerve_leiden_v2",
)

# Persist the cluster column back onto the h5ad on disk so downstream rules
# can reuse it.
adata.write_h5ad(snakemake.input.h5ad, compression=H5AD_COMPRESSION)

n_cells = adata.n_obs
samples = sorted(adata.obs["sample_id"].astype(str).unique().tolist())
n_samples = len(samples)
log_transformation(log, "nerve_batch_qc_v2",
                   f"{n_cells} cells across {n_samples} samples")

# ----------------------------------------------------------------------------
# Per-cluster purity (shared helper)
# ----------------------------------------------------------------------------
purity = compute_cluster_purity(
    adata,
    cluster_col="nerve_leiden_v2",
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
n_clusters = int(len(purity_df))
n_pass = int(pass_overall.sum())
n_fail_dominant = int((~pass_dominant).sum())
n_fail_diversity = int((~pass_diversity).sum())
median_dom = float(np.median(purity_df["dominant_sample_fraction"].to_numpy()))
median_norm_entropy = float(np.median(purity_df["normalised_entropy"].to_numpy()))

log_transformation(
    log, "nerve_batch_qc_v2",
    f"v2 verdict: {n_pass}/{n_clusters} clusters PASS "
    f"(dominant<{DOMINANT_FRACTION_MAX:g}, >={MIN_CONTRIBUTING_SAMPLES} samples). "
    f"Median dominant fraction = {median_dom:.3f}; "
    f"median normalised entropy = {median_norm_entropy:.3f} of {expected_uniform:.3f} bits.",
    status="SUCCESS" if n_pass == n_clusters else "WARNING",
)
if n_fail_dominant or n_fail_diversity:
    failed = purity_df.loc[~pass_overall, ["cluster", "dominant_sample", "dominant_sample_fraction", "n_contributing_samples"]]
    log_transformation(log, "nerve_batch_qc_v2",
                       f"Failing clusters:\n{failed.to_string(index=False)}",
                       status="WARNING")

# ----------------------------------------------------------------------------
# UMAP figures — same layout as v1 batch_qc
# ----------------------------------------------------------------------------
umap = adata.obsm["X_umap"]
sample_arr = adata.obs["sample_id"].astype(str).to_numpy()
cluster_arr = adata.obs["nerve_leiden_v2"].astype(str).to_numpy()
clusters = sorted(np.unique(cluster_arr).tolist(), key=lambda c: int(c))

cmap_samples = plt.get_cmap("tab20", n_samples)
sample_to_color = {s: cmap_samples(i) for i, s in enumerate(samples)}
sample_short = {s: s[:8] for s in samples}

fig_main, axes_main = plt.subplots(1, 2, figsize=(18, 8))
ax_left, ax_right = axes_main
for sid in samples:
    mask = sample_arr == sid
    ax_left.scatter(umap[mask, 0], umap[mask, 1], s=1.5,
                    c=[sample_to_color[sid]], alpha=0.55, linewidths=0,
                    label=sample_short[sid])
ax_left.set_xlabel("UMAP 1"); ax_left.set_ylabel("UMAP 2")
ax_left.set_title("v2 nerve-cell UMAP (X_scANVI) — coloured by sample")
ax_left.legend(markerscale=4, fontsize=7, bbox_to_anchor=(1.02, 1),
               loc="upper left", title="sample", title_fontsize=8)

cmap_clusters = plt.get_cmap("tab20", max(20, len(clusters)))
for ci, cl in enumerate(clusters):
    mask = cluster_arr == cl
    ax_right.scatter(umap[mask, 0], umap[mask, 1], s=1.5,
                     c=[cmap_clusters(ci % 20)], alpha=0.55, linewidths=0)
ax_right.set_xlabel("UMAP 1"); ax_right.set_ylabel("UMAP 2")
ax_right.set_title(f"v2 nerve-cell UMAP — coloured by Leiden (r={LEIDEN_RESOLUTION})")
fig_main.tight_layout()
fig_main.savefig(snakemake.output.umap_main, dpi=140, bbox_inches="tight")
plt.close(fig_main)
verify_artifact(snakemake.output.umap_main, min_size_bytes=10_000)

# Small multiples per sample.
ncols = 5
nrows = int(np.ceil(n_samples / ncols))
fig_sm, axes_sm = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows),
                                sharex=True, sharey=True)
axes_flat = axes_sm.ravel()
xlim = (umap[:, 0].min() - 0.5, umap[:, 0].max() + 0.5)
ylim = (umap[:, 1].min() - 0.5, umap[:, 1].max() + 0.5)
for i, sid in enumerate(samples):
    ax = axes_flat[i]
    this_sample = sample_arr == sid
    other = ~this_sample
    ax.scatter(umap[other, 0], umap[other, 1], s=0.6, c="#dadada", alpha=0.4, linewidths=0)
    ax.scatter(umap[this_sample, 0], umap[this_sample, 1], s=1.6,
               c=[sample_to_color[sid]], alpha=0.85, linewidths=0)
    ax.set_title(f"{sample_short[sid]}  (n={int(this_sample.sum())})", fontsize=8)
    ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xticks([]); ax.set_yticks([])
for j in range(n_samples, len(axes_flat)):
    axes_flat[j].set_axis_off()
fig_sm.suptitle(
    "v2 per-sample contribution to X_scANVI UMAP (others greyed)",
    fontsize=11, y=1.005,
)
fig_sm.tight_layout()
fig_sm.savefig(snakemake.output.umap_per_sample, dpi=130, bbox_inches="tight")
plt.close(fig_sm)
verify_artifact(snakemake.output.umap_per_sample, min_size_bytes=10_000)

# ----------------------------------------------------------------------------
# Provenance
# ----------------------------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.purity,
    rule_name="nerve_batch_qc_v2",
    input_paths=[snakemake.input.h5ad],
    tool_versions={
        "anndata": ad.__version__,
        "scanpy": sc.__version__,
        "matplotlib": matplotlib.__version__,
        "pandas": pd.__version__,
    },
    parameters={
        "n_cells": int(n_cells),
        "n_clusters": n_clusters,
        "n_samples": int(n_samples),
        "leiden_resolution": LEIDEN_RESOLUTION,
        "dominant_fraction_max": DOMINANT_FRACTION_MAX,
        "min_contributing_fraction": MIN_CONTRIBUTING_FRACTION,
        "min_contributing_samples": MIN_CONTRIBUTING_SAMPLES,
        "n_clusters_pass": n_pass,
        "n_clusters_fail_dominant": n_fail_dominant,
        "n_clusters_fail_diversity": n_fail_diversity,
        "median_dominant_fraction": round(median_dom, 4),
        "median_normalised_entropy": round(median_norm_entropy, 4),
        "pass_fraction": round(n_pass / max(n_clusters, 1), 4),
        "random_seed": SEED,
    },
    description=(
        "Per-cluster sample-purity verdict on the v2 scANVI nerve-cell latent "
        "space. Identical purity thresholds to v1.0.0 nerve_batch_qc so the "
        "v1 vs v2 comparison is apples-to-apples."
    ),
    ontology_operation="operation:3186",  # EDAM: Statistical calculation
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log, "nerve_batch_qc_v2",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.purity,
        snakemake.output.umap_main,
        snakemake.output.umap_per_sample,
    ],
)
