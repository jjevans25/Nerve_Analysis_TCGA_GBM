"""Sweep ``leiden_resolution`` to test whether coarser clusters reduce batch-QC failures.

Question answered: at v1.0.0 (resolution 1.0) Leiden produced 24 clusters and 6/24
failed batch purity, all in the small-cluster tail — a diagnostic signature of
over-splitting. This sweep re-cuts the *same* scVI kNN graph at multiple
resolutions and reports purity + silhouette per resolution, so the researcher
can pick a coarser cut without re-training scVI.

Read-only with respect to ``nerve_cells.h5ad`` — operates on an in-memory copy.
Inputs: ``data/processed/nerve_cells.h5ad`` (already carries ``X_scVI`` and the kNN graph).
Outputs: sweep CSV + 3-panel PNG + provenance JSON.
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
import sklearn
from sklearn.metrics import silhouette_score

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


def _crash_hook(exc_type, exc_value, exc_tb):
    """Persist any uncaught exception into snakemake.log[0] before exiting."""
    import traceback
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as f:
        f.write("=== UNCAUGHT EXCEPTION ===\n")
        traceback.print_exception(exc_type, exc_value, exc_tb, file=f)
        f.write("===========================\n")
    sys.__excepthook__(exc_type, exc_value, exc_tb)


sys.excepthook = _crash_hook

seed = int(snakemake.params.random_seed)
os.environ["PYTHONHASHSEED"] = str(seed)
np.random.seed(seed)

resolutions: list[float] = [float(r) for r in snakemake.params.resolutions]
dominant_max = float(snakemake.params.dominant_fraction_max)
min_contributing_fraction = float(snakemake.params.min_contributing_fraction)
min_contributing_samples = int(snakemake.params.min_contributing_samples)
silhouette_sample_size = int(getattr(snakemake.params, "silhouette_sample_size", 5000))

log_transformation(log, "nerve_leiden_resolution_sweep", f"Loading {snakemake.input.h5ad}")
adata = ad.read_h5ad(snakemake.input.h5ad)
if adata.n_obs == 0 or "X_scVI" not in adata.obsm:
    raise RuntimeError("[FAIR-ALERT] nerve_cells.h5ad is empty or missing obsm['X_scVI']")
if "neighbors" not in adata.uns:
    raise RuntimeError("[FAIR-ALERT] nerve_cells.h5ad has no neighbors graph; cannot run Leiden")
if "sample_id" not in adata.obs:
    raise RuntimeError("[FAIR-ALERT] nerve_cells.h5ad has no obs['sample_id']")

n_cells = adata.n_obs
X_scvi = np.asarray(adata.obsm["X_scVI"])
log_transformation(
    log, "nerve_leiden_resolution_sweep",
    f"{n_cells} cells; sweeping {len(resolutions)} resolutions: {resolutions}",
)

# Pre-sample a fixed subset for silhouette so comparisons across resolutions are consistent.
rng = np.random.default_rng(seed)
if n_cells > silhouette_sample_size:
    sil_idx = rng.choice(n_cells, size=silhouette_sample_size, replace=False)
else:
    sil_idx = np.arange(n_cells)
X_sil = X_scvi[sil_idx]

rows: list[dict] = []
for r in resolutions:
    key = f"leiden_r{r:g}"
    log_transformation(log, "nerve_leiden_resolution_sweep", f"Leiden @ resolution={r}")
    # Reuses adata.uns['neighbors'] / .obsp['connectivities'] — no re-embedding.
    sc.tl.leiden(adata, resolution=r, random_state=seed, key_added=key)

    purity = compute_cluster_purity(
        adata,
        cluster_col=key,
        batch_col="sample_id",
        dominant_max=dominant_max,
        min_contributing_fraction=min_contributing_fraction,
        min_contributing_samples=min_contributing_samples,
    )
    df = purity.df
    n_clusters = int(len(df))
    cluster_sizes = df["n_cells"].to_numpy()

    labels_full = adata.obs[key].astype(str).to_numpy()
    labels_sil = labels_full[sil_idx]
    n_unique_sil = len(np.unique(labels_sil))
    if 2 <= n_unique_sil < len(labels_sil):
        sil = float(silhouette_score(X_sil, labels_sil, random_state=seed))
    else:
        sil = float("nan")

    rows.append(
        {
            "resolution": r,
            "n_clusters": n_clusters,
            "n_clusters_pass": int(df["pass_overall"].sum()),
            "n_clusters_fail_dominant": int((~df["pass_dominant"]).sum()),
            "n_clusters_fail_diversity": int((~df["pass_diversity"]).sum()),
            "pass_fraction": round(float(df["pass_overall"].mean()), 4),
            "median_dominant_fraction": round(float(df["dominant_sample_fraction"].median()), 4),
            "median_normalised_entropy": round(float(df["normalised_entropy"].median()), 4),
            "mean_cluster_size": round(float(cluster_sizes.mean()), 2),
            "min_cluster_size": int(cluster_sizes.min()),
            "n_clusters_under_500_cells": int((cluster_sizes < 500).sum()),
            "silhouette_scvi": round(sil, 4) if not np.isnan(sil) else float("nan"),
        }
    )

sweep_df = pd.DataFrame(rows).sort_values("resolution").reset_index(drop=True)
sweep_df.to_csv(snakemake.output.csv, index=False)
verify_artifact(snakemake.output.csv, min_size_bytes=64)
log_transformation(
    log, "nerve_leiden_resolution_sweep",
    f"Sweep table:\n{sweep_df.to_string(index=False)}",
)

# Sanity check: row at r=1.0 should reproduce v1.0.0 (24 clusters, 18 pass).
if 1.0 in resolutions:
    row_1 = sweep_df.loc[sweep_df["resolution"] == 1.0].iloc[0]
    status = "SUCCESS" if (row_1["n_clusters"] == 24 and row_1["n_clusters_pass"] == 18) else "WARNING"
    log_transformation(
        log, "nerve_leiden_resolution_sweep",
        f"Sanity @ r=1.0: n_clusters={int(row_1['n_clusters'])} n_pass={int(row_1['n_clusters_pass'])} "
        f"(expected 24 / 18 from v1.0.0).",
        status=status,
    )

# ---------------------------------------------------------------------------
# Three-panel figure
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
ax_a, ax_b, ax_c = axes

# (a) cluster counts vs resolution
ax_a.plot(sweep_df["resolution"], sweep_df["n_clusters"], "o-", color="#1f77b4", label="n_clusters")
ax_a.plot(sweep_df["resolution"], sweep_df["n_clusters_pass"], "s--", color="#2ca02c", label="n_pass")
ax_a.set_xlabel("leiden_resolution")
ax_a.set_ylabel("count")
ax_a.set_title("(a) Cluster count and pass count")
ax_a.legend(fontsize=9)
ax_a.grid(alpha=0.3)

# (b) purity metrics
ax_b.plot(sweep_df["resolution"], sweep_df["median_dominant_fraction"], "o-", color="#d62728",
          label="median dominant fraction")
ax_b.plot(sweep_df["resolution"], sweep_df["median_normalised_entropy"], "s--", color="#9467bd",
          label="median normalised entropy")
ax_b.axhline(dominant_max, color="#d62728", ls=":", lw=1, alpha=0.6, label=f"dominance threshold = {dominant_max}")
ax_b.set_xlabel("leiden_resolution")
ax_b.set_ylabel("value (lower dom = more diverse)")
ax_b.set_title("(b) Purity metrics")
ax_b.legend(fontsize=8)
ax_b.grid(alpha=0.3)

# (c) silhouette
ax_c.plot(sweep_df["resolution"], sweep_df["silhouette_scvi"], "o-", color="#ff7f0e")
ax_c.set_xlabel("leiden_resolution")
ax_c.set_ylabel("silhouette (on X_scVI, n=5,000 sampled)")
ax_c.set_title("(c) Cluster separability in scVI latent space")
ax_c.grid(alpha=0.3)

fig.suptitle(
    "Leiden resolution sweep — nerve-cell subspace (read-only; scVI integration unchanged)",
    fontsize=11, y=1.02,
)
fig.tight_layout()
fig.savefig(snakemake.output.figure, dpi=140, bbox_inches="tight")
plt.close(fig)
verify_artifact(snakemake.output.figure, min_size_bytes=10_000)

# ---------------------------------------------------------------------------
# FAIR provenance
# ---------------------------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.csv,
    rule_name="nerve_leiden_resolution_sweep",
    input_paths=[snakemake.input.h5ad],
    tool_versions={
        "anndata": ad.__version__,
        "scanpy": sc.__version__,
        "scikit-learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "pandas": pd.__version__,
        "numpy": np.__version__,
    },
    parameters={
        "n_cells": int(n_cells),
        "resolutions": resolutions,
        "random_seed": seed,
        "dominant_fraction_max": dominant_max,
        "min_contributing_fraction": min_contributing_fraction,
        "min_contributing_samples": min_contributing_samples,
        "silhouette_sample_size": int(silhouette_sample_size),
        "n_resolutions": len(resolutions),
    },
    description=(
        "Sweep of Leiden resolution on the nerve-cell scVI kNN graph, evaluating "
        "per-cluster batch purity (dominance + Shannon diversity) and "
        "scVI-space silhouette to assess whether a coarser cut reduces "
        "patient-dominated clusters without collapsing biological structure."
    ),
    ontology_operation="operation:3432",  # EDAM: Clustering
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log, "nerve_leiden_resolution_sweep",
    "Complete",
    status="SUCCESS",
    artifact_paths=[snakemake.output.csv, snakemake.output.figure, snakemake.output.provenance],
)
