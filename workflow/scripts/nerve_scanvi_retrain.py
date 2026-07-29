"""Re-train scVI -> scANVI on the nerve-cell subset with canonical cell-type labels.

Question answered: the v1.0.0 X_scVI carries residual patient-private structure
that survives every Leiden resolution (sweep, separate plan). scANVI adds a
classification head over biological labels (neuron / OPC / oligo / astrocyte,
marker-derived), pushing the latent space toward biology and away from
patient identity. ``sample_id`` remains ``batch_key`` — the new signal is the
``labels_key`` head.

This script trains a NEW model from scratch on the nerve subset only. The
global v1.0.0 ``results/models/scvi_model`` is not touched.
"""

import os
import sys
import warnings
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import (
    H5AD_COMPRESSION,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]
SEED = int(snakemake.params.random_seed)
torch.manual_seed(SEED)
os.environ["PYTHONHASHSEED"] = str(SEED)

device = str(snakemake.params.device)
precision = str(snakemake.params.precision)
if device == "mps":
    if not torch.backends.mps.is_available():
        warnings.warn("[MPS-ALERT] MPS requested but not available — falling back to CPU.")
        device = "cpu"
    else:
        torch.set_default_dtype(torch.float32)
        os.environ.setdefault("PYTORCH_MPS_HIGH_WATERMARK_RATIO", "0.0")

log_transformation(log, "nerve_scanvi_retrain",
                   f"Device: {device} | Precision: {precision} | Seed: {SEED}")

# scvi-tools import happens after the env-var setup so the lib sees them.
import scvi  # noqa: E402

scvi.settings.seed = SEED

# ----------------------------------------------------------------------------
# Load counts + labels
# ----------------------------------------------------------------------------
adata = ad.read_h5ad(snakemake.input.counts_h5ad)
log_transformation(log, "nerve_scanvi_retrain",
                   f"{adata.n_obs} cells x {adata.n_vars} genes")

for col in (snakemake.params.batch_key, snakemake.params.labels_key):
    if col not in adata.obs.columns:
        raise RuntimeError(f"[FAIR-ALERT] required obs column '{col}' missing")

label_counts = adata.obs[snakemake.params.labels_key].value_counts()
log_transformation(log, "nerve_scanvi_retrain",
                   f"Label distribution:\n{label_counts.to_string()}")

# ----------------------------------------------------------------------------
# Stage 1: baseline scVI
# ----------------------------------------------------------------------------
scvi.model.SCVI.setup_anndata(
    adata,
    layer=None,
    batch_key=snakemake.params.batch_key,
)
scvi_model = scvi.model.SCVI(
    adata,
    n_latent=int(snakemake.params.n_latent),
    n_layers=int(snakemake.params.n_layers),
)
log_transformation(log, "nerve_scanvi_retrain",
                   f"Training baseline scVI (max_epochs={snakemake.params.scvi_max_epochs})")
scvi_model.train(
    accelerator=device,
    devices=1,
    max_epochs=int(snakemake.params.scvi_max_epochs),
    early_stopping=True,
)
scvi_history = {k: list(v.iloc[:, 0]) for k, v in scvi_model.history.items()}

# Save the baseline so a scANVI-stage crash doesn't lose the long pretrain.
baseline_dir = Path(snakemake.output.model_dir).parent / "nerve_scvi_baseline"
scvi_model.save(str(baseline_dir), overwrite=True)
log_transformation(log, "nerve_scanvi_retrain",
                   f"Baseline scVI checkpointed at {baseline_dir}")

# ----------------------------------------------------------------------------
# Stage 2: scANVI fine-tune with cell_type labels
# ----------------------------------------------------------------------------
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    labels_key=str(snakemake.params.labels_key),
    unlabeled_category=str(snakemake.params.unlabeled_category),
)
log_transformation(log, "nerve_scanvi_retrain",
                   f"Training scANVI (max_epochs={snakemake.params.scanvi_max_epochs}, "
                   f"n_samples_per_label={snakemake.params.n_samples_per_label})")
scanvi_model.train(
    accelerator=device,
    devices=1,
    max_epochs=int(snakemake.params.scanvi_max_epochs),
    n_samples_per_label=int(snakemake.params.n_samples_per_label),
    early_stopping=True,
)
scanvi_history = {k: list(v.iloc[:, 0]) for k, v in scanvi_model.history.items()}

# ----------------------------------------------------------------------------
# Extract latent, predicted labels, kNN, UMAP
# ----------------------------------------------------------------------------
X_scanvi = scanvi_model.get_latent_representation()
if not np.all(np.isfinite(X_scanvi)):
    raise RuntimeError("[FAIR-ALERT] X_scANVI contains NaN/inf values")
adata.obsm["X_scANVI"] = X_scanvi.astype(np.float32)

adata.obs["scanvi_predicted_celltype"] = scanvi_model.predict()

# Held-out classification accuracy on the labelled subset.
labels_true = adata.obs[snakemake.params.labels_key].astype(str)
labels_pred = adata.obs["scanvi_predicted_celltype"].astype(str)
labelled_mask = labels_true != str(snakemake.params.unlabeled_category)
if labelled_mask.sum() > 0:
    acc = float((labels_true[labelled_mask] == labels_pred[labelled_mask]).mean())
else:
    acc = float("nan")
log_transformation(log, "nerve_scanvi_retrain",
                   f"Classifier accuracy on labelled cells: {acc:.4f}")

log_transformation(log, "nerve_scanvi_retrain", "Computing neighbors + UMAP on X_scANVI")
sc.pp.neighbors(adata, use_rep="X_scANVI", random_state=SEED)
sc.tl.umap(adata, random_state=SEED)

# ----------------------------------------------------------------------------
# Save model + latent h5ad + training curves
# ----------------------------------------------------------------------------
scanvi_model.save(snakemake.output.model_dir, overwrite=True)

adata.write_h5ad(snakemake.output.latent_h5ad, compression=H5AD_COMPRESSION)
verify_artifact(snakemake.output.latent_h5ad, min_size_bytes=1_000_000)

# Training-curves figure (two panels: scVI ELBO and scANVI losses).
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
for k, v in scvi_history.items():
    if "loss" in k.lower() or "elbo" in k.lower():
        ax.plot(range(len(v)), v, label=k)
ax.set_xlabel("epoch")
ax.set_ylabel("value")
ax.set_title("Stage 1 — baseline scVI")
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

ax = axes[1]
for k, v in scanvi_history.items():
    if any(s in k.lower() for s in ("loss", "elbo", "classification", "accuracy")):
        ax.plot(range(len(v)), v, label=k)
ax.set_xlabel("epoch")
ax.set_ylabel("value")
ax.set_title("Stage 2 — scANVI fine-tune")
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

fig.suptitle("scANVI re-train training curves (nerve subset)", y=1.02)
fig.tight_layout()
fig.savefig(snakemake.output.curves, dpi=140, bbox_inches="tight")
plt.close(fig)
verify_artifact(snakemake.output.curves, min_size_bytes=10_000)

# ----------------------------------------------------------------------------
# Provenance
# ----------------------------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.latent_h5ad,
    rule_name="nerve_scanvi_retrain",
    input_paths=[snakemake.input.counts_h5ad],
    tool_versions={
        "scvi-tools": scvi.__version__,
        "torch": torch.__version__,
        "anndata": ad.__version__,
        "scanpy": sc.__version__,
    },
    parameters={
        "device": device,
        "precision": precision,
        "n_latent": int(snakemake.params.n_latent),
        "n_layers": int(snakemake.params.n_layers),
        "batch_key": str(snakemake.params.batch_key),
        "labels_key": str(snakemake.params.labels_key),
        "unlabeled_category": str(snakemake.params.unlabeled_category),
        "scvi_max_epochs": int(snakemake.params.scvi_max_epochs),
        "scanvi_max_epochs": int(snakemake.params.scanvi_max_epochs),
        "n_samples_per_label": int(snakemake.params.n_samples_per_label),
        "random_seed": SEED,
        "classifier_accuracy_labelled": round(acc, 4) if not np.isnan(acc) else None,
        "label_distribution": {str(k): int(v) for k, v in label_counts.items()},
    },
    description=(
        "scANVI latent representation of the nerve-cell subset. Baseline scVI "
        "trained with sample_id as batch_key; scANVI fine-tune anchored on "
        "marker-derived cell_type labels (neuron/OPC/oligo/astrocyte + Unknown)."
    ),
    ontology_operation="data:3917",  # EDAM: Gene expression matrix
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log, "nerve_scanvi_retrain",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.model_dir,
        snakemake.output.latent_h5ad,
        snakemake.output.curves,
    ],
)
