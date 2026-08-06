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

import gc
import hashlib
import json
import os
import shutil
import sys
import warnings
from pathlib import Path

# numba defaults to an OpenMP threading layer, which collides with the libomp
# that torch/MPS has already initialised in this process. On 2026-08-02
# sc.pp.neighbors (numba-jitted pynndescent) segfaulted twice at 377k cells —
# SIGSEGV in __kmp_launch_worker, once at 32.8 GB peak RSS and again at 16.1 GB,
# which rules out memory pressure and identifies a runtime conflict. "workqueue"
# is numba's own non-OpenMP pool and is documented as always safe. This must be
# set before numba is imported, i.e. before scanpy pulls in pynndescent/umap.
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
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
# Checkpoint identity guard
# ----------------------------------------------------------------------------
# These checkpoints live at UNDECLARED paths so Snakemake cannot delete them when
# a job fails — that is what makes them crash-safe. The same property makes them
# dangerous: they also survive a legitimate change to the input, and nothing here
# used to check. On 2026-08-05 a baseline trained on the pre-fix 377,343-cell
# nerve subset was loaded against the corrected 60,034-cell subset and the rule
# died ~29 min into a 10 h run.
#
# So: fingerprint the training problem, store it beside the checkpoint, and reuse
# only on an exact match. Crash-safe *and* wrong-safe.


def _fingerprint() -> dict:
    """Identity of the training problem — anything here changing invalidates a checkpoint."""
    roster = hashlib.sha256("\n".join(map(str, adata.obs_names)).encode()).hexdigest()
    genes = hashlib.sha256("\n".join(map(str, adata.var_names)).encode()).hexdigest()
    labels = adata.obs[str(snakemake.params.labels_key)].astype(str)
    return {
        "n_obs":       int(adata.n_obs),
        "n_vars":      int(adata.n_vars),
        "obs_sha256":  roster,
        "var_sha256":  genes,
        "batch_key":   str(snakemake.params.batch_key),
        "labels_key":  str(snakemake.params.labels_key),
        "label_set":   sorted(labels.unique().tolist()),
        "n_latent":    int(snakemake.params.n_latent),
        "n_layers":    int(snakemake.params.n_layers),
        "seed":        SEED,
    }


FINGERPRINT = _fingerprint()
_FP_NAME = "input_fingerprint.json"


def _reusable(ckpt_dir: Path, what: str) -> bool:
    """Is this checkpoint from the same training problem we are solving now?"""
    if not ckpt_dir.exists():
        return False
    fp_path = ckpt_dir / _FP_NAME
    if not fp_path.exists():
        log_transformation(log, "nerve_scanvi_retrain",
                           f"[FAIR-ALERT] {what} checkpoint at {ckpt_dir} predates fingerprinting "
                           "— cannot prove it matches this input. Retraining from scratch.",
                           status="WARNING")
        return False
    with open(fp_path) as fh:
        stored = json.load(fh)
    if stored == FINGERPRINT:
        return True
    diffs = [k for k in FINGERPRINT if stored.get(k) != FINGERPRINT[k]]
    log_transformation(log, "nerve_scanvi_retrain",
                       f"[FAIR-ALERT] {what} checkpoint at {ckpt_dir} was trained on a DIFFERENT "
                       f"input — mismatched fields: {diffs} "
                       f"(stored n_obs={stored.get('n_obs')}, current n_obs={FINGERPRINT['n_obs']}). "
                       "Discarding it and retraining from scratch.",
                       status="WARNING")
    shutil.rmtree(ckpt_dir)
    return False


def _stamp(ckpt_dir: Path) -> None:
    """Record which training problem this checkpoint belongs to."""
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    with open(ckpt_dir / _FP_NAME, "w") as fh:
        json.dump(FINGERPRINT, fh, indent=2)


# ----------------------------------------------------------------------------
# Stage 1: baseline scVI
# ----------------------------------------------------------------------------
scvi.model.SCVI.setup_anndata(
    adata,
    layer=None,
    batch_key=snakemake.params.batch_key,
)
# The baseline is checkpointed so a stage-2 crash doesn't lose the ~4.5 h
# pretrain — but until 2026-08-02 nothing ever *read* it back, so the protection
# the comment promised was never actually realised and every retry paid the full
# pretrain again. Reuse it when present. Note this path is deliberately NOT a
# declared rule output: Snakemake removes a failed job's declared outputs ("since
# they might be corrupted"), which is exactly how the stage-2 model was lost on
# the 2026-08-02 14:44 failure while this checkpoint survived.
baseline_dir = Path(snakemake.output.model_dir).parent / "nerve_scvi_baseline"

if _reusable(baseline_dir, "baseline scVI"):
    scvi_model = scvi.model.SCVI.load(str(baseline_dir), adata=adata)
    log_transformation(log, "nerve_scanvi_retrain",
                       f"Reusing baseline scVI checkpoint at {baseline_dir} "
                       "(skipping stage-1 pretrain)")
else:
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
    scvi_model.save(str(baseline_dir), overwrite=True)
    _stamp(baseline_dir)
    log_transformation(log, "nerve_scanvi_retrain",
                       f"Baseline scVI checkpointed at {baseline_dir}")

# history_ is persisted inside model.pt, so the curves figure works either way.
scvi_history = {k: list(v.iloc[:, 0]) for k, v in scvi_model.history.items()}

# ----------------------------------------------------------------------------
# Stage 2: scANVI fine-tune with cell_type labels
# ----------------------------------------------------------------------------
scanvi_sidecar_in = Path(snakemake.output.model_dir).parent / "nerve_scanvi_stage2"
if _reusable(scanvi_sidecar_in, "stage-2 scANVI"):
    scanvi_model = scvi.model.SCANVI.load(str(scanvi_sidecar_in), adata=adata)
    log_transformation(log, "nerve_scanvi_retrain",
                       f"Reusing stage-2 scANVI checkpoint at {scanvi_sidecar_in} "
                       "(skipping scANVI fine-tune)")
else:
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

# Persist the trained model *before* the neighbors/UMAP step, because a crash in
# that purely visual step must not cost the ~2.5 h stage-2 train (it did on
# 2026-08-01). Write to an undeclared sidecar as well: on failure Snakemake
# deletes the declared output, so the sidecar is what actually survives and what
# the reuse branch below can pick up.
scanvi_dir = Path(snakemake.output.model_dir)
scanvi_sidecar = scanvi_dir.parent / "nerve_scanvi_stage2"
scanvi_model.save(str(scanvi_sidecar), overwrite=True)
_stamp(scanvi_sidecar)
scanvi_model.save(str(scanvi_dir), overwrite=True)
log_transformation(log, "nerve_scanvi_retrain",
                   f"scANVI model saved at {scanvi_dir} (+ crash-safe sidecar at "
                   f"{scanvi_sidecar}) before neighbors/UMAP")

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

# Release the model and the MPS allocator pool before the kNN graph. The scANVI
# weights and torch's cached MPS blocks are dead weight from here on, and
# sc.pp.neighbors needs headroom on a 36 GB machine: at 377k cells the previous
# run exhausted memory here and OpenMP segfaulted spawning a worker thread.
del scanvi_model
gc.collect()
if device == "mps":
    torch.mps.empty_cache()

log_transformation(log, "nerve_scanvi_retrain", "Computing neighbors + UMAP on X_scANVI")
sc.pp.neighbors(adata, use_rep="X_scANVI", random_state=SEED)
sc.tl.umap(adata, random_state=SEED)

# ----------------------------------------------------------------------------
# Save latent h5ad + training curves (model already saved above)
# ----------------------------------------------------------------------------
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
