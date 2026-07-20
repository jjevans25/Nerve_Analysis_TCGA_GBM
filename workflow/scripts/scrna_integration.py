"""Train scVI VAE on QC-filtered AnnData objects; MPS-accelerated on Apple Silicon."""

import os
import sys
import warnings

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch

sys.path.insert(0, "workflow/scripts")
from counts_utils import recover_counts_from_log1p
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

# --- Reproducibility -----------------------------------------------------------
SEED = snakemake.params.random_seed
torch.manual_seed(SEED)
os.environ["PYTHONHASHSEED"] = str(SEED)

# --- MPS configuration --------------------------------------------------------
device = snakemake.params.device          # "mps" | "cpu" | "cuda"
precision = snakemake.params.precision    # "float32"

if device == "mps":
    if not torch.backends.mps.is_available():
        warnings.warn("[MPS-ALERT] MPS requested but not available — falling back to CPU.")
        device = "cpu"
    else:
        # MPS float16 is often slower than float32; enforce float32
        torch.set_default_dtype(torch.float32)
        # Prevent memory high-watermark OOM kills during long VAE training
        os.environ.setdefault("PYTORCH_MPS_HIGH_WATERMARK_RATIO", "0.0")

log_transformation(snakemake.log[0], "scrna_integration",
                   f"Device: {device} | Precision: {precision} | Seed: {SEED}")

# --- Load and concatenate QC-filtered samples ---------------------------------
# scvi/torch are imported here, after the MPS env setup above, so the
# PYTORCH_MPS_HIGH_WATERMARK_RATIO / float32 settings take effect before import.
import scvi  # noqa: E402

adatas = []
for f in snakemake.input.h5ads:
    a = ad.read_h5ad(f)
    # Sparse .X makes the outer-join pad missing genes with 0, not NaN. The
    # baseline arm's .X is dense with heterogeneous per-sample gene sets, which
    # NaN-padded the latent and crashed scVI; sparsifying removes that failure
    # mode structurally. Census cohorts are already sparse (no-op here).
    if not sp.issparse(a.X):
        a.X = sp.csr_matrix(a.X)
    adatas.append(a)

# join="outer": never drop a gene merely because one sample lacks it.
# fill_value=0: belt-and-suspenders so no arm can ever NaN-pad missing genes.
adata = ad.concat(adatas, join="outer", merge="same", fill_value=0)
adata.obs_names_make_unique()
n_genes_concat = adata.n_vars

# Gene filtering happens ONCE, here, across the whole cohort — never per sample.
sc.pp.filter_genes(adata, min_cells=snakemake.params.min_cells)

# scVI requires raw integer counts. The baseline arm's .X is Seurat SCT log1p —
# recover integer counts via expm1; census cohorts already ship integer counts
# in .X. Explicit per-arm flag, never an X.max() heuristic. .X itself is left
# untouched (log1p for the baseline) so every downstream consumer that reads .X
# keeps its published behavior.
if snakemake.params.counts_from_log1p:
    adata.layers["counts"] = recover_counts_from_log1p(adata.X)
    adata.uns["counts_recovery"] = "expm1_round_int32"
else:
    adata.layers["counts"] = sp.csr_matrix(adata.X)
    adata.uns["counts_recovery"] = "identity"

# Fail fast if the scVI input is not clean non-negative integers.
_counts_data = adata.layers["counts"].data
if _counts_data.size and (
    np.isnan(_counts_data).any()
    or (_counts_data < 0).any()
    or not np.all(_counts_data == _counts_data.astype(np.int32))
):
    raise RuntimeError("[FAIR-ALERT] counts layer is not clean non-negative integers")

log_transformation(snakemake.log[0], "scrna_integration",
                   f"Concatenated {len(adatas)} samples → {adata.n_obs} cells; "
                   f"genes {n_genes_concat} → {adata.n_vars} "
                   f"(global filter_genes(min_cells={snakemake.params.min_cells})); "
                   f"scVI counts layer via {adata.uns['counts_recovery']}")

# --- scVI setup ---------------------------------------------------------------
scvi.settings.seed = SEED
scvi.settings.num_threads = snakemake.resources.threads

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",          # raw integer counts (recovered for the baseline arm)
    batch_key=snakemake.params.batch_key,
)

model = scvi.model.SCVI(
    adata,
    n_latent=snakemake.params.n_latent,
    n_layers=snakemake.params.n_layers,
)

# Train on MPS via PyTorch Lightning accelerator
model.train(
    accelerator=device,
    devices=1,
    max_epochs=400,
    early_stopping=True,
)

# --- Save model and latent representation -------------------------------------
model.save(snakemake.output.model_dir, overwrite=True)

adata.obsm["X_scVI"] = model.get_latent_representation()
adata.write_h5ad(snakemake.output.latent_h5ad)

verify_artifact(snakemake.output.latent_h5ad, min_size_bytes=1024)

# --- FAIR provenance ----------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.latent_h5ad,
    rule_name="scrna_integration",
    input_paths=snakemake.input.h5ads,
    tool_versions={
        "scvi-tools": scvi.__version__,
        "torch": torch.__version__,
        "anndata": ad.__version__,
    },
    parameters={
        "device": device,
        "precision": precision,
        "n_latent": snakemake.params.n_latent,
        "n_layers": snakemake.params.n_layers,
        "batch_key": snakemake.params.batch_key,
        "random_seed": SEED,
        "min_cells": snakemake.params.min_cells,
        "n_genes_concat": n_genes_concat,
        "n_genes_kept": int(adata.n_vars),
        "counts_from_log1p": bool(snakemake.params.counts_from_log1p),
        "counts_recovery": adata.uns["counts_recovery"],
    },
    description="scVI VAE latent representation of batch-corrected scRNA-seq data",
    ontology_operation="data:3917",  # EDAM: Gene expression matrix
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(snakemake.log[0], "scrna_integration",
                   f"Model saved to {snakemake.output.model_dir}", status="SUCCESS",
                   artifact_paths=[snakemake.output.latent_h5ad])
