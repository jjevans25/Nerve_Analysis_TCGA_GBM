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
from fair_utils import H5AD_COMPRESSION, log_transformation, stamp_artifact, verify_artifact, write_provenance

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

# Release the per-sample objects immediately. Without this the full input set
# (~12 GB for the 170-sample Census cohort, ~21 GB for the 17-sample baseline)
# stays resident alongside the concatenated copy for the rest of the script —
# which is what forced the cohort to be halved on this 36 GB machine.
n_samples_concat = len(adatas)
adatas.clear()
del adatas

adata.obs_names_make_unique()
n_genes_concat = adata.n_vars

# Gene filtering happens ONCE, here, across the whole cohort — never per sample.
sc.pp.filter_genes(adata, min_cells=snakemake.params.min_cells)

# Guard the int32 headroom. scipy promotes CSR index arrays to int64 once nnz
# exceeds 2**31, which silently doubles the matrix footprint (13 GB -> 26 GB at
# full-cohort scale). If this ever trips, reduce the gene space (HVG selection)
# rather than the cell count.
if sp.issparse(adata.X) and adata.X.indices.dtype != np.int32:
    log_transformation(snakemake.log[0], "scrna_integration",
        f"[FAIR-ALERT] CSR indices promoted to {adata.X.indices.dtype} "
        f"(nnz={adata.X.nnz}); matrix footprint has doubled", status="WARNING")

# scVI requires raw integer counts. The baseline arm's .X is Seurat SCT log1p —
# recover integer counts via expm1; census cohorts already ship integer counts
# in .X. Explicit per-arm flag, never an X.max() heuristic. .X itself is left
# untouched (log1p for the baseline) so every downstream consumer that reads .X
# keeps its published behavior.
#
# When .X already IS counts (census arm) we do NOT materialize a "counts" layer:
# it would be a bit-identical duplicate of the matrix (~10 GB for the Census
# cohort), held through the whole of training and written into the output h5ad.
# scVI reads .X directly when layer=None, so the copy buys nothing.
if snakemake.params.counts_from_log1p:
    adata.layers["counts"] = recover_counts_from_log1p(adata.X)
    adata.uns["counts_recovery"] = "expm1_round_int32"
    scvi_layer = "counts"
    _counts_matrix = adata.layers["counts"]
else:
    adata.uns["counts_recovery"] = "identity"
    scvi_layer = None
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    _counts_matrix = adata.X
adata.uns["scvi_counts_source"] = scvi_layer or "X"

# Fail fast if the scVI input is not clean non-negative integers.
_counts_data = _counts_matrix.data
if _counts_data.size and (
    np.isnan(_counts_data).any()
    or (_counts_data < 0).any()
    or not np.all(_counts_data == _counts_data.astype(np.int32))
):
    raise RuntimeError("[FAIR-ALERT] counts input is not clean non-negative integers")
del _counts_matrix, _counts_data

log_transformation(snakemake.log[0], "scrna_integration",
                   f"Concatenated {n_samples_concat} samples → {adata.n_obs} cells; "
                   f"genes {n_genes_concat} → {adata.n_vars} "
                   f"(global filter_genes(min_cells={snakemake.params.min_cells})); "
                   f"scVI counts via {adata.uns['scvi_counts_source']} "
                   f"({adata.uns['counts_recovery']})")

# --- scVI setup ---------------------------------------------------------------
scvi.settings.seed = SEED
scvi.settings.num_threads = snakemake.resources.threads

scvi.model.SCVI.setup_anndata(
    adata,
    # "counts" layer for the baseline arm (recovered via expm1); None for cohorts
    # whose .X already holds raw integer counts — see the counts block above.
    layer=scvi_layer,
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
adata.write_h5ad(snakemake.output.latent_h5ad, compression=H5AD_COMPRESSION)

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
