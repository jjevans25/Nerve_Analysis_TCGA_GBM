"""Train scVI VAE on QC-filtered AnnData objects; MPS-accelerated on Apple Silicon."""

import os
import sys
import warnings

import anndata as ad
import torch

sys.path.insert(0, "workflow/scripts")
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
import scvi
from pytorch_lightning import Trainer

adatas = [ad.read_h5ad(f) for f in snakemake.input.h5ads]
adata = ad.concat(adatas, merge="same")
adata.obs_names_make_unique()

log_transformation(snakemake.log[0], "scrna_integration",
                   f"Concatenated {len(adatas)} samples → {adata.n_obs} cells, {adata.n_vars} genes")

# --- scVI setup ---------------------------------------------------------------
scvi.settings.seed = SEED
scvi.settings.num_threads = snakemake.resources.threads

scvi.model.SCVI.setup_anndata(
    adata,
    layer=None,              # expects raw counts in .X
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
    },
    description="scVI VAE latent representation of batch-corrected scRNA-seq data",
    ontology_operation="data:3917",  # EDAM: Gene expression matrix
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(snakemake.log[0], "scrna_integration",
                   f"Model saved to {snakemake.output.model_dir}", status="SUCCESS",
                   artifact_paths=[snakemake.output.latent_h5ad])
