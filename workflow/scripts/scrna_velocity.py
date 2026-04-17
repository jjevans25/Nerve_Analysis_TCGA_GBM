"""RNA velocity estimation with scVelo on MPS-integrated latent space."""

import sys

import anndata as ad
import matplotlib.pyplot as plt
import scvelo as scv

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
log_transformation(log, "scrna_velocity", "Loading integrated latent AnnData")

adata = ad.read_h5ad(snakemake.input.latent_h5ad)

# --- Merge spliced/unspliced counts from loom files ---------------------------
loom_adatas = [scv.read(f, cache=True) for f in snakemake.input.loom_files]
loom = ad.concat(loom_adatas, merge="same")
loom.obs_names_make_unique()

scv.utils.merge(adata, loom)
log_transformation(log, "scrna_velocity",
                   f"Merged loom data: {adata.n_obs} cells with spliced/unspliced counts")

# --- Preprocessing for velocity -----------------------------------------------
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# --- Velocity estimation (dynamical model) ------------------------------------
scv.tl.recover_dynamics(adata, n_jobs=snakemake.resources.threads)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

log_transformation(log, "scrna_velocity", "Velocity graph computed")

# --- Save outputs -------------------------------------------------------------
adata.write_h5ad(snakemake.output.velocity_h5ad)
verify_artifact(snakemake.output.velocity_h5ad, min_size_bytes=1024)

fig, ax = plt.subplots(figsize=(8, 6))
scv.pl.velocity_embedding_stream(adata, basis="X_scVI", ax=ax, show=False)
fig.savefig(snakemake.output.velocity_plot, dpi=150, bbox_inches="tight")
plt.close(fig)

# --- FAIR provenance ----------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.velocity_h5ad,
    rule_name="scrna_velocity",
    input_paths=[snakemake.input.latent_h5ad, *snakemake.input.loom_files],
    tool_versions={"scvelo": scv.__version__, "anndata": ad.__version__},
    parameters={"velocity_mode": "dynamical", "n_top_genes": 2000, "n_pcs": 30},
    description="RNA velocity estimates (dynamical model) on scVI latent space",
    ontology_operation="operation:0430",  # EDAM: Gene expression analysis
)
write_provenance(prov, "provenance")

log_transformation(log, "scrna_velocity", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.velocity_h5ad, snakemake.output.velocity_plot])
