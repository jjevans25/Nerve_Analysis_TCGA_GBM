# workflow/rules/integration.smk
# scVI-based batch integration and latent space learning.

import os

rule scrna_integration:
    """Train scVI VAE on QC-filtered samples; records seed, latent dims, and loss in provenance."""
    input:
        h5ads = expand(
            os.path.join(config["dirs"]["data_processed"], "{sample}_qc.h5ad"),
            sample=config.get("samples", []),
        ),
    output:
        model_dir   = directory(os.path.join(config["dirs"]["models"], "scvi_model")),
        latent_h5ad = os.path.join(config["dirs"]["data_processed"], "integrated_latent.h5ad"),
        provenance  = os.path.join(config["dirs"]["provenance"],     "scvi_integration_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "scvi_integration.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["gpu_threads"],
    params:
        device      = config["hardware"]["device"],
        precision   = config["hardware"]["precision"],
        n_latent    = config["scrna"]["n_latent"],
        n_layers    = config["scrna"]["n_layers"],
        batch_key   = config["scrna"]["batch_key"],
        random_seed = config["scrna"]["random_seed"],
        min_cells   = config["scrna"]["min_cells"],
        # Baseline .X is Seurat SCT log1p (from loom_to_h5ad) — a structural
        # property of the source, not tunable. scVI reads recovered counts.
        counts_from_log1p = True,
    script:
        "../scripts/scrna_integration.py"


rule scrna_velocity:
    """Estimate RNA velocity on integrated data using scVelo."""
    input:
        latent_h5ad = os.path.join(config["dirs"]["data_processed"], "integrated_latent.h5ad"),
        loom_files  = expand(
            os.path.join(config["dirs"]["data_raw"], "{sample}.loom"),
            sample=config.get("samples", []),
        ),
    output:
        velocity_h5ad = os.path.join(config["dirs"]["data_processed"], "velocity.h5ad"),
        velocity_plot = os.path.join(config["dirs"]["figures"],        "velocity_streamplot.png"),
        provenance    = os.path.join(config["dirs"]["provenance"],     "velocity_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "scrna_velocity.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    script:
        "../scripts/scrna_velocity.py"
