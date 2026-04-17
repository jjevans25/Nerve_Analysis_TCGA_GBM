# workflow/rules/annotation.smk
# Cell-type annotation (marker scoring + scVI latent KNN) and CNV-based malignancy labeling.

import os


rule scrna_annotate:
    """Annotate cell types using scVI latent space clustering and canonical marker gene scoring."""
    input:
        latent_h5ad = os.path.join(config["dirs"]["data_processed"], "integrated_latent.h5ad"),
        model_dir   = os.path.join(config["dirs"]["models"], "scvi_model"),
    output:
        h5ad       = os.path.join(config["dirs"]["data_processed"], "annotated.h5ad"),
        summary    = os.path.join(config["dirs"]["tables"],         "annotation_summary.csv"),
        provenance = os.path.join(config["dirs"]["provenance"],     "annotation_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "scrna_annotate.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        leiden_resolution = config["nerve_cells"]["leiden_resolution"],
        markers           = config["nerve_cells"]["markers"],
        random_seed       = config["scrna"]["random_seed"],
        census_version    = config["databases"]["cellxgene_census_version"],
    script:
        "../scripts/scrna_annotate.py"


rule scrna_malignancy:
    """CNV-based tumor/normal separation; labels each cell with is_malignant flag."""
    input:
        h5ad = os.path.join(config["dirs"]["data_processed"], "annotated.h5ad"),
    output:
        h5ad       = os.path.join(config["dirs"]["data_processed"], "malignancy_labeled.h5ad"),
        cnv_plot   = os.path.join(config["dirs"]["figures"],        "cnv_heatmap.png"),
        provenance = os.path.join(config["dirs"]["provenance"],     "malignancy_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "scrna_malignancy.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        random_seed    = config["scrna"]["random_seed"],
        census_version = config["databases"]["cellxgene_census_version"],
    script:
        "../scripts/scrna_malignancy.py"
