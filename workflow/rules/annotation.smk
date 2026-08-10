# workflow/rules/annotation.smk
# Cell-type annotation (marker scoring + scVI latent KNN) and CNV-based malignancy labeling.

import os


rule scrna_annotate:
    """Annotate cell types using scVI latent space clustering and canonical marker gene scoring."""
    input:
        latent_h5ad = os.path.join(config["dirs"]["data_processed"], "integrated_latent.h5ad"),
        model_dir   = os.path.join(config["dirs"]["models"], "scvi_model"),
        symbol_map  = config["gene_symbol_map"]["cache_tsv"],
    output:
        h5ad       = os.path.join(config["dirs"]["data_processed"], "annotated.h5ad"),
        summary    = os.path.join(config["dirs"]["tables"],         "annotation_summary.csv"),
        cluster_scores = os.path.join(config["dirs"]["tables"],      "annotation_cluster_scores.csv"),
        provenance = os.path.join(config["dirs"]["provenance"],     "annotation_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "scrna_annotate.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        leiden_resolution  = ANNOTATE_LEIDEN_RESOLUTION,
        markers            = config["nerve_cells"]["markers"],
        annotation_markers = ANNOTATION_MARKERS,
        ambiguous_margin   = ANNOTATE_AMBIGUOUS_MARGIN,
        panel_compartment  = ANNOTATE_PANEL_COMPARTMENT,
        random_seed        = config["scrna"]["random_seed"],
        census_version     = config["databases"]["cellxgene_census_version"],
    script:
        "../scripts/scrna_annotate.py"


rule scrna_malignancy:
    """CNV-based tumor/normal separation; labels each cell with is_malignant flag."""
    input:
        h5ad = os.path.join(config["dirs"]["data_processed"], "annotated.h5ad"),
        gene_positions = os.path.join(config["gene_positions"]["download_dir"],
                                      config["gene_positions"]["positions_filename"]),
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
        cnv_chunk_size = config["scrna"]["cnv_chunk_size"],
        cnv_window     = config["scrna"]["cnv_window"],
        cnv_clip       = config["scrna"]["cnv_clip"],
        cnv_threshold_sd = config["scrna"]["cnv_threshold_sd"],
        reference_labels = config["scrna"]["cnv_reference_labels"],
        reference_confidence_quantile = config["scrna"]["cnv_reference_confidence_quantile"],
        min_reference_cells = config["scrna"]["cnv_min_reference_cells"],
        min_genes_placed_fraction = config["scrna"]["cnv_min_genes_placed_fraction"],
        cnv_exclusion_sd = config["scrna"]["cnv_exclusion_sd"],
        cnv_gain_contigs = config["scrna"]["cnv_gain_contigs"],
        cnv_loss_contigs = config["scrna"]["cnv_loss_contigs"],
    script:
        "../scripts/scrna_malignancy.py"
