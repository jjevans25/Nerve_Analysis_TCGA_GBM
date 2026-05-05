# workflow/rules/nerve_cells.smk
# Nerve cell subset extraction and heterogeneity characterization.

import os


rule nerve_cell_subset:
    """Extract non-malignant neural/glial cells; re-cluster in nerve-cell subspace."""
    input:
        h5ad     = os.path.join(config["dirs"]["data_processed"],  "malignancy_labeled.h5ad"),
        clinical = os.path.join(config["dirs"]["data_external"],   "gdc_clinical.tsv"),
    output:
        h5ad        = os.path.join(config["dirs"]["data_processed"], "nerve_cells.h5ad"),
        umap        = os.path.join(config["dirs"]["figures"],        "nerve_cells_umap.png"),
        composition = os.path.join(config["dirs"]["tables"],         "nerve_cluster_composition.csv"),
        provenance  = os.path.join(config["dirs"]["provenance"],     "nerve_subset_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_cell_subset.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        leiden_resolution = config["nerve_cells"]["leiden_resolution"],
        cell_types        = config["nerve_cells"]["cell_types"],
        markers           = config["nerve_cells"]["markers"],
        random_seed       = config["scrna"]["random_seed"],
        n_top_genes       = config["scrna"]["n_top_genes"],
    script:
        "../scripts/nerve_cell_subset.py"


rule nerve_cell_heterogeneity:
    """Differential expression, GSEA, marker dot plot, and per-sample cluster abundance."""
    input:
        h5ad       = os.path.join(config["dirs"]["data_processed"], "nerve_cells.h5ad"),
        symbol_map = config["gene_symbol_map"]["cache_tsv"],
    output:
        markers    = os.path.join(config["dirs"]["tables"],  "nerve_cluster_markers.csv"),
        enrichment = os.path.join(config["dirs"]["tables"],  "nerve_enrichment.csv"),
        dotplot    = os.path.join(config["dirs"]["figures"], "nerve_dotplot.png"),
        abundance  = os.path.join(config["dirs"]["figures"], "nerve_abundance_heatmap.png"),
        provenance = os.path.join(config["dirs"]["provenance"], "heterogeneity_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_cell_heterogeneity.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        markers     = config["nerve_cells"]["markers"],
        random_seed = config["scrna"]["random_seed"],
    script:
        "../scripts/nerve_cell_heterogeneity.py"


rule nerve_cluster_annotations:
    """Per-cluster biological labels: dominant cell type + top markers + canonical module scores."""
    input:
        h5ad       = os.path.join(config["dirs"]["data_processed"], "nerve_cells.h5ad"),
        markers    = os.path.join(config["dirs"]["tables"],         "nerve_cluster_markers.csv"),
        symbol_map = config["gene_symbol_map"]["cache_tsv"],
    output:
        annotations = os.path.join(config["dirs"]["tables"],     "nerve_cluster_annotations.csv"),
        provenance  = os.path.join(config["dirs"]["provenance"], "nerve_cluster_annotations_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_cluster_annotations.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        markers     = config["nerve_cells"]["markers"],
        random_seed = config["scrna"]["random_seed"],
    script:
        "../scripts/nerve_cluster_annotations.py"


rule nerve_clinical_association:
    """Test cluster-abundance differences across clinical covariates and age correlation."""
    input:
        h5ad     = os.path.join(config["dirs"]["data_processed"], "nerve_cells.h5ad"),
        clinical = os.path.join(config["dirs"]["data_external"],  "gdc_clinical.tsv"),
    output:
        stats      = os.path.join(config["dirs"]["tables"],     "nerve_clinical_association.csv"),
        heatmap    = os.path.join(config["dirs"]["figures"],    "nerve_clinical_pvalue_heatmap.png"),
        boxplots   = os.path.join(config["dirs"]["figures"],    "nerve_clinical_boxplots.png"),
        provenance = os.path.join(config["dirs"]["provenance"], "nerve_clinical_association_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_clinical_association.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        random_seed = config["scrna"]["random_seed"],
    script:
        "../scripts/nerve_clinical_association.py"


rule nerve_tumor_interaction:
    """Tumor-nerve cell-cell communication via LIANA+ consensus ligand-receptor inference."""
    input:
        malig = os.path.join(config["dirs"]["data_processed"], "malignancy_labeled.h5ad"),
        nerve = os.path.join(config["dirs"]["data_processed"], "nerve_cells.h5ad"),
    output:
        lr_table   = os.path.join(config["dirs"]["tables"],     "nerve_tumor_interactions.csv"),
        top_pairs  = os.path.join(config["dirs"]["tables"],     "nerve_tumor_top_pairs.csv"),
        heatmap    = os.path.join(config["dirs"]["figures"],    "nerve_tumor_sig_heatmap.png"),
        dotplot    = os.path.join(config["dirs"]["figures"],    "nerve_tumor_dotplot.png"),
        provenance = os.path.join(config["dirs"]["provenance"], "nerve_tumor_interaction_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_tumor_interaction.log"),
    conda:
        "../envs/scrna.yaml",
    threads: config["resources"]["default_threads"],
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
    params:
        random_seed = config["scrna"]["random_seed"],
    script:
        "../scripts/nerve_tumor_interaction.py"
