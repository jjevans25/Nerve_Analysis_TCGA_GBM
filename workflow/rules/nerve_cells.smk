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
