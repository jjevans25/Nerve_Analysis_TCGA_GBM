# workflow/rules/immune.smk
# Tumor-associated immune compartment subclustering + three-way
# nerve-tumor-immune cell-cell communication.

import os


rule immune_cell_subset:
    """Extract non-malignant immune cells; re-cluster in the immune subspace on X_scVI."""
    input:
        h5ad = os.path.join(config["dirs"]["data_processed"], "malignancy_labeled.h5ad"),
    output:
        h5ad        = os.path.join(config["dirs"]["data_processed"], "immune_cells.h5ad"),
        umap        = os.path.join(config["dirs"]["figures"],        "immune_cells_umap.png"),
        composition = os.path.join(config["dirs"]["tables"],         "immune_cluster_composition.csv"),
        purity      = os.path.join(config["dirs"]["tables"],         "immune_cluster_sample_purity.csv"),
        provenance  = os.path.join(config["dirs"]["provenance"],     "immune_cell_subset_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "immune_cell_subset.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        source_label              = config["immune_cells"]["source_label"],
        leiden_resolution         = config["immune_cells"]["leiden_resolution"],
        n_top_genes               = config["scrna"]["n_top_genes"],
        random_seed               = config["scrna"]["random_seed"],
        dominant_fraction_max     = config["immune_cells"]["batch_qc"]["dominant_fraction_max"],
        min_contributing_fraction = config["immune_cells"]["batch_qc"]["min_contributing_fraction"],
        min_contributing_samples  = config["immune_cells"]["batch_qc"]["min_contributing_samples"],
    script:
        "../scripts/immune_cell_subset.py"


rule immune_cluster_annotations:
    """Assign immune subtype labels (microglia/TAM/T/NK/DC) per immune Leiden cluster."""
    input:
        h5ad = os.path.join(config["dirs"]["data_processed"], "immune_cells.h5ad"),
    output:
        h5ad          = os.path.join(config["dirs"]["data_processed"], "immune_cells_labeled.h5ad"),
        annotations   = os.path.join(config["dirs"]["tables"],         "immune_cluster_annotations.csv"),
        subtype_purity= os.path.join(config["dirs"]["tables"],         "immune_subtype_sample_purity.csv"),
        provenance    = os.path.join(config["dirs"]["provenance"],     "immune_cluster_annotations_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "immune_cluster_annotations.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        subtype_markers           = config["immune_cells"]["subtype_markers"],
        random_seed               = config["scrna"]["random_seed"],
        dominant_fraction_max     = config["immune_cells"]["batch_qc"]["dominant_fraction_max"],
        min_contributing_fraction = config["immune_cells"]["batch_qc"]["min_contributing_fraction"],
        min_contributing_samples  = config["immune_cells"]["batch_qc"]["min_contributing_samples"],
    script:
        "../scripts/immune_cluster_annotations.py"


if not BASELINE_PINNED:
    # PINNED (baseline.pinned) — regenerating this would rewrite frozen v1.3.0
    # artifacts and/or needs the deleted, unreproducible nerve_cells.h5ad.
    rule nerve_tumor_immune_interaction:
        """Three-way nerve-tumor-immune LR communication via LIANA+ consensus rank."""
        input:
            malig  = os.path.join(config["dirs"]["data_processed"], "malignancy_labeled.h5ad"),
            nerve  = os.path.join(config["dirs"]["data_processed"], "nerve_cells.h5ad"),
            immune = os.path.join(config["dirs"]["data_processed"], "immune_cells_labeled.h5ad"),
        output:
            lr_table   = os.path.join(config["dirs"]["tables"],     "nerve_tumor_immune_interactions.csv"),
            top_pairs  = os.path.join(config["dirs"]["tables"],     "nerve_tumor_immune_top_pairs.csv"),
            heatmap    = os.path.join(config["dirs"]["figures"],    "nerve_tumor_immune_sig_heatmap.png"),
            dotplot    = os.path.join(config["dirs"]["figures"],    "nerve_tumor_immune_dotplot.png"),
            provenance = os.path.join(config["dirs"]["provenance"], "nerve_tumor_immune_interaction_provenance.json"),
        log:
            os.path.join(config["dirs"]["logs"], "nerve_tumor_immune_interaction.log"),
        conda:
            "../envs/scrna.yaml",
        threads: config["resources"]["default_threads"],
        resources:
            mem_mb = 64000,
        params:
            random_seed = config["scrna"]["random_seed"],
            # Reference .X is Seurat SCT log1p already — do not re-normalize.
            normalize_counts = False,
        script:
            "../scripts/nerve_tumor_immune_interaction.py"
