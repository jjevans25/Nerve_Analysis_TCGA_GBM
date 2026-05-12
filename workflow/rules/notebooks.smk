# workflow/rules/notebooks.smk
# Snakemake rules for Marimo reactive notebook execution.
# Notebooks run in batch (export) mode to produce FAIR HTML artifacts.


rule explore_gbm_notebook:
    """Export interactive GBM exploration notebook to HTML artifact."""
    input:
        notebook  = "notebooks/01_explore_gbm_data.py",
        manifest  = config["loom_manifest"],
    output:
        html       = os.path.join(config["dirs"]["figures"],    "01_explore_gbm_data.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "explore_gbm_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "explore_gbm_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    script:
        "../scripts/run_notebook_export.py"


rule nerve_enrichment_notebook:
    """Export nerve-cell GSEA enrichment explorer to HTML (consumes nerve_enrichment_with_qc.csv)."""
    input:
        notebook   = "notebooks/02_nerve_enrichment_explorer.py",
        enrichment = os.path.join(config["dirs"]["tables"], "nerve_enrichment_with_qc.csv"),
        markers    = os.path.join(config["dirs"]["tables"], "nerve_cluster_markers_with_qc.csv"),
    output:
        html       = os.path.join(config["dirs"]["figures"],    "02_nerve_enrichment_explorer.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "nerve_enrichment_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_enrichment_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    script:
        "../scripts/run_notebook_export.py"
