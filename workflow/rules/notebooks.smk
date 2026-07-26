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


rule nerve_tumor_exploration_notebook:
    """Export nerve-tumor LIANA LR interaction explorer to HTML (consumes nerve_tumor_*_with_qc.csv)."""
    input:
        notebook     = "notebooks/nerve_tumor_exploration.py",
        interactions = os.path.join(config["dirs"]["tables"], "nerve_tumor_interactions_with_qc.csv"),
        top_pairs    = os.path.join(config["dirs"]["tables"], "nerve_tumor_top_pairs_with_qc.csv"),
    output:
        html       = os.path.join(config["dirs"]["figures"],    "nerve_tumor_exploration.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "nerve_tumor_exploration_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_tumor_exploration_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    script:
        "../scripts/run_notebook_export.py"


rule tme_nerve_immune_notebook:
    """Export the TME context explorer to HTML (compartment census, patient
    composition, cell-type-labelled LR browser, curated lead targets).

    Reference cohort only: the CELLxGENE Census LIANA tables were computed on raw
    UMI counts rather than log1p data, so they are excluded until re-run. See
    markdowns/blocker_census_liana_raw_counts.md.

    This rule's inputs are all pinned v1.3.0 artifacts. `--allowed-rules` is no
    longer needed to protect them: with `baseline.pinned` set, the rules that
    would rebuild them are not defined (see common.smk BASELINE_PINNED), so no
    invocation can schedule the nerve cascade. `--rerun-triggers mtime` is still
    worth passing to avoid unrelated full-pipeline re-runs.
    """
    input:
        notebook           = "notebooks/04_tme_nerve_immune_explorer.py",
        interactions       = os.path.join(config["dirs"]["tables"], "nerve_tumor_immune_interactions_with_qc.csv"),
        top_pairs          = os.path.join(config["dirs"]["tables"], "nerve_tumor_immune_top_pairs_with_qc.csv"),
        nerve_annotations  = os.path.join(config["dirs"]["tables"], "nerve_cluster_annotations.csv"),
        immune_annotations = os.path.join(config["dirs"]["tables"], "immune_cluster_annotations.csv"),
        immune_composition = os.path.join(config["dirs"]["tables"], "immune_cluster_composition.csv"),
        nerve_purity       = os.path.join(config["dirs"]["tables"], "nerve_cluster_sample_purity.csv"),
        immune_purity      = os.path.join(config["dirs"]["tables"], "immune_subtype_sample_purity.csv"),
        annotation_summary = os.path.join(config["dirs"]["tables"], "annotation_summary.csv"),
        lead_targets       = os.path.join(config["dirs"]["tables"], "nerve_crosstalk_lead_targets.csv"),
        clinical           = os.path.join(config["dirs"]["tables"], "nerve_clinical_association.csv"),
        lr_provenance      = os.path.join(config["dirs"]["provenance"], "nerve_tumor_immune_interaction_provenance.json"),
    output:
        html       = os.path.join(config["dirs"]["figures"],    "04_tme_nerve_immune_explorer.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "tme_nerve_immune_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "tme_nerve_immune_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    script:
        "../scripts/run_notebook_export.py"


rule nerve_tumor_immune_notebook:
    """Export three-way nerve-tumor-immune LR explorer to HTML (consumes nerve_tumor_immune_*_with_qc.csv)."""
    input:
        notebook     = "notebooks/03_nerve_tumor_immune_explorer.py",
        interactions = os.path.join(config["dirs"]["tables"], "nerve_tumor_immune_interactions_with_qc.csv"),
        top_pairs    = os.path.join(config["dirs"]["tables"], "nerve_tumor_immune_top_pairs_with_qc.csv"),
    output:
        html       = os.path.join(config["dirs"]["figures"],    "03_nerve_tumor_immune_explorer.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "nerve_tumor_immune_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_tumor_immune_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    script:
        "../scripts/run_notebook_export.py"
