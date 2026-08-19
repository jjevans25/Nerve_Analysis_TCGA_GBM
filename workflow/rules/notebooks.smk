# workflow/rules/notebooks.smk
# Snakemake rules for Marimo reactive notebook execution.
# Notebooks run in batch (export) mode to produce FAIR HTML artifacts.



# ── Reference-cohort notebook rules REMOVED 2026-08-19 ────────────────────────
# explore_gbm_notebook, nerve_enrichment_notebook, nerve_tumor_exploration_notebook,
# tme_nerve_immune_notebook and nerve_tumor_immune_notebook are gone, and their
# notebooks now live in notebooks/archive/.
#
# All five read the pinned v1.3.0 reference, whose nerve compartment was built by the
# logic removed from nerve_cell_subset.py on 2026-08-06 (59% malignant / 11% neural).
# Their output cannot be corrected: the reference is pinned and unreproducible, and its
# v1.1.0 frozen-barcode insulator bypasses the compartment fix entirely, so re-running
# them would faithfully reproduce the defect.
#
# The rules were deleted rather than repointed at notebooks/archive/ so the pipeline has
# no path to regenerating void HTML. The notebooks remain readable in-tree and in git.
#
# Their Census successors are below.
# ─────────────────────────────────────────────────────────────────────────────


rule ds_census_cohort_qc_notebook:
    """Export the cohort QC and composition explorer to HTML.

    Census successor to the archived explore_gbm_notebook. The notebook globs the 170
    per-donor `*_qc_metrics.csv` / `*_gene_presence.csv` files rather than declaring
    them: naming 340 inputs would make the rule unreadable, and `directory()` on the
    tables dir would make every unrelated table a rerun trigger.

    `annotation_summary.csv` is the honest dependency edge — it is produced by
    ds_scrna_annotate, which sits downstream of ds_scrna_qc for every sample, so it
    cannot exist before the globbed files do.
    """
    input:
        notebook   = "notebooks/01_census_cohort_qc.py",
        annotation = os.path.join(config["dirs"]["tables"], "{dataset}", "annotation_summary.csv"),
    output:
        html       = os.path.join(config["dirs"]["figures"],    "{dataset}", "01_census_cohort_qc.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "{dataset}", "census_cohort_qc_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_census_cohort_qc_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = 8000,
        threads = 1,
    params:
        env = lambda wc: {"GBM_DATASET": wc.dataset},
    script:
        "../scripts/run_notebook_export.py"


rule ds_census_nerve_enrichment_notebook:
    """Export the nerve-cluster GSEA enrichment explorer to HTML.

    Census successor to the archived nerve_enrichment_notebook. Same question over a
    compartment that is 95.4% neural rather than 11%.
    """
    input:
        notebook    = "notebooks/02_census_nerve_enrichment.py",
        enrichment  = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_enrichment_with_qc.csv"),
        markers     = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_cluster_markers_with_qc.csv"),
        annotations = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_cluster_annotations.csv"),
    output:
        html       = os.path.join(config["dirs"]["figures"],    "{dataset}", "02_census_nerve_enrichment.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "{dataset}", "census_nerve_enrichment_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_census_nerve_enrichment_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = 4000,
        threads = 1,
    params:
        env = lambda wc: {"GBM_DATASET": wc.dataset},
    script:
        "../scripts/run_notebook_export.py"


rule ds_census_compartment_audit_notebook:
    """Export the compartment integrity audit (Test Oracle) explorer to HTML.

    No reference-cohort ancestor and cannot have one: the audit cross-tabulates each
    compartment against the CELLxGENE Census author annotation, and the 17-sample
    TCGA reference carries no such annotation. This notebook is the evidence that
    distinguishes the Census cohort from the archived one.
    """
    input:
        notebook       = "notebooks/06_census_compartment_audit.py",
        gates          = os.path.join(config["dirs"]["tables"], "{dataset}", "compartment_audit_gates.csv"),
        audit          = os.path.join(config["dirs"]["tables"], "{dataset}", "compartment_audit.csv"),
        confusion      = os.path.join(config["dirs"]["tables"], "{dataset}", "malignancy_confusion.csv"),
        cluster_audit  = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_compartment_cluster_audit.csv"),
        purity         = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_cluster_sample_purity.csv"),
        purity_v2      = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_cluster_sample_purity_v2.csv"),
        celltype_labels = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_celltype_label_summary.csv"),
    output:
        html       = os.path.join(config["dirs"]["figures"],    "{dataset}", "06_census_compartment_audit.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "{dataset}", "census_compartment_audit_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_census_compartment_audit_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = 4000,
        threads = 1,
    params:
        env = lambda wc: {"GBM_DATASET": wc.dataset},
    script:
        "../scripts/run_notebook_export.py"

rule ds_census_nerve_immune_notebook:
    """Export the replication-cohort nerve x immune explorer to HTML.

    Cohort-namespaced twin of tme_nerve_immune_notebook. Separate notebook rather
    than a cohort switch because the cohorts differ in what exists: no clinical
    metadata here, no per-cohort curated target list, 169 donors instead of 17.

    Requires the 2026-07-26 LIANA normalization fix — before it this cohort's
    tables were computed on raw UMI counts and were invalid. The notebook checks
    for the defect's signature at load time and refuses to vouch for the data.
    """
    input:
        notebook           = "notebooks/05_census_nerve_immune_explorer.py",
        interactions       = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_tumor_immune_interactions_with_qc.csv"),
        top_pairs          = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_tumor_immune_top_pairs_with_qc.csv"),
        nerve_annotations  = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_cluster_annotations.csv"),
        immune_annotations = os.path.join(config["dirs"]["tables"], "{dataset}", "immune_cluster_annotations.csv"),
        nerve_purity       = os.path.join(config["dirs"]["tables"], "{dataset}", "nerve_cluster_sample_purity.csv"),
        immune_cl_purity   = os.path.join(config["dirs"]["tables"], "{dataset}", "immune_cluster_sample_purity.csv"),
        immune_purity      = os.path.join(config["dirs"]["tables"], "{dataset}", "immune_subtype_sample_purity.csv"),
        annotation_summary = os.path.join(config["dirs"]["tables"], "{dataset}", "annotation_summary.csv"),
        lr_provenance      = os.path.join(config["dirs"]["provenance"], "{dataset}", "nerve_tumor_immune_interaction_provenance.json"),
        # Panel E's source. Un-wildcarded: it spans both Census arms. Declared so
        # rebuilding the shortlist re-renders the notebook — Panel E asserts the table
        # against a live recomputation, so a stale render hides a regression.
        lead_axes_postfix  = os.path.join(config["dirs"]["tables"], "nerve_immune_lead_axes_postfix.csv"),
        # Removed 2026-08-19 with Panel F:
        #   concordance / shared_pairs — ds_cohort_concordance still runs and still
        #     writes them; this notebook simply no longer reads them.
        #   lead_targets (nerve_crosstalk_lead_targets.csv) — the withdrawn 2026-07-15
        #     reference shortlist. It was declared but never opened by this notebook
        #     (only notebook 04 reads it), so it was a spurious dependency tying the
        #     Census notebook to the pinned reference cohort.
    output:
        html       = os.path.join(config["dirs"]["figures"],    "{dataset}", "05_census_nerve_immune_explorer.html"),
        provenance = os.path.join(config["dirs"]["provenance"], "{dataset}", "census_nerve_immune_notebook_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_census_nerve_immune_notebook.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    params:
        # marimo has no argv passthrough; the notebook reads GBM_DATASET.
        env = lambda wc: {"GBM_DATASET": wc.dataset},
    script:
        "../scripts/run_notebook_export.py"

