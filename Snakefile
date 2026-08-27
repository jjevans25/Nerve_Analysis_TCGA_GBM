# =============================================================
# Snakefile — Base Pipeline
# Autonomous Agentic Biomedical Research Environment
# Apple M4 Max | Snakemake 9.x | FAIR4RS compliant
#
# Usage:
#   snakemake --use-conda --cores all
#   snakemake --use-conda --cores all --report results/snakemake_report.html
#   snakemake --list   # show all available rules
# =============================================================

# -------------------------------------------------------------
# Configuration
# -------------------------------------------------------------
configfile: "config/config.yaml"

# -------------------------------------------------------------
# Modular rule imports (common must be first)
# -------------------------------------------------------------
include: "workflow/rules/common.smk"
include: "workflow/rules/fair.smk"
include: "workflow/rules/ingest.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/integration.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/nerve_cells.smk"
include: "workflow/rules/immune.smk"
include: "workflow/rules/proteomics.smk"
include: "workflow/rules/notebooks.smk"
include: "workflow/rules/datasets.smk"
include: "workflow/rules/leads.smk"

# Ensure required directories exist before any rule runs
from pathlib import Path
for _d in config["dirs"].values():
    Path(_d).mkdir(parents=True, exist_ok=True)

# Convenience aliases
SAMPLES    = config.get("samples", [])
MS_SAMPLES = config.get("ms_samples", [])
DATASETS   = config.get("datasets", {})   # replication cohorts (cohort-namespaced track)

# Utility rules that run locally (no cluster submission)
localrules: all, clean_logs, show_dag, list_artifacts

# -------------------------------------------------------------
# Default target: all terminal artifacts + FAIR validation
# When no samples are defined, only the FAIR report runs.
# -------------------------------------------------------------
rule all:
    input:
        p(config["dirs"]["results"], "fair_validation_report.json"),
        *([p(config["dirs"]["tables"],         "qc_summary.csv")]                   if SAMPLES    else []),
        # `pinned_target`, not a bare path: this is a v1.3.0 reference artifact like the
        # nerve_* targets below, and it was the one root-level `.h5ad` still demanded
        # unconditionally. While `baseline.pinned` is true its producer (scrna_integration)
        # IS still defined, so a missing file here schedules an ~11 h scVI retrain of the
        # very cohort the pin exists to protect — and cascades back through scrna_qc to the
        # GDC looms. Guarding it lets the 67 GB reference chain be reclaimed from disk
        # without `rule all` trying to rebuild it. See CHANGELOG 2026-08-27.
        *(pinned_target(config["dirs"]["data_processed"], "integrated_latent.h5ad")  if SAMPLES    else []),
        *([p(config["dirs"]["data_external"],  "gdc_clinical.tsv")]                 if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "annotation_summary.csv")]           if SAMPLES    else []),
        *([p(config["dirs"]["figures"],        "cnv_heatmap.png")]                  if SAMPLES    else []),
        *(pinned_target(config["dirs"]["figures"],        "nerve_cells_umap.png")             if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_cluster_markers.csv")        if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_enrichment.csv")             if SAMPLES    else []),
        *(pinned_target(config["dirs"]["figures"],        "nerve_dotplot.png")                if SAMPLES    else []),
        *(pinned_target(config["dirs"]["figures"],        "nerve_abundance_heatmap.png")      if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_cluster_sample_purity.csv")  if SAMPLES    else []),
        *(pinned_target(config["dirs"]["figures"],        "nerve_cells_umap_by_sample.png")   if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_leiden_resolution_sweep.csv") if SAMPLES   else []),
        *(pinned_target(config["dirs"]["figures"],        "nerve_leiden_resolution_sweep.png") if SAMPLES   else []),
        *([p(config["dirs"]["tables"],         "nerve_celltype_label_summary.csv")] if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_cluster_sample_purity_v2.csv")] if SAMPLES  else []),
        *([p(config["dirs"]["figures"],        "nerve_cells_umap_by_sample_v2.png")] if SAMPLES   else []),
        *([p(config["dirs"]["figures"],        "nerve_scanvi_training_curves.png")] if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_cluster_annotations.csv")    if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_clinical_association.csv")   if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_tumor_interactions.csv")     if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_enrichment_with_qc.csv")     if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_cluster_markers_with_qc.csv") if SAMPLES   else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_tumor_interactions_with_qc.csv") if SAMPLES else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_tumor_top_pairs_with_qc.csv") if SAMPLES   else []),
        *([p(config["dirs"]["figures"],        "immune_cells_umap.png")]            if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "immune_cluster_annotations.csv")]   if SAMPLES    else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_tumor_immune_interactions.csv") if SAMPLES else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_tumor_immune_interactions_with_qc.csv") if SAMPLES else []),
        *(pinned_target(config["dirs"]["tables"],         "nerve_tumor_immune_top_pairs_with_qc.csv") if SAMPLES else []),
        *([p(config["dirs"]["tables"],         "protein_quant_matrix.csv")]         if MS_SAMPLES else []),
        # ── Reference-cohort notebooks: DEMOTED 2026-08-19, deliberately not built ──
        # 01/02/03/04 and nerve_tumor_exploration all read the pinned v1.3.0 reference,
        # whose nerve compartment was built by the logic `nerve_cell_subset.py` removed
        # on 2026-08-06 (it measured 59% malignant / 11% neural). Every nerve claim they
        # render is void, and they cannot be corrected: `data/processed/nerve_cells.h5ad`
        # is gone and v1.3.0 is pinned and unreproducible.
        #
        # Their rules are kept in workflow/rules/notebooks.smk and can be invoked
        # explicitly, so the code stays reusable if a v1.4.0 baseline is ever rebuilt
        # through the corrected pipeline. What is removed is the automatic rebuild —
        # leaving them here regenerated superseded HTMLs into results/figures/ next to
        # notebook 05's, where they looked equally current.
        #
        #   scripts/run_snakemake.sh results/figures/04_tme_nerve_immune_explorer.html \
        #     --use-conda --cores 1 --allowed-rules tme_nerve_immune_notebook
        # Replication cohorts: terminal concordance target pulls each dataset's
        # full cohort-namespaced chain (Stage A→D). Reference outputs above are
        # untouched. Gated on the `datasets:` config block being present.
        *[p(config["dirs"]["tables"], d, "cohort_concordance_summary.json") for d in DATASETS],
        # Test Oracle: compartment masks vs Census author annotation. Gates the
        # scientific validity of every interaction table above.
        *[p(config["dirs"]["tables"], d, "compartment_audit_gates.csv") for d in DATASETS],
        # scANVI-v2 nerve branch (full parity; side-branch not pulled by concordance).
        *[p(config["dirs"]["tables"],  d, "nerve_cluster_sample_purity_v2.csv") for d in DATASETS],
        *[p(config["dirs"]["tables"],  d, "nerve_celltype_label_summary.csv")   for d in DATASETS],
        *[p(config["dirs"]["figures"], d, "nerve_scanvi_training_curves.png")    for d in DATASETS],
        # Replication-cohort explorers. These five (per arm) are the notebook set —
        # the reference-cohort notebooks were archived on 2026-08-19, see above.
        *[p(config["dirs"]["figures"], d, "01_census_cohort_qc.html")           for d in DATASETS],
        *[p(config["dirs"]["figures"], d, "02_census_nerve_enrichment.html")    for d in DATASETS],
        *[p(config["dirs"]["figures"], d, "05_census_nerve_immune_explorer.html") for d in DATASETS],
        *[p(config["dirs"]["figures"], d, "06_census_compartment_audit.html")   for d in DATASETS],
        # Cross-arm lead shortlist. Un-wildcarded: spans both census arms, so it
        # sits at the tables root rather than in either arm's namespace. Gated on
        # both arms being configured, since it intersects them.
        *([p(config["dirs"]["tables"], "nerve_immune_lead_axes_postfix.csv")]
          if {config["lead_axes"]["full_arm"], config["lead_axes"]["capped_arm"]} <= set(DATASETS)
          else []),


# -------------------------------------------------------------
# Utility rules
# -------------------------------------------------------------

rule clean_logs:
    """Remove all log files to start a fresh run."""
    output:
        touch(p(config["dirs"]["logs"], "clean_logs.done")),
    log:
        p(config["dirs"]["logs"], "clean_logs.log"),
    conda:
        "workflow/envs/base.yaml",
    params:
        log_dir = lambda w, output: str(Path(output[0]).parent),
    shell:
        "find {params.log_dir} -name '*.log' -not -name 'clean_logs.log' -delete 2> {log}"


rule show_dag:
    """Render the pipeline DAG to a PNG for visual inspection."""
    output:
        dag = p("results", "pipeline_dag.png"),
    log:
        p(config["dirs"]["logs"], "show_dag.log"),
    conda:
        "workflow/envs/base.yaml",
    shell:
        "snakemake --dag 2> {log} | dot -Tpng > {output.dag}"


rule list_artifacts:
    """Print all provenance records currently registered."""
    output:
        touch(p(config["dirs"]["logs"], "list_artifacts.done")),
    log:
        p(config["dirs"]["logs"], "list_artifacts.log"),
    conda:
        "workflow/envs/base.yaml",
    params:
        provenance_dir = config["dirs"]["provenance"],
    script:
        "workflow/scripts/list_artifacts.py"
