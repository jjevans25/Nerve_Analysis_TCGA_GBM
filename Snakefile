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
include: "workflow/rules/proteomics.smk"
include: "workflow/rules/notebooks.smk"

# Ensure required directories exist before any rule runs
from pathlib import Path
for _d in config["dirs"].values():
    Path(_d).mkdir(parents=True, exist_ok=True)

# Convenience aliases
SAMPLES    = config.get("samples", [])
MS_SAMPLES = config.get("ms_samples", [])

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
        *([p(config["dirs"]["data_processed"], "integrated_latent.h5ad")]           if SAMPLES    else []),
        *([p(config["dirs"]["data_external"],  "gdc_clinical.tsv")]                 if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "annotation_summary.csv")]           if SAMPLES    else []),
        *([p(config["dirs"]["figures"],        "cnv_heatmap.png")]                  if SAMPLES    else []),
        *([p(config["dirs"]["figures"],        "nerve_cells_umap.png")]             if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_cluster_markers.csv")]        if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_enrichment.csv")]             if SAMPLES    else []),
        *([p(config["dirs"]["figures"],        "nerve_dotplot.png")]                if SAMPLES    else []),
        *([p(config["dirs"]["figures"],        "nerve_abundance_heatmap.png")]      if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_cluster_sample_purity.csv")]  if SAMPLES    else []),
        *([p(config["dirs"]["figures"],        "nerve_cells_umap_by_sample.png")]   if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_cluster_annotations.csv")]    if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_clinical_association.csv")]   if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_tumor_interactions.csv")]     if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_enrichment_with_qc.csv")]     if SAMPLES    else []),
        *([p(config["dirs"]["tables"],         "nerve_cluster_markers_with_qc.csv")] if SAMPLES   else []),
        *([p(config["dirs"]["tables"],         "nerve_tumor_interactions_with_qc.csv")] if SAMPLES else []),
        *([p(config["dirs"]["tables"],         "nerve_tumor_top_pairs_with_qc.csv")] if SAMPLES   else []),
        *([p(config["dirs"]["tables"],         "protein_quant_matrix.csv")]         if MS_SAMPLES else []),
        *([p(config["dirs"]["figures"],        "01_explore_gbm_data.html")]         if SAMPLES    else []),
        *([p(config["dirs"]["figures"],        "02_nerve_enrichment_explorer.html")] if SAMPLES    else []),


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
