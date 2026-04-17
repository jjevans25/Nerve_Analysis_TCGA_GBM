# workflow/rules/qc.smk
# Quality control rules for scRNA-seq samples.
# Each rule logs QC thresholds and registers provenance before filtering.

import os

rule scrna_qc:
    """Run per-sample QC filtering; logs thresholds before dropping any cells."""
    input:
        h5ad = os.path.join(config["dirs"]["data_processed"], "{sample}.h5ad"),
    output:
        h5ad       = os.path.join(config["dirs"]["data_processed"], "{sample}_qc.h5ad"),
        qc_metrics = os.path.join(config["dirs"]["tables"],         "{sample}_qc_metrics.csv"),
        provenance = os.path.join(config["dirs"]["provenance"],     "{sample}_qc_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{sample}_qc.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        min_genes    = config["scrna"]["min_genes"],
        max_genes    = config["scrna"]["max_genes"],
        min_cells    = config["scrna"]["min_cells"],
        max_pct_mito = config["scrna"]["max_pct_mito"],
    script:
        "../scripts/scrna_qc.py"


rule scrna_qc_report:
    """Aggregate per-sample QC metrics into a single summary table."""
    input:
        metrics = expand(
            os.path.join(config["dirs"]["tables"], "{sample}_qc_metrics.csv"),
            sample=config.get("samples", []),
        ),
    output:
        summary = os.path.join(config["dirs"]["tables"], "qc_summary.csv"),
    log:
        os.path.join(config["dirs"]["logs"], "qc_summary.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 4000,
        threads = 1,
    script:
        "../scripts/scrna_qc_report.py"
