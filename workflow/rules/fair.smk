# workflow/rules/fair.smk
# FAIR compliance rules: provenance indexing, metadata validation, report generation.

import os

rule fair_validate_metadata:
    """Verify that every processed artifact has an associated provenance JSON."""
    input:
        provenance_dir = config["dirs"]["provenance"],
    output:
        validation_report = os.path.join(config["dirs"]["results"], "fair_validation_report.json"),
    log:
        os.path.join(config["dirs"]["logs"], "fair_validate_metadata.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    script:
        "../scripts/fair_validate_metadata.py"


rule fair_snakemake_report:
    """Generate a full Snakemake provenance report for the current pipeline run."""
    output:
        report = os.path.join(config["dirs"]["results"], "snakemake_report.html"),
    log:
        os.path.join(config["dirs"]["logs"], "snakemake_report.log"),
    conda:
        "../envs/notebooks.yaml",
    shell:
        "snakemake --report {output.report} 2> {log}"
