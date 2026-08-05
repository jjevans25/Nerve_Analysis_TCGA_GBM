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


rule download_gene_positions:
    """Fetch Ensembl gene coordinates so CNV inference can order genes along the genome.

    Cohort-independent (coordinates are a property of the assembly, not the
    samples), so both Census arms and the reference cohort share one copy —
    same pattern as the MSigDB download.
    """
    output:
        gtf       = os.path.join(config["gene_positions"]["download_dir"],
                                 config["gene_positions"]["gtf_filename"]),
        positions = os.path.join(config["gene_positions"]["download_dir"],
                                 config["gene_positions"]["positions_filename"]),
        manifest  = os.path.join(config["gene_positions"]["download_dir"],
                                 "gene_positions_manifest.json"),
        provenance = os.path.join(config["dirs"]["provenance"], "gene_positions_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "download_gene_positions.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 4000,
        threads = 1,
    params:
        release         = config["gene_positions"]["release"],
        url             = config["gene_positions"]["url"],
        sha256          = config["gene_positions"]["sha256"],
        request_timeout = config["gene_positions"]["request_timeout"],
        max_retries     = config["gene_positions"]["max_retries"],
    script:
        "../scripts/download_gene_positions.py"


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


rule freeze_baseline_provenance:
    """Bundle every per-rule provenance JSON into a single tracked baseline file.

    Run on demand (not auto-triggered): `snakemake --use-conda freeze_baseline_provenance`.
    The output is the only file inside `provenance/` that gets committed to git
    (see `.gitignore` exception for `provenance/baseline_*.json`).
    """
    output:
        bundle = os.path.join(config["dirs"]["provenance"],
                              f"baseline_{config['baseline']['version']}.json"),
    log:
        os.path.join(config["dirs"]["logs"], "freeze_baseline_provenance.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = 1000,
        threads = 1,
    params:
        version        = config["baseline"]["version"],
        # Wrapped in a function to deactivate Snakemake's automatic wildcard
        # expansion: the summary text may contain literal braces (e.g. cluster
        # sets like {0-14, 16-27} or size tuples {773,692,378,11}) that would
        # otherwise be misparsed as {wildcard} placeholders.
        summary        = lambda wildcards: config["baseline"]["summary"],
        provenance_dir = config["dirs"]["provenance"],
    script:
        "../scripts/freeze_baseline_provenance.py"


rule freeze_pinned_reference:
    """Hash every pinned v1.3.0 artifact into a drift-detection manifest.

    Run on demand after any deliberate change to `baseline.pinned_artifacts`:
    `snakemake --use-conda freeze_pinned_reference`. Not referenced by `rule all` —
    re-freezing must be an explicit act, otherwise drift would be silently
    absorbed into a new manifest.
    """
    output:
        manifest = os.path.join(config["dirs"]["provenance"],
                                f"pinned_reference_{config['baseline']['version']}.json"),
    log:
        os.path.join(config["dirs"]["logs"], "freeze_pinned_reference.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    params:
        version   = config["baseline"]["version"],
        artifacts = config["baseline"].get("pinned_artifacts", []),
    script:
        "../scripts/freeze_pinned_reference.py"


rule verify_pinned_reference:
    """Re-hash the pinned v1.3.0 artifacts and fail on any drift.

    The only drift detector for these files: while `baseline.pinned` is true their
    producing rules are not defined, so Snakemake will never notice a modified or
    deleted pinned table. Run on demand:
    `snakemake --use-conda verify_pinned_reference`.
    """
    input:
        manifest = os.path.join(config["dirs"]["provenance"],
                                f"pinned_reference_{config['baseline']['version']}.json"),
    output:
        report = os.path.join(config["dirs"]["results"], "pinned_reference_verification.json"),
    log:
        os.path.join(config["dirs"]["logs"], "verify_pinned_reference.log"),
    conda:
        "../envs/notebooks.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    params:
        artifacts = config["baseline"].get("pinned_artifacts", []),
    script:
        "../scripts/verify_pinned_reference.py"
