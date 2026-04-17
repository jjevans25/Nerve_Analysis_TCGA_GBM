# workflow/rules/proteomics.smk
# AlphaPept-based mass spectrometry processing.
# Runs in an isolated conda env separate from scRNA-seq tools.

import os

rule proteomics_search:
    """Run AlphaPept spectral matching and peptide identification on raw MS files."""
    input:
        raw_ms = os.path.join(config["dirs"]["data_raw"],      "{sample}.raw"),
        fasta  = os.path.join(config["dirs"]["data_external"], "proteome.fasta"),
    output:
        results    = os.path.join(config["dirs"]["data_processed"], "{sample}_ms_results.csv"),
        provenance = os.path.join(config["dirs"]["provenance"],     "{sample}_proteomics_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{sample}_proteomics_search.log"),
    conda:
        "../envs/proteomics.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        fdr                = config["proteomics"]["fdr_threshold"],
        min_peptides       = config["proteomics"]["min_peptides"],
        enzyme             = config["proteomics"]["enzyme"],
        missed_cleavages   = config["proteomics"]["missed_cleavages"],
        mass_tolerance_ppm = config["proteomics"]["mass_tolerance_ppm"],
    script:
        "../scripts/proteomics_search.py"


rule proteomics_quantification:
    """Aggregate per-sample peptide tables into a protein-level quantification matrix."""
    input:
        results = expand(
            os.path.join(config["dirs"]["data_processed"], "{sample}_ms_results.csv"),
            sample=config.get("ms_samples", []),
        ),
    output:
        quant_matrix = os.path.join(config["dirs"]["tables"],    "protein_quant_matrix.csv"),
        provenance   = os.path.join(config["dirs"]["provenance"], "proteomics_quant_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "proteomics_quantification.log"),
    conda:
        "../envs/proteomics.yaml",
    resources:
        mem_mb  = 8000,
        threads = config["resources"]["default_threads"],
    script:
        "../scripts/proteomics_quantification.py"
