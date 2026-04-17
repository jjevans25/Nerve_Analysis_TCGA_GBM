# workflow/rules/ingest.smk
# Loom-to-AnnData staging and GDC clinical metadata fetch.

import glob
import os


def _find_loom(wildcards) -> str:
    """Resolve the GDC loom file path for a given sample UUID."""
    pattern = os.path.join(config["loom_dir"], wildcards.sample, "*.loom")
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No loom file found for sample {wildcards.sample}")
    return matches[0]


rule loom_to_h5ad:
    """Convert GDC loom to AnnData h5ad; assigns batch metadata and emits gene presence report."""
    wildcard_constraints:
        sample = r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
    input:
        loom = _find_loom,
    output:
        h5ad          = os.path.join(config["dirs"]["data_processed"], "{sample}.h5ad"),
        gene_presence = os.path.join(config["dirs"]["tables"],         "{sample}_gene_presence.csv"),
        provenance    = os.path.join(config["dirs"]["provenance"],     "{sample}_ingest_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{sample}_ingest.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 2,
    params:
        sample_id     = lambda wc: wc.sample,
        nerve_markers = config["nerve_cells"]["markers"],
    script:
        "../scripts/loom_to_h5ad.py"


rule gdc_clinical_fetch:
    """Fetch clinical metadata from the GDC public REST API for all samples."""
    output:
        clinical   = os.path.join(config["dirs"]["data_external"], "gdc_clinical.tsv"),
        provenance = os.path.join(config["dirs"]["provenance"],    "gdc_clinical_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "gdc_clinical_fetch.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    params:
        sample_ids = config.get("samples", []),
        api_base   = config["gdc_api"]["base_url"],
        api_fields = config["gdc_api"]["fields"],
    script:
        "../scripts/gdc_clinical_fetch.py"
