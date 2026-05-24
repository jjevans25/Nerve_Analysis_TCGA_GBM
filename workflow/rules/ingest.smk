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


def _all_looms(_wildcards=None) -> list[str]:
    """Resolve every sample's loom path — used to build the gene-symbol map once."""
    return [_find_loom(type("W", (), {"sample": s})) for s in config.get("samples", [])]


rule build_gene_symbol_map:
    """Query MyGene.info once for all loom Ensembl IDs; write a TSV cache reused by every sample."""
    input:
        looms = _all_looms,
    output:
        cache      = config["gene_symbol_map"]["cache_tsv"],
        provenance = os.path.join(config["dirs"]["provenance"], "gene_symbol_map_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "build_gene_symbol_map.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 4000,
        threads = 1,
    params:
        ensembl_release = config["databases"]["ensembl_release"],
        chunk_size      = config["gene_symbol_map"]["chunk_size"],
        request_timeout = config["gene_symbol_map"]["request_timeout"],
    script:
        "../scripts/build_gene_symbol_map.py"


rule loom_to_h5ad:
    """Convert GDC loom to AnnData h5ad; assigns batch metadata + MyGene symbol map.

    Nerve-cell marker annotation (is_nerve_marker + gene_presence) was moved
    to scrna_qc in v1.2.0 so panel-config changes do not trigger sample-level
    reruns. See markdowns/failing_cluster_diagnosis.md v1.1.0 supplement.
    """
    wildcard_constraints:
        sample = r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
    input:
        loom       = _find_loom,
        symbol_map = config["gene_symbol_map"]["cache_tsv"],
    output:
        h5ad       = os.path.join(config["dirs"]["data_processed"], "{sample}.h5ad"),
        provenance = os.path.join(config["dirs"]["provenance"],     "{sample}_ingest_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{sample}_ingest.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 2,
    params:
        sample_id = lambda wc: wc.sample,
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
    retries: 2
    resources:
        mem_mb  = 2000,
        threads = 1,
    params:
        sample_ids      = config.get("samples", []),
        api_base        = config["gdc_api"]["base_url"],
        api_fields      = config["gdc_api"]["fields"],
        request_timeout = config["gdc_api"]["request_timeout"],
        max_retries     = config["gdc_api"]["max_retries"],
    script:
        "../scripts/gdc_clinical_fetch.py"
