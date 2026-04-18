"""AlphaPept spectral search and peptide ID. Runs Numba JIT on CPU/ARM64 — MPS not supported by Numba."""

import sys

import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
log_transformation(log, "proteomics_search",
    f"AlphaPept search on {snakemake.input.raw_ms} | "
    "Note: Numba JIT on ARM64 CPU (MPS not supported by Numba)")

# AlphaPept API
from alphapept.interface import run_complete_workflow
from alphapept.settings import load_settings

settings = load_settings()
settings["experiment"]["fasta_paths"] = [snakemake.input.fasta]
settings["experiment"]["file_paths"]  = [snakemake.input.raw_ms]

# Peptide search parameters from config
settings["search"]["peptide_fdr"]        = snakemake.params.fdr
settings["search"]["min_pep_length"]     = 7
settings["search"]["missed_cleavages"]   = snakemake.params.missed_cleavages
settings["search"]["enzyme"]             = snakemake.params.enzyme
settings["search"]["mass_tolerance_ppm"] = snakemake.params.mass_tolerance_ppm

log_transformation(log, "proteomics_search",
    f"FDR={snakemake.params.fdr}, enzyme={snakemake.params.enzyme}, "
    f"missed_cleavages={snakemake.params.missed_cleavages}")

results = run_complete_workflow(settings)

# Extract peptide-level results to CSV
peptide_df = results.get("peptide_table", pd.DataFrame())
if peptide_df.empty:
    raise ValueError("[FAIR-ALERT] AlphaPept returned an empty peptide table — check raw MS file and FASTA.")

peptide_df.to_csv(snakemake.output.results, index=False)
verify_artifact(snakemake.output.results, min_size_bytes=512)

# --- FAIR provenance ----------------------------------------------------------
import alphapept
prov = stamp_artifact(
    output_path=snakemake.output.results,
    rule_name="proteomics_search",
    input_paths=[snakemake.input.raw_ms, snakemake.input.fasta],
    tool_versions={"alphapept": alphapept.__version__},
    parameters={
        "fdr":                snakemake.params.fdr,
        "enzyme":             snakemake.params.enzyme,
        "missed_cleavages":   snakemake.params.missed_cleavages,
        "mass_tolerance_ppm": snakemake.params.mass_tolerance_ppm,
        "acceleration":       "Numba JIT (ARM64 CPU)",
    },
    description="AlphaPept peptide identification from raw MS data",
    ontology_operation="operation:3767",  # EDAM: Protein identification
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "proteomics_search",
    f"Identified {len(peptide_df)} peptides", status="SUCCESS",
    artifact_paths=[snakemake.output.results])
