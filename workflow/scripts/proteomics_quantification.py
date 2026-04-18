"""Aggregate per-sample peptide tables into a protein-level quantification matrix."""

import sys

import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
log_transformation(log, "proteomics_quantification",
                   f"Aggregating {len(snakemake.input.results)} sample result files")

frames = [pd.read_csv(f) for f in snakemake.input.results]
combined = pd.concat(frames, ignore_index=True)

# Pivot to protein × sample matrix using MaxLFQ intensity if present
intensity_col = "intensity" if "intensity" in combined.columns else combined.columns[-1]
protein_col   = "protein" if "protein" in combined.columns else combined.columns[0]

quant_matrix = combined.pivot_table(
    index=protein_col,
    columns="sample",
    values=intensity_col,
    aggfunc="sum",
    fill_value=0,
)

quant_matrix.to_csv(snakemake.output.quant_matrix)
verify_artifact(snakemake.output.quant_matrix, min_size_bytes=256)

log_transformation(log, "proteomics_quantification",
    f"Quantification matrix: {quant_matrix.shape[0]} proteins × {quant_matrix.shape[1]} samples")

# --- FAIR provenance ----------------------------------------------------------
import alphapept
prov = stamp_artifact(
    output_path=snakemake.output.quant_matrix,
    rule_name="proteomics_quantification",
    input_paths=list(snakemake.input.results),
    tool_versions={"alphapept": alphapept.__version__, "pandas": pd.__version__},
    parameters={"aggregation": "sum", "intensity_column": intensity_col},
    description="Protein-level quantification matrix across all MS samples",
    ontology_operation="operation:3630",  # EDAM: Protein quantification
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "proteomics_quantification", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.quant_matrix])
