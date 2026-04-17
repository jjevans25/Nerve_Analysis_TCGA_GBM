"""Aggregate per-sample QC metrics into a single summary table."""

import sys
import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation

frames = [pd.read_csv(m) for m in snakemake.input.metrics]
summary = pd.concat(frames, ignore_index=True)
summary.to_csv(snakemake.output.summary, index=False)

log_transformation(
    log_path=snakemake.log[0],
    rule_name="scrna_qc_report",
    message=f"Aggregated {len(frames)} QC tables.",
    artifact_paths=[snakemake.output.summary],
)
