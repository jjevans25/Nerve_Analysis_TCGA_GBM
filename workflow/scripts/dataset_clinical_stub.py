"""Minimal clinical-metadata stub for a replication cohort (Stage A helper).

The reference cohort feeds nerve_cell_subset / nerve_clinical_association a GDC
clinical TSV keyed on ``file_uuid`` and joined onto ``obs['sample_id']``. External
cohorts (e.g. the CELLxGENE Census pull) have no GDC clinical record, so this
writes a schema-compatible stub — one row per sample, clinical columns blank —
that satisfies the join without inventing data. nerve_cell_subset.py fills the
resulting NaNs with the 'unknown' sentinel, so downstream clinical association
simply reports 'unknown' rather than crashing.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
samples = list(snakemake.params.samples)

CLINICAL_COLS = ["primary_diagnosis", "tumor_grade", "prior_malignancy",
                 "tissue_type", "gender", "race", "age_at_index"]

df = pd.DataFrame({"file_uuid": [str(s) for s in samples]})
for c in CLINICAL_COLS:
    df[c] = ""  # blank → NaN on reread → 'unknown' sentinel in nerve_cell_subset

Path(snakemake.output.clinical).parent.mkdir(parents=True, exist_ok=True)
df.to_csv(snakemake.output.clinical, sep="\t", index=False)
verify_artifact(snakemake.output.clinical, min_size_bytes=8)
log_transformation(log, "dataset_clinical_stub",
    f"Wrote clinical stub for {len(samples)} samples (columns blank → 'unknown')")

prov = stamp_artifact(
    output_path=snakemake.output.clinical,
    rule_name="dataset_clinical_stub",
    input_paths=[],
    tool_versions={"pandas": pd.__version__},
    parameters={"dataset": snakemake.params.dataset, "n_samples": len(samples)},
    description="Schema-compatible blank clinical stub for an external replication cohort",
    ontology_operation="operation:0004",
)
write_provenance(prov, snakemake.output.provenance)
log_transformation(log, "dataset_clinical_stub", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.clinical])
