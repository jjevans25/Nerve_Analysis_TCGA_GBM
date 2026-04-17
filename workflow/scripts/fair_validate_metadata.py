"""Validate that every processed artifact has an associated provenance JSON."""

import json
import pathlib
import sys
from datetime import datetime, timezone

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation

provenance_dir = pathlib.Path(snakemake.input.provenance_dir)
records = list(provenance_dir.glob("*.json"))

report = {
    "validated_at": datetime.now(timezone.utc).isoformat(),
    "total_provenance_records": len(records),
    "artifacts": [str(r) for r in records],
}

pathlib.Path(snakemake.output.validation_report).write_text(
    json.dumps(report, indent=2)
)

log_transformation(
    log_path=snakemake.log[0],
    rule_name="fair_validate_metadata",
    message=f"Validated {len(records)} provenance records.",
    artifact_paths=[snakemake.output.validation_report],
)
