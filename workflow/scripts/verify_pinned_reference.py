"""Detect drift in the pinned v1.3.0 reference artifacts.

Answers "has anything in the frozen reference cohort changed since it was
pinned?". Because the producing rules are not defined while ``baseline.pinned``
is true, Snakemake itself will never notice a modified or deleted pinned table —
this rule is the only check.

Re-hashes every artifact in ``provenance/pinned_reference_<version>.json`` and
raises on any mismatch, missing file, or artifact that is configured in
``baseline.pinned_artifacts`` but absent from the manifest.

Run on demand: ``snakemake --use-conda verify_pinned_reference``.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "workflow/scripts")
from fair_utils import file_sha256, log_transformation, verify_artifact  # noqa: E402

log = snakemake.log[0]  # type: ignore[name-defined]
manifest_path = Path(snakemake.input.manifest)  # type: ignore[name-defined]
configured: list[str] = list(snakemake.params.artifacts)  # type: ignore[name-defined]
out_path = Path(snakemake.output.report)  # type: ignore[name-defined]

manifest = json.loads(manifest_path.read_text())
entries = manifest["artifacts"]

log_transformation(
    log, "verify_pinned_reference",
    f"Verifying {len(entries)} artifacts against {manifest_path} "
    f"(frozen {manifest.get('frozen_at_utc')})",
)

ok: list[str] = []
modified: list[dict] = []
missing: list[str] = []

for entry in entries:
    path = Path(entry["path"])
    if not path.exists():
        missing.append(entry["path"])
        continue
    actual = file_sha256(path)
    if actual != entry["sha256"]:
        modified.append({
            "path": entry["path"],
            "expected_sha256": entry["sha256"],
            "actual_sha256": actual,
            "expected_size": entry["size_bytes"],
            "actual_size": path.stat().st_size,
        })
    else:
        ok.append(entry["path"])

# An artifact added to config but never re-frozen is unprotected, not "clean".
manifest_paths = {e["path"] for e in entries}
unmanifested = sorted(set(configured) - manifest_paths)

passed = not (modified or missing or unmanifested)
report = {
    "checked_at_utc": datetime.now(timezone.utc).isoformat(),
    "manifest": str(manifest_path),
    "baseline_version": manifest.get("baseline_version"),
    "manifest_frozen_at_utc": manifest.get("frozen_at_utc"),
    "pass": passed,
    "n_checked": len(entries),
    "n_unchanged": len(ok),
    "n_modified": len(modified),
    "n_missing": len(missing),
    "n_configured_but_unmanifested": len(unmanifested),
    "modified": modified,
    "missing": missing,
    "configured_but_unmanifested": unmanifested,
}

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(json.dumps(report, indent=2) + "\n")
verify_artifact(out_path, min_size_bytes=50)

if passed:
    log_transformation(
        log, "verify_pinned_reference",
        f"PASS — all {len(ok)} pinned artifacts unchanged",
        artifact_paths=[str(out_path)],
    )
else:
    detail = (
        f"{len(modified)} modified, {len(missing)} missing, "
        f"{len(unmanifested)} configured but not in manifest"
    )
    # Snakemake deletes a failed job's outputs, so `report` will not survive this
    # raise. The log is preserved on failure — embed the full finding there so the
    # diagnostic is not lost with the file it was written to.
    log_transformation(
        log, "verify_pinned_reference",
        f"[FAIR-ALERT] Pinned reference drift detected — {detail}. "
        f"Full report (the output file is removed by Snakemake on failure): "
        + json.dumps(report, indent=2),
        status="FAILURE",
        artifact_paths=[str(out_path)],
    )
    raise RuntimeError(
        f"Pinned v1.3.0 reference has drifted ({detail}). See {out_path}. "
        "These artifacts are unreproducible — restore from git/backup rather "
        "than re-running the pipeline."
    )
