"""Freeze a versioned baseline by bundling all per-rule provenance JSONs into a
single tracked artifact.

Per the project's `.gitignore`, individual `provenance/*.json` files are
gitignored — they're regenerated every time a rule runs. To make a baseline
result set citable + reproducible across collaborators, we collect every
existing provenance record into one bundle (`provenance/baseline_<version>.json`)
that IS tracked (via a `!provenance/baseline_*.json` exception in
`.gitignore`).

Output schema:
{
  "baseline_version":  "v1.0.0",
  "frozen_at_utc":     "2026-...",
  "git_commit":        "abc123...",
  "git_branch":        "main",
  "n_provenance_files": 49,
  "summary":           "...",
  "rules":             [ {rule, provenance_file, sha256, ...}, ... ],
}
"""

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "workflow/scripts")
from fair_utils import file_sha256, log_transformation


def _git(cmd: list[str]) -> str:
    """Run a `git` command and return stripped stdout (empty string on failure)."""
    try:
        return subprocess.check_output(
            ["git", *cmd], stderr=subprocess.DEVNULL, text=True
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


log = snakemake.log[0]
version: str = str(snakemake.params.version)
provenance_dir = Path(snakemake.params.provenance_dir)
output_path = Path(snakemake.output.bundle)
summary: str = str(snakemake.params.summary)

log_transformation(
    log,
    "freeze_baseline_provenance",
    f"Scanning {provenance_dir} for provenance records",
)

# Collect every provenance record except the bundle itself
prov_files = sorted(
    p for p in provenance_dir.glob("*.json")
    if not p.name.startswith("baseline_")
)
log_transformation(
    log,
    "freeze_baseline_provenance",
    f"Found {len(prov_files)} provenance files to bundle",
)

records: list[dict] = []
for p in prov_files:
    try:
        prov = json.loads(p.read_text())
    except json.JSONDecodeError as exc:
        log_transformation(
            log,
            "freeze_baseline_provenance",
            f"WARNING: could not parse {p.name}: {exc}",
            status="WARNING",
        )
        continue
    records.append(
        {
            "provenance_file": str(p.relative_to(provenance_dir.parent)),
            "provenance_sha256": file_sha256(p),
            "rule": prov.get("snakemake_rule") or prov.get("rule"),
            "artifact_id": prov.get("artifact_id"),
            "artifact_path": prov.get("path"),
            "artifact_sha256": prov.get("sha256"),
            "created_at": prov.get("created_at"),
            "ontology_operation": prov.get("ontology_operation"),
            "tool_versions": prov.get("tool_versions"),
        }
    )

bundle = {
    "baseline_version": version,
    "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
    "summary": summary,
    "git_commit": _git(["rev-parse", "HEAD"]) or None,
    "git_branch": _git(["rev-parse", "--abbrev-ref", "HEAD"]) or None,
    "git_remote_url": _git(["config", "--get", "remote.origin.url"]) or None,
    "n_provenance_files": len(records),
    "rules": records,
}

output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(bundle, indent=2, sort_keys=False))

log_transformation(
    log,
    "freeze_baseline_provenance",
    f"Wrote baseline bundle: {output_path} "
    f"({len(records)} rule records, git_commit={bundle['git_commit'][:12] if bundle['git_commit'] else 'none'})",
    status="SUCCESS",
    artifact_paths=[str(output_path)],
)
