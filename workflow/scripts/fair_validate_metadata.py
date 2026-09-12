"""Validate FAIR provenance: every artifact declared, and produced by the pinned tools.

Biological question addressed: none directly — this is the audit gate that decides
whether any of the cohort's biological claims are reusable. A result whose producing
environment is unknown cannot be reproduced, so this rule is what converts the
project's "Reusable" claim from prose into something that fails a build.

Three checks, in increasing strength:

1. **Coverage.** Count every provenance record, recursing into the per-dataset
   namespaces. The previous version globbed ``provenance/*.json`` only and therefore
   never saw the ~650 records under ``provenance/<dataset>/``, i.e. it reported on the
   retired reference cohort while the live Census arms went unaudited.

2. **Declared-pin agreement (Defect 1).** Each record's ``tool_versions`` is compared
   against the ``pkg==version`` pins in ``workflow/envs/*.yaml``. This is the standing
   check for the conda-env shadowing defect: before 2026-08-05 the declared envs were
   not the envs that executed, and the symptom visible from outside was exactly this —
   a recorded version that no declared pin allows. A conflict fails the build unless it
   is in ``ACCEPTED_VERSION_EXCEPTIONS`` below, so a *new* drift can never be silent.

3. **Enforcement-era accounting.** Records created before the conda-enforcement fix are
   counted and listed, so it stays visible which artifacts predate it without having to
   re-derive that fact by hand each time.

Malformed ``tool_versions`` entries (a module path where a version belongs) are reported
as warnings rather than failures: they make the record less reusable, but they do not
imply the wrong environment ran.
"""

# NOTE: no `from __future__ import annotations` here. Snakemake prepends its own
# preamble to every `script:` file, so a __future__ import is no longer the first
# statement and raises SyntaxError at job start. Python 3.12 evaluates `X | None`
# and builtin generics natively, so it buys nothing here anyway.
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation  # noqa: E402

# Commit a375811, 2026-08-05 19:46:05 -0400 — `scripts/run_snakemake.sh`, the point
# from which `--use-conda` actually enforces workflow/envs/*.yaml for `script:` rules.
# See markdowns/task_conda_env_enforcement.md.
CONDA_ENFORCEMENT_FIX = datetime(2026, 8, 5, 23, 46, 5, tzinfo=timezone.utc)

# A recorded version that no declared pin allows, which is nevertheless ACCEPTED.
# Each entry is a decision with a reason, not a suppression: keyed by
# (provenance path suffix, package) so it cannot accidentally widen to other records.
ACCEPTED_VERSION_EXCEPTIONS: dict[tuple[str, str], str] = {
    ("gbm_cellxgene_56c4912d_full/scvi_integration_provenance.json", "torch"): (
        "torch 2.11.0 (venv) vs pinned 2.12.0. Produced 2026-07-30, before the "
        "conda-enforcement fix. This is the ONLY load-bearing artifact in either "
        "Census arm that predates it — every rule downstream (annotate, malignancy, "
        "both subsets, both interaction rules, compartment audit, scANVI, "
        "concordance) was rebuilt 2026-08-06 under the enforced env. Re-deriving it "
        "means an ~11 h scVI retrain that renumbers Leiden clusters and cascades "
        "through every table, so it is a researcher decision, not a scheduling one. "
        "Tracked in markdowns/task_conda_env_enforcement.md."
    ),
    ("nerve_clinical_association_provenance.json", "scanpy"): (
        "Records scanpy=0.12.10, which is anndata's version, not scanpy's: the script "
        "carried `\"scanpy\": ad.__version__` as a deliberate \"stand-in for the env\". "
        "The source bug is fixed (the key is gone), but this artifact belongs to the "
        "v1.3.0 reference cohort, which is pinned and unreproducible — "
        "data/processed/nerve_cells.h5ad was deleted — so the record cannot be "
        "regenerated. It is accepted as a known-bad historical record, NOT as evidence "
        "the right environment ran. See [[baseline-reference-pinned-do-not-regenerate]]."
    ),
    ("gene_positions_provenance.json", "requests"): (
        "requests 2.33.1 vs pinned 2.32.3. `download_gene_positions` fetches a static "
        "genomic coordinate table; the artifact is verified by content, not by the "
        "HTTP client that retrieved it."
    ),
}

# Historical FAIR freezes. CLAUDE.md forbids rewriting these, so they are also not
# held to current pins — they record what was true at freeze time, by design.
FROZEN_RECORD_PATTERN = re.compile(r"baseline_v\d+\.\d+\.\d+\.json$")

# Per-sample rules: 170 of each per arm. Excluded from the pre-fix listing only, so the
# listing stays readable; they are still counted and still pin-checked.
PER_SAMPLE_PATTERN = re.compile(r"_(ingest|qc)_provenance\.json$")

VERSION_RE = re.compile(r"^\d+(\.\d+)*([a-zA-Z0-9.\-+]*)$")


def declared_pins(env_dir: Path) -> dict[str, set[str]]:
    """Parse `pkg==version` pins from every environment spec the workflow declares."""
    pins: dict[str, set[str]] = {}
    for spec in sorted(env_dir.glob("*.yaml")):
        for line in spec.read_text().splitlines():
            match = re.match(r"\s*-\s*([A-Za-z0-9_.\-]+)==([0-9][^\s#]*)", line)
            if match:
                name = match.group(1).lower().replace("-", "_")
                pins.setdefault(name, set()).add(match.group(2))
    return pins


def created_at(record: dict[str, object]) -> datetime | None:
    """Read a provenance record's UTC creation time, tolerating a naive timestamp."""
    raw = record.get("created_at") or record.get("timestamp")
    if not isinstance(raw, str):
        return None
    try:
        stamp = datetime.fromisoformat(raw.replace("Z", "+00:00"))
    except ValueError:
        return None
    return stamp if stamp.tzinfo else stamp.replace(tzinfo=timezone.utc)


def main() -> None:
    prov_dir = Path(snakemake.input.provenance_dir)  # noqa: F821
    env_dir = Path("workflow/envs")
    pins = declared_pins(env_dir)

    records = sorted(prov_dir.glob("**/*.json"))
    conflicts: list[dict[str, object]] = []
    accepted: list[dict[str, object]] = []
    malformed: list[dict[str, object]] = []
    unreadable: list[str] = []
    pre_fix: list[dict[str, str]] = []
    n_with_versions = 0

    for path in records:
        rel = str(path.relative_to(prov_dir))
        try:
            record = json.loads(path.read_text())
        except (json.JSONDecodeError, OSError) as exc:
            unreadable.append(f"{rel}: {exc}")
            continue
        if not isinstance(record, dict):
            unreadable.append(f"{rel}: top level is not an object")
            continue

        stamp = created_at(record)
        if (stamp and stamp < CONDA_ENFORCEMENT_FIX
                and not PER_SAMPLE_PATTERN.search(rel)
                and not FROZEN_RECORD_PATTERN.search(rel)):
            pre_fix.append({"record": rel, "created_at": stamp.isoformat()})

        if FROZEN_RECORD_PATTERN.search(rel):
            continue

        versions = record.get("tool_versions")
        if not isinstance(versions, dict) or not versions:
            continue
        n_with_versions += 1

        for tool, observed in versions.items():
            name = str(tool).lower().replace("-", "_")
            observed = str(observed)
            if name not in pins:
                continue
            if not VERSION_RE.match(observed):
                malformed.append({"record": rel, "tool": tool, "observed": observed})
                continue
            if observed in pins[name]:
                continue
            finding: dict[str, object] = {
                "record": rel,
                "tool": tool,
                "observed": observed,
                "declared": sorted(pins[name]),
            }
            reason = ACCEPTED_VERSION_EXCEPTIONS.get((rel, str(tool)))
            if reason:
                accepted.append({**finding, "accepted_because": reason})
            else:
                conflicts.append(finding)

    report = {
        "validated_at": datetime.now(timezone.utc).isoformat(),
        "pass": not conflicts and not unreadable,
        "total_provenance_records": len(records),
        "records_with_tool_versions": n_with_versions,
        "declared_pins": len(pins),
        "conda_enforcement_fix": CONDA_ENFORCEMENT_FIX.isoformat(),
        "n_records_predating_enforcement_fix": len(pre_fix),
        "records_predating_enforcement_fix": sorted(
            pre_fix, key=lambda r: r["created_at"]),
        "n_version_conflicts": len(conflicts),
        "version_conflicts": conflicts,
        "n_accepted_exceptions": len(accepted),
        "accepted_exceptions": accepted,
        "n_malformed_versions": len(malformed),
        "malformed_versions": malformed,
        "unreadable_records": unreadable,
    }

    out = Path(snakemake.output.validation_report)  # noqa: F821
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))

    summary = (
        f"{len(records)} provenance records "
        f"({n_with_versions} carrying tool_versions) against {len(pins)} declared pins: "
        f"{len(conflicts)} conflict(s), {len(accepted)} accepted exception(s), "
        f"{len(malformed)} malformed version string(s), "
        f"{len(pre_fix)} record(s) predating the conda-enforcement fix"
    )
    log_transformation(
        log_path=snakemake.log[0],  # noqa: F821
        rule_name="fair_validate_metadata",
        message=summary,
        artifact_paths=[str(out)],
    )

    if unreadable:
        raise RuntimeError(
            "[FAIR-ALERT] unreadable provenance record(s): " + "; ".join(unreadable))
    if conflicts:
        detail = "; ".join(
            f"{c['record']}: {c['tool']}={c['observed']} not in {c['declared']}"
            for c in conflicts)
        raise RuntimeError(
            "[FAIR-ALERT] artifact(s) recorded a tool version no declared environment "
            "pin allows, i.e. the environment that ran was not the environment "
            f"declared: {detail}. Either re-derive under the enforced env (always via "
            "scripts/run_snakemake.sh) or, if the drift is knowingly accepted, add it "
            "to ACCEPTED_VERSION_EXCEPTIONS with the reason."
        )


main()
