"""Record a content hash for every pinned v1.3.0 reference artifact.

Answers the biological-provenance question "is the frozen reference cohort I am
reading today byte-identical to the one the published results were derived
from?". The v1.3.0 nerve cascade cannot be re-run (``nerve_cells.h5ad`` was
deleted and is unreproducible — see CHANGELOG 2026-07-21), so the surviving
tables are the only record of that analysis and drift in them is silent.

Deliberately does NOT reuse ``provenance/baseline_v1.3.0.json``: that bundle is
dated 2026-05-25, several tables were legitimately regenerated on 2026-07-11/21,
and it stores ``artifact_sha256: null`` for ``annotate_cluster_qc``. Comparing
against it would report false drift.

Run on demand: ``snakemake --use-conda freeze_pinned_reference``.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, "workflow/scripts")
from fair_utils import file_sha256, log_transformation, verify_artifact  # noqa: E402

log = snakemake.log[0]  # type: ignore[name-defined]
version: str = snakemake.params.version  # type: ignore[name-defined]
artifacts: list[str] = list(snakemake.params.artifacts)  # type: ignore[name-defined]
out_path = Path(snakemake.output.manifest)  # type: ignore[name-defined]


def describe(rel_path: str) -> dict:
    """Hash one pinned artifact and capture the fields needed to detect drift."""
    path = Path(rel_path)
    stat = path.stat()
    return {
        "path": rel_path,
        "sha256": file_sha256(path),
        "size_bytes": stat.st_size,
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
    }


log_transformation(
    log, "freeze_pinned_reference",
    f"Hashing {len(artifacts)} pinned {version} artifacts",
)

missing = [a for a in artifacts if not Path(a).exists()]
if missing:
    # Hard-fail rather than emit a partial manifest: a manifest that silently
    # omits artifacts would later "verify" clean while they are still gone.
    log_transformation(
        log, "freeze_pinned_reference",
        f"[FAIR-ALERT] {len(missing)} configured artifacts do not exist: {missing}",
        status="FAILURE",
    )
    raise FileNotFoundError(
        f"{len(missing)} artifacts in baseline.pinned_artifacts are missing: {missing}"
    )

entries = [describe(a) for a in sorted(artifacts)]

manifest = {
    "baseline_version": version,
    "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
    "description": (
        "Content hashes of the pinned v1.3.0 reference artifacts. The rules that "
        "produce these are not defined while baseline.pinned is true, so this is "
        "the only drift detector for them."
    ),
    "n_artifacts": len(entries),
    "artifacts": entries,
}

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(json.dumps(manifest, indent=2) + "\n")
verify_artifact(out_path, min_size_bytes=100)

log_transformation(
    log, "freeze_pinned_reference",
    f"Wrote manifest for {len(entries)} artifacts to {out_path}",
    artifact_paths=[str(out_path)],
)
