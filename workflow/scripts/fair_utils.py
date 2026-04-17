"""FAIR compliance utilities: provenance logging, artifact registration, metadata stamping."""

from __future__ import annotations

import hashlib
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def artifact_id(path: Path | str) -> str:
    """Return a stable UUID5 for a file path, seeded on its absolute path."""
    return str(uuid.uuid5(uuid.NAMESPACE_URL, str(Path(path).resolve())))


def file_sha256(path: Path | str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def stamp_artifact(
    output_path: Path | str,
    rule_name: str,
    input_paths: list[Path | str],
    tool_versions: dict[str, str],
    parameters: dict[str, Any],
    description: str,
    ontology_operation: str = "EDAM",
) -> dict:
    """Build a FAIR provenance record for a Snakemake output artifact."""
    output_path = Path(output_path)
    provenance = {
        "artifact_id": artifact_id(output_path),
        "path": str(output_path.resolve()),
        "sha256": file_sha256(output_path) if output_path.exists() else None,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "description": description,
        "snakemake_rule": rule_name,
        "ontology_operation": ontology_operation,
        "inputs": [
            {
                "path": str(Path(p).resolve()),
                "sha256": file_sha256(p) if Path(p).exists() else None,
                "artifact_id": artifact_id(p),
            }
            for p in input_paths
        ],
        "tool_versions": tool_versions,
        "parameters": parameters,
    }
    return provenance


def write_provenance(
    provenance: dict,
    provenance_dir: Path | str,
) -> Path:
    """Write provenance record to the provenance directory as JSON."""
    provenance_dir = Path(provenance_dir)
    provenance_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    rule = provenance.get("snakemake_rule", "unknown")
    out_file = provenance_dir / f"{rule}_{ts}_{provenance['artifact_id'][:8]}.json"
    with open(out_file, "w") as f:
        json.dump(provenance, f, indent=2)
    return out_file


def log_transformation(
    log_path: Path | str,
    rule_name: str,
    message: str,
    status: str = "SUCCESS",
    artifact_paths: list[str] | None = None,
) -> None:
    """Append a transformation record to the Snakemake log file."""
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "rule": rule_name,
        "status": status,
        "message": message,
        "artifacts": artifact_paths or [],
    }
    with open(log_path, "a") as f:
        f.write(json.dumps(entry) + "\n")


def verify_artifact(path: Path | str, min_size_bytes: int = 1) -> None:
    """Goal-backward check: raise if artifact is missing or below minimum size."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"[FAIR-ALERT] Expected artifact not found: {p}")
    if p.stat().st_size < min_size_bytes:
        raise ValueError(f"[FAIR-ALERT] Artifact exists but is suspiciously small ({p.stat().st_size} bytes): {p}")
