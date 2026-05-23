"""FAIR compliance utilities: provenance logging, artifact registration, metadata stamping."""

from __future__ import annotations

import hashlib
import json
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


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
    output_path: Path | str,
) -> Path:
    """Write provenance record to the Snakemake-declared output path."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(provenance, f, indent=2)
    return output_path


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


@dataclass(frozen=True)
class PurityResult:
    """Per-cluster batch purity result + cohort-level constants needed for reporting."""
    df: pd.DataFrame
    n_samples: int
    expected_uniform_entropy: float


def compute_cluster_purity(
    adata,
    cluster_col: str,
    batch_col: str,
    *,
    dominant_max: float,
    min_contributing_fraction: float,
    min_contributing_samples: int,
) -> PurityResult:
    """Per-cluster batch purity (dominance + Shannon diversity).

    Biology question: did batch correction produce biologically-driven clusters,
    or do clusters reflect patient identity? Returns a verdict per cluster.

    Output df columns: ``cluster, n_cells, dominant_sample,
    dominant_sample_fraction, n_contributing_samples, shannon_entropy_bits,
    normalised_entropy, pass_dominant, pass_diversity, pass_overall``.
    Clusters are sorted by integer value of the label (matches scanpy Leiden output).
    """
    samples = sorted(adata.obs[batch_col].astype(str).unique().tolist())
    n_samples = len(samples)
    clusters = sorted(
        adata.obs[cluster_col].astype(str).unique().tolist(),
        key=lambda c: int(c),
    )

    counts = (
        adata.obs.assign(
            **{
                cluster_col: adata.obs[cluster_col].astype(str),
                batch_col: adata.obs[batch_col].astype(str),
            }
        )
        .groupby([cluster_col, batch_col], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(index=clusters, columns=samples, fill_value=0)
    )
    totals = counts.sum(axis=1)
    proportions = counts.div(totals.replace(0, np.nan), axis=0).fillna(0.0)

    dominant_fraction = proportions.max(axis=1)
    dominant_sample = proportions.idxmax(axis=1)
    n_contributing_samples = (proportions >= min_contributing_fraction).sum(axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        log_p = np.where(proportions > 0, np.log2(proportions), 0.0)
    shannon = -(proportions.values * log_p).sum(axis=1)
    expected_uniform = float(np.log2(n_samples)) if n_samples > 1 else 0.0
    normalised_entropy = (shannon / expected_uniform) if expected_uniform > 0 else np.zeros_like(shannon)

    pass_dominant = (dominant_fraction < dominant_max).values
    pass_diversity = (n_contributing_samples >= min_contributing_samples).values
    pass_overall = pass_dominant & pass_diversity

    df = pd.DataFrame(
        {
            "cluster": clusters,
            "n_cells": totals.values.astype(int),
            "dominant_sample": dominant_sample.values,
            "dominant_sample_fraction": dominant_fraction.values.round(4),
            "n_contributing_samples": n_contributing_samples.values.astype(int),
            "shannon_entropy_bits": np.round(shannon, 4),
            "normalised_entropy": np.round(normalised_entropy, 4),
            "pass_dominant": pass_dominant,
            "pass_diversity": pass_diversity,
            "pass_overall": pass_overall,
        }
    )
    return PurityResult(df=df, n_samples=n_samples, expected_uniform_entropy=expected_uniform)
