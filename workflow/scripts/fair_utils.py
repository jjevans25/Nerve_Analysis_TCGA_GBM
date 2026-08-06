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


#: Compression applied to every pipeline .h5ad. AnnData defaults to none, which
#: left data/processed at 143 GB of uncompressed HDF5. Single-cell count matrices
#: are highly compressible (the raw Census pull is 7.1 GB gzipped for a 39.4 GB
#: payload), and read cost is negligible next to the disk saved.
H5AD_COMPRESSION = "gzip"


@dataclass(frozen=True)
class PurityResult:
    """Per-cluster batch purity result + cohort-level constants needed for reporting."""
    df: pd.DataFrame
    n_samples: int
    expected_uniform_entropy: float


def cap_cells_per_group(
    labels,
    max_per_group: int,
    seed: int,
):
    """Cap each interaction group to `max_per_group` cells; returns a boolean mask.

    Biology question this protects: LIANA scores a ligand-receptor pair from each
    group's MEAN expression and expression proportion, with a permutation null
    built by shuffling group labels. Neither quantity needs every cell — a group
    mean over 15,000 cells already has a negligible standard error.

    Why it is necessary here: `rank_aggregate` builds a zero-centred `scaled`
    layer, and zero-centring densifies. At 952,087 cells x 24,135 genes that is
    ~92 GB, which on a 36 GB machine drove the process into swap thrashing —
    measured at ~1 second of CPU per 20 seconds of wall clock, i.e. never
    finishing. Capping bounds that layer directly.

    Why it also improves the statistics: the compartments are wildly unbalanced
    (immune subtypes of 100k+ cells against nerve clusters of ~1k). A permutation
    null over such groups is dominated by the largest ones. Capping balances the
    design as a side effect.

    Sampling is per group, without replacement, from a seeded generator, so the
    selection is reproducible. Groups at or below the cap are kept whole.
    """
    import numpy as _np
    import pandas as _pd

    # Work positionally, not on the index: obs_names can repeat after a subset,
    # and a label-based lookup would then silently select the wrong cells.
    codes, _ = _pd.factorize(_pd.Series(labels).astype(str).to_numpy())
    rng = _np.random.default_rng(seed)
    keep = _np.zeros(codes.shape[0], dtype=bool)
    order = _np.argsort(codes, kind="stable")
    bounds = _np.flatnonzero(_np.diff(codes[order])) + 1
    for group_pos in _np.split(order, bounds):
        if group_pos.size > max_per_group:
            group_pos = rng.choice(group_pos, size=max_per_group, replace=False)
        keep[group_pos] = True
    return keep


NERVE_GROUP_PREFIX = "nerve_"


def nerve_group_key(label: str) -> str:
    """Cluster key behind a nerve LIANA group label: 'nerve_c24' -> '24', 'nerve_neuron' -> 'neuron'.

    The nerve compartment is no longer only per-Leiden-cluster groups. Neurons are
    carried as one pooled group (`nerve_neuron`) because ~4.3k of them across 170
    donors is a group, not a cluster set. Every consumer that used to assume the
    `nerve_c{N}` shape must go through this instead of slicing the prefix off by
    hand — that assumption crashed the three-way interaction rule on 2026-08-06.
    """
    s = str(label)
    if not s.startswith(NERVE_GROUP_PREFIX):
        return s
    rest = s[len(NERVE_GROUP_PREFIX):]
    # 'c24' -> '24', but leave a named group like 'neuron' alone.
    if rest.startswith("c") and rest[1:].isdigit():
        return rest[1:]
    return rest


def nerve_group_sort_key(label: str) -> tuple[int, int, str]:
    """Order nerve groups: numbered clusters ascending, then named groups alphabetically."""
    key = nerve_group_key(label)
    return (0, int(key), "") if key.isdigit() else (1, 0, key)


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

    # Sort integer-like cluster ids numerically (scanpy Leiden output); fall back
    # to lexical order for non-numeric labels (e.g. immune subtype names). Numeric
    # labels are unaffected — this preserves the original nerve/immune-Leiden order.
    def _cluster_sort_key(c: str) -> tuple[int, float, str]:
        return (0, int(c), "") if c.isdigit() else (1, 0.0, c)

    clusters = sorted(
        adata.obs[cluster_col].astype(str).unique().tolist(),
        key=_cluster_sort_key,
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
