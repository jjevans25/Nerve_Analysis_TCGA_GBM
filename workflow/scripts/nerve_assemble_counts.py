"""Reassemble raw-count nerve-cell AnnData for scANVI re-training.

Question answered: scvi-tools requires raw integer counts in ``.X``, but
``data/processed/nerve_cells.h5ad`` ships log1p-transformed floats and has no
``layers['counts']`` slot. This script reconstructs the missing counts by
expm1+rounding the per-sample QC files' log1p .X (verified to be log1p of raw
integer counts, no library-size normalization), subsetting to the nerve obs,
and reindexing var to match the v1.0.0 nerve gene set. Pure data marshalling
— no biology re-derived. v1.0.0 artifacts are not touched.
"""

import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "workflow/scripts")
from counts_utils import recover_counts_from_log1p as _recover_counts
from fair_utils import (
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]


def _crash_hook(exc_type, exc_value, exc_tb):
    """Persist any uncaught exception into snakemake.log[0] before exiting."""
    import traceback
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as f:
        f.write("=== UNCAUGHT EXCEPTION ===\n")
        traceback.print_exception(exc_type, exc_value, exc_tb, file=f)
        f.write("===========================\n")
    sys.__excepthook__(exc_type, exc_value, exc_tb)


sys.excepthook = _crash_hook

os.environ["PYTHONHASHSEED"] = str(int(snakemake.params.random_seed))

per_sample_h5ads: list[str] = list(snakemake.input.qc_h5ads)
nerve_ref_path: str = str(snakemake.input.nerve_ref)
out_path = Path(snakemake.output.h5ad)

log_transformation(
    log, "nerve_assemble_counts",
    f"Loading nerve obs/var reference from {nerve_ref_path}",
)
nerve_ref = ad.read_h5ad(nerve_ref_path, backed="r")
nerve_obs_names = set(nerve_ref.obs_names)
nerve_var_index = nerve_ref.var.index.copy()
# v1.0.0 nerve_cells.h5ad obs["total_counts"] is the sum of log1p(.X), not raw
# counts (calculate_qc_metrics ran on log1p data). Use Seurat's nCount_SCT
# which holds the per-cell sum of the SCT-corrected integer counts that
# log1p was applied to — this is the upper bound we recover.
nerve_obs_nCount_SCT = nerve_ref.obs["nCount_SCT"].copy()
n_obs_expected = nerve_ref.n_obs
n_vars_expected = nerve_ref.n_vars
log_transformation(
    log, "nerve_assemble_counts",
    f"Reference: {n_obs_expected} cells, {n_vars_expected} vars to align to",
)


# Count recovery now lives in counts_utils.recover_counts_from_log1p (imported
# above as _recover_counts) so scrna_integration.py and this script share one
# implementation. Body is unchanged from the previous local definition.

pieces: list[ad.AnnData] = []
n_cells_collected = 0
for qc_path in per_sample_h5ads:
    a = ad.read_h5ad(qc_path)
    sample_id = Path(qc_path).stem.replace("_qc", "")
    keep_mask = np.array([n in nerve_obs_names for n in a.obs_names])
    n_keep = int(keep_mask.sum())
    if n_keep == 0:
        log_transformation(log, "nerve_assemble_counts",
                           f"  {sample_id}: 0 nerve cells (skipping)", status="WARNING")
        continue
    a = a[keep_mask].copy()
    # Recover counts from log1p .X.
    a.X = _recover_counts(a.X)
    # Reindex var to the v1.0.0 nerve gene set (all 20,420 expected to be present).
    missing = [g for g in nerve_var_index if g not in a.var.index]
    if missing:
        raise RuntimeError(
            f"[FAIR-ALERT] {sample_id}: missing {len(missing)} nerve vars "
            f"(first 5: {missing[:5]})"
        )
    a = a[:, nerve_var_index].copy()
    pieces.append(a)
    n_cells_collected += n_keep
    log_transformation(log, "nerve_assemble_counts",
                       f"  {sample_id}: kept {n_keep} nerve cells")

log_transformation(
    log, "nerve_assemble_counts",
    f"Concatenating {len(pieces)} sample pieces -> {n_cells_collected} cells",
)
adata = ad.concat(pieces, axis=0, join="outer", index_unique=None, merge="same")
adata.var = nerve_ref.var.copy()  # carry forward the v1.0.0 var schema verbatim

# Re-order obs to match the v1.0.0 nerve_cells.h5ad ordering so downstream
# rules can rely on positional cell correspondence.
adata = adata[nerve_ref.obs_names].copy()

# ----------------------------------------------------------------------------
# Verification gates
# ----------------------------------------------------------------------------
if adata.n_obs != n_obs_expected:
    raise RuntimeError(
        f"[FAIR-ALERT] n_obs mismatch: got {adata.n_obs}, expected {n_obs_expected}"
    )
if not (adata.var.index == nerve_var_index).all():
    raise RuntimeError("[FAIR-ALERT] var.index does not match v1.0.0 nerve var set")

X = adata.X
data_arr = X.data if sp.issparse(X) else np.asarray(X).ravel()
nonzero = data_arr[data_arr != 0]
if not np.all(nonzero == nonzero.astype(np.int32)):
    raise RuntimeError("[FAIR-ALERT] recovered .X is not integer-valued")
if (nonzero < 0).any():
    raise RuntimeError("[FAIR-ALERT] recovered .X has negative entries")

# Compare per-cell sum of recovered SCT counts (over the 20,420 nerve var
# subset) to Seurat's nCount_SCT (the per-cell sum over the FULL gene set
# pre-subset). Reduced-set total must be <= full-set total (we dropped genes,
# not added). Allow 1.01x rounding slack — recovered values come from
# expm1+round, so a per-cell discrepancy of a few counts is expected.
hvg_total = np.asarray(X.sum(axis=1)).ravel()
slack = 1.01
violation = hvg_total > (nerve_obs_nCount_SCT.values * slack)
if violation.any():
    n_bad = int(violation.sum())
    examples = np.where(violation)[0][:3]
    ex = [(int(i), float(hvg_total[i]), float(nerve_obs_nCount_SCT.values[i]))
          for i in examples]
    raise RuntimeError(
        f"[FAIR-ALERT] {n_bad} cells have reduced-set count > 1.01x nCount_SCT — "
        f"gene-subset mismatch or recovery error. Examples (idx, hvg_total, nCount_SCT): {ex}"
    )
recovery_ratio = float(np.median(hvg_total / np.maximum(nerve_obs_nCount_SCT.values, 1)))
if recovery_ratio < 0.50:
    raise RuntimeError(
        f"[FAIR-ALERT] median recovered-counts / nCount_SCT = {recovery_ratio:.3f} "
        f"— too much signal lost in the var reindex; investigate gene set."
    )
log_transformation(
    log, "nerve_assemble_counts",
    f"Verification: n_obs={adata.n_obs}; integer counts confirmed; "
    f"median recovered-count / nCount_SCT = {recovery_ratio:.3f} (expected ~0.85-1.0)",
)

# Stamp the new obs columns we need downstream.
adata.obs["sample_id"] = adata.obs["sample_id"].astype(str)
adata.obs["batch"] = adata.obs["batch"].astype(str)
adata.obs["hvg_total_counts"] = hvg_total.astype(np.int32)

out_path.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out_path)
verify_artifact(out_path, min_size_bytes=1_000_000)

prov = stamp_artifact(
    output_path=out_path,
    rule_name="nerve_assemble_counts",
    input_paths=[nerve_ref_path, *per_sample_h5ads],
    tool_versions={
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "scipy": sp.__name__ + "@" + getattr(sp, "__version__", "n/a"),
    },
    parameters={
        "n_cells": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "n_samples": len(pieces),
        "recovery_method": "expm1_round_int32",
    },
    description=(
        "Raw-count nerve-cell AnnData reconstructed from per-sample log1p QC "
        "files (no library-size normalization detected upstream); reindexed "
        "to v1.0.0 nerve gene set. Input for scANVI re-training."
    ),
    ontology_operation="operation:3431",  # EDAM: Deposition (data assembly)
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log, "nerve_assemble_counts",
    "Complete",
    status="SUCCESS",
    artifact_paths=[str(out_path), str(snakemake.output.provenance)],
)
