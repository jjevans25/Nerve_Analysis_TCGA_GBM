"""Cross-cohort concordance of the three-way nerve–tumor–immune LR result (Stage D).

Answers the project's replication question: do the tumor→immune→nerve
ligand–receptor interactions found in the reference TCGA-GBM cohort reproduce in
an independent GBM cohort? Reads both cohorts' ``nerve_tumor_immune_interactions
_with_qc.csv`` and computes, on the batch-QC-passing significant pairs:
  - Jaccard overlap of significant LR pairs (+ shared/unique pair table),
  - Spearman rank correlation of interaction magnitude on the shared pairs,
  - a comparison figure and a compact summary JSON.
Provenance-stamped like every other rule.
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
P = snakemake.params
dataset = P.dataset
pval_max = float(getattr(P, "pval_max", 0.05))

KEY_COLS = ["source_compartment", "target_compartment", "ligand_complex", "receptor_complex"]


def _load(path: str, label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    log_transformation(log, "cohort_concordance", f"{label}: {len(df)} rows from {path}")
    return df


def _significant(df: pd.DataFrame) -> pd.DataFrame:
    """Batch-QC-passing, p<=pval_max rows (defensive to missing columns)."""
    m = pd.Series(True, index=df.index)
    if "batch_qc_pass" in df.columns:
        m &= df["batch_qc_pass"].astype(str).str.lower().isin(["true", "1", "1.0"])
    if "cellphone_pvals" in df.columns:
        m &= pd.to_numeric(df["cellphone_pvals"], errors="coerce").fillna(1.0) <= pval_max
    return df[m].copy()


def _pair_key(df: pd.DataFrame) -> pd.Series:
    have = [c for c in KEY_COLS if c in df.columns]
    return df[have].astype(str).agg("|".join, axis=1)


def _magnitude(df: pd.DataFrame) -> pd.Series:
    """Higher = stronger interaction (for rank correlation)."""
    if "lrscore" in df.columns:
        return pd.to_numeric(df["lrscore"], errors="coerce")
    if "magnitude_rank" in df.columns:  # lower rank = stronger → negate
        return -pd.to_numeric(df["magnitude_rank"], errors="coerce")
    return pd.to_numeric(df.get("lr_means", pd.Series(np.nan, index=df.index)), errors="coerce")


ref = _load(snakemake.input.reference, "reference")
new = _load(snakemake.input.dataset, f"dataset[{dataset}]")

ref_sig, new_sig = _significant(ref), _significant(new)
ref_sig["_key"], new_sig["_key"] = _pair_key(ref_sig), _pair_key(new_sig)
ref_sig["_mag"], new_sig["_mag"] = _magnitude(ref_sig), _magnitude(new_sig)

# Collapse to best (max magnitude) row per pair key per cohort
ref_best = ref_sig.sort_values("_mag").drop_duplicates("_key", keep="last").set_index("_key")
new_best = new_sig.sort_values("_mag").drop_duplicates("_key", keep="last").set_index("_key")

ref_keys, new_keys = set(ref_best.index), set(new_best.index)
shared = sorted(ref_keys & new_keys)
union = ref_keys | new_keys
jaccard = (len(shared) / len(union)) if union else 0.0

# Spearman on shared-pair magnitudes
if len(shared) >= 3:
    rho, pval = spearmanr(ref_best.loc[shared, "_mag"].values,
                          new_best.loc[shared, "_mag"].values)
    rho, pval = float(rho), float(pval)
else:
    rho, pval = float("nan"), float("nan")

# Shared / unique pair table
rows = []
for k in sorted(union):
    rows.append({
        "pair_key": k,
        "in_reference": k in ref_keys,
        "in_dataset": k in new_keys,
        "reference_magnitude": float(ref_best.loc[k, "_mag"]) if k in ref_keys else np.nan,
        "dataset_magnitude": float(new_best.loc[k, "_mag"]) if k in new_keys else np.nan,
    })
pairs_df = pd.DataFrame(rows).sort_values(
    ["in_reference", "in_dataset"], ascending=False)
Path(snakemake.output.shared_pairs).parent.mkdir(parents=True, exist_ok=True)
pairs_df.to_csv(snakemake.output.shared_pairs, index=False)

# --- Figure: shared-pair magnitude scatter + overlap bar ----------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
if shared:
    axes[0].scatter(ref_best.loc[shared, "_mag"], new_best.loc[shared, "_mag"],
                    s=14, alpha=0.6, edgecolor="none")
    axes[0].set_xlabel("Reference cohort magnitude")
    axes[0].set_ylabel(f"{dataset} magnitude")
    axes[0].set_title(f"Shared significant LR pairs (n={len(shared)})\n"
                      f"Spearman ρ={rho:.3f} (p={pval:.2g})")
else:
    axes[0].text(0.5, 0.5, "No shared significant pairs", ha="center", va="center")
    axes[0].set_axis_off()

axes[1].bar(["reference-only", "shared", f"{dataset}-only"],
            [len(ref_keys - new_keys), len(shared), len(new_keys - ref_keys)],
            color=["#4C72B0", "#55A868", "#C44E52"])
axes[1].set_ylabel("significant LR pairs")
axes[1].set_title(f"Overlap (Jaccard={jaccard:.3f})")
fig.tight_layout()
fig.savefig(snakemake.output.figure, dpi=150, bbox_inches="tight")
plt.close(fig)

summary = {
    "dataset": dataset,
    "reference_table": str(snakemake.input.reference),
    "dataset_table": str(snakemake.input.dataset),
    "pval_max": pval_max,
    "pair_key_columns": [c for c in KEY_COLS if c in ref.columns],
    "n_reference_significant_pairs": len(ref_keys),
    "n_dataset_significant_pairs": len(new_keys),
    "n_shared_pairs": len(shared),
    "n_union_pairs": len(union),
    "jaccard_overlap": round(jaccard, 4),
    "spearman_rho_shared": None if np.isnan(rho) else round(rho, 4),
    "spearman_pvalue_shared": None if np.isnan(pval) else round(pval, 6),
}
Path(snakemake.output.summary).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.summary, "w") as f:
    json.dump(summary, f, indent=2)
verify_artifact(snakemake.output.summary, min_size_bytes=16)
log_transformation(log, "cohort_concordance",
    f"Jaccard={jaccard:.3f}, shared={len(shared)}, Spearman ρ={rho:.3f}")

prov = stamp_artifact(
    output_path=snakemake.output.summary,
    rule_name="cohort_concordance",
    input_paths=[snakemake.input.reference, snakemake.input.dataset],
    tool_versions={"pandas": pd.__version__, "numpy": np.__version__},
    parameters=summary,
    description="Cross-cohort concordance of three-way nerve-tumor-immune LR interactions",
    ontology_operation="operation:3501",  # EDAM: Enrichment / comparison
)
write_provenance(prov, snakemake.output.provenance)
log_transformation(log, "cohort_concordance", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.summary, snakemake.output.figure])
