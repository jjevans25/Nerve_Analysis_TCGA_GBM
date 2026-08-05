"""Test Oracle: do the nerve / tumor / immune compartments contain the cells they claim to?

Every downstream ligand-receptor conclusion in this project is a statement about
signalling *between compartments*. That claim is only meaningful if each mask
actually holds the cell population it is named after. This rule cross-tabulates
every compartment mask against the CELLxGENE Census author annotation
(`obs['cell_type']`, harmonized to the CL ontology), which is the only external
ground truth available for these cohorts.

It reads existing artifacts only — no re-run — so it can record a failing
baseline before any fix, and then gate the fix afterwards. See
markdowns/plan_compartment_integrity_fix.md.
"""

import json
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

ARM: str = snakemake.params.dataset
CLASS_MAP: dict[str, list[str]] = snakemake.params.census_class_map
IMMUNE_CLASSES: list[str] = list(snakemake.params.immune_classes)
GATES: dict[str, float] = snakemake.params.gates
PER_ARM: dict[str, object] = snakemake.params.per_arm or {}
ENFORCE: bool = bool(snakemake.params.enforce)

MALIGNANT_CLASS = "malignant"
UNMAPPED = "unmapped"


def _class_lookup() -> dict[str, str]:
    """Invert the configured class map to cell_type -> coarse class."""
    lookup: dict[str, str] = {}
    for klass, members in CLASS_MAP.items():
        for member in members:
            lookup[str(member)] = str(klass)
    return lookup


def load_obs(path: str, columns: list[str]) -> pd.DataFrame:
    """Read only the obs frame of an artifact — the audit never touches .X."""
    adata = ad.read_h5ad(path, backed="r")
    missing = [c for c in columns if c not in adata.obs.columns]
    if missing:
        raise KeyError(
            f"[FAIR-ALERT] {path} is missing required obs columns {missing}. "
            "The compartment audit cannot run without the Census author "
            "annotation; if `cell_type` is absent the artifact was written by a "
            "rule that overwrote it (defect D6)."
        )
    obs = adata.obs[columns].copy()
    del adata
    return obs


def classify(obs: pd.DataFrame, lookup: dict[str, str]) -> pd.DataFrame:
    """Attach the coarse Census class (neural / malignant / myeloid / …) to each cell."""
    obs = obs.copy()
    obs["census_cell_type"] = obs["cell_type"].astype(str)
    obs["census_class"] = obs["census_cell_type"].map(lookup).fillna(UNMAPPED)
    return obs


def composition(obs: pd.DataFrame, compartment: str) -> pd.DataFrame:
    """What is this compartment actually made of, by Census annotation?"""
    n_total = len(obs)
    rows: list[dict[str, object]] = []
    for level, column in (("class", "census_class"), ("cell_type", "census_cell_type")):
        counts = obs[column].value_counts()
        for key, n in counts.items():
            rows.append({
                "arm":                   ARM,
                "compartment":           compartment,
                "level":                 level,
                "key":                   str(key),
                "n_cells":               int(n),
                "share_of_compartment":  float(n) / n_total if n_total else 0.0,
            })
    rows.append({
        "arm":                  ARM,
        "compartment":          compartment,
        "level":                "total",
        "key":                  "n_cells",
        "n_cells":              int(n_total),
        "share_of_compartment": 1.0 if n_total else 0.0,
    })
    return pd.DataFrame(rows)


def class_fraction(obs: pd.DataFrame, klass: str) -> float:
    """Fraction of a compartment that the Census assigns to one coarse class."""
    if not len(obs):
        return 0.0
    return float((obs["census_class"] == klass).sum()) / len(obs)


def nerve_cluster_audit(nerve: pd.DataFrame) -> pd.DataFrame:
    """Per nerve_leiden cluster: is this cluster neural, or is it something else entirely?

    The S1PR1 finding turned on exactly this table — nerve_c24 is 92.9%
    endothelial, so a receptor called "nerve-side" was really vascular.
    """
    classes = sorted(set(CLASS_MAP.keys()) | {UNMAPPED})
    rows: list[dict[str, object]] = []
    for cluster, block in nerve.groupby("nerve_leiden", observed=True):
        n = len(block)
        top_type = block["census_cell_type"].value_counts()
        row: dict[str, object] = {
            "arm":                    ARM,
            "nerve_leiden":           str(cluster),
            "n_cells":                int(n),
            "dominant_census_type":   str(top_type.index[0]) if len(top_type) else "none",
            "dominant_share":         float(top_type.iloc[0]) / n if n and len(top_type) else 0.0,
            "endothelial_fraction":   float(
                (block["census_cell_type"] == "endothelial cell").sum()
            ) / n if n else 0.0,
        }
        for klass in classes:
            row[f"frac_{klass}"] = class_fraction(block, klass)
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.sort_values("n_cells", ascending=False).reset_index(drop=True)


def malignancy_confusion(malig_obs: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float]]:
    """Does the CNV `is_malignant` flag agree with the Census 'malignant cell' call?

    Precision answers "is the tumor compartment actually tumor"; recall answers
    "how much of the tumor did we find". Both feed hard gates on the CNV rebuild.
    """
    truth = (malig_obs["census_cell_type"] == "malignant cell").to_numpy()
    pred = malig_obs["is_malignant"].astype(bool).to_numpy()

    tp = int(np.sum(pred & truth))
    fp = int(np.sum(pred & ~truth))
    fn = int(np.sum(~pred & truth))
    tn = int(np.sum(~pred & ~truth))

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    rows: list[dict[str, object]] = [
        {"arm": ARM, "metric": "true_positive", "value": tp},
        {"arm": ARM, "metric": "false_positive", "value": fp},
        {"arm": ARM, "metric": "false_negative", "value": fn},
        {"arm": ARM, "metric": "true_negative", "value": tn},
        {"arm": ARM, "metric": "precision", "value": precision},
        {"arm": ARM, "metric": "recall", "value": recall},
        {"arm": ARM, "metric": "f1", "value": f1},
        {"arm": ARM, "metric": "n_predicted_malignant", "value": tp + fp},
        {"arm": ARM, "metric": "n_census_malignant", "value": tp + fn},
    ]

    # What the false positives actually are — this is how "the tumor compartment
    # is ~41% myeloid" becomes reproducible rather than prose.
    fp_types = malig_obs.loc[pred & ~truth, "census_cell_type"].value_counts()
    for cell_type, n in fp_types.items():
        rows.append({
            "arm":    ARM,
            "metric": f"false_positive_is::{cell_type}",
            "value":  int(n),
        })

    stats = {"precision": precision, "recall": recall, "f1": f1}
    return pd.DataFrame(rows), stats


def evaluate_gates(observed: dict[str, float]) -> pd.DataFrame:
    """Score every configured gate; the rule fails on any breach when enforcing."""
    specs: list[tuple[str, str, str, float | None]] = [
        ("nerve_neural_fraction",       "nerve_neural_fraction",       ">=", GATES.get("nerve_neural_fraction_min")),
        ("tumor_malignant_fraction",    "tumor_malignant_fraction",    ">=", GATES.get("tumor_malignant_fraction_min")),
        ("malignancy_recall",           "malignancy_recall",           ">=", GATES.get("malignancy_recall_min")),
        ("malignancy_precision",        "malignancy_precision",        ">=", GATES.get("malignancy_precision_min")),
        ("immune_purity",               "immune_purity",               ">=", GATES.get("immune_purity_min")),
        ("max_nerve_cluster_endothelial_fraction",
         "max_nerve_cluster_endothelial_fraction", "<=",
         GATES.get("max_nerve_cluster_endothelial_fraction")),
        ("immune_size_fraction_of_baseline", "immune_size_fraction_of_baseline", ">=",
         GATES.get("immune_size_fraction_of_baseline_min")),
        ("nerve_compartment_n_min",     "nerve_compartment_n",         ">=", PER_ARM.get("nerve_n_min")),
        ("nerve_compartment_n_max",     "nerve_compartment_n",         "<=", PER_ARM.get("nerve_n_max")),
        ("neuron_group_n_min",          "neuron_group_n",              ">=", PER_ARM.get("neuron_n_min")),
    ]

    rows: list[dict[str, object]] = []
    for name, observed_key, comparator, threshold in specs:
        value = observed.get(observed_key)
        if threshold is None or value is None:
            verdict = "NOT_SET"
        elif comparator == ">=":
            verdict = "PASS" if value >= threshold else "FAIL"
        else:
            verdict = "PASS" if value <= threshold else "FAIL"
        rows.append({
            "arm":        ARM,
            "gate":       name,
            "observed":   value,
            "comparator": comparator,
            "threshold":  threshold,
            "verdict":    verdict,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Load — obs only, never .X
# ---------------------------------------------------------------------------
log_transformation(log, "compartment_audit", f"Auditing arm '{ARM}' against Census cell_type")

lookup = _class_lookup()

malig_obs = classify(
    load_obs(snakemake.input.malig, ["cell_type", "cell_type_predicted", "is_malignant"]), lookup
)
nerve_obs = classify(load_obs(snakemake.input.nerve, ["cell_type", "nerve_leiden"]), lookup)
immune_obs = classify(load_obs(snakemake.input.immune, ["cell_type", "immune_subtype"]), lookup)

# The nerve compartment is split into glia and neurons (see nerve_cell_subset);
# audit them separately, because a 3.4k-cell neuron group can be entirely wrong
# without moving the combined figure at all.
_nerve_full = ad.read_h5ad(snakemake.input.nerve, backed="r")
if "nerve_subcompartment" in _nerve_full.obs.columns:
    nerve_obs["nerve_subcompartment"] = _nerve_full.obs["nerve_subcompartment"].astype(str).values
else:
    nerve_obs["nerve_subcompartment"] = "glia"
del _nerve_full

unmapped = sorted(set(malig_obs.loc[malig_obs["census_class"] == UNMAPPED, "census_cell_type"]))
if unmapped:
    log_transformation(log, "compartment_audit",
        f"[FAIR-ALERT] {len(unmapped)} Census cell_type values are absent from "
        f"compartment_audit.census_class_map and are counted as '{UNMAPPED}': {unmapped}",
        status="WARNING")

tumor_obs = malig_obs[malig_obs["is_malignant"].astype(bool)].copy()

log_transformation(log, "compartment_audit",
    f"Compartment sizes — cohort: {len(malig_obs)}, nerve: {len(nerve_obs)}, "
    f"tumor (is_malignant): {len(tumor_obs)}, immune: {len(immune_obs)}")

# ---------------------------------------------------------------------------
# Composition + confusion
# ---------------------------------------------------------------------------
_frames = [
    composition(malig_obs, "cohort"),
    composition(nerve_obs, "nerve"),
    composition(tumor_obs, "tumor"),
    composition(immune_obs, "immune"),
]
for _sub in sorted(nerve_obs["nerve_subcompartment"].unique()):
    _frames.append(
        composition(nerve_obs[nerve_obs["nerve_subcompartment"] == _sub], f"nerve_{_sub}")
    )
audit = pd.concat(_frames, ignore_index=True)

cluster_audit = nerve_cluster_audit(nerve_obs)
confusion, malig_stats = malignancy_confusion(malig_obs)

immune_purity = sum(class_fraction(immune_obs, k) for k in IMMUNE_CLASSES)
immune_baseline = PER_ARM.get("immune_baseline_n")

max_endo = float(cluster_audit["endothelial_fraction"].max()) if len(cluster_audit) else 0.0
immune_size_frac = (
    float(len(immune_obs)) / float(immune_baseline) if immune_baseline else None
)

observed: dict[str, float] = {
    "nerve_neural_fraction":                  class_fraction(nerve_obs, "neural"),
    "tumor_malignant_fraction":               class_fraction(tumor_obs, MALIGNANT_CLASS),
    "malignancy_recall":                      malig_stats["recall"],
    "malignancy_precision":                   malig_stats["precision"],
    "immune_purity":                          immune_purity,
    "max_nerve_cluster_endothelial_fraction": max_endo,
    "nerve_compartment_n":                    float(len(nerve_obs)),
    "immune_size_fraction_of_baseline":       immune_size_frac,
    "neuron_group_n": float(
        (nerve_obs["nerve_subcompartment"] == "neuron").sum()
    ),
}

gates = evaluate_gates(observed)

log_transformation(log, "compartment_audit",
    "Headline metrics:\n" + json.dumps(
        {k: (round(v, 4) if isinstance(v, float) else v) for k, v in observed.items()},
        indent=2,
    ))
log_transformation(log, "compartment_audit", "Gate results:\n" + gates.to_string(index=False))

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
_outputs = (
    (audit,         snakemake.output.audit,         "compartment_audit.csv"),
    (cluster_audit, snakemake.output.cluster_audit, "nerve_compartment_cluster_audit.csv"),
    (confusion,     snakemake.output.confusion,     "malignancy_confusion.csv"),
    (gates,         snakemake.output.gates,         "compartment_audit_gates.csv"),
)
for frame, path, _ in _outputs:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    verify_artifact(path, min_size_bytes=32)

# Durable sidecar copy at an UNDECLARED path.
#
# Snakemake deletes the declared outputs of a failed job. For an ordinary rule
# that is correct — a half-written artifact is worse than none. For this rule it
# is exactly backwards: when a gate fails, the cross-tabs explaining WHY are the
# whole point, and they were just deleted along with the failure. (Learned the
# hard way twice in this project; see the CHANGELOG note on undeclared sidecars.)
sidecar = Path(snakemake.params.sidecar_dir) / ARM
sidecar.mkdir(parents=True, exist_ok=True)
for frame, _, name in _outputs:
    frame.to_csv(sidecar / name, index=False)
log_transformation(log, "compartment_audit",
    f"Durable copies written to {sidecar} (undeclared — survives a failed gate)")

prov = stamp_artifact(
    output_path=snakemake.output.audit,
    rule_name="compartment_audit",
    input_paths=[snakemake.input.malig, snakemake.input.nerve, snakemake.input.immune],
    tool_versions={"anndata": ad.__version__, "pandas": pd.__version__, "numpy": np.__version__},
    parameters={
        "arm":               ARM,
        "enforce":           ENFORCE,
        "observed":          observed,
        "gates":             gates.to_dict(orient="records"),
        "unmapped_cell_types": unmapped,
        "oracle":            "CELLxGENE Census author annotation (obs['cell_type'], CL ontology)",
    },
    description="Compartment-vs-Census cross-tabulation with configurable integrity gates",
    ontology_operation="operation:2428",  # EDAM: Validation
)
write_provenance(prov, snakemake.output.provenance)

failed = gates[gates["verdict"] == "FAIL"]
if len(failed):
    message = (
        f"[FAIR-ALERT] {len(failed)} compartment-integrity gate(s) failed for arm "
        f"'{ARM}':\n{failed.to_string(index=False)}"
    )
    if ENFORCE:
        log_transformation(log, "compartment_audit", message, status="ERROR")
        raise RuntimeError(
            message + "\n\nThe compartment masks do not contain the cells they claim to. "
            "Downstream ligand-receptor results would be between-compartment claims about "
            "within-compartment signalling. See markdowns/plan_compartment_integrity_fix.md."
        )
    log_transformation(log, "compartment_audit",
        message + "\n(compartment_audit.enforce is false — recording baseline, not failing)",
        status="WARNING")

log_transformation(log, "compartment_audit", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.audit, snakemake.output.cluster_audit,
                                   snakemake.output.confusion, snakemake.output.gates])
