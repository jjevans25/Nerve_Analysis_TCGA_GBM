"""Cross-arm reproducible nerve-immune signalling axes, tiered by druggability.

Answers: which ligand-receptor axes between the nerve compartment and the immune
compartment reproduce in *both* census arms after the compartment-integrity fix,
and which of them are pharmacologically reachable? Cross-arm agreement replaces
the reference-cohort comparison, which is unusable post-fix because the reference
nerve compartment was built by the removed logic.

Reimplements the out-of-band generator exported to
``claude_science/code_export/{01_build_corrected_axis_tables,
02_chembl_druggability_pipeline,03_lead_targets_report_postfix}.py``, which had no
producing rule and whose output could neither be rebuilt nor restored from git
(``results/`` is gitignored). Verified byte-identical against the 2026-08-07
artifact it replaces.

The table is deliberately a concatenation of two independently-built halves, and
they do NOT share a selection rule:

  oligodendrocyte half (128 rows)   neuron half (55 rows)
  ------------------------------    ---------------------------------
  rows from the FULL arm            rows from the CAPPED arm
  batch_qc_pass required            NO QC filter (pooled neurons fail
                                      donor QC in the full arm)
  all nerve clusters                nerve_cluster == 'nerve_neuron' only
  sorted ['best_mag','min_pval']    sorted ['best_mag'] only
  rank_full  dense 1..128           rank_capped dense 1..55
  capped_best_mag populated         capped_best_mag never created -> NaN

Consequence, preserved here on purpose: the column named ``best_mag`` means
"full arm, QC-passing" for 128 rows and "capped arm, no QC, neuron-only" for 55,
and ``n_rows`` is arm-inconsistent the same way. Notebook 05 Panel E reports this
as the 128/128-vs-128/167 discrepancy. Changing it is a separate, reviewable
change -- this script's contract is exact reproduction.

The four drug-annotation columns (tier_v2, agents_flagged, glioma_trials,
withdrawn) cannot be recomputed: they came from ChEMBL and ClinicalTrials.gov
calls made through a Claude Science MCP connector whose response caches did not
survive, and one input (a prior-session artifact) is permanently lost. They are
supplied instead by a committed, version-pinned snapshot joined on ``axis``,
which is lossless because all four are a pure function of the axis string
(verified: zero disagreement across the 39 axes appearing on both nerve sides).
See ``reference/drug_annotation/MANIFEST.json`` for the provenance and limits of
that snapshot, and ``rule refresh_drug_annotation`` to mint a fresh one.
"""

# NOTE: no `from __future__ import annotations` here. Snakemake's `script:`
# directive prepends its own preamble to this file, so a __future__ import is no
# longer the first statement and raises SyntaxError at runtime. Same constraint
# as download_gene_positions.py.

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    log_transformation,
    nerve_group_key,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

RULE = "nerve_immune_lead_axes"

log = snakemake.log[0]  # type: ignore[name-defined]
P = snakemake.params  # type: ignore[name-defined]

PVAL_MAX = float(P.pval_max)
MAG_MAX = float(P.magnitude_rank_max)
NEURON_GROUP = str(P.neuron_group)
SIDE_OLIGO = str(P.compartment_side_oligo)
SIDE_NEURON = str(P.compartment_side_neuron)

# The generator's boolean coercion: these columns arrive as a mix of real bools
# and the strings pandas wrote on a previous round-trip.
BOOL_MAP = {True: True, "True": True, False: False, "False": False}

DRUG_COLS = ["tier_v2", "agents_flagged", "glioma_trials", "withdrawn"]

# Column order of the emitted file. `qc_pass` is listed by the generator but
# exists in neither half, so the `if c in df.columns` filter below drops it and
# the file carries 15 columns, not 16. Kept here so the lineage stays legible.
COLS = ["rank_full", "rank_capped", "axis", "tier_v2", "interfaces", "nerve_side",
        "immune", "best_mag", "capped_best_mag", "min_pval", "n_rows", "qc_pass",
        "agents_flagged", "glioma_trials", "withdrawn"]


def _prep(df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    """Attach nerve cell-type labels and collapse reciprocal LR rows to undirected axes."""
    df = df.copy()
    for col in ("batch_qc_pass", "nerve_batch_qc_pass", "immune_batch_qc_pass"):
        if col in df.columns:
            df[col] = df[col].map(BOOL_MAP)
    df["batch_qc_pass"] = df["batch_qc_pass"].fillna(True).astype(bool)
    df["nerve_cluster"] = df["nerve_cluster"].fillna("").astype(str)
    df["immune_subtype"] = df["immune_subtype"].fillna("").astype(str)

    # nerve_group_key, not a prefix slice: the compartment carries both
    # per-Leiden groups (nerve_c24) and one pooled neuron group (nerve_neuron).
    df["cluster_key"] = df["nerve_cluster"].map(nerve_group_key).astype(str).str.strip()

    ann = ann.copy()
    ann["cluster_key"] = ann["cluster"].astype(str)
    ann["nerve_cell_type"] = ann["label"].astype(str).str.split("|", expand=True)[1].str.strip()
    df = df.merge(ann[["cluster_key", "nerve_cell_type"]], on="cluster_key", how="left")
    df.loc[df["cluster_key"] == "neuron", "nerve_cell_type"] = "neuron(pooled)"
    df["nerve_cell_type"] = df["nerve_cell_type"].fillna("")

    # Undirected axis: lexicographic sort makes APP|CD74 the same axis whichever
    # side was the ligand (measured 109 ligand-first / 70 receptor-first).
    df["axis"] = ["|".join(sorted([str(lig), str(rec)]))
                  for lig, rec in zip(df.ligand_complex, df.receptor_complex)]
    return df


def _axes_table(df: pd.DataFrame, qc_only: bool = True) -> pd.DataFrame:
    """Best-of aggregation per undirected axis across the three compartment interfaces."""
    sig = df[(df.cellphone_pvals <= PVAL_MAX) & (df.magnitude_rank <= MAG_MAX)]
    if qc_only:
        sig = sig[sig.batch_qc_pass]
    grouped = sig.groupby("axis").agg(
        best_mag=("magnitude_rank", "min"),
        best_spec=("specificity_rank", "min"),
        best_lrscore=("lrscore", "max"),
        min_pval=("cellphone_pvals", "min"),
        n_rows=("axis", "size"),
        interfaces=("compartment_pair", lambda x: ",".join(sorted(set(x)))),
        directions=("direction", lambda x: ",".join(sorted(set(x)))),
        nerve_types=("nerve_cell_type", lambda x: ",".join(sorted({v for v in x if v}))),
        immune=("immune_subtype", lambda x: ",".join(sorted({v for v in x if v}))),
    ).reset_index()
    # Sort key is load-bearing: ranks below are positional indices into this order.
    return grouped.sort_values(["best_mag", "min_pval"]).reset_index(drop=True)


def _pooled_neuron_axes(df: pd.DataFrame) -> pd.DataFrame:
    """Same aggregation restricted to the pooled neuron group, with NO batch-QC filter.

    The pooled neuron group passes donor QC in the capped arm (40% dominant donor,
    10 donors) but fails in the full arm (50.3%, 7 donors). Filtering on QC here
    would empty the neuron side rather than caveat it.
    """
    sig = df[(df.nerve_cluster == NEURON_GROUP)
             & (df.cellphone_pvals <= PVAL_MAX)
             & (df.magnitude_rank <= MAG_MAX)]
    grouped = sig.groupby("axis").agg(
        best_mag=("magnitude_rank", "min"),
        min_pval=("cellphone_pvals", "min"),
        n_rows=("axis", "size"),
        interfaces=("compartment_pair", lambda x: ",".join(sorted(set(x)))),
        immune=("immune_subtype", lambda x: ",".join(sorted({v for v in x if v}))),
    ).reset_index()
    # Sorted on best_mag ONLY -- not ['best_mag','min_pval'] like the oligo half.
    # pandas' stable sort means adding min_pval here would shift neuron ranks.
    return grouped.sort_values("best_mag")


def _join_drug_annotation(df: pd.DataFrame, drug: pd.DataFrame) -> pd.DataFrame:
    """Attach the pinned ChEMBL/trials snapshot; every axis must be covered."""
    lookup = drug.set_index("axis")
    missing = sorted(set(df.axis) - set(lookup.index))
    if missing:
        raise ValueError(
            f"[FAIR-ALERT] {len(missing)} axes have no row in the pinned drug "
            f"annotation snapshot and would silently emit empty tiers: "
            f"{missing[:10]}{'...' if len(missing) > 10 else ''}. The upstream "
            f"interaction tables have changed relative to the snapshot -- refresh it "
            f"via `rule refresh_drug_annotation` rather than shipping blank columns."
        )
    for col in DRUG_COLS:
        df[col] = df.axis.map(lookup[col])
    return df


def _select(df: pd.DataFrame, side: str, sort_rank: str) -> pd.DataFrame:
    """Project to the emitted column order and sort within the compartment side."""
    df = df.copy()
    df["compartment_side"] = side
    keep = [c for c in COLS if c in df.columns] + ["compartment_side"]
    return df[keep].sort_values(["tier_v2", sort_rank])


# --- Load both census arms ----------------------------------------------------
prepped, ax, pn = {}, {}, {}
for arm_key, lr_path, ann_path in (
    ("full", snakemake.input.full_lr, snakemake.input.full_ann),  # type: ignore[name-defined]
    ("capped", snakemake.input.capped_lr, snakemake.input.capped_ann),  # type: ignore[name-defined]
):
    lr = pd.read_csv(lr_path, low_memory=False)
    ann = pd.read_csv(ann_path)
    prepped[arm_key] = _prep(lr, ann)
    ax[arm_key] = _axes_table(prepped[arm_key])
    pn[arm_key] = _pooled_neuron_axes(prepped[arm_key])
    log_transformation(log, RULE,
                       f"{arm_key} arm: {len(lr)} LR rows -> {len(ax[arm_key])} QC-passing axes, "
                       f"{len(pn[arm_key])} pooled-neuron axes")

# --- Oligodendrocyte half: cross-arm shared, rows from the FULL arm ------------
shared = set(ax["full"].axis) & set(ax["capped"].axis)
capped_by_axis = ax["capped"].set_index("axis")
lead = ax["full"][ax["full"].axis.isin(shared)].copy().reset_index(drop=True)
lead["rank_full"] = lead.index + 1
lead["rank_capped"] = lead.axis.map({a: i + 1 for i, a in enumerate(ax["capped"].axis)})
lead["capped_best_mag"] = lead.axis.map(capped_by_axis.best_mag)
lead["nerve_side"] = lead.nerve_types.replace("", "—")

# --- Neuron half: cross-arm shared, rows from the CAPPED arm -------------------
neuron_shared = set(pn["full"].axis) & set(pn["capped"].axis)
neu = pn["capped"][pn["capped"].axis.isin(neuron_shared)].copy().reset_index(drop=True)
neu["rank_capped"] = neu.index + 1
neu["rank_full"] = neu.axis.map({a: i + 1 for i, a in enumerate(pn["full"].axis)})
neu["nerve_side"] = "neuron(pooled)"

# --- Pinned drug annotation ---------------------------------------------------
drug = pd.read_csv(snakemake.input.drug)  # type: ignore[name-defined]
lead = _join_drug_annotation(lead, drug)
neu = _join_drug_annotation(neu, drug)

# --- Emit ---------------------------------------------------------------------
leads = pd.concat(
    [_select(lead, SIDE_OLIGO, "rank_full"), _select(neu, SIDE_NEURON, "rank_capped")],
    ignore_index=True,
)
Path(snakemake.output.table).parent.mkdir(parents=True, exist_ok=True)  # type: ignore[name-defined]
leads.to_csv(snakemake.output.table, index=False)  # type: ignore[name-defined]
verify_artifact(snakemake.output.table, min_size_bytes=1024)  # type: ignore[name-defined]

tier_counts = leads.tier_v2.value_counts().to_dict()
summary = {
    "n_axes": int(len(leads)),
    "n_oligodendrocyte_side": int((leads.compartment_side == SIDE_OLIGO).sum()),
    "n_neuron_side": int((leads.compartment_side == SIDE_NEURON).sum()),
    "n_distinct_axes": int(leads.axis.nunique()),
    "pval_max": PVAL_MAX,
    "magnitude_rank_max": MAG_MAX,
    "neuron_group": NEURON_GROUP,
    "tier_v2_counts": {str(k): int(v) for k, v in tier_counts.items()},
    "drug_annotation_snapshot": str(snakemake.input.drug),  # type: ignore[name-defined]
    "best_mag_is_arm_inconsistent": True,
    "best_mag_note": (
        "best_mag/n_rows are full-arm QC-passing for the oligodendrocyte side and "
        "capped-arm non-QC pooled-neuron for the neuron side; capped_best_mag is "
        "populated on the oligodendrocyte side only. Reproduced from the 2026-08-07 "
        "generator on purpose."
    ),
}
log_transformation(log, RULE,
                   f"{len(leads)} axes ({summary['n_oligodendrocyte_side']} oligodendrocyte-side, "
                   f"{summary['n_neuron_side']} neuron-side), tiers={tier_counts}")

prov = stamp_artifact(
    output_path=snakemake.output.table,  # type: ignore[name-defined]
    rule_name=RULE,
    input_paths=[
        snakemake.input.full_lr,  # type: ignore[name-defined]
        snakemake.input.capped_lr,  # type: ignore[name-defined]
        snakemake.input.full_ann,  # type: ignore[name-defined]
        snakemake.input.capped_ann,  # type: ignore[name-defined]
        snakemake.input.drug,  # type: ignore[name-defined]
    ],
    tool_versions={"pandas": pd.__version__},
    parameters=summary,
    description=("Cross-arm reproducible nerve-immune signalling axes tiered by "
                 "druggability (post-compartment-fix lead shortlist)"),
    ontology_operation="operation:3501",  # EDAM: Enrichment / comparison
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]
log_transformation(log, RULE, "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.table])  # type: ignore[name-defined]
