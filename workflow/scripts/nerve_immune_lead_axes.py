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

The generator emitted ``best_mag`` and ``n_rows`` straight off whichever arm each
half was built from, so a single column meant "full arm, QC-passing" for 128 rows
and "capped arm, no QC, neuron-only" for the other 55 -- with nothing in the row
to say which. Notebook 05 Panel E surfaced this as the 128/128-vs-128/167
discrepancy. Both columns are now replaced by arm-labelled ones --
``full_best_mag``, ``capped_best_mag``, ``full_n_rows``, ``capped_n_rows`` -- each
fully populated, because every axis in either half is a cross-arm intersection and
so always has a counterpart in the other arm.

WHAT THIS DOES NOT FIX, and must not be read as fixing: the *selection rule* still
differs between the halves. The oligodendrocyte side's numbers come from a
QC-passing, all-cluster aggregation and the neuron side's from a no-QC,
pooled-neuron-only one. That difference is what defines the two halves, it is
labelled by ``compartment_side``, and it is inherent rather than incidental.
Comparing ``full_best_mag`` across the two halves is still comparing two different
filters -- the arm is now explicit, the filter is not.

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

# Column order of the emitted file.
#
# The generator emitted `best_mag` and `n_rows`, which were arm-inconsistent: full-arm
# for the oligodendrocyte half and capped-arm for the neuron half, with nothing in the
# row saying which. Both are replaced by arm-labelled columns, all four fully populated
# (every axis in each half is a cross-arm intersection, so its counterpart always
# exists). The generator's `qc_pass` is dropped — it existed in neither half.
#
# `selection_arm` / `qc_filter_applied` were added 2026-09-11. The two halves are
# built under genuinely different selection rules (see the module docstring), and
# until now `compartment_side` was the only thing in the row that implied which --
# it takes reading the generator to learn that "oligodendrocyte (QC-passing)" means
# full-arm-with-QC while "neuron(pooled)" means capped-arm-without. A reader
# filtering the CSV or the notebook table could not see it at all. These two make
# the row state its own provenance.
COLS = ["rank_full", "rank_capped", "axis", "tier_v2", "interfaces", "nerve_side",
        "immune", "full_best_mag", "capped_best_mag", "full_n_rows", "capped_n_rows",
        "min_pval", "agents_flagged", "glioma_trials", "withdrawn",
        "selection_arm", "qc_filter_applied"]


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


def _select(df: pd.DataFrame, side: str, sort_rank: str,
            selection_arm: str, qc_filter_applied: bool) -> pd.DataFrame:
    """Project to the emitted column order and sort within the compartment side.

    `selection_arm` / `qc_filter_applied` record which arm supplied this half's rows
    and whether batch QC gated them, so the asymmetry between the halves is readable
    from the row rather than only from this generator.
    """
    df = df.copy()
    df["compartment_side"] = side
    df["selection_arm"] = selection_arm
    df["qc_filter_applied"] = qc_filter_applied
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
# Rows come from the full arm, so its aggregates are already on the frame; the capped
# arm's are mapped across. Every shared axis is present in both by construction.
lead["full_best_mag"] = lead["best_mag"]
lead["full_n_rows"] = lead["n_rows"]
lead["capped_best_mag"] = lead.axis.map(capped_by_axis.best_mag)
lead["capped_n_rows"] = lead.axis.map(capped_by_axis.n_rows)
lead["nerve_side"] = lead.nerve_types.replace("", "—")

# --- Neuron half: cross-arm shared, rows from the CAPPED arm -------------------
neuron_shared = set(pn["full"].axis) & set(pn["capped"].axis)
full_neuron_by_axis = pn["full"].set_index("axis")
neu = pn["capped"][pn["capped"].axis.isin(neuron_shared)].copy().reset_index(drop=True)
neu["rank_capped"] = neu.index + 1
neu["rank_full"] = neu.axis.map({a: i + 1 for i, a in enumerate(pn["full"].axis)})
# Mirror image of the oligodendrocyte half: rows come from the CAPPED arm here, so it
# is the full arm's aggregates that are mapped across.
neu["capped_best_mag"] = neu["best_mag"]
neu["capped_n_rows"] = neu["n_rows"]
neu["full_best_mag"] = neu.axis.map(full_neuron_by_axis.best_mag)
neu["full_n_rows"] = neu.axis.map(full_neuron_by_axis.n_rows)
neu["nerve_side"] = "neuron(pooled)"

# --- Pinned drug annotation ---------------------------------------------------
drug = pd.read_csv(snakemake.input.drug)  # type: ignore[name-defined]
lead = _join_drug_annotation(lead, drug)
neu = _join_drug_annotation(neu, drug)

# --- Emit ---------------------------------------------------------------------
leads = pd.concat(
    [_select(lead, SIDE_OLIGO, "rank_full", "full", True),
     _select(neu, SIDE_NEURON, "rank_capped", "capped", False)],
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
    "magnitude_columns_are_arm_labelled": True,
    "magnitude_columns_note": (
        "full_best_mag/capped_best_mag/full_n_rows/capped_n_rows replace the "
        "generator's arm-inconsistent best_mag/n_rows and are fully populated. The "
        "SELECTION RULE still differs by compartment_side (oligodendrocyte = "
        "QC-passing all-cluster; neuron = no-QC pooled-neuron), so these columns are "
        "comparable within a compartment side, not across the two."
    ),
    "selection_rule_is_row_level": True,
    "selection_rule_note": (
        "selection_arm and qc_filter_applied state, per row, which arm supplied it "
        "and whether batch QC gated its selection. Added 2026-09-11: the asymmetry "
        "was previously legible only from compartment_side plus the generator "
        "source, so any consumer reading the CSV or the notebook table on its own "
        "could not tell that the two halves are selected under different rules."
    ),
    "selection_rule_counts": {
        f"{arm}/qc={qc}": int(n) for (arm, qc), n in
        leads.groupby(["selection_arm", "qc_filter_applied"]).size().items()
    },
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
