"""Three-way nerve-tumor-immune cell-cell communication (LIANA+ consensus rank).

Extends the tumor-nerve interaction analysis (``nerve_tumor_interaction.py``) by
bringing the tumor-associated immune compartment into the picture. Surfaces
candidate ligand-receptor axes across all three cross-compartment interfaces:

    nerve  <->  tumor      (malignant ↔ nerve Leiden clusters)
    nerve  <->  immune     (nerve clusters ↔ immune subtypes)
    immune <->  tumor      (immune subtypes ↔ malignant)

Approach
--------
1. Assemble three compartments into one AnnData:
   - **tumor**  : ``is_malignant == True`` cells from ``malignancy_labeled.h5ad``
     (single ``malignant`` group).
   - **nerve**  : 27 nerve Leiden clusters from ``nerve_cells.h5ad``
     (per-cluster groups ``nerve_c{N}``).
   - **immune** : subtype groups from ``immune_cells_labeled.h5ad``
     (``immune_{subtype}`` — microglia / tam / t_cell / nk_cell / dendritic).
2. Use HGNC symbols as var index (LIANA consensus resource is HGNC-keyed).
3. Run ``liana.mt.rank_aggregate`` restricting the source/target search space to
   **cross-compartment** directional pairings only (within-compartment pairs,
   e.g. nerve_c1↔nerve_c2, are excluded — they are not the question here).
4. Persist the full LR table (tagged with ``compartment_pair`` and ``direction``),
   the top-N magnitude-ranked pairs per directional pairing, a compartment-pair
   significance heatmap, and a LIANA dotplot of the top consensus interactions.
"""

import os
import sys

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

matplotlib.use("Agg")

import liana as li  # noqa: E402

sys.path.insert(0, "workflow/scripts")
from counts_utils import LOG1P_MAX_PLAUSIBLE, is_log1p_scale  # noqa: E402
from fair_utils import (  # noqa: E402
    nerve_group_sort_key,  # noqa: E402
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

# Tunables --------------------------------------------------------------------
N_PERMS = 1000
EXPR_PROP = 0.10
TOP_N_PER_PAIRING = 10
MAGNITUDE_RANK_SIG = 0.05
RESOURCE_NAME = "consensus"
NORMALIZE_TARGET_SUM = 1e4    # counts-per-10k, applied only when normalize_counts

log = snakemake.log[0]  # type: ignore[name-defined]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)  # type: ignore[name-defined]


# --- Load and label the three compartments -----------------------------------
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Loading {snakemake.input.malig}, {snakemake.input.nerve}, {snakemake.input.immune}",  # type: ignore[name-defined]
)
# EXPRESSION comes from the shared parent; the compartment files supply only
# LABELS. Two reasons, both load-bearing:
#
# 1. CORRECTNESS. The three files are on DIFFERENT SCALES.
#    malignancy_labeled.h5ad carries raw counts, while nerve_cell_subset and
#    immune_cell_subset each ran normalize_total + log1p before writing. The old
#    code concatenated them and then normalized the WHOLE object, so tumor cells
#    were transformed once and nerve/immune cells twice. Every cross-compartment
#    ligand-receptor comparison was therefore between differently-transformed
#    data — silently, since the combined X.max() is dominated by the raw tumor
#    cells and so looks like counts to any scale check.
#
# 2. MEMORY. All three compartments are subsets of the same parent, so
#    concatenating them duplicated the payload: parent + three subsets + the
#    concat output. Once the immune compartment grew to 593k cells (T/NK/B cells
#    joined it), the combined object reached 952k of the cohort's 1,006,344
#    cells and the rule was SIGKILLed at 30+ GB. Reading the parent backed and
#    materialising the labelled subset once costs ~16 GB instead.
parent = ad.read_h5ad(snakemake.input.malig, backed="r")  # type: ignore[name-defined]
nerve_obs = ad.read_h5ad(snakemake.input.nerve, backed="r").obs  # type: ignore[name-defined]
immune_obs = ad.read_h5ad(snakemake.input.immune, backed="r").obs  # type: ignore[name-defined]

# --- mint the group label for each compartment ------------------------------
# Neurons are carried as ONE group, glia per cluster. Splitting ~4.3k neurons
# spread across 170 donors into per-cluster groups would give LIANA a handful of
# cells per group and turn sampling noise into "interactions"; pooling glia
# would throw away real cluster structure. The `nerve_` prefix is kept on both
# so the compartment key stays comparable with the pinned v1.3.0 reference —
# renaming it would break every concordance pair.
if "nerve_subcompartment" in nerve_obs.columns:
    nerve_labels = pd.Series(
        np.where(
            nerve_obs["nerve_subcompartment"].astype(str).eq("neuron").to_numpy(),
            "nerve_neuron",
            "nerve_c" + nerve_obs["nerve_leiden"].astype(str),
        ),
        index=nerve_obs.index,
    )
else:
    nerve_labels = "nerve_c" + nerve_obs["nerve_leiden"].astype(str)
immune_labels = "immune_" + immune_obs["immune_subtype"].astype(str)

tumor_idx = parent.obs_names[parent.obs["is_malignant"].astype(bool).to_numpy()]
nerve_idx = nerve_labels.index
immune_idx = immune_labels.index

# The three masks must partition their cells: nerve and immune both exclude
# is_malignant upstream, so any overlap means a mask changed meaning.
for a_name, a_idx, b_name, b_idx in (
    ("tumor", tumor_idx, "nerve", nerve_idx),
    ("tumor", tumor_idx, "immune", immune_idx),
    ("nerve", nerve_idx, "immune", immune_idx),
):
    overlap = a_idx.intersection(b_idx)
    if len(overlap):
        raise RuntimeError(
            f"[FAIR-ALERT] {a_name} and {b_name} compartments share {len(overlap)} "
            f"cells (e.g. {list(overlap[:5])}). Compartments must be disjoint or a "
            "cell would appear on both sides of its own interaction."
        )

cell_label = pd.Series(pd.NA, index=parent.obs_names, dtype=object)
compartment = pd.Series(pd.NA, index=parent.obs_names, dtype=object)
cell_label.loc[tumor_idx] = "malignant"
compartment.loc[tumor_idx] = "tumor"
cell_label.loc[nerve_idx] = nerve_labels.reindex(nerve_idx).to_numpy()
compartment.loc[nerve_idx] = "nerve"
cell_label.loc[immune_idx] = immune_labels.reindex(immune_idx).to_numpy()
compartment.loc[immune_idx] = "immune"

keep_cells = cell_label.notna().to_numpy()

# Gene selection is done on the backed view too, so the payload is materialised
# exactly once. Reindex var by HGNC symbol (LIANA's consensus resource is
# HGNC-keyed), dropping unmapped and duplicate symbols.
if "gene_symbol" not in parent.var.columns:
    raise KeyError("gene_symbol missing from var — cannot run LIANA on Ensembl IDs")
_sym = parent.var["gene_symbol"]
keep_genes = (_sym.notna() & ~_sym.astype(str).duplicated(keep="first")).to_numpy()

n_tumor_cells = int(len(tumor_idx))
n_nerve_cells = int(len(nerve_idx))
n_immune_cells = int(len(immune_idx))

log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Sourcing expression from the parent ({parent.n_obs:,} cells): keeping "
    f"{int(keep_cells.sum()):,} labelled cells x {int(keep_genes.sum()):,} "
    f"symbol-mapped genes",
)

combined = parent[keep_cells, keep_genes].to_memory()
del parent
combined.var = combined.var.copy()
combined.var.index = combined.var["gene_symbol"].astype(str).values
combined.var.index.name = "gene_symbol"
combined.obs["cell_label"] = cell_label[keep_cells].to_numpy()
combined.obs["compartment"] = compartment[keep_cells].to_numpy()
combined.obs_names_make_unique()

nerve_group_count = int(pd.Series(nerve_labels).nunique())

log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Compartment sizes — tumor: {n_tumor_cells:,}, nerve: {n_nerve_cells:,} "
    f"({nerve_group_count} groups), immune: {n_immune_cells:,} "
    f"({int(immune_labels.nunique())} subtypes)",
)
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Combined AnnData: {combined.n_obs:,} cells x {combined.n_vars:,} genes "
    "(materialised once from the parent, so all three compartments sit on one "
    "consistent raw-count scale)",
)

# --- Normalization: LIANA requires log1p input -------------------------------
# Cohort-aware, mirroring the `counts_from_log1p` flag used by scrna_integration
# and nerve_assemble_counts. The reference cohort's .X arrives as Seurat SCT
# log1p and must be left alone; the Census cohort carries raw UMIs end to end
# (correct for scVI, which wants counts) and must be normalized here.
#
# This previously only LOGGED an assumption, which let the Census cohort run to
# completion on raw counts and emit silently invalid results — all 71,189 rows
# with an empty specificity_rank and 13,492 with lr_logfc = inf. The guard below
# now makes that failure mode impossible.
normalize_counts = bool(snakemake.params.normalize_counts)  # type: ignore[name-defined]
x_max_before = float(combined.X.max())

if normalize_counts:
    sc.pp.normalize_total(combined, target_sum=NORMALIZE_TARGET_SUM)
    sc.pp.log1p(combined)
    log_transformation(
        log,
        "nerve_tumor_immune_interaction",
        f"Normalized for LIANA: normalize_total(target_sum={NORMALIZE_TARGET_SUM:g}) "
        f"+ log1p. X.max() {x_max_before:.3f} -> {float(combined.X.max()):.3f}",
    )
else:
    log_transformation(
        log,
        "nerve_tumor_immune_interaction",
        f"normalize_counts=False — using .X as supplied. X.max() = {x_max_before:.3f}",
    )

if not is_log1p_scale(combined.X):
    log_transformation(
        log,
        "nerve_tumor_immune_interaction",
        f"[FAIR-ALERT] X.max() = {float(combined.X.max()):.3f} exceeds the log1p "
        f"plausibility bound ({LOG1P_MAX_PLAUSIBLE}) — refusing to run LIANA.",
        status="FAILURE",
    )
    raise ValueError(
        f"Expression matrix is not on a log1p scale (max = {float(combined.X.max()):.3f}, "
        f"bound = {LOG1P_MAX_PLAUSIBLE}). LIANA's scoring assumes log1p input; running on "
        f"raw counts yields empty specificity ranks and infinite logFCs rather than an "
        f"error. Set `normalize_counts: true` for this cohort. "
        f"See markdowns/blocker_census_liana_raw_counts.md."
    )

# --- Build cross-compartment directional pair list ---------------------------
label_to_compartment = (
    combined.obs[["cell_label", "compartment"]]
    .drop_duplicates()
    .set_index("cell_label")["compartment"]
    .astype(str)
    .to_dict()
)
tumor_groups = sorted(g for g, c in label_to_compartment.items() if c == "tumor")
nerve_groups = sorted(
    (g for g, c in label_to_compartment.items() if c == "nerve"),
    key=nerve_group_sort_key,
)
immune_groups = sorted(g for g, c in label_to_compartment.items() if c == "immune")

pairs_records: list[dict[str, str]] = []
for a_groups, b_groups in [
    (tumor_groups, nerve_groups),
    (tumor_groups, immune_groups),
    (nerve_groups, immune_groups),
]:
    for g1 in a_groups:
        for g2 in b_groups:
            pairs_records.append({"source": g1, "target": g2})
            pairs_records.append({"source": g2, "target": g1})
groupby_pairs = pd.DataFrame(pairs_records)
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Cross-compartment search space: {len(groupby_pairs)} directional pairings "
    f"({len(tumor_groups)} tumor × {len(nerve_groups)} nerve × {len(immune_groups)} immune)",
)

# --- Run LIANA consensus rank ------------------------------------------------
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Running li.mt.rank_aggregate (n_perms={N_PERMS}, expr_prop={EXPR_PROP}, "
    f"resource={RESOURCE_NAME})",
)
li.mt.rank_aggregate(
    combined,
    groupby="cell_label",
    resource_name=RESOURCE_NAME,
    expr_prop=EXPR_PROP,
    groupby_pairs=groupby_pairs,
    n_perms=N_PERMS,
    seed=int(snakemake.params.random_seed),  # type: ignore[name-defined]
    use_raw=False,
    n_jobs=int(snakemake.threads),  # type: ignore[name-defined]
    verbose=True,
    inplace=True,
    key_added="liana_res",
)

lr_full = combined.uns["liana_res"].copy()
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"LIANA returned {len(lr_full)} (source,target,ligand,receptor) rows",
)

# --- Annotate compartments, direction, and cluster/subtype -------------------
lr_full["source_compartment"] = lr_full["source"].map(label_to_compartment)
lr_full["target_compartment"] = lr_full["target"].map(label_to_compartment)

# Defensive: keep only cross-compartment rows (LIANA respects groupby_pairs).
lr_full = lr_full[
    lr_full["source_compartment"] != lr_full["target_compartment"]
].copy()

lr_full["compartment_pair"] = [
    "-".join(sorted((s, t)))
    for s, t in zip(lr_full["source_compartment"], lr_full["target_compartment"])
]
lr_full["direction"] = (
    lr_full["source_compartment"] + "_to_" + lr_full["target_compartment"]
)


def _nerve_cluster(row: pd.Series) -> str:
    # "nerve_" not "nerve_c": the pooled neuron group is `nerve_neuron`, and
    # matching on the old prefix silently blanked its nerve_cluster column.
    for side in (row["source"], row["target"]):
        if str(side).startswith("nerve_"):
            return str(side)
    return ""


def _immune_subtype(row: pd.Series) -> str:
    for side in (row["source"], row["target"]):
        if str(side).startswith("immune_"):
            return str(side).removeprefix("immune_")
    return ""


lr_full["nerve_cluster"] = lr_full.apply(_nerve_cluster, axis=1)
lr_full["immune_subtype"] = lr_full.apply(_immune_subtype, axis=1)

lr_full = lr_full.sort_values(
    by=["compartment_pair", "direction", "magnitude_rank"],
    ascending=[True, True, True],
).reset_index(drop=True)
lr_full.to_csv(snakemake.output.lr_table, index=False)  # type: ignore[name-defined]
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Wrote {len(lr_full)} three-way LR rows to {snakemake.output.lr_table}; "  # type: ignore[name-defined]
    f"compartment_pair counts: {lr_full['compartment_pair'].value_counts().to_dict()}",
)

# --- Top-N per (compartment_pair, direction, nerve_cluster, immune_subtype) ---
top_per = (
    lr_full.sort_values("magnitude_rank")
    .groupby(
        ["compartment_pair", "direction", "nerve_cluster", "immune_subtype"],
        observed=True,
        group_keys=False,
    )
    .head(TOP_N_PER_PAIRING)
    .reset_index(drop=True)
)
top_per.to_csv(snakemake.output.top_pairs, index=False)  # type: ignore[name-defined]
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Top {TOP_N_PER_PAIRING} pairs/pairing → {len(top_per)} rows",
)

# --- Significance heatmap: # sig LR pairs per compartment_pair × direction ----
sig = lr_full[lr_full["magnitude_rank"] < MAGNITUDE_RANK_SIG]
sig_counts = (
    sig.groupby(["compartment_pair", "direction"], observed=True)
    .size()
    .unstack("direction", fill_value=0)
    .astype(int)
)
fig, ax = plt.subplots(
    figsize=(max(5, 0.8 * sig_counts.shape[1] + 3), max(3, 0.6 * sig_counts.shape[0] + 2))
)
sns.heatmap(
    sig_counts,
    ax=ax,
    cmap="rocket_r",
    linewidths=0.3,
    linecolor="white",
    annot=True,
    fmt="d",
    cbar_kws={"label": f"# LR pairs (magnitude_rank < {MAGNITUDE_RANK_SIG})"},
)
ax.set_title(
    "Nerve↔tumor↔immune cell-cell communication\n"
    f"Significant LR pairs per compartment pair × direction\n"
    f"(LIANA+ consensus, {RESOURCE_NAME}, n_perms={N_PERMS})"
)
ax.set_xlabel("Direction (source → target)")
ax.set_ylabel("Compartment pair")
fig.tight_layout()
fig.savefig(snakemake.output.heatmap, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
plt.close(fig)
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Wrote significance heatmap with {int(sig_counts.values.sum())} total sig pairs",
)

# --- LIANA dotplot of top consensus interactions -----------------------------
n_pairs_for_plot = min(25, len(lr_full))
lr_global_top = lr_full.sort_values("magnitude_rank").head(n_pairs_for_plot).copy()
try:
    fig = li.pl.dotplot(
        liana_res=lr_global_top,
        colour="magnitude_rank",
        size="specificity_rank",
        inverse_colour=True,
        inverse_size=True,
        top_n=n_pairs_for_plot,
        orderby="magnitude_rank",
        orderby_ascending=True,
        figure_size=(max(8, 0.5 * combined.obs["cell_label"].nunique()), 0.4 * n_pairs_for_plot + 2),
    )
    fig.save(snakemake.output.dotplot, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
except Exception as exc:  # plotnine API differences across versions
    log_transformation(
        log,
        "nerve_tumor_immune_interaction",
        f"WARNING: LIANA dotplot failed ({exc}); writing fallback bar chart",
        status="WARNING",
    )
    fig, ax = plt.subplots(figsize=(9, 0.35 * n_pairs_for_plot + 2))
    plot_df = lr_global_top.iloc[::-1]
    label = (
        plot_df["source"].astype(str)
        + " → "
        + plot_df["target"].astype(str)
        + " | "
        + plot_df["ligand_complex"].astype(str)
        + "→"
        + plot_df["receptor_complex"].astype(str)
    )
    ax.barh(label.values, -np.log10(plot_df["magnitude_rank"].clip(lower=1e-6).values))
    ax.set_xlabel("-log10(magnitude_rank)")
    ax.set_title(f"Top {n_pairs_for_plot} nerve↔tumor↔immune LR pairs (LIANA+ consensus)")
    fig.tight_layout()
    fig.savefig(snakemake.output.dotplot, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
    plt.close(fig)
log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    f"Wrote dotplot with top-{n_pairs_for_plot} interactions",
)

# --- FAIR provenance ---------------------------------------------------------
verify_artifact(snakemake.output.lr_table)  # type: ignore[name-defined]
prov = stamp_artifact(
    output_path=snakemake.output.lr_table,  # type: ignore[name-defined]
    rule_name="nerve_tumor_immune_interaction",
    input_paths=[snakemake.input.malig, snakemake.input.nerve, snakemake.input.immune],  # type: ignore[name-defined]
    tool_versions={
        "liana": li.__version__,
        "anndata": ad.__version__,
        "pandas": pd.__version__,
    },
    parameters={
        "n_perms": N_PERMS,
        "expr_prop": EXPR_PROP,
        "resource_name": RESOURCE_NAME,
        "magnitude_rank_sig": MAGNITUDE_RANK_SIG,
        "top_n_per_pairing": TOP_N_PER_PAIRING,
        "n_tumor_cells": int(n_tumor_cells),
        "n_nerve_cells": int(n_nerve_cells),
        "n_immune_cells": int(n_immune_cells),
        "n_nerve_clusters": len(nerve_groups),
        "n_immune_subtypes": len(immune_groups),
        "n_directional_pairings": int(len(groupby_pairs)),
        "n_lr_rows_total": int(len(lr_full)),
        "n_lr_rows_sig": int((lr_full["magnitude_rank"] < MAGNITUDE_RANK_SIG).sum()),
        "compartment_pair_counts": lr_full["compartment_pair"].value_counts().to_dict(),
    },
    description="Three-way nerve-tumor-immune ligand-receptor interaction inference (LIANA+ consensus rank)",
    ontology_operation="operation:3501",  # EDAM: Enrichment analysis (closest available)
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]

log_transformation(
    log,
    "nerve_tumor_immune_interaction",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.lr_table,  # type: ignore[name-defined]
        snakemake.output.top_pairs,  # type: ignore[name-defined]
        snakemake.output.heatmap,  # type: ignore[name-defined]
        snakemake.output.dotplot,  # type: ignore[name-defined]
        snakemake.output.provenance,  # type: ignore[name-defined]
    ],
)
