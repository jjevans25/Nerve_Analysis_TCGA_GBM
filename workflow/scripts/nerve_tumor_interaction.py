"""Tumor-nerve cell-cell communication analysis (LIANA+ consensus rank).

Implements §4 / Recommended Next Analysis #4 of
``markdowns/next_steps_interpretation.md``: surface candidate cancer-neuron
and cancer-glia ligand-receptor axes between the malignant compartment and
each non-malignant nerve-cell Leiden cluster.

Approach
--------
1. Pull `is_malignant == True` cells from ``malignancy_labeled.h5ad``
   (single "malignant" group) and concatenate with all 36 nerve-cell
   Leiden clusters from ``nerve_cells.h5ad`` (per-cluster groups
   ``nerve_c{N}``).
2. Use HGNC symbols as the var index (LIANA's consensus resource is
   HGNC-keyed); drop genes with no symbol.
3. Run ``liana.mt.rank_aggregate`` with the consensus L-R resource and
   restrict the source/target search space to malignant↔nerve pairings
   only (37 groups would be 1,332 unrestricted combinations; restricting
   to malignant↔nerve cuts that to 72 directional pairings — both
   sender→receiver directions retained).
4. Persist the full LR table, the top-10 magnitude-ranked pairs per
   nerve cluster (each direction), a LIANA dotplot of the top consensus
   interactions, and a heatmap of significant LR-pair counts per
   directional pairing.
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
TOP_N_PER_CLUSTER = 10
MAGNITUDE_RANK_SIG = 0.05
RESOURCE_NAME = "consensus"
NORMALIZE_TARGET_SUM = 1e4    # counts-per-10k, applied only when normalize_counts

log = snakemake.log[0]  # type: ignore[name-defined]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)  # type: ignore[name-defined]

# --- Load and combine inputs -------------------------------------------------
log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Loading {snakemake.input.malig} and {snakemake.input.nerve}",  # type: ignore[name-defined]
)
# EXPRESSION comes from the shared parent; nerve_cells.h5ad supplies only LABELS.
#
# CORRECTNESS: the two files are on DIFFERENT SCALES. malignancy_labeled.h5ad
# carries raw counts, while nerve_cell_subset ran normalize_total + log1p before
# writing. The old code concatenated them and then normalized the WHOLE object,
# so malignant cells were transformed once and nerve cells twice — every
# tumor-vs-nerve comparison was between differently-transformed data. It was
# silent because the combined X.max() is dominated by the raw malignant cells and
# so looks like counts to any scale check.
#
# MEMORY: both compartments are subsets of the same parent, so concatenating them
# duplicated the payload. Reading the parent backed and materialising the
# labelled subset once is also what keeps the three-way twin of this rule inside
# 36 GB (it was SIGKILLed at 30+ GB on 2026-08-06).
parent = ad.read_h5ad(snakemake.input.malig, backed="r")  # type: ignore[name-defined]
nerve_obs = ad.read_h5ad(snakemake.input.nerve, backed="r").obs  # type: ignore[name-defined]

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

malig_idx = parent.obs_names[parent.obs["is_malignant"].astype(bool).to_numpy()]
nerve_idx = nerve_labels.index
overlap = malig_idx.intersection(nerve_idx)
if len(overlap):
    raise RuntimeError(
        f"[FAIR-ALERT] tumor and nerve compartments share {len(overlap)} cells "
        f"(e.g. {list(overlap[:5])}). They must be disjoint or a cell would appear "
        "on both sides of its own interaction."
    )

cell_label = pd.Series(pd.NA, index=parent.obs_names, dtype=object)
src = pd.Series(pd.NA, index=parent.obs_names, dtype=object)
cell_label.loc[malig_idx] = "malignant"
src.loc[malig_idx] = "malignant"
cell_label.loc[nerve_idx] = nerve_labels.reindex(nerve_idx).to_numpy()
src.loc[nerve_idx] = "nerve"
keep_cells = cell_label.notna().to_numpy()

# Gene selection on the backed view too, so the payload is materialised exactly
# once. LIANA's consensus resource is HGNC-keyed.
if "gene_symbol" not in parent.var.columns:
    raise KeyError("gene_symbol missing from var — cannot run LIANA on Ensembl IDs")
_sym = parent.var["gene_symbol"]
keep_genes = (_sym.notna() & ~_sym.astype(str).duplicated(keep="first")).to_numpy()

n_malignant_cells = int(len(malig_idx))
n_nerve_cells = int(len(nerve_idx))

log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Malignant subset: {n_malignant_cells:,} cells; nerve subset: "
    f"{n_nerve_cells:,} cells ({int(pd.Series(nerve_labels).nunique())} groups)",
)

combined = parent[keep_cells, keep_genes].to_memory()
del parent
combined.var = combined.var.copy()
combined.var.index = combined.var["gene_symbol"].astype(str).values
combined.var.index.name = "gene_symbol"
combined.obs["cell_label"] = cell_label[keep_cells].to_numpy()
combined.obs["src_h5ad"] = src[keep_cells].to_numpy()
combined.obs_names_make_unique()
log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Combined AnnData: {combined.n_obs:,} cells x {combined.n_vars:,} genes "
    "(materialised once from the parent, so both compartments sit on one "
    "consistent raw-count scale)",
)

# --- Normalization: LIANA requires log1p input -------------------------------
# Cohort-aware, mirroring the `counts_from_log1p` flag used by scrna_integration
# and nerve_assemble_counts. The reference cohort's .X arrives as Seurat SCT
# log1p and must be left alone; the Census cohort carries raw UMIs end to end
# (correct for scVI, which wants counts) and must be normalized here.
#
# This previously only LOGGED an assumption, which let the Census cohort run to
# completion on raw counts and emit silently invalid results. The guard below
# now makes that failure mode impossible.
normalize_counts = bool(snakemake.params.normalize_counts)  # type: ignore[name-defined]
x_max_before = float(combined.X.max())

if normalize_counts:
    sc.pp.normalize_total(combined, target_sum=NORMALIZE_TARGET_SUM)
    sc.pp.log1p(combined)
    log_transformation(
        log,
        "nerve_tumor_interaction",
        f"Normalized for LIANA: normalize_total(target_sum={NORMALIZE_TARGET_SUM:g}) "
        f"+ log1p. X.max() {x_max_before:.3f} -> {float(combined.X.max()):.3f}",
    )
else:
    log_transformation(
        log,
        "nerve_tumor_interaction",
        f"normalize_counts=False — using .X as supplied. X.max() = {x_max_before:.3f}",
    )

if not is_log1p_scale(combined.X):
    log_transformation(
        log,
        "nerve_tumor_interaction",
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

# --- Restrict groupby pairs to malignant↔nerve only --------------------------
nerve_groups = sorted(
    combined.obs.loc[combined.obs["cell_label"] != "malignant", "cell_label"].unique(),
    key=nerve_group_sort_key,
)
pairs_records = []
for g in nerve_groups:
    pairs_records.append({"source": "malignant", "target": g})
    pairs_records.append({"source": g, "target": "malignant"})
groupby_pairs = pd.DataFrame(pairs_records)
log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Restricting search to {len(groupby_pairs)} malignant↔nerve directional pairs "
    f"(over {len(nerve_groups)} nerve clusters)",
)

# --- Run LIANA consensus rank ------------------------------------------------
log_transformation(
    log,
    "nerve_tumor_interaction",
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
    "nerve_tumor_interaction",
    f"LIANA returned {len(lr_full)} (source,target,ligand,receptor) rows",
)

# Defensive filter: only retain malignant↔nerve directions (LIANA respects
# groupby_pairs, but persist the contract explicitly).
nerve_set = set(nerve_groups)
lr_full = lr_full[
    ((lr_full["source"] == "malignant") & (lr_full["target"].isin(nerve_set)))
    | ((lr_full["target"] == "malignant") & (lr_full["source"].isin(nerve_set)))
].copy()

# Add a direction column for downstream slicing.
lr_full["direction"] = np.where(
    lr_full["source"] == "malignant",
    "malignant_to_nerve",
    "nerve_to_malignant",
)
lr_full["nerve_cluster"] = np.where(
    lr_full["direction"] == "malignant_to_nerve",
    lr_full["target"],
    lr_full["source"],
)

# Sort for stable output.
lr_full = lr_full.sort_values(
    by=["nerve_cluster", "direction", "magnitude_rank"],
    ascending=[True, True, True],
).reset_index(drop=True)

lr_full.to_csv(snakemake.output.lr_table, index=False)  # type: ignore[name-defined]
log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Wrote {len(lr_full)} tumor-nerve LR rows to {snakemake.output.lr_table}",  # type: ignore[name-defined]
)

# --- Top-N per (nerve_cluster, direction) ------------------------------------
top_per = (
    lr_full.sort_values("magnitude_rank")
    .groupby(["nerve_cluster", "direction"], observed=True, group_keys=False)
    .head(TOP_N_PER_CLUSTER)
    .reset_index(drop=True)
)
top_per.to_csv(snakemake.output.top_pairs, index=False)  # type: ignore[name-defined]
log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Top {TOP_N_PER_CLUSTER} pairs/cluster/direction → {len(top_per)} rows",
)

# --- Significance heatmap (#sig LR pairs per directional pairing) ------------
sig = lr_full[lr_full["magnitude_rank"] < MAGNITUDE_RANK_SIG]
sig_counts = (
    sig.groupby(["direction", "nerve_cluster"], observed=True)
    .size()
    .unstack("direction", fill_value=0)
    .reindex(nerve_groups)
    .fillna(0)
    .astype(int)
)
for col in ["malignant_to_nerve", "nerve_to_malignant"]:
    if col not in sig_counts.columns:
        sig_counts[col] = 0
sig_counts = sig_counts[["malignant_to_nerve", "nerve_to_malignant"]]

fig, ax = plt.subplots(figsize=(5, max(6, 0.28 * len(sig_counts))))
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
    f"Tumor↔nerve cell-cell communication\n"
    f"Significant LR pairs per directional pairing\n"
    f"(LIANA+ consensus, {RESOURCE_NAME}, n_perms={N_PERMS})"
)
ax.set_xlabel("Direction")
ax.set_ylabel("Nerve-cell cluster")
fig.tight_layout()
fig.savefig(snakemake.output.heatmap, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
plt.close(fig)
log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Wrote significance heatmap with {int(sig_counts.sum().sum())} total sig pairs across "
    f"{len(nerve_groups)} clusters",
)

# --- LIANA dotplot of top consensus interactions -----------------------------
# Plot top global magnitude-ranked LR pairs across all nerve clusters in one
# figure (LIANA's built-in dotplot expects the in-place uns dict).
n_pairs_for_plot = min(25, len(lr_full))
# Globally rank by magnitude for the dotplot (the on-disk CSV is sorted by
# cluster for readability; the plot wants a single cohort-wide top list).
lr_global_top = lr_full.sort_values("magnitude_rank").head(n_pairs_for_plot).copy()
try:
    plot_df = lr_global_top
    fig = li.pl.dotplot(
        liana_res=plot_df,
        colour="magnitude_rank",
        size="specificity_rank",
        inverse_colour=True,
        inverse_size=True,
        top_n=n_pairs_for_plot,
        orderby="magnitude_rank",
        orderby_ascending=True,
        figure_size=(max(8, 0.4 * len(nerve_groups)), 0.4 * n_pairs_for_plot + 2),
    )
    fig.save(snakemake.output.dotplot, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
except Exception as exc:  # plotnine API differences across versions
    log_transformation(
        log,
        "nerve_tumor_interaction",
        f"WARNING: LIANA dotplot failed ({exc}); writing fallback bar chart",
        status="WARNING",
    )
    fig, ax = plt.subplots(figsize=(8, 0.35 * n_pairs_for_plot + 2))
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
    ax.set_title(f"Top {n_pairs_for_plot} tumor↔nerve LR pairs (LIANA+ consensus)")
    fig.tight_layout()
    fig.savefig(snakemake.output.dotplot, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
    plt.close(fig)

log_transformation(
    log,
    "nerve_tumor_interaction",
    f"Wrote dotplot with top-{n_pairs_for_plot} interactions",
)

# --- FAIR provenance ---------------------------------------------------------
verify_artifact(snakemake.output.lr_table)  # type: ignore[name-defined]
prov = stamp_artifact(
    output_path=snakemake.output.lr_table,  # type: ignore[name-defined]
    rule_name="nerve_tumor_interaction",
    input_paths=[snakemake.input.malig, snakemake.input.nerve],  # type: ignore[name-defined]
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
        "top_n_per_cluster": TOP_N_PER_CLUSTER,
        "n_malignant_cells": int(n_malignant_cells),
        "n_nerve_cells": int(n_nerve_cells),
        "n_nerve_clusters": int(len(nerve_groups)),
        "n_lr_rows_total": int(len(lr_full)),
        "n_lr_rows_sig": int((lr_full["magnitude_rank"] < MAGNITUDE_RANK_SIG).sum()),
    },
    description="Tumor-nerve ligand-receptor interaction inference (LIANA+ consensus rank)",
    ontology_operation="operation:3501",  # EDAM: Enrichment analysis (closest available)
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]

log_transformation(
    log,
    "nerve_tumor_interaction",
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
