"""Evidence pack for the cl15 'ependymal' label sanity check (v1.2.0 Issue 1).

Mirrors scripts/diagnose_failing_clusters.py: gathers per-cluster diagnostic
evidence into a JSON sidecar + figures for hand-synthesis. Does NOT pick a
hypothesis verdict (H1/H2/H3/H4) — that is researcher judgment.

Diagnostics produced (each indexed by cluster in the JSON output):
    1. score_* distributions for the ependymal-dominant clusters (cl11/15/19)
       and astrocyte-dominant controls (cl1/9/14).
    2. Ependymal marker-panel raw expression (% cells expressing, mean
       log-expression) per cluster. Missing panel members are flagged.
    3. score_ependymal − max(other 4 scores) margin distribution for cl15
       cells (from the labeled counts file's 5-score panel).
    4. Sub-Leiden on cl15 cells alone (resolution=0.5, on X_scVI latent).
       Reports sub-cluster sizes, per-sub-cluster score means, and top
       differentially-expressed markers per sub-cluster.
    5. Per-sample composition of cl15 + cross-reference to
       nerve_cluster_sample_purity.csv.

Inputs (read-only, all already on disk from v1.1.0 cut):
    data/processed/nerve_cells.h5ad
    data/processed/nerve_cells_counts_labeled.h5ad
    results/tables/nerve_cluster_markers.csv
    results/tables/nerve_cluster_sample_purity.csv
    config/config.yaml

Outputs:
    results/tables/cl15_ependymal_diagnosis.json
    results/figures/cl15_score_boxplots.png
    results/figures/cl15_ependymal_panel_dotplot.png
    results/figures/cl15_argmax_margin_histogram.png
    results/figures/cl15_patient_composition.png

Not a Snakemake rule — top-level diagnostic, run once.
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse

ROOT = Path(__file__).resolve().parent.parent

EPENDYMAL_DOMINANT = ["11", "15", "19"]
ASTROCYTE_CONTROLS = ["1", "9", "14"]
ALL_FOCUS = EPENDYMAL_DOMINANT + ASTROCYTE_CONTROLS

FIG_DIR = ROOT / "results/figures"
TBL_DIR = ROOT / "results/tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TBL_DIR.mkdir(parents=True, exist_ok=True)


def _parse_args() -> argparse.Namespace:
    """Parse the optional version tag used to suffix all output artifacts.

    Lets the same evidence pack run against successive pipeline cuts without
    clobbering the previous cut's JSON/figures (FAIR: preserve prior artifacts
    for side-by-side comparison). Empty tag reproduces the original v1.1.0
    filenames.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--version-tag",
        default="",
        help="Suffix appended to output filenames, e.g. 'v1_2_0' -> "
        "cl15_ependymal_diagnosis_v1_2_0.json. Default '' = v1.1.0 names.",
    )
    return p.parse_args()


ARGS = _parse_args()
SUFFIX = f"_{ARGS.version_tag}" if ARGS.version_tag else ""


def tbl_path(stem: str, ext: str = "json") -> Path:
    """Versioned table path under results/tables."""
    return TBL_DIR / f"{stem}{SUFFIX}.{ext}"


def fig_path(stem: str) -> Path:
    """Versioned figure path under results/figures."""
    return FIG_DIR / f"{stem}{SUFFIX}.png"


# ---------- inputs -----------------------------------------------------------
def load_config_ependymal_panel() -> list[str]:
    cfg = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    return list(cfg["nerve_cells"]["markers"]["ependymal"])


print("=" * 76)
print("Loading inputs")
print("=" * 76)

EPENDYMAL_PANEL = load_config_ependymal_panel()
print(f"  ependymal panel (config): {EPENDYMAL_PANEL}")

# Backed read of the main nerve AnnData. We materialise only the obs columns
# we need + the panel-gene expression columns to keep RAM bounded.
adata_nerve_path = ROOT / "data/processed/nerve_cells.h5ad"
adata_label_path = ROOT / "data/processed/nerve_cells_counts_labeled.h5ad"

adata_nerve = ad.read_h5ad(adata_nerve_path, backed="r")
adata_label = ad.read_h5ad(adata_label_path, backed="r")
print(f"  nerve  : {adata_nerve.shape}")
print(f"  label  : {adata_label.shape}")

assert (adata_nerve.obs_names == adata_label.obs_names).all(), (
    "Cell index mismatch between nerve_cells.h5ad and "
    "nerve_cells_counts_labeled.h5ad — diagnostic margin step requires "
    "aligned indices."
)

SCORE_COLS_FULL = [c for c in adata_nerve.obs.columns if c.startswith("score_")]
SCORE_COLS_LABEL = [c for c in adata_label.obs.columns if c.startswith("score_")]
print(f"  nerve score_*   ({len(SCORE_COLS_FULL)}): {SCORE_COLS_FULL}")
print(f"  label score_*   ({len(SCORE_COLS_LABEL)}): {SCORE_COLS_LABEL}")

obs_nerve = adata_nerve.obs[
    ["nerve_leiden", "sample_id", "cell_type_predicted",
     "cell_type_confidence", "pct_counts_mt", "n_genes_by_counts"]
    + SCORE_COLS_FULL
].copy()
obs_nerve["nerve_leiden"] = obs_nerve["nerve_leiden"].astype(str)
obs_nerve["sample_id"]    = obs_nerve["sample_id"].astype(str)

obs_label = adata_label.obs[
    ["cell_type", "cell_type_score_max"] + SCORE_COLS_LABEL
].copy()

# var index → gene_symbol map (uppercased).
gene_symbols = adata_nerve.var["gene_symbol"].astype(str).str.upper()
sym_to_idx: dict[str, int] = {}
for i, s in enumerate(gene_symbols.tolist()):
    if s and s != "NAN" and s not in sym_to_idx:
        sym_to_idx[s] = i

panel_present: dict[str, int] = {}
panel_missing: list[str] = []
for g in EPENDYMAL_PANEL:
    idx = sym_to_idx.get(g.upper())
    if idx is None:
        panel_missing.append(g)
    else:
        panel_present[g] = idx
print(f"  panel present in var: {list(panel_present.keys())}")
print(f"  panel MISSING from var (filtered upstream): {panel_missing}")


# ---------- diagnostic 1: score_* distributions per cluster ------------------
print("\n" + "=" * 76)
print("Diagnostic 1: score_* distributions for focus clusters")
print("=" * 76)


def score_summary(obs_df: pd.DataFrame, cluster: str) -> dict:
    """Per-cluster mean/median/quartile for each score_* column."""
    sub = obs_df[obs_df["nerve_leiden"] == cluster]
    out = {"n_cells": int(len(sub))}
    for c in SCORE_COLS_FULL:
        v = sub[c].to_numpy(dtype=float)
        out[c] = {
            "mean": float(np.nanmean(v)),
            "median": float(np.nanmedian(v)),
            "q25": float(np.nanpercentile(v, 25)),
            "q75": float(np.nanpercentile(v, 75)),
        }
    # also majority cell_type_predicted in cluster
    out["cell_type_predicted_top"] = (
        sub["cell_type_predicted"].astype(str).value_counts(normalize=True)
        .round(3).head(5).to_dict()
    )
    out["cell_type_confidence_median"] = float(
        np.nanmedian(sub["cell_type_confidence"].to_numpy(dtype=float))
    )
    return out


score_dist = {cl: score_summary(obs_nerve, cl) for cl in ALL_FOCUS}
for cl in ALL_FOCUS:
    print(f"  cl{cl:>2}  n={score_dist[cl]['n_cells']:>5}  "
          f"top_predicted={score_dist[cl]['cell_type_predicted_top']}")

# Figure 1: per-cluster boxplots, one panel per cluster, score_* on x-axis.
fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharey=True)
short_labels = [c.replace("score_", "") for c in SCORE_COLS_FULL]
for ax, cl in zip(axes.flat, ALL_FOCUS):
    sub = obs_nerve[obs_nerve["nerve_leiden"] == cl]
    data = [sub[c].to_numpy(dtype=float) for c in SCORE_COLS_FULL]
    ax.boxplot(data, tick_labels=short_labels, showfliers=False)
    kind = "ependymal-dom" if cl in EPENDYMAL_DOMINANT else "astro-control"
    ax.set_title(f"cl{cl} ({kind}, n={len(sub)})")
    ax.tick_params(axis="x", rotation=70)
    ax.axhline(0, color="grey", lw=0.5)
fig.suptitle("score_* distributions per cluster (cl15 sanity check)")
fig.tight_layout()
fig1_path = fig_path("cl15_score_boxplots")
fig.savefig(fig1_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"  -> {fig1_path}")


# ---------- diagnostic 2: ependymal panel expression -------------------------
print("\n" + "=" * 76)
print("Diagnostic 2: ependymal panel raw expression per cluster")
print("=" * 76)


def panel_expression_for_clusters(
    cluster_ids: list[str],
    panel_idx: dict[str, int],
    backed_adata: ad.AnnData,
    obs_df: pd.DataFrame,
) -> dict:
    """For each cluster, % cells expressing and mean log-expression of each
    panel gene. Uses backed .X (log1p-normalised on disk)."""
    out: dict = {}
    for cl in cluster_ids:
        cell_mask = (obs_df["nerve_leiden"] == cl).to_numpy()
        if cell_mask.sum() == 0:
            out[cl] = {"n_cells": 0}
            continue
        cell_pos = np.where(cell_mask)[0]
        per_gene: dict = {}
        for gene, vidx in panel_idx.items():
            # Backed slice along cells × single gene → dense vector.
            col = backed_adata.X[cell_pos, vidx]
            if sparse.issparse(col):
                col = col.toarray().ravel()
            else:
                col = np.asarray(col).ravel()
            per_gene[gene] = {
                "pct_expressing": float(np.mean(col > 0) * 100.0),
                "mean_log_expr": float(np.mean(col)),
                "p90_log_expr": float(np.percentile(col, 90)),
            }
        out[cl] = {"n_cells": int(cell_mask.sum()), "genes": per_gene}
    return out


# Cohort baseline: stratify "any expression" against a non-target cluster pool.
panel_expr = panel_expression_for_clusters(
    ALL_FOCUS, panel_present, adata_nerve, obs_nerve
)
for cl in ALL_FOCUS:
    print(f"\n  cl{cl} (n={panel_expr[cl]['n_cells']}):")
    for g in panel_present:
        s = panel_expr[cl]["genes"][g]
        print(f"    {g:<8} pct={s['pct_expressing']:>5.1f}%  "
              f"mean={s['mean_log_expr']:>5.3f}  p90={s['p90_log_expr']:>5.3f}")

# Figure 2: dotplot — rows = clusters, cols = panel genes.
genes_in_order = list(panel_present.keys())
pct_mat = np.array(
    [[panel_expr[cl]["genes"][g]["pct_expressing"] for g in genes_in_order]
     for cl in ALL_FOCUS]
)
mean_mat = np.array(
    [[panel_expr[cl]["genes"][g]["mean_log_expr"] for g in genes_in_order]
     for cl in ALL_FOCUS]
)
fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(genes_in_order) + 3),
                                0.6 * len(ALL_FOCUS) + 2))
# dot size = pct_expressing, color = mean_log_expr
for i, cl in enumerate(ALL_FOCUS):
    for j, g in enumerate(genes_in_order):
        size = max(20.0, pct_mat[i, j] * 12.0)
        ax.scatter(j, i, s=size, c=[mean_mat[i, j]],
                   cmap="viridis", vmin=0.0, vmax=max(mean_mat.max(), 1e-3),
                   edgecolors="black", linewidths=0.4)
ax.set_xticks(range(len(genes_in_order)))
ax.set_xticklabels(genes_in_order, rotation=45, ha="right")
ax.set_yticks(range(len(ALL_FOCUS)))
ax.set_yticklabels([f"cl{cl}" for cl in ALL_FOCUS])
ax.set_title("Ependymal panel expression (dot=%cells, color=mean log1p)")
if panel_missing:
    ax.text(0, -1.5,
            f"Panel genes missing from var: {', '.join(panel_missing)}",
            color="darkred", fontsize=9)
fig.tight_layout()
fig2_path = fig_path("cl15_ependymal_panel_dotplot")
fig.savefig(fig2_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"\n  -> {fig2_path}")


# ---------- diagnostic 3: argmax margin distribution for cl15 ----------------
print("\n" + "=" * 76)
print("Diagnostic 3: ependymal argmax-margin distribution for cl15")
print("=" * 76)

cl15_mask = (obs_nerve["nerve_leiden"] == "15").to_numpy()
n_cl15 = int(cl15_mask.sum())
print(f"  cl15 n_cells = {n_cl15}")

# Use the 5-score panel from the labeled file (astrocyte / ependymal / neuron /
# oligodendrocyte / opc — the panel used by nerve_celltype_labels.py).
EPENDYMAL_SCORE = "score_ependymal"
other_scores = [c for c in SCORE_COLS_LABEL if c != EPENDYMAL_SCORE]
print(f"  ependymal-vs-others using label scores: {other_scores}")

label_scores_cl15 = obs_label.loc[cl15_mask, SCORE_COLS_LABEL].to_numpy(dtype=float)
ep_score = label_scores_cl15[:, SCORE_COLS_LABEL.index(EPENDYMAL_SCORE)]
other_idx = [SCORE_COLS_LABEL.index(c) for c in other_scores]
runner_up = label_scores_cl15[:, other_idx].max(axis=1)
margin = ep_score - runner_up

margin_summary = {
    "n_cells": n_cl15,
    "mean": float(np.mean(margin)),
    "median": float(np.median(margin)),
    "p10": float(np.percentile(margin, 10)),
    "p25": float(np.percentile(margin, 25)),
    "p75": float(np.percentile(margin, 75)),
    "p90": float(np.percentile(margin, 90)),
    "frac_margin_below_0":   float(np.mean(margin < 0.0)),
    "frac_margin_below_0p02": float(np.mean(margin < 0.02)),
    "frac_margin_below_0p05": float(np.mean(margin < 0.05)),
    "frac_margin_below_0p10": float(np.mean(margin < 0.10)),
    "labeled_as_ependymal_in_obs": float(
        (obs_label.loc[cl15_mask, "cell_type"].astype(str) == "ependymal").mean()
    ),
}
print(f"  margin median = {margin_summary['median']:.4f}")
print(f"  fraction with margin < 0.05 = "
      f"{margin_summary['frac_margin_below_0p05']:.3f}")
print(f"  fraction labeled 'ependymal' in obs.cell_type = "
      f"{margin_summary['labeled_as_ependymal_in_obs']:.3f}")

fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(margin, bins=60, color="#4472c4", edgecolor="black", linewidth=0.3)
for thr, label in [(0.0, "0"), (0.05, "0.05"), (0.10, "0.10")]:
    ax.axvline(thr, color="red", linestyle="--", linewidth=0.8,
               label=f"margin = {label}")
ax.set_xlabel("score_ependymal − max(score_astro, neuron, oligo, opc)")
ax.set_ylabel("cl15 cell count")
ax.set_title(f"cl15 (n={n_cl15}) ependymal argmax-margin distribution")
ax.legend()
fig.tight_layout()
fig3_path = fig_path("cl15_argmax_margin_histogram")
fig.savefig(fig3_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"  -> {fig3_path}")


# ---------- diagnostic 4: sub-Leiden on cl15 ---------------------------------
print("\n" + "=" * 76)
print("Diagnostic 4: sub-Leiden split on cl15 (resolution=0.5, on X_scVI)")
print("=" * 76)

# Load cl15 cells into memory (1854 cells × full gene set), but only the
# scVI latent + needed obs. We use latent for neighbours and the raw .X (log1p)
# for downstream DE.
cl15_pos = np.where(cl15_mask)[0]
X_scvi = adata_nerve.obsm["X_scVI"][cl15_pos]
print(f"  X_scVI subset shape: {X_scvi.shape}")

# Materialise the X submatrix for cl15 cells (sparse → keep sparse).
X_cl15 = adata_nerve.X[cl15_pos, :]
if not sparse.issparse(X_cl15):
    X_cl15 = sparse.csr_matrix(X_cl15)

sub = ad.AnnData(
    X=X_cl15,
    obs=adata_nerve.obs.iloc[cl15_pos].copy(),
    var=adata_nerve.var.copy(),
    obsm={"X_scVI": X_scvi},
)
sub.obs_names = adata_nerve.obs_names[cl15_pos]

sc.pp.neighbors(sub, use_rep="X_scVI", random_state=0)
sc.tl.leiden(sub, resolution=0.5, key_added="cl15_subleiden", random_state=0,
             flavor="igraph", n_iterations=2, directed=False)
sub_sizes = sub.obs["cl15_subleiden"].astype(str).value_counts().sort_index()
print(f"  sub-cluster sizes: {sub_sizes.to_dict()}")

# Per-sub-cluster score means (use the full 11-score panel from nerve obs).
sub_scores: dict = {}
for s_id, n in sub_sizes.items():
    s_mask = (sub.obs["cl15_subleiden"].astype(str) == s_id).to_numpy()
    sub_scores[s_id] = {
        "n_cells": int(n),
        "score_means": {
            c: float(np.nanmean(sub.obs.loc[s_mask, c].to_numpy(dtype=float)))
            for c in SCORE_COLS_FULL
        },
    }

# DE markers per sub-cluster (top 10 by score). Only meaningful if ≥2 splits.
sub_markers: dict = {}
if sub_sizes.shape[0] >= 2:
    sc.tl.rank_genes_groups(sub, "cl15_subleiden", method="wilcoxon",
                            n_genes=15, use_raw=False)
    rgg = sub.uns["rank_genes_groups"]
    for s_id in sub_sizes.index:
        names = list(rgg["names"][s_id])
        # Map ensembl → symbol
        symbols = []
        for ens in names:
            v_loc = sub.var_names.get_loc(ens) if ens in sub.var_names else None
            sym = (sub.var["gene_symbol"].iloc[v_loc]
                   if v_loc is not None else None)
            symbols.append(str(sym) if sym else ens)
        scores_s = [float(s) for s in rgg["scores"][s_id]]
        sub_markers[s_id] = {
            "top_genes": symbols,
            "top_scores": scores_s,
        }
        print(f"  sub-cluster {s_id} (n={sub_sizes[s_id]}) top: "
              f"{symbols[:8]}")
else:
    print("  Only 1 sub-cluster found at resolution=0.5; no DE computed.")

# Patient composition per sub-cluster
sub_patient: dict = {}
for s_id in sub_sizes.index:
    s_mask = (sub.obs["cl15_subleiden"].astype(str) == s_id).to_numpy()
    sids = sub.obs.loc[s_mask, "sample_id"].astype(str)
    vc = sids.value_counts(normalize=True).round(3)
    sub_patient[s_id] = vc.head(5).to_dict()

sub_leiden_summary = {
    "n_sub_clusters": int(sub_sizes.shape[0]),
    "sizes": sub_sizes.to_dict(),
    "score_means": sub_scores,
    "top_markers": sub_markers,
    "patient_composition": sub_patient,
}


# ---------- diagnostic 5: cl15 patient composition + purity row -------------
print("\n" + "=" * 76)
print("Diagnostic 5: cl15 patient composition")
print("=" * 76)

cl15_samples = (
    obs_nerve.loc[cl15_mask, "sample_id"]
    .astype(str).value_counts(normalize=True).round(4)
)
print("  top 10 cl15 sample fractions:")
for s, f in cl15_samples.head(10).items():
    print(f"    {s[:8]}…  {f:.4f}")

purity = pd.read_csv(ROOT / "results/tables/nerve_cluster_sample_purity.csv")
purity["cluster"] = purity["cluster"].astype(str)
purity_cl15 = purity[purity["cluster"] == "15"].iloc[0].to_dict()
print(f"\n  purity row for cl15: {purity_cl15}")

fig, ax = plt.subplots(figsize=(10, 5))
top10 = cl15_samples.head(10)
ax.bar(range(len(top10)), top10.values, color="#4472c4", edgecolor="black")
ax.set_xticks(range(len(top10)))
ax.set_xticklabels([s[:8] for s in top10.index], rotation=45, ha="right")
ax.set_ylabel("fraction of cl15 cells")
ax.set_title(f"cl15 patient composition "
             f"(dominant_fraction={purity_cl15['dominant_sample_fraction']:.3f}, "
             f"n_contrib={purity_cl15['n_contributing_samples']})")
fig.tight_layout()
fig5_path = fig_path("cl15_patient_composition")
fig.savefig(fig5_path, dpi=130, bbox_inches="tight")
plt.close(fig)
print(f"  -> {fig5_path}")


# ---------- per-hypothesis evidence summary ---------------------------------
# Frame evidence as 'consistent / inconsistent' per the plan's H1–H4. Picking
# a winner is researcher judgment; this section is a structured pointer only.
print("\n" + "=" * 76)
print("Evidence-vs-hypothesis pointer (researcher reads from here)")
print("=" * 76)

cl15_panel_ciliated = [g for g in ("DNAH7", "DNAH9", "DNAH11", "CFAP54", "RSPH1")
                       if g in panel_present]
cl15_panel_tf = [g for g in ("FOXJ1", "RFX3") if g in panel_present]

def avg_panel(genes: list[str]) -> float:
    if not genes:
        return float("nan")
    return float(np.mean([
        panel_expr["15"]["genes"][g]["pct_expressing"] for g in genes
    ]))

cl15_tf_pct = avg_panel(cl15_panel_tf)
cl15_dynein_pct = avg_panel(cl15_panel_ciliated)
print(f"  cl15 avg %expressing FOXJ1/RFX3-style TFs : {cl15_tf_pct:.1f}%")
print(f"  cl15 avg %expressing DNAH/CFAP/RSPH1     : {cl15_dynein_pct:.1f}%")
print(f"  cl15 ependymal-argmax-margin median       : "
      f"{margin_summary['median']:.4f}")
print(f"  cl15 fraction margin < 0.05               : "
      f"{margin_summary['frac_margin_below_0p05']:.3f}")
print(f"  cl15 sub-Leiden n_sub_clusters @res=0.5   : "
      f"{sub_leiden_summary['n_sub_clusters']}")
print(f"  cl15 dominant sample fraction             : "
      f"{purity_cl15['dominant_sample_fraction']:.3f}")
print()
print("  NOTE: This is evidence, not a verdict. H1/H2/H3/H4 selection is")
print("  a researcher call — see markdowns/DO_THIS_NEXT_v1.2.0_open_issues.md")
print("  §'Issue 1' hypothesis table.")


# ---------- save JSON --------------------------------------------------------
out = {
    "config_ependymal_panel": EPENDYMAL_PANEL,
    "panel_present_in_var": list(panel_present.keys()),
    "panel_missing_from_var": panel_missing,
    "focus_clusters": {
        "ependymal_dominant": EPENDYMAL_DOMINANT,
        "astrocyte_controls": ASTROCYTE_CONTROLS,
    },
    "score_cols_nerve":  SCORE_COLS_FULL,
    "score_cols_label":  SCORE_COLS_LABEL,
    "diagnostic_1_score_distributions": score_dist,
    "diagnostic_2_panel_expression":    panel_expr,
    "diagnostic_3_argmax_margin_cl15":  margin_summary,
    "diagnostic_4_subleiden_cl15":      sub_leiden_summary,
    "diagnostic_5_patient_composition": {
        "cl15_top_samples": cl15_samples.head(10).to_dict(),
        "purity_row": purity_cl15,
    },
    "evidence_pointer": {
        "cl15_pct_expr_FOXJ1_RFX3_avg": cl15_tf_pct,
        "cl15_pct_expr_dynein_cfap_rsph_avg": cl15_dynein_pct,
        "cl15_margin_median": margin_summary["median"],
        "cl15_margin_frac_below_0p05": margin_summary["frac_margin_below_0p05"],
        "cl15_n_subclusters_res0p5": sub_leiden_summary["n_sub_clusters"],
        "cl15_dominant_sample_fraction":
            float(purity_cl15["dominant_sample_fraction"]),
    },
}

out_path = tbl_path("cl15_ependymal_diagnosis")
out_path.write_text(json.dumps(out, indent=2, default=str))
print(f"\nWrote structured evidence: {out_path}")

adata_nerve.file.close()
adata_label.file.close()
