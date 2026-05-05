"""Test whether nerve-cell cluster abundance is associated with clinical covariates.

Implements §4.3 / Recommended Next Analysis #3 of
``markdowns/next_steps_interpretation.md``:

  1. *Which nerve-cell clusters differ in abundance between diagnosis groups?*
     → Mann-Whitney U on per-sample cluster proportions for each available
     binary categorical covariate (``tissue_type``, ``gender``, ``race``;
     ``primary_diagnosis`` / ``tumor_grade`` / ``prior_malignancy`` are
     constants in the current TCGA-GBM TSV and emit a "skipped" diagnostic
     row).
  2. *Is there a nerve-cell state whose proportion correlates with
     ``age_at_index``?*
     → Spearman correlation per cluster against the per-sample
     ``age_at_index``. The current TSV records "unknown" for every sample,
     so this branch emits a "skipped" diagnostic row but the code path is
     ready when the clinical TSV is enriched.

Outputs
-------
``results/tables/nerve_clinical_association.csv``
    One row per (covariate, cluster) — full results plus diagnostic rows.

``results/figures/nerve_clinical_pvalue_heatmap.png``
    -log10(padj) heatmap, 36 clusters × testable covariates.

``results/figures/nerve_clinical_boxplots.png``
    Per-sample proportion boxplots split by group, for clusters that are
    significant (padj < 0.10) in any covariate; falls back to the top-12
    smallest-padj clusters when nothing crosses 0.10.
"""

import os
import sys
import warnings

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

# Covariate config -------------------------------------------------------------
CATEGORICAL_COVARIATES = ["tissue_type", "gender", "race", "primary_diagnosis",
                          "tumor_grade", "prior_malignancy"]
CONTINUOUS_COVARIATES = ["age_at_index"]
SIG_THRESHOLD = 0.10
TOP_N_FALLBACK = 12

log = snakemake.log[0]  # type: ignore[name-defined]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)  # type: ignore[name-defined]


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Return BH-adjusted p-values; NaNs in input remain NaN in output."""
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan)
    mask = ~np.isnan(pvals)
    p = pvals[mask]
    n = p.size
    if n == 0:
        return out
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out_masked = np.empty(n)
    out_masked[order] = adj
    out[mask] = out_masked
    return out


def coerce_continuous(s: pd.Series) -> pd.Series:
    """Cast a Series of strings (with 'unknown'/'Not Reported') to numeric."""
    if s.dtype.name == "category":
        s = s.astype(object)
    return pd.to_numeric(s, errors="coerce")


# --- Load AnnData and clinical TSV --------------------------------------------
log_transformation(
    log,
    "nerve_clinical_association",
    f"Loading {snakemake.input.h5ad} and {snakemake.input.clinical}",  # type: ignore[name-defined]
)
adata = ad.read_h5ad(snakemake.input.h5ad)  # type: ignore[name-defined]
clinical = pd.read_csv(snakemake.input.clinical, sep="\t")  # type: ignore[name-defined]

# Per-sample covariate frame (one row per sample_id). Prefer the values
# already attached to ``adata.obs`` (set upstream by ``nerve_cell_subset``),
# fall back to the TSV if needed.
sample_covars = (
    adata.obs.groupby("sample_id", observed=True)[
        [c for c in CATEGORICAL_COVARIATES + CONTINUOUS_COVARIATES if c in adata.obs.columns]
    ]
    .first()
    .copy()
)
# Pull missing columns from the TSV (joined on sample_id == file_uuid).
missing = [c for c in CATEGORICAL_COVARIATES + CONTINUOUS_COVARIATES if c not in sample_covars.columns]
if missing and "file_uuid" in clinical.columns:
    extra = clinical.set_index("file_uuid")[
        [c for c in missing if c in clinical.columns]
    ]
    sample_covars = sample_covars.join(extra, how="left")

# Coerce continuous covariates to numeric (handles "unknown" → NaN).
for c in CONTINUOUS_COVARIATES:
    if c in sample_covars.columns:
        sample_covars[c] = coerce_continuous(sample_covars[c])

# --- Per-sample cluster proportions -------------------------------------------
counts = (
    adata.obs.groupby(["sample_id", "nerve_leiden"], observed=True)
    .size()
    .unstack(fill_value=0)
)
# Defensive: cells_per_sample sums must be > 0 here (nerve_cells.h5ad already
# filters); guard anyway.
totals = counts.sum(axis=1).replace(0, np.nan)
proportions = counts.div(totals, axis=0).fillna(0.0)
clusters = list(proportions.columns)
log_transformation(
    log,
    "nerve_clinical_association",
    f"Computed proportions for {len(proportions)} samples × {len(clusters)} clusters",
)

# Align covariates to the proportions index.
sample_covars = sample_covars.reindex(proportions.index)

# --- Run tests ----------------------------------------------------------------
rows: list[dict] = []
testable_cats: list[str] = []

for cov in CATEGORICAL_COVARIATES:
    if cov not in sample_covars.columns:
        rows.append({
            "covariate": cov, "test": "n/a", "n_groups": 0,
            "group_labels": "", "n_samples": 0, "cluster": "",
            "statistic": np.nan, "pvalue": np.nan, "padj": np.nan,
            "note": "covariate not present in inputs — skipped",
        })
        continue

    series = sample_covars[cov].astype(object)
    # Treat sentinel "unknown" / "Not Reported" / NaN as missing.
    sentinel = {"unknown", "Unknown", "Not Reported", "not reported", "", "nan"}
    valid_mask = series.notna() & ~series.astype(str).isin(sentinel)
    series_clean = series[valid_mask]
    groups = series_clean.dropna().unique().tolist()

    if len(groups) < 2:
        rows.append({
            "covariate": cov, "test": "skipped", "n_groups": len(groups),
            "group_labels": ",".join(map(str, groups)), "n_samples": int(valid_mask.sum()),
            "cluster": "", "statistic": np.nan, "pvalue": np.nan, "padj": np.nan,
            "note": f"only {len(groups)} usable group(s) — insufficient variability",
        })
        log_transformation(
            log, "nerve_clinical_association",
            f"[FAIR-ALERT] {cov}: only {len(groups)} group(s) — skipped",
            status="WARNING",
        )
        continue

    if len(groups) > 2:
        test_name = "kruskal"
    else:
        test_name = "mannwhitneyu"
    testable_cats.append(cov)

    raw_pvals: list[float] = []
    raw_stats: list[float] = []
    for cluster in clusters:
        per_group = [
            proportions.loc[series_clean[series_clean == g].index, cluster].values
            for g in groups
        ]
        try:
            if test_name == "mannwhitneyu":
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    stat, pval = stats.mannwhitneyu(
                        per_group[0], per_group[1], alternative="two-sided"
                    )
            else:
                stat, pval = stats.kruskal(*per_group)
        except ValueError:
            stat, pval = np.nan, np.nan
        raw_stats.append(float(stat) if not np.isnan(stat) else np.nan)
        raw_pvals.append(float(pval) if not np.isnan(pval) else np.nan)

    padj = benjamini_hochberg(np.array(raw_pvals))
    for cluster, s, p, pa in zip(clusters, raw_stats, raw_pvals, padj):
        rows.append({
            "covariate": cov, "test": test_name, "n_groups": len(groups),
            "group_labels": ",".join(map(str, groups)), "n_samples": int(valid_mask.sum()),
            "cluster": str(cluster), "statistic": s, "pvalue": p, "padj": pa,
            "note": "",
        })

testable_cont: list[str] = []
for cov in CONTINUOUS_COVARIATES:
    if cov not in sample_covars.columns:
        rows.append({
            "covariate": cov, "test": "n/a", "n_groups": 0,
            "group_labels": "", "n_samples": 0, "cluster": "",
            "statistic": np.nan, "pvalue": np.nan, "padj": np.nan,
            "note": "covariate not present in inputs — skipped",
        })
        continue
    series = sample_covars[cov]
    valid_mask = series.notna()
    n_valid = int(valid_mask.sum())
    if n_valid < 3:
        rows.append({
            "covariate": cov, "test": "skipped", "n_groups": 0,
            "group_labels": "", "n_samples": n_valid, "cluster": "",
            "statistic": np.nan, "pvalue": np.nan, "padj": np.nan,
            "note": f"only {n_valid} numeric values available — insufficient for Spearman",
        })
        log_transformation(
            log, "nerve_clinical_association",
            f"[FAIR-ALERT] {cov}: only {n_valid} numeric values — Spearman skipped",
            status="WARNING",
        )
        continue
    testable_cont.append(cov)
    raw_pvals, raw_stats = [], []
    x = series[valid_mask].astype(float).values
    for cluster in clusters:
        y = proportions.loc[series[valid_mask].index, cluster].values
        try:
            res = stats.spearmanr(x, y)
            raw_stats.append(float(res.statistic))
            raw_pvals.append(float(res.pvalue))
        except ValueError:
            raw_stats.append(np.nan)
            raw_pvals.append(np.nan)
    padj = benjamini_hochberg(np.array(raw_pvals))
    for cluster, s, p, pa in zip(clusters, raw_stats, raw_pvals, padj):
        rows.append({
            "covariate": cov, "test": "spearman", "n_groups": 0,
            "group_labels": "", "n_samples": n_valid, "cluster": str(cluster),
            "statistic": s, "pvalue": p, "padj": pa, "note": "",
        })

stats_df = pd.DataFrame(rows)
stats_df.to_csv(snakemake.output.stats, index=False)  # type: ignore[name-defined]
log_transformation(
    log, "nerve_clinical_association",
    f"Wrote {len(stats_df)} stats rows; testable categorical={testable_cats}; "
    f"testable continuous={testable_cont}",
)

# --- Heatmap of -log10(padj) across testable covariates -----------------------
testable_all = testable_cats + testable_cont
if testable_all:
    pivot = (
        stats_df[stats_df["covariate"].isin(testable_all)]
        .pivot_table(index="cluster", columns="covariate", values="padj")
        .reindex([str(c) for c in clusters])
    )

    def _cluster_sort_key(c: str) -> tuple[int, str]:
        try:
            return (0, f"{int(c):06d}")
        except (TypeError, ValueError):
            return (1, str(c))

    pivot = pivot.loc[sorted(pivot.index, key=_cluster_sort_key)]
    neglog = -np.log10(pivot.clip(lower=1e-300))
    fig, ax = plt.subplots(
        figsize=(max(4, 1.4 * len(testable_all)), max(6, 0.25 * len(pivot)))
    )
    sns.heatmap(
        neglog, ax=ax, cmap="rocket_r", linewidths=0.3, linecolor="white",
        cbar_kws={"label": "-log10(padj)"}, vmin=0,
    )
    ax.set_title("Nerve-cell cluster abundance vs. clinical covariates\n(BH-adjusted)")
    ax.set_xlabel("Covariate")
    ax.set_ylabel("Nerve-cell cluster (leiden)")
    fig.tight_layout()
    fig.savefig(snakemake.output.heatmap, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
    plt.close(fig)
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(
        0.5, 0.5,
        "[FAIR-ALERT] No testable clinical covariates available",
        ha="center", va="center", transform=ax.transAxes, fontsize=11,
    )
    ax.set_axis_off()
    fig.savefig(snakemake.output.heatmap, dpi=100, bbox_inches="tight")  # type: ignore[name-defined]
    plt.close(fig)

# --- Boxplots for clusters with signal ----------------------------------------
sig_clusters: list[str] = []
if testable_cats:
    cat_results = stats_df[stats_df["covariate"].isin(testable_cats)]
    sig = cat_results[cat_results["padj"] < SIG_THRESHOLD]["cluster"].unique().tolist()
    if sig:
        sig_clusters = sig
    else:
        # Fallback: top N smallest padj across testable categorical covariates.
        ranked = (
            cat_results.dropna(subset=["padj"])
            .sort_values("padj")
            .drop_duplicates("cluster")
            .head(TOP_N_FALLBACK)["cluster"]
            .tolist()
        )
        sig_clusters = ranked

if testable_cats and sig_clusters:
    nrows = len(sig_clusters)
    ncols = len(testable_cats)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(max(3.0, 2.6 * ncols), max(2.5, 2.0 * nrows)),
        squeeze=False,
    )
    for i, cluster in enumerate(sig_clusters):
        for j, cov in enumerate(testable_cats):
            ax = axes[i, j]
            series = sample_covars[cov].astype(object)
            sentinel = {"unknown", "Unknown", "Not Reported", "not reported", "", "nan"}
            valid = series.notna() & ~series.astype(str).isin(sentinel)
            df_plot = pd.DataFrame({
                "group": series[valid].astype(str).values,
                "proportion": proportions.loc[series[valid].index, cluster].values,
            })
            sns.boxplot(data=df_plot, x="group", y="proportion", ax=ax,
                        color="lightgray", showfliers=False)
            sns.stripplot(data=df_plot, x="group", y="proportion", ax=ax,
                          color="black", size=4, jitter=True, alpha=0.7)
            row_match = stats_df[
                (stats_df["covariate"] == cov) & (stats_df["cluster"] == str(cluster))
            ]
            padj_str = (
                f"padj={row_match['padj'].iloc[0]:.2g}"
                if len(row_match) and not np.isnan(row_match["padj"].iloc[0])
                else "padj=n/a"
            )
            if i == 0:
                ax.set_title(f"{cov}\n{padj_str}", fontsize=9)
            else:
                ax.set_title(padj_str, fontsize=9)
            ax.set_xlabel("")
            ax.set_ylabel(f"c{cluster}" if j == 0 else "")
            ax.tick_params(axis="x", labelsize=8, rotation=30)
            ax.tick_params(axis="y", labelsize=8)
    fig.suptitle(
        f"Per-sample nerve-cluster proportions by covariate "
        f"({'sig (padj<{:.2f})'.format(SIG_THRESHOLD) if any(stats_df['padj'].dropna() < SIG_THRESHOLD) else f'top-{TOP_N_FALLBACK} ranked by padj'})",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(snakemake.output.boxplots, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
    plt.close(fig)
else:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(
        0.5, 0.5,
        "[FAIR-ALERT] No testable categorical covariates — boxplot skipped",
        ha="center", va="center", transform=ax.transAxes, fontsize=11,
    )
    ax.set_axis_off()
    fig.savefig(snakemake.output.boxplots, dpi=100, bbox_inches="tight")  # type: ignore[name-defined]
    plt.close(fig)

log_transformation(
    log, "nerve_clinical_association",
    f"Boxplots: rendered {len(sig_clusters)} clusters × {len(testable_cats)} covariates",
)

# --- FAIR provenance -----------------------------------------------------------
verify_artifact(snakemake.output.stats)  # type: ignore[name-defined]
prov = stamp_artifact(
    output_path=snakemake.output.stats,  # type: ignore[name-defined]
    rule_name="nerve_clinical_association",
    input_paths=[snakemake.input.h5ad, snakemake.input.clinical],  # type: ignore[name-defined]
    tool_versions={
        "scanpy": ad.__version__,  # anndata version stand-in for the env
        "anndata": ad.__version__,
        "pandas": pd.__version__,
        "scipy": stats.__name__ + "@" + getattr(__import__("scipy"), "__version__", "?"),
    },
    parameters={
        "n_clusters": int(len(clusters)),
        "n_samples": int(proportions.shape[0]),
        "categorical_covariates_tested": testable_cats,
        "categorical_covariates_skipped": [
            c for c in CATEGORICAL_COVARIATES if c not in testable_cats
        ],
        "continuous_covariates_tested": testable_cont,
        "continuous_covariates_skipped": [
            c for c in CONTINUOUS_COVARIATES if c not in testable_cont
        ],
        "sig_threshold": SIG_THRESHOLD,
        "top_n_fallback": TOP_N_FALLBACK,
        "n_sig_clusters": int(len(sig_clusters)),
    },
    description="Per-cluster proportion vs. clinical covariate association tests",
    ontology_operation="operation:2238",  # EDAM: Statistical inference
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]

log_transformation(
    log, "nerve_clinical_association",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.stats,  # type: ignore[name-defined]
        snakemake.output.heatmap,  # type: ignore[name-defined]
        snakemake.output.boxplots,  # type: ignore[name-defined]
        snakemake.output.provenance,  # type: ignore[name-defined]
    ],
)
