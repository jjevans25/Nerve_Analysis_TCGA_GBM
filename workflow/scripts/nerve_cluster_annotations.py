"""Per-cluster biological labels for the nerve-cell Leiden clusters.

Combines:
  - Dominant ``cell_type_predicted`` from the upstream annotation step.
  - Top 5 marker symbols from ``nerve_cluster_markers.csv`` (Wilcoxon).
  - Canonical nerve-cell module scores via ``scanpy.tl.score_genes``
    (neuron / excitatory_neuron / inhibitory_neuron / opc / oligodendrocyte / astrocyte).

Implements §4.2 of ``markdowns/next_steps_interpretation.md``: produces
``results/tables/nerve_cluster_annotations.csv`` with the columns
``cluster, label, top_markers, interpretation`` plus auxiliary per-module
score columns to support downstream clinical-association work.
"""

import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

LOW_SCORE_THRESHOLD = 0.05
TOP_N_MARKERS = 5

log = snakemake.log[0]  # type: ignore[name-defined]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)  # type: ignore[name-defined]

log_transformation(
    log,
    "nerve_cluster_annotations",
    f"Loading {snakemake.input.h5ad} and {snakemake.input.markers}",  # type: ignore[name-defined]
)
adata = ad.read_h5ad(snakemake.input.h5ad)  # type: ignore[name-defined]
markers_df = pd.read_csv(snakemake.input.markers)  # type: ignore[name-defined]

# --- Build HGNC-symbol ↔ Ensembl-ID lookup (mirrors nerve_cell_heterogeneity:36-69) ---
if "gene_symbol" not in adata.var.columns:
    log_transformation(
        log,
        "nerve_cluster_annotations",
        "gene_symbol missing from var — re-joining from MyGene cache",
        status="WARNING",
    )
    _symbol_map = pd.read_csv(snakemake.input.symbol_map, sep="\t", dtype=str)  # type: ignore[name-defined]
    if "ensembl_id" not in adata.var.columns:
        adata.var["ensembl_id"] = adata.var_names
    _var_joined = adata.var.merge(
        _symbol_map[["ensembl_id", "gene_symbol", "chromosome"]],
        on="ensembl_id",
        how="left",
    )
    _var_joined.index = adata.var.index
    adata.var = _var_joined

symbol_to_ensembl: dict[str, str] = (
    adata.var.dropna(subset=["gene_symbol"])
    .reset_index()
    .drop_duplicates(subset="gene_symbol", keep="first")
    .set_index("gene_symbol")["index"]
    .to_dict()
)
ensembl_to_symbol: dict[str, str] = adata.var.dropna(subset=["gene_symbol"])[
    "gene_symbol"
].to_dict()
log_transformation(
    log,
    "nerve_cluster_annotations",
    f"Built symbol↔Ensembl lookup: {len(symbol_to_ensembl)} unique symbols "
    f"({adata.var['gene_symbol'].notna().sum()}/{adata.n_vars} genes mapped)",
)

# --- Empty-data guard (mirrors nerve_cell_heterogeneity:78-113) ----------------
n_clusters = adata.obs["nerve_leiden"].nunique() if "nerve_leiden" in adata.obs else 0
if adata.n_obs == 0 or n_clusters == 0:
    log_transformation(
        log,
        "nerve_cluster_annotations",
        "[FAIR-ALERT] Empty nerve-cell AnnData — writing placeholder annotation CSV.",
        status="WARNING",
    )
    pd.DataFrame(
        columns=["cluster", "label", "top_markers", "interpretation"]
    ).to_csv(snakemake.output.annotations, index=False)  # type: ignore[name-defined]
    prov_empty = stamp_artifact(
        output_path=snakemake.output.annotations,  # type: ignore[name-defined]
        rule_name="nerve_cluster_annotations",
        input_paths=[snakemake.input.h5ad, snakemake.input.markers],  # type: ignore[name-defined]
        tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
        parameters={"n_clusters": 0, "n_cells": 0, "note": "no nerve cells found"},
        description="Placeholder: no nerve cells available for annotation",
        ontology_operation="operation:3431",
    )
    write_provenance(prov_empty, snakemake.output.provenance)  # type: ignore[name-defined]
    raise SystemExit(0)

log_transformation(
    log,
    "nerve_cluster_annotations",
    f"{adata.n_obs} cells across {n_clusters} nerve-cell clusters",
)

# --- Step A: dominant cell_type_predicted per cluster --------------------------
def _mode_or_na(s: pd.Series) -> str:
    counts = s.dropna().value_counts()
    return counts.index[0] if len(counts) else "unknown"


cluster_size = adata.obs.groupby("nerve_leiden", observed=True).size().rename("n_cells")
dominant_ct = (
    adata.obs.groupby("nerve_leiden", observed=True)["cell_type_predicted"]
    .agg(_mode_or_na)
    .rename("dominant_cell_type")
)
dominant_pct = (
    adata.obs.groupby("nerve_leiden", observed=True)["cell_type_predicted"]
    .agg(lambda s: (s == _mode_or_na(s)).mean() * 100.0 if len(s.dropna()) else 0.0)
    .rename("dominant_pct")
)

# --- Step B: top 5 marker symbols per cluster ----------------------------------
def _pick_symbol(row: pd.Series) -> str:
    sym = row.get("gene_symbol")
    if isinstance(sym, str) and sym.strip():
        return sym
    return str(row["names"])


if len(markers_df) and "gene_symbol" not in markers_df.columns:
    markers_df["gene_symbol"] = markers_df["names"].map(ensembl_to_symbol)

top_markers_by_cluster: dict[str, list[str]] = {}
for cluster, grp in markers_df.groupby("cluster"):
    ranked = grp.sort_values("scores", ascending=False).head(TOP_N_MARKERS)
    top_markers_by_cluster[str(cluster)] = [
        _pick_symbol(r) for _, r in ranked.iterrows()
    ]

# --- Step C: canonical nerve-cell module scores --------------------------------
# Use the same normalization fallback as nerve_cell_heterogeneity:118-126.
if adata.raw is not None:
    adata_de = adata.raw.to_adata()
    adata_de.obs = adata.obs.copy()
else:
    adata_de = adata.copy()
    if not ("log1p" in adata_de.uns or adata_de.X.max() < 20):
        sc.pp.normalize_total(adata_de, target_sum=1e4)
        sc.pp.log1p(adata_de)

# Ensure var_names of adata_de match Ensembl IDs (raw.to_adata preserves them).
module_score_columns: list[str] = []
modules_used: dict[str, list[str]] = {}
for module, symbols in snakemake.params.markers.items():  # type: ignore[name-defined]
    ensembl_list = [
        symbol_to_ensembl[s]
        for s in symbols
        if s in symbol_to_ensembl and symbol_to_ensembl[s] in adata_de.var_names
    ]
    score_col = f"score_{module}"
    if not ensembl_list:
        log_transformation(
            log,
            "nerve_cluster_annotations",
            f"WARNING: no genes from module '{module}' present in dataset — score set to NaN",
            status="WARNING",
        )
        adata_de.obs[score_col] = np.nan
    else:
        sc.tl.score_genes(
            adata_de,
            gene_list=ensembl_list,
            score_name=score_col,
            random_state=snakemake.params.random_seed,  # type: ignore[name-defined]
            use_raw=False,
        )
        modules_used[module] = ensembl_list
    module_score_columns.append(score_col)

cluster_scores = (
    adata_de.obs.groupby("nerve_leiden", observed=True)[module_score_columns]
    .mean()
)
log_transformation(
    log,
    "nerve_cluster_annotations",
    f"Computed module scores for {len(modules_used)}/{len(snakemake.params.markers)} modules",  # type: ignore[name-defined]
)

# --- Step D: assemble per-cluster annotation rows ------------------------------
score_only = cluster_scores  # already cluster × module
best_module = score_only.idxmax(axis=1).str.replace("score_", "", regex=False)
best_score = score_only.max(axis=1)

rows: list[dict] = []
for cluster_label, dominant in dominant_ct.items():
    cluster_str = str(cluster_label)
    top5 = top_markers_by_cluster.get(cluster_str, [])
    pct = float(dominant_pct.loc[cluster_label])
    n = int(cluster_size.loc[cluster_label])
    bm = str(best_module.loc[cluster_label]) if cluster_label in best_module.index else "n/a"
    bs = float(best_score.loc[cluster_label]) if cluster_label in best_score.index else float("nan")

    label = f"c{cluster_str} | {dominant} | {bm}"
    top_str = ", ".join(top5) if top5 else "(no significant markers)"
    interp = (
        f"Dominant {dominant} ({pct:.0f}% of {n} cells); "
        f"top markers {top_str}; "
        f"highest canonical score = {bm} ({bs:.2f})."
    )
    if not np.isnan(bs) and bs < LOW_SCORE_THRESHOLD:
        interp += " Low canonical-marker signal — possible mixed/uncharacterized state."

    row = {
        "cluster": cluster_str,
        "label": label,
        "top_markers": top_str,
        "interpretation": interp,
    }
    for col in module_score_columns:
        row[col] = (
            float(cluster_scores.loc[cluster_label, col])
            if cluster_label in cluster_scores.index
            else float("nan")
        )
    rows.append(row)

annotations_df = pd.DataFrame(rows)


def _cluster_sort_key(c: str) -> tuple[int, str]:
    try:
        return (0, f"{int(c):06d}")
    except (TypeError, ValueError):
        return (1, str(c))


annotations_df = annotations_df.sort_values(
    by="cluster", key=lambda s: s.map(_cluster_sort_key)
).reset_index(drop=True)

annotations_df.to_csv(snakemake.output.annotations, index=False)  # type: ignore[name-defined]
log_transformation(
    log,
    "nerve_cluster_annotations",
    f"Wrote {len(annotations_df)} cluster annotations to {snakemake.output.annotations}",  # type: ignore[name-defined]
)

# --- FAIR provenance -----------------------------------------------------------
verify_artifact(snakemake.output.annotations)  # type: ignore[name-defined]
prov = stamp_artifact(
    output_path=snakemake.output.annotations,  # type: ignore[name-defined]
    rule_name="nerve_cluster_annotations",
    input_paths=[snakemake.input.h5ad, snakemake.input.markers, snakemake.input.symbol_map],  # type: ignore[name-defined]
    tool_versions={
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "pandas": pd.__version__,
    },
    parameters={
        "n_clusters": int(n_clusters),
        "n_cells": int(adata.n_obs),
        "top_n_markers": TOP_N_MARKERS,
        "low_score_threshold": LOW_SCORE_THRESHOLD,
        "modules_scored": list(modules_used.keys()),
        "modules_skipped": [
            m for m in snakemake.params.markers.keys() if m not in modules_used  # type: ignore[name-defined]
        ],
    },
    description="Per-cluster biological labels combining dominant cell type, top markers, and canonical module scores",
    ontology_operation="operation:3431",  # EDAM: Cell type annotation
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]

log_transformation(
    log,
    "nerve_cluster_annotations",
    "Complete",
    status="SUCCESS",
    artifact_paths=[snakemake.output.annotations, snakemake.output.provenance],  # type: ignore[name-defined]
)
