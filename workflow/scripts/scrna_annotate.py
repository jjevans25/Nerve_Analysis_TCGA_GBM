"""Cell-type annotation using scVI latent clustering and canonical marker gene scoring."""

import os
import sys

# numba's OpenMP pool deadlocks/SIGSEGVs against torch's libomp once both are
# loaded. Must be set BEFORE anything imports scanpy (which pulls numba via
# umap-learn). See markdowns/ and the 2026-08-02 run notes.
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, "workflow/scripts")
from counts_utils import (
    denormalize_log1p_inplace,
    is_log1p_scale,
    library_sizes,
    normalize_log1p_inplace,
)
from fair_utils import H5AD_COMPRESSION, log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

log_transformation(log, "scrna_annotate", f"Loading {snakemake.input.latent_h5ad}")
adata = ad.read_h5ad(snakemake.input.latent_h5ad)

# --- Neighbors + UMAP from scVI latent representation -----------------------
sc.settings.seed = snakemake.params.random_seed
sc.pp.neighbors(adata, use_rep="X_scVI", random_state=snakemake.params.random_seed)
sc.tl.umap(adata, random_state=snakemake.params.random_seed)

log_transformation(log, "scrna_annotate",
    f"Computed neighbors and UMAP from X_scVI ({adata.obsm['X_scVI'].shape[1]} dims)")

# --- Leiden clustering -------------------------------------------------------
res = snakemake.params.leiden_resolution
sc.tl.leiden(adata, resolution=res, random_state=snakemake.params.random_seed,
             key_added="leiden", flavor="igraph", directed=False, n_iterations=2)
n_clusters = adata.obs["leiden"].nunique()
log_transformation(log, "scrna_annotate",
    f"Leiden clustering (resolution={res}) → {n_clusters} clusters")

# --- Canonical marker gene sets for scoring ----------------------------------
# Two sources, merged. `nerve_cells.markers` holds the neural panels and is also
# read by scrna_qc (is_nerve_marker) across every per-sample job, so it is left
# untouched here — editing it would invalidate 170 QC artifacts per arm and
# cascade into the scVI train. Non-neural panels, and any annotation-only
# override of a neural panel, live in `annotation_markers` and win on conflict.
MARKER_SETS: dict[str, list[str]] = {
    label: list(genes)
    for label, genes in snakemake.params.markers.items()
}
_overrides: list[str] = []
for label, genes in (snakemake.params.annotation_markers or {}).items():
    if label in MARKER_SETS and list(genes) != MARKER_SETS[label]:
        _overrides.append(label)
    MARKER_SETS[label] = list(genes)

log_transformation(log, "scrna_annotate",
    f"Marker panels: {len(MARKER_SETS)} total "
    f"({len(snakemake.params.markers)} neural from nerve_cells.markers, "
    f"{len(snakemake.params.annotation_markers or {})} from annotation_markers"
    + (f"; annotation-only overrides of neural panels: {_overrides}" if _overrides else "")
    + ")")

# --- Build HGNC-symbol → Ensembl-ID lookup ----------------------------------
# Markers above are HGNC symbols; adata.var_names are versioned Ensembl IDs.
# The loom_to_h5ad rule attaches gene_symbol via the MyGene.info cache, but
# ad.concat(merge="same") in scrna_integration silently drops the column when
# NaN handling diverges across samples. Re-join from the same cache here.
if "gene_symbol" not in adata.var.columns:
    log_transformation(log, "scrna_annotate",
        "gene_symbol missing from integrated var — re-joining from MyGene cache",
        status="WARNING")
    _symbol_map = pd.read_csv(snakemake.input.symbol_map, sep="\t", dtype=str)
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
log_transformation(log, "scrna_annotate",
    f"Built symbol→Ensembl lookup: {len(symbol_to_ensembl)} unique symbols "
    f"({adata.var['gene_symbol'].notna().sum()}/{adata.n_vars} genes mapped)")

# --- D1 fix: score on a log1p scale, not on raw UMI counts -------------------
# sc.tl.score_genes assumes log-normalized input and does not check. The Census
# cohorts carry raw integer UMIs in .X (correct for scVI upstream), so scoring
# them directly compared genes on wildly different absolute scales: panel
# mean_confidence spanned 0.13-24.2 here versus 0.11-0.73 in the reference arm.
# See markdowns/blocker_census_annotation_scoring.md.
#
# A normalized *copy* is not affordable (16.5 GB on top of 16.5 GB, on a 36 GB
# machine), so .X is normalized in place, scored, and restored to counts before
# anything is written. Downstream consumers all require raw counts.
if not sp.isspmatrix_csr(adata.X):
    adata.X = sp.csr_matrix(adata.X)

_was_log1p = is_log1p_scale(adata.X)
_x_max_before = float(adata.X.max())
_x_sum_before = float(adata.X.data.sum())
_lib = library_sizes(adata.X)

if _was_log1p:
    log_transformation(log, "scrna_annotate",
        f"Input .X is already on a log1p scale (max={_x_max_before:.2f}) — scoring as-is",
        status="WARNING")
else:
    normalize_log1p_inplace(adata.X, _lib, target_sum=1e4)
    if not is_log1p_scale(adata.X):
        raise RuntimeError(
            f"[FAIR-ALERT] .X still exceeds the log1p plausibility bound after "
            f"normalization (max={float(adata.X.max()):.2f}). Refusing to score — "
            "this is the exact silent-wrong-answer mode that invalidated the "
            "Census annotation. See markdowns/blocker_census_annotation_scoring.md."
        )
    log_transformation(log, "scrna_annotate",
        f"Normalized .X in place for scoring: raw max {_x_max_before:.0f} → "
        f"log1p max {float(adata.X.max()):.2f} (target_sum=1e4)")

scored_sets: list[str] = []
for label, genes in MARKER_SETS.items():
    available_ensembl = [symbol_to_ensembl[s] for s in genes if s in symbol_to_ensembl]
    available_ensembl = [e for e in available_ensembl if e in adata.var_names]
    unmapped = [s for s in genes if s not in symbol_to_ensembl]
    if len(available_ensembl) >= 1:
        sc.tl.score_genes(adata, gene_list=available_ensembl, score_name=f"score_{label}",
                          random_state=snakemake.params.random_seed)
        scored_sets.append(label)
        log_transformation(log, "scrna_annotate",
            f"Scored '{label}': {len(available_ensembl)}/{len(genes)} markers mapped"
            + (f"; unmapped: {unmapped}" if unmapped else ""))
    else:
        log_transformation(log, "scrna_annotate",
            f"WARNING: no markers for '{label}' mapped to dataset (symbols tried: {genes})",
            status="WARNING")

log_transformation(log, "scrna_annotate",
    f"Scored {len(scored_sets)} cell-type gene sets: {scored_sets}")

# Restore raw counts before anything else touches .X.
if not _was_log1p:
    denormalize_log1p_inplace(adata.X, _lib, target_sum=1e4)
    _x_max_after = float(adata.X.max())
    _x_sum_after = float(adata.X.data.sum())
    _drift = abs(_x_sum_after - _x_sum_before) / max(_x_sum_before, 1.0)
    if _drift > 1e-4 or abs(_x_max_after - _x_max_before) > 1.0:
        raise RuntimeError(
            f"[FAIR-ALERT] counts round-trip did not restore .X: sum drift {_drift:.2e}, "
            f"max {_x_max_before:.0f} → {_x_max_after:.0f}. Refusing to write an "
            "artifact whose .X is neither counts nor log1p."
        )
    log_transformation(log, "scrna_annotate",
        f"Restored raw counts in .X (max={_x_max_after:.0f}, relative sum drift {_drift:.2e})")

# --- D2 fix: z-score panels before argmax, and allow 'ambiguous' -------------
# Raw score_genes values are not comparable across panels: a panel built from
# highly-expressed genes wins the argmax on scale alone, regardless of biology.
# That is how 322k cells became "astrocyte" against 304 in the author annotation.
# Z-scoring each panel across cells puts every panel on the same footing, so the
# argmax compares *relative enrichment* rather than absolute expression.
score_cols = [f"score_{s}" for s in scored_sets]
z_cols = [f"z_{s}" for s in scored_sets]

for label in scored_sets:
    v = adata.obs[f"score_{label}"].to_numpy(dtype=np.float64)
    sd = v.std()
    adata.obs[f"z_{label}"] = (
        ((v - v.mean()) / sd) if sd > 0 else np.zeros_like(v)
    ).astype(np.float32)

margin_floor = float(snakemake.params.ambiguous_margin)
panel_compartment: dict[str, str] = dict(snakemake.params.panel_compartment or {})
label_map: dict[str, str] = {}
cluster_report: list[dict[str, object]] = []

if z_cols:
    cluster_means = (
        adata.obs[["leiden"] + z_cols]
        .groupby("leiden", observed=True)[z_cols]
        .mean()
    )
    for cluster in cluster_means.index:
        row = cluster_means.loc[cluster].sort_values(ascending=False)
        top, second = row.index[0], (row.index[1] if len(row) > 1 else None)
        margin = float(row.iloc[0] - row.iloc[1]) if len(row) > 1 else float("inf")
        top_panel = top.replace("z_", "")
        second_panel = second.replace("z_", "") if second else ""
        # A cluster whose best and second-best panels are indistinguishable is
        # not evidence for the winner — it is evidence that no panel fits. But
        # only a tie that crosses compartments threatens compartment integrity;
        # a tie between two immune lineages lands in the same compartment either
        # way, and discarding it would lose a correctly placed cell.
        top_comp = panel_compartment.get(top_panel, top_panel)
        second_comp = panel_compartment.get(second_panel, second_panel)
        cross_compartment = bool(second_panel) and top_comp != second_comp
        called = "ambiguous" if (margin < margin_floor and cross_compartment) else top_panel
        label_map[cluster] = called
        cluster_report.append({
            "leiden":            str(cluster),
            "n_cells":           int((adata.obs["leiden"] == cluster).sum()),
            "label":             called,
            "top_panel":         top_panel,
            "top_compartment":   top_comp,
            "top_z":             float(row.iloc[0]),
            "second_panel":      second_panel,
            "second_compartment": second_comp,
            "second_z":          float(row.iloc[1]) if len(row) > 1 else float("nan"),
            "margin":            margin,
            "cross_compartment_tie": cross_compartment and margin < margin_floor,
        })
else:
    log_transformation(log, "scrna_annotate",
        "WARNING: no marker gene sets scored — all clusters labeled 'unscored'", status="WARNING")
    for cluster in adata.obs["leiden"].unique():
        label_map[cluster] = "unscored"

adata.obs["cell_type_predicted"] = adata.obs["leiden"].map(label_map).astype(str)

n_ambiguous = int((adata.obs["cell_type_predicted"] == "ambiguous").sum())
if n_ambiguous:
    log_transformation(log, "scrna_annotate",
        f"{n_ambiguous} cells ({100 * n_ambiguous / adata.n_obs:.1f}%) in clusters whose "
        f"top-minus-second z margin fell below {margin_floor} → labeled 'ambiguous'",
        status="WARNING")

if cluster_report:
    pd.DataFrame(cluster_report).sort_values("n_cells", ascending=False).to_csv(
        snakemake.output.cluster_scores, index=False)
else:
    pd.DataFrame(columns=["leiden", "n_cells", "label", "top_panel", "top_z",
                          "second_panel", "second_z", "margin"]).to_csv(
        snakemake.output.cluster_scores, index=False)

# Confidence: per-cell margin between top and second-best z-score.
if len(z_cols) >= 2:
    score_matrix = adata.obs[z_cols].values
    sorted_scores = np.sort(score_matrix, axis=1)
    adata.obs["cell_type_confidence"] = (sorted_scores[:, -1] - sorted_scores[:, -2]).astype(np.float32)
else:
    adata.obs["cell_type_confidence"] = np.float32(1.0)

# Record schema in uns for FAIR interoperability
adata.uns["annotation_schema"] = {
    "method":            "leiden_marker_scoring_zscored",
    "leiden_resolution": res,
    "marker_sets":       {k: v for k, v in MARKER_SETS.items() if k in scored_sets},
    "scoring_scale":     "log1p(normalize_total(1e4)); .X restored to raw counts after scoring",
    "assignment":        "per-cluster argmax over cell-z-scored panel scores",
    "ambiguous_margin":  margin_floor,
    "ontology":          "Allen Brain Atlas cell type taxonomy",
    "edam_operation":    "operation:3432",  # EDAM: Clustering
}
adata.uns["normalization"] = {
    "X":            "raw UMI counts (integer)",
    "scoring_only": "normalize_total(target_sum=1e4) + log1p, applied in place and reverted",
}

log_transformation(log, "scrna_annotate",
    f"Cell type distribution: {adata.obs['cell_type_predicted'].value_counts().to_dict()}")

# --- Annotation summary ------------------------------------------------------
summary = (
    adata.obs.groupby("cell_type_predicted")
    .agg(n_cells=("cell_type_predicted", "count"),
         mean_confidence=("cell_type_confidence", "mean"))
    .reset_index()
    .sort_values("n_cells", ascending=False)
)
summary.to_csv(snakemake.output.summary, index=False)

# --- Write output ------------------------------------------------------------
adata.write_h5ad(snakemake.output.h5ad, compression=H5AD_COMPRESSION)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="scrna_annotate",
    input_paths=[snakemake.input.latent_h5ad],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "leiden_resolution": res,
        "n_clusters":        n_clusters,
        "scored_sets":       scored_sets,
        "n_cells":           adata.n_obs,
        "scoring_scale":     "log1p_1e4_reverted",
        "ambiguous_margin":  margin_floor,
        "n_ambiguous_cells": n_ambiguous,
        "marker_sets":       {k: v for k, v in MARKER_SETS.items() if k in scored_sets},
    },
    description="Cell-type annotated AnnData with Leiden clusters, UMAP, and marker scores",
    ontology_operation="operation:3432",
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "scrna_annotate", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.summary,
                                   snakemake.output.cluster_scores])
