"""Per-sample QC filtering; logs all thresholds before dropping any cells.

Also annotates is_nerve_marker on var and writes the per-sample gene_presence
report. Both were moved here from loom_to_h5ad in v1.2.0 so panel-config
changes do not trigger sample-level reruns.
"""

import sys
import warnings

import anndata as ad
import pandas as pd
import scanpy as sc

sys.path.insert(0, "workflow/scripts")
from fair_utils import H5AD_COMPRESSION, log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
sample_id: str = snakemake.params.sample_id
nerve_markers: dict = snakemake.params.markers

log_transformation(log, "scrna_qc", f"Loading {snakemake.input.h5ad}")

adata = ad.read_h5ad(snakemake.input.h5ad)
n_cells_raw = adata.n_obs

# --- Compute QC metrics before filtering (FAIR: log thresholds first) ---------
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

params = snakemake.params
log_transformation(log, "scrna_qc",
    f"QC thresholds — min_genes={params.min_genes}, max_genes={params.max_genes}, "
    f"max_pct_mito={params.max_pct_mito}% (gene filtering deferred to scrna_integration)")

# Save QC metrics before any filtering
adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(
    snakemake.output.qc_metrics
)

# --- Nerve-marker annotation + per-sample gene-presence report ----------------
# gene_presence.csv reflects pre-filter membership (matches v1.1.0 semantics).
# Genes are no longer filtered in this rule; see the note above the cell filters.
if "gene_symbol" in adata.var.columns:
    all_nerve_markers: list[str] = [
        g for genes in nerve_markers.values() for g in genes
    ]
    adata.var["is_nerve_marker"] = (
        adata.var["gene_symbol"].isin(all_nerve_markers).fillna(False)
    )
    n_present = int(adata.var["is_nerve_marker"].sum())
    log_transformation(log, "scrna_qc",
        f"sample={sample_id}: {n_present} nerve-marker genes present pre-filter")

    symbol_set = set(adata.var["gene_symbol"].dropna().tolist())
    rows = []
    for cell_type, markers_list in nerve_markers.items():
        for gene in markers_list:
            rows.append({
                "cell_type": cell_type,
                "gene":      gene,
                "present":   gene in symbol_set,
                "sample_id": sample_id,
            })
    presence_df = pd.DataFrame(rows)
    presence_df.to_csv(snakemake.output.gene_presence, index=False)

    missing = presence_df[~presence_df["present"]]["gene"].tolist()
    if missing:
        warnings.warn(
            f"[FAIR-ALERT] {len(missing)} nerve-cell markers absent from "
            f"{sample_id}: {missing}"
        )
        log_transformation(log, "scrna_qc",
            f"WARNING: {len(missing)} nerve markers missing: {missing}",
            status="WARNING")
else:
    log_transformation(log, "scrna_qc",
        "[FAIR-ALERT] gene_symbol column missing from var; cannot annotate "
        "is_nerve_marker / gene_presence", status="WARNING")
    adata.var["is_nerve_marker"] = False
    pd.DataFrame(columns=["cell_type", "gene", "present", "sample_id"]).to_csv(
        snakemake.output.gene_presence, index=False
    )
    all_nerve_markers = []

# --- Filter cells -------------------------------------------------------------
# Genes are deliberately NOT filtered per sample. Filtering here and then
# intersecting in scrna_integration deletes any gene absent from a single sample
# — that collapsed the 170-sample cohort to 1,519 marker-free genes. Genes are
# filtered once, globally, after concatenation. See
# markdowns/diagnosis_rerun_gene_intersection.md
sc.pp.filter_cells(adata, min_genes=params.min_genes)
sc.pp.filter_cells(adata, max_genes=params.max_genes)
adata = adata[adata.obs.pct_counts_mt < params.max_pct_mito].copy()

n_cells_kept = adata.n_obs
log_transformation(log, "scrna_qc",
    f"Filtered: {n_cells_raw} → {n_cells_kept} cells retained "
    f"({n_cells_raw - n_cells_kept} removed)")

# Record normalization step in uns for downstream reference
adata.uns["normalization"] = "raw_counts_pre_qc"

adata.write_h5ad(snakemake.output.h5ad, compression=H5AD_COMPRESSION)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

# --- FAIR provenance ----------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="scrna_qc",
    input_paths=[snakemake.input.h5ad],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "min_genes":           params.min_genes,
        "max_genes":           params.max_genes,
        "max_pct_mito":        params.max_pct_mito,
        "cells_raw":           n_cells_raw,
        "cells_kept":          n_cells_kept,
        "nerve_markers_found": int(adata.var["is_nerve_marker"].sum()),
        "nerve_markers_total": len(all_nerve_markers),
    },
    description="Per-sample QC-filtered AnnData with is_nerve_marker annotation; raw counts preserved in .X",
    ontology_operation="operation:3435",  # EDAM: Gene expression data filtering
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "scrna_qc", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.qc_metrics])
