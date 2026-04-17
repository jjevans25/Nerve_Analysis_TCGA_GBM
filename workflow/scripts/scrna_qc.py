"""Per-sample QC filtering; logs all thresholds before dropping any cells."""

import sys

import anndata as ad
import scanpy as sc

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
log_transformation(log, "scrna_qc", f"Loading {snakemake.input.h5ad}")

adata = ad.read_h5ad(snakemake.input.h5ad)
n_cells_raw = adata.n_obs

# --- Compute QC metrics before filtering (FAIR: log thresholds first) ---------
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

params = snakemake.params
log_transformation(log, "scrna_qc",
    f"QC thresholds — min_genes={params.min_genes}, max_genes={params.max_genes}, "
    f"min_cells={params.min_cells}, max_pct_mito={params.max_pct_mito}%")

# Save QC metrics before any filtering
adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].to_csv(
    snakemake.output.qc_metrics
)

# --- Filter cells and genes ---------------------------------------------------
sc.pp.filter_cells(adata, min_genes=params.min_genes)
sc.pp.filter_cells(adata, max_genes=params.max_genes)
sc.pp.filter_genes(adata, min_cells=params.min_cells)
adata = adata[adata.obs.pct_counts_mt < params.max_pct_mito].copy()

n_cells_kept = adata.n_obs
log_transformation(log, "scrna_qc",
    f"Filtered: {n_cells_raw} → {n_cells_kept} cells retained "
    f"({n_cells_raw - n_cells_kept} removed)")

# Record normalization step in uns for downstream reference
adata.uns["normalization"] = "raw_counts_pre_qc"

adata.write_h5ad(snakemake.output.h5ad)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

# --- FAIR provenance ----------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="scrna_qc",
    input_paths=[snakemake.input.h5ad],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "min_genes":    params.min_genes,
        "max_genes":    params.max_genes,
        "min_cells":    params.min_cells,
        "max_pct_mito": params.max_pct_mito,
        "cells_raw":    n_cells_raw,
        "cells_kept":   n_cells_kept,
    },
    description="Per-sample QC-filtered AnnData; raw counts preserved in .X",
    ontology_operation="operation:3435",  # EDAM: Gene expression data filtering
)
write_provenance(prov, "provenance")

log_transformation(log, "scrna_qc", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.qc_metrics])
