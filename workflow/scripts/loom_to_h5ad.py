"""Convert a GDC loom file to AnnData h5ad; assigns batch metadata + symbol map.

Nerve-cell marker annotation (is_nerve_marker var-column + gene_presence
report) was moved to scrna_qc in v1.2.0 so panel-config changes do not
trigger sample-level reruns. See markdowns/failing_cluster_diagnosis.md
v1.1.0 supplement.
"""

import sys
from pathlib import Path

import anndata as ad
import loompy
import numpy as np
import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
sample_id = snakemake.params.sample_id

log_transformation(log, "loom_to_h5ad", f"Loading loom: {snakemake.input.loom}")

# --- Load loom ----------------------------------------------------------------
with loompy.connect(str(snakemake.input.loom), mode="r") as ds:
    n_genes, n_cells = ds.shape

    # Transpose: loom is genes×cells; AnnData expects cells×genes
    matrix = ds[:, :].T.astype(np.float32)

    gene_names = [str(g) for g in ds.ra.get("Gene", [f"Gene_{i}" for i in range(n_genes)])]
    cell_ids   = [str(c) for c in ds.ca.get("CellID", [f"Cell_{i}" for i in range(n_cells)])]

    # Collect any extra column attributes as obs metadata
    extra_ca: dict[str, list] = {}
    for key in ds.ca.keys():
        if key != "CellID":
            try:
                extra_ca[key] = list(ds.ca[key])
            except Exception:
                pass

log_transformation(log, "loom_to_h5ad",
    f"Loaded {n_cells} cells × {n_genes} genes from {Path(snakemake.input.loom).name}")

# --- Build AnnData ------------------------------------------------------------
adata = ad.AnnData(X=matrix)
adata.var_names = gene_names  # canonical ID is versioned Ensembl (per config.fair.ontology_gene)
adata.obs_names = [f"{sample_id}_{c}" for c in cell_ids]

# Batch / sample metadata (required for scVI batch correction)
adata.obs["batch"]     = sample_id
adata.obs["sample_id"] = sample_id

for key, vals in extra_ca.items():
    if len(vals) == n_cells:
        adata.obs[key] = vals

# --- Join MyGene.info symbol/chromosome cache --------------------------------
# GDC looms store versioned Ensembl IDs (ENSG00000136492.9). Marker scoring needs
# HGNC symbols, and MT detection needs chromosome. The cache is built once by
# build_gene_symbol_map and contains every Ensembl ID in the loom gene universe.
symbol_map = pd.read_csv(snakemake.input.symbol_map, sep="\t", dtype=str)
adata.var["ensembl_id"]            = adata.var_names
adata.var["ensembl_id_no_version"] = adata.var_names.str.split(".").str[0]

var_joined = adata.var.merge(
    symbol_map[["ensembl_id", "gene_symbol", "chromosome"]],
    on="ensembl_id",
    how="left",
)
var_joined.index = adata.var.index
adata.var = var_joined

n_mapped = adata.var["gene_symbol"].notna().sum()
log_transformation(log, "loom_to_h5ad",
    f"Joined symbol map: {n_mapped}/{adata.n_vars} genes mapped to HGNC symbol "
    f"({100 * n_mapped / adata.n_vars:.1f}%)")

# Flag mitochondrial genes by chromosome (authoritative) — the old MT- prefix
# check never matched because var_names are Ensembl IDs, not symbols.
adata.var["mt"] = adata.var["chromosome"].fillna("").astype(str).str.upper() == "MT"
log_transformation(log, "loom_to_h5ad",
    f"Flagged {int(adata.var['mt'].sum())} mitochondrial genes (chromosome == 'MT')")

log_transformation(log, "loom_to_h5ad", f"batch='{sample_id}' assigned")

# --- Write h5ad ---------------------------------------------------------------
adata.write_h5ad(snakemake.output.h5ad)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

# --- FAIR provenance ----------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="loom_to_h5ad",
    input_paths=[snakemake.input.loom, snakemake.input.symbol_map],
    tool_versions={"anndata": ad.__version__, "loompy": loompy.__version__},
    parameters={
        "sample_id":           sample_id,
        "n_cells":             n_cells,
        "n_genes":             n_genes,
        "n_genes_mapped":      int(n_mapped),
        "n_mt_genes":          int(adata.var["mt"].sum()),
    },
    description="GDC loom converted to AnnData h5ad with batch metadata and MyGene symbol map",
    ontology_operation="operation:2409",  # EDAM: Format conversion
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "loom_to_h5ad", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad])
