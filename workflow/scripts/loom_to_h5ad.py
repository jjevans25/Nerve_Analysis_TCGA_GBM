"""Convert a GDC loom file to AnnData h5ad; assigns batch metadata and emits gene presence report."""

import sys
import warnings
from pathlib import Path

import anndata as ad
import loompy
import numpy as np
import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
sample_id = snakemake.params.sample_id
nerve_markers: dict = snakemake.params.nerve_markers

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
adata.var_names = gene_names
adata.obs_names = [f"{sample_id}_{c}" for c in cell_ids]

# Batch / sample metadata (required for scVI batch correction)
adata.obs["batch"]     = sample_id
adata.obs["sample_id"] = sample_id

for key, vals in extra_ca.items():
    if len(vals) == n_cells:
        adata.obs[key] = vals

# Flag mitochondrial genes for downstream QC
adata.var["mt"] = pd.Index(gene_names).str.upper().str.startswith("MT-")

# Mark gene set membership for nerve-cell markers
all_nerve_markers: list[str] = [
    g for genes in nerve_markers.values() for g in genes
]
adata.var["is_nerve_marker"] = adata.var_names.isin(all_nerve_markers)

log_transformation(log, "loom_to_h5ad",
    f"batch='{sample_id}' assigned; {adata.var['is_nerve_marker'].sum()} nerve markers present")

# --- Gene presence report -----------------------------------------------------
gene_set = set(adata.var_names)
rows = []
for cell_type, markers in nerve_markers.items():
    for gene in markers:
        rows.append({
            "cell_type":  cell_type,
            "gene":       gene,
            "present":    gene in gene_set,
            "sample_id":  sample_id,
        })

presence_df = pd.DataFrame(rows)
presence_df.to_csv(snakemake.output.gene_presence, index=False)

missing = presence_df[~presence_df["present"]]["gene"].tolist()
if missing:
    warnings.warn(
        f"[FAIR-ALERT] {len(missing)} nerve-cell markers absent from {sample_id}: {missing}"
    )
    log_transformation(log, "loom_to_h5ad",
        f"WARNING: {len(missing)} nerve markers missing: {missing}", status="WARNING")

# --- Write h5ad ---------------------------------------------------------------
adata.write_h5ad(snakemake.output.h5ad)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

# --- FAIR provenance ----------------------------------------------------------
prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="loom_to_h5ad",
    input_paths=[snakemake.input.loom],
    tool_versions={"anndata": ad.__version__, "loompy": loompy.__version__},
    parameters={
        "sample_id":           sample_id,
        "n_cells":             n_cells,
        "n_genes":             n_genes,
        "nerve_markers_found": int(adata.var["is_nerve_marker"].sum()),
        "nerve_markers_total": len(all_nerve_markers),
    },
    description="GDC loom converted to AnnData h5ad with batch metadata and nerve-marker annotation",
    ontology_operation="operation:2409",  # EDAM: Format conversion
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "loom_to_h5ad", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.gene_presence])
