"""Generalized per-sample ingest for replication cohorts (Stage A).

Dispatches on the dataset entry's ``source`` and emits one canonical
``data/processed/<dataset>/<sample>.h5ad`` per sample, matching the contract the
loom path (loom_to_h5ad.py) produces for the reference cohort:
  - ``.X`` = raw counts (float32)
  - ``.var`` carries ``ensembl_id`` (canonical index) + ``gene_symbol`` + ``mt``
  - ``.obs`` carries ``batch`` / ``sample_id`` / ``dataset``
so every downstream reused script (scrna_qc, scrna_integration, scrna_annotate,
…) works unchanged.

Currently implements the ``h5ad`` branch (CELLxGENE Census pull), validated by
the 2026-07-16 Stage-A smoke test. Other sources raise NotImplementedError with
the extension point documented, per the plan's dispatch design.
"""

import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
P = snakemake.params

dataset: str = P.dataset
sample_id: str = P.sample_id
source: str = P.source
raw_file: str = P.raw_file
filt: dict = dict(P.filter or {})
sample_key: str = P.sample_key
cap = P.subsample_per_donor          # int or None
seed: int = P.random_seed

os.environ["PYTHONHASHSEED"] = str(seed)


def _ingest_h5ad() -> ad.AnnData:
    """h5ad source (e.g. CELLxGENE Census pull): backed read → filter to this
    sample → optional per-sample cap → canonical var/obs. Never loads the full
    matrix; only this sample's rows come into memory."""
    log_transformation(log, "ingest_dataset",
        f"[{dataset}/{sample_id}] backed read {raw_file} (source=h5ad)")
    backed = ad.read_h5ad(raw_file, backed="r")
    obs = backed.obs

    mask = np.ones(obs.shape[0], dtype=bool)
    for col, val in filt.items():
        mask &= (obs[col].astype(str) == str(val)).to_numpy()
    mask &= (obs[sample_key].astype(str) == str(sample_id)).to_numpy()

    n_sel = int(mask.sum())
    if n_sel == 0:
        raise ValueError(f"[FAIR-ALERT] no cells for {sample_key}={sample_id} "
                         f"under filter {filt} in {raw_file}")
    sub = backed[mask].to_memory()
    log_transformation(log, "ingest_dataset",
        f"[{dataset}/{sample_id}] selected {n_sel} cells")

    if cap and sub.n_obs > int(cap):
        sc.pp.subsample(sub, n_obs=int(cap), random_state=seed)
        log_transformation(log, "ingest_dataset",
            f"[{dataset}/{sample_id}] per-sample cap {cap} applied → {sub.n_obs} cells (seed={seed})")
    return sub


DISPATCH = {"h5ad": _ingest_h5ad}
if source not in DISPATCH:
    raise NotImplementedError(
        f"ingest_dataset source '{source}' not implemented yet. "
        f"Available: {sorted(DISPATCH)}. Extend DISPATCH for gdc_loom/mtx/cellxgene_census."
    )

adata = DISPATCH[source]()

# --- Force raw-count float32 .X ------------------------------------------------
X = adata.X
data = X.data if hasattr(X, "data") else np.asarray(X).ravel()
frac_int = float(np.mean(data[:500000] == np.round(data[:500000]))) if data.size else 1.0
if frac_int < 0.999:
    log_transformation(log, "ingest_dataset",
        f"[FAIR-ALERT] .X does not look like raw counts (frac_integer={frac_int:.3f})",
        status="WARNING")
if X.dtype != np.float32:
    adata.X = X.astype(np.float32)

# --- Canonical var: ensembl_id (index) + gene_symbol + mt ----------------------
gene_id_type = P.gene_id_type
if {"feature_id", "feature_name"} <= set(adata.var.columns):
    adata.var["ensembl_id"] = adata.var["feature_id"].astype(str).values
    adata.var["gene_symbol"] = adata.var["feature_name"].astype(str).values
elif "gene_symbol" not in adata.var.columns:
    # h5ad without Census var columns: assume var_names already hold the id type
    adata.var["ensembl_id"] = adata.var_names.astype(str)
    adata.var["gene_symbol"] = adata.var_names.astype(str)

adata.var_names = adata.var["ensembl_id"].values  # canonical index = Ensembl (matches loom path)
adata.var_names_make_unique()
adata.var["mt"] = adata.var["gene_symbol"].astype(str).str.upper().str.startswith("MT-")
n_mapped = int(adata.var["gene_symbol"].replace({"": np.nan, "nan": np.nan}).notna().sum())
log_transformation(log, "ingest_dataset",
    f"[{dataset}/{sample_id}] gene_symbol mapped {n_mapped}/{adata.n_vars} "
    f"({100 * n_mapped / adata.n_vars:.1f}%); {int(adata.var['mt'].sum())} MT genes "
    f"(gene_id_type={gene_id_type})")

# --- Canonical obs: dataset / batch / sample_id --------------------------------
adata.obs["dataset"] = dataset
adata.obs["batch"] = str(sample_id)
adata.obs["sample_id"] = str(sample_id)
adata.obs_names = [f"{sample_id}_{c}" for c in adata.obs_names]
adata.uns["normalization"] = "raw_counts_pre_qc"

# --- Write + provenance --------------------------------------------------------
Path(snakemake.output.h5ad).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(snakemake.output.h5ad)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="ingest_dataset",
    input_paths=[raw_file],
    tool_versions={"anndata": ad.__version__, "scanpy": sc.__version__},
    parameters={
        "dataset": dataset, "sample_id": sample_id, "source": source,
        "filter": filt, "sample_key": sample_key, "subsample_per_donor": cap,
        "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "n_genes_mapped": n_mapped, "gene_id_type": gene_id_type, "random_seed": seed,
    },
    description=f"Replication-cohort sample ingested to canonical h5ad ({source} source)",
    ontology_operation="operation:2409",  # EDAM: Format conversion
)
write_provenance(prov, snakemake.output.provenance)
log_transformation(log, "ingest_dataset", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad])
