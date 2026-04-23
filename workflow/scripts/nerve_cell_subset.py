"""Extract non-malignant neural/glial cells; re-cluster in nerve-cell subspace."""

import os
import sys

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

log_transformation(log, "nerve_cell_subset", f"Loading {snakemake.input.h5ad}")
adata = ad.read_h5ad(snakemake.input.h5ad)
n_total = adata.n_obs

# ---------------------------------------------------------------------------
# Step 1: Filter to non-malignant nerve cells
# ---------------------------------------------------------------------------
cell_types: list[str] = snakemake.params.cell_types

# Normalize type strings for comparison: lowercase, strip, collapse _↔space.
# scrna_annotate emits labels like "excitatory_neuron"; config lists them as
# "excitatory neuron" — both should match.
def _norm(s: str) -> str:
    return s.lower().strip().replace("_", " ")

adata.obs["_ctype_norm"] = adata.obs["cell_type_predicted"].astype(str).map(_norm)
type_targets = {_norm(t) for t in cell_types}

# Also accept partial matches for flexibility (e.g., "excitatory neuron" matches "neuron")
def _is_nerve(label: str) -> bool:
    label = _norm(label)
    if label in type_targets:
        return True
    # partial: scored as "neuron" covers both excitatory and inhibitory
    return any(label in t or t in label for t in type_targets)

nerve_mask = (
    adata.obs["_ctype_norm"].apply(_is_nerve)
    & (~adata.obs["is_malignant"])
)
adata_nerve = adata[nerve_mask].copy()
n_nerve = adata_nerve.n_obs

log_transformation(log, "nerve_cell_subset",
    f"Retained {n_nerve}/{n_total} non-malignant nerve cells "
    f"({100 * n_nerve / n_total:.1f}%)")
log_transformation(log, "nerve_cell_subset",
    f"Cell-type breakdown: {adata_nerve.obs['cell_type_predicted'].value_counts().to_dict()}")

if n_nerve < 2:
    log_transformation(log, "nerve_cell_subset",
        f"[FAIR-ALERT] 0 nerve cells selected — all cells labeled 'unscored' or malignant. "
        "This occurs with 1000-gene subsampled loom files that lack canonical markers. "
        "Writing placeholder outputs so the pipeline can complete.",
        status="WARNING")
    # Write placeholder outputs so downstream rules and rule-all targets exist
    adata_nerve = adata[:0].copy()  # empty AnnData, correct var schema
    adata_nerve.obs["nerve_leiden"] = pd.Series(dtype="category")
    adata_nerve.write_h5ad(snakemake.output.h5ad)
    # Placeholder UMAP
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "[FAIR-ALERT] No nerve cells found in subsampled data",
            ha="center", va="center", transform=ax.transAxes, fontsize=11)
    ax.set_axis_off()
    fig.savefig(snakemake.output.umap, dpi=100, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(columns=["nerve_leiden", "cell_type_predicted", "sample_id", "n_cells"]).to_csv(
        snakemake.output.composition, index=False
    )
    prov = stamp_artifact(
        output_path=snakemake.output.h5ad,
        rule_name="nerve_cell_subset",
        input_paths=[snakemake.input.h5ad, snakemake.input.clinical],
        tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
        parameters={"n_cells_total": n_total, "n_cells_nerve": 0, "note": "no nerve cells found"},
        description="Placeholder: no non-malignant nerve cells found in subsampled data",
        ontology_operation="operation:3432",
    )
    write_provenance(prov, snakemake.output.provenance)
    log_transformation(log, "nerve_cell_subset", "Placeholder outputs written", status="WARNING",
                       artifact_paths=[snakemake.output.h5ad, snakemake.output.umap])
    raise SystemExit(0)

if n_nerve < 20:
    log_transformation(log, "nerve_cell_subset",
        f"WARNING: only {n_nerve} nerve cells — downstream analysis may be underpowered",
        status="WARNING")

# ---------------------------------------------------------------------------
# Step 2: Join clinical metadata
# ---------------------------------------------------------------------------
clinical_df = pd.read_csv(snakemake.input.clinical, sep="\t", dtype=str)
clinical_df = clinical_df.set_index("file_uuid")

# sample_id in obs corresponds to file UUID
_clinical_cols = ["primary_diagnosis", "tumor_grade", "prior_malignancy",
                  "tissue_type", "gender", "race", "age_at_index"]
adata_nerve.obs = adata_nerve.obs.join(
    clinical_df[_clinical_cols],
    on="sample_id",
    how="left",
)
# h5py cannot serialize object columns mixing strings and float NaN. Coerce the
# joined clinical metadata to plain strings with "unknown" sentinel (e.g.,
# prior_malignancy and age_at_index are empty in the source TSV).
for _col in _clinical_cols:
    adata_nerve.obs[_col] = adata_nerve.obs[_col].fillna("unknown").astype(str)

_n_with_dx = (adata_nerve.obs["primary_diagnosis"] != "unknown").sum()
log_transformation(log, "nerve_cell_subset",
    f"Joined clinical metadata; {_n_with_dx} cells have primary_diagnosis "
    f"(clinical columns coerced to str with 'unknown' fill)")

# ---------------------------------------------------------------------------
# Step 3: Re-cluster in nerve-cell subspace
# ---------------------------------------------------------------------------
# Re-select HVGs within this subset
sc.pp.normalize_total(adata_nerve, target_sum=1e4)
sc.pp.log1p(adata_nerve)
sc.pp.highly_variable_genes(
    adata_nerve,
    n_top_genes=min(snakemake.params.n_top_genes, adata_nerve.n_vars),
    subset=False,
)
adata_nerve.raw = adata_nerve

# PCA on nerve-cell HVGs
sc.tl.pca(adata_nerve, n_comps=min(30, n_nerve - 1, adata_nerve.n_vars - 1),
          random_state=snakemake.params.random_seed)
sc.pp.neighbors(adata_nerve, use_rep="X_pca", random_state=snakemake.params.random_seed)
sc.tl.umap(adata_nerve, random_state=snakemake.params.random_seed)
sc.tl.leiden(adata_nerve,
             resolution=snakemake.params.leiden_resolution,
             random_state=snakemake.params.random_seed,
             key_added="nerve_leiden")

n_clusters = adata_nerve.obs["nerve_leiden"].nunique()
log_transformation(log, "nerve_cell_subset",
    f"Re-clustered → {n_clusters} nerve-cell leiden clusters "
    f"(resolution={snakemake.params.leiden_resolution})")

# ---------------------------------------------------------------------------
# Step 4: UMAP plot colored by nerve cluster and cell type
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sc.pl.umap(adata_nerve, color="nerve_leiden", ax=axes[0], show=False,
           title="Nerve Cells — Leiden Clusters")
sc.pl.umap(adata_nerve, color="cell_type_predicted", ax=axes[1], show=False,
           title="Nerve Cells — Predicted Cell Type")

fig.suptitle("Non-malignant Nerve Cell Subspace", fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(snakemake.output.umap, dpi=150, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Step 5: Cluster composition table
# ---------------------------------------------------------------------------
comp = (
    adata_nerve.obs
    .groupby(["nerve_leiden", "cell_type_predicted", "sample_id"])
    .size()
    .reset_index(name="n_cells")
)
comp.to_csv(snakemake.output.composition, index=False)

# --- Write output ------------------------------------------------------------
adata_nerve.write_h5ad(snakemake.output.h5ad)
verify_artifact(snakemake.output.h5ad, min_size_bytes=512)

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="nerve_cell_subset",
    input_paths=[snakemake.input.h5ad, snakemake.input.clinical],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "cell_types_selected":  cell_types,
        "leiden_resolution":    snakemake.params.leiden_resolution,
        "n_cells_total":        n_total,
        "n_cells_nerve":        n_nerve,
        "n_leiden_clusters":    n_clusters,
        "n_top_genes":          snakemake.params.n_top_genes,
    },
    description="Non-malignant nerve-cell AnnData subspace with clinical metadata and re-clustering",
    ontology_operation="operation:3432",
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "nerve_cell_subset", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.umap,
                                   snakemake.output.composition])
