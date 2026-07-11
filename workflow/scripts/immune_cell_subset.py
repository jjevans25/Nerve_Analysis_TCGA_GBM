"""Extract the tumor-associated immune compartment; re-cluster in immune subspace.

Biological question: what myeloid/lymphoid subpopulations make up the GBM
immune infiltrate, and how do they partition when re-clustered away from the
neural/glial and malignant cells?

The whole-dataset annotation (`scrna_annotate`) assigns a single coarse
``microglia`` argmax label to the entire immune blob because its Leiden
resolution merges myeloid and lymphoid cells. Here we subset those cells and
re-cluster on the existing scVI batch-corrected latent (``X_scVI`` — no model
retraining) so that microglia / TAM / T / NK / DC subtypes can separate. Mirrors
the ``nerve_cell_subset`` preprocessing so the immune compartment is on equal
footing with the nerve compartment in the downstream three-way LIANA analysis.
"""

import os
import sys

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    compute_cluster_purity,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]  # type: ignore[name-defined]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)  # type: ignore[name-defined]

log_transformation(log, "immune_cell_subset", f"Loading {snakemake.input.h5ad}")  # type: ignore[name-defined]
adata = ad.read_h5ad(snakemake.input.h5ad)  # type: ignore[name-defined]
n_total = adata.n_obs

# ---------------------------------------------------------------------------
# Step 1: Filter to non-malignant immune cells
# ---------------------------------------------------------------------------
source_label: str = snakemake.params.source_label  # type: ignore[name-defined]


def _norm(s: str) -> str:
    return s.lower().strip().replace("_", " ")


adata.obs["_ctype_norm"] = adata.obs["cell_type_predicted"].astype(str).map(_norm)
immune_mask = (adata.obs["_ctype_norm"] == _norm(source_label)) & (
    ~adata.obs["is_malignant"].astype(bool)
)
adata_immune = adata[immune_mask].copy()
del adata_immune.obs["_ctype_norm"]
n_immune = adata_immune.n_obs

log_transformation(
    log,
    "immune_cell_subset",
    f"Retained {n_immune}/{n_total} non-malignant immune cells "
    f"(cell_type_predicted == '{source_label}'; {100 * n_immune / n_total:.1f}%)",
)

if n_immune < 20:
    log_transformation(
        log,
        "immune_cell_subset",
        f"[FAIR-ALERT] only {n_immune} immune cells selected — the "
        f"'{source_label}' population is missing or vanishingly small. "
        "Downstream subtype resolution and interaction inference will be "
        "underpowered.",
        status="WARNING",
    )

# ---------------------------------------------------------------------------
# Step 2: Re-cluster in the immune subspace (clustering on X_scVI)
# ---------------------------------------------------------------------------
# Mirror nerve_cell_subset preprocessing: re-select HVGs + PCA as a reference
# embedding, but build the kNN graph on the scVI batch-corrected latent so the
# clusters reflect cell state, not patient identity.
sc.pp.normalize_total(adata_immune, target_sum=1e4)
sc.pp.log1p(adata_immune)
sc.pp.highly_variable_genes(
    adata_immune,
    n_top_genes=min(snakemake.params.n_top_genes, adata_immune.n_vars),  # type: ignore[name-defined]
    subset=False,
)
adata_immune.raw = adata_immune

sc.tl.pca(
    adata_immune,
    n_comps=min(30, n_immune - 1, adata_immune.n_vars - 1),
    random_state=snakemake.params.random_seed,  # type: ignore[name-defined]
)
sc.pp.neighbors(
    adata_immune, use_rep="X_scVI", random_state=snakemake.params.random_seed  # type: ignore[name-defined]
)
sc.tl.umap(adata_immune, random_state=snakemake.params.random_seed)  # type: ignore[name-defined]
sc.tl.leiden(
    adata_immune,
    resolution=snakemake.params.leiden_resolution,  # type: ignore[name-defined]
    random_state=snakemake.params.random_seed,  # type: ignore[name-defined]
    key_added="immune_leiden",
    flavor="igraph",
    directed=False,
    n_iterations=2,
)
n_clusters = adata_immune.obs["immune_leiden"].nunique()
log_transformation(
    log,
    "immune_cell_subset",
    f"Re-clustered → {n_clusters} immune Leiden clusters "
    f"(resolution={snakemake.params.leiden_resolution})",  # type: ignore[name-defined]
)

# ---------------------------------------------------------------------------
# Step 3: UMAP plot colored by immune cluster
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(
    adata_immune,
    color="immune_leiden",
    ax=ax,
    show=False,
    title=f"Tumor-associated immune subspace — {n_clusters} Leiden clusters",
)
fig.tight_layout()
fig.savefig(snakemake.output.umap, dpi=150, bbox_inches="tight")  # type: ignore[name-defined]
plt.close(fig)

# ---------------------------------------------------------------------------
# Step 4: Cluster composition + per-cluster batch purity
# ---------------------------------------------------------------------------
comp = (
    adata_immune.obs.groupby(
        ["immune_leiden", "sample_id"], observed=True
    )
    .size()
    .reset_index(name="n_cells")
)
comp.to_csv(snakemake.output.composition, index=False)  # type: ignore[name-defined]

purity = compute_cluster_purity(
    adata_immune,
    cluster_col="immune_leiden",
    batch_col="sample_id",
    dominant_max=snakemake.params.dominant_fraction_max,  # type: ignore[name-defined]
    min_contributing_fraction=snakemake.params.min_contributing_fraction,  # type: ignore[name-defined]
    min_contributing_samples=snakemake.params.min_contributing_samples,  # type: ignore[name-defined]
)
purity.df.to_csv(snakemake.output.purity, index=False)  # type: ignore[name-defined]
n_pass = int(purity.df["pass_overall"].sum())
log_transformation(
    log,
    "immune_cell_subset",
    f"Batch purity: {n_pass}/{len(purity.df)} clusters PASS "
    f"(dominant-fraction < {snakemake.params.dominant_fraction_max})",  # type: ignore[name-defined]
)

# ---------------------------------------------------------------------------
# Step 5: Write outputs + provenance
# ---------------------------------------------------------------------------
adata_immune.write_h5ad(snakemake.output.h5ad)  # type: ignore[name-defined]
verify_artifact(snakemake.output.h5ad, min_size_bytes=512)  # type: ignore[name-defined]

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,  # type: ignore[name-defined]
    rule_name="immune_cell_subset",
    input_paths=[snakemake.input.h5ad],  # type: ignore[name-defined]
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "source_label": source_label,
        "leiden_resolution": snakemake.params.leiden_resolution,  # type: ignore[name-defined]
        "n_cells_total": n_total,
        "n_cells_immune": n_immune,
        "n_leiden_clusters": int(n_clusters),
        "n_clusters_pass_qc": n_pass,
        "n_top_genes": snakemake.params.n_top_genes,  # type: ignore[name-defined]
    },
    description="Non-malignant tumor-associated immune AnnData subspace, re-clustered on X_scVI",
    ontology_operation="operation:3432",  # EDAM: Clustering
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]

log_transformation(
    log,
    "immune_cell_subset",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.h5ad,  # type: ignore[name-defined]
        snakemake.output.umap,  # type: ignore[name-defined]
        snakemake.output.composition,  # type: ignore[name-defined]
        snakemake.output.purity,  # type: ignore[name-defined]
    ],
)
