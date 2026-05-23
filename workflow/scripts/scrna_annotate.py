"""Cell-type annotation using scVI latent clustering and canonical marker gene scoring."""

import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

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
MARKER_SETS: dict[str, list[str]] = {
    "neuron":            snakemake.params.markers.get("neuron", []),
    "excitatory_neuron": snakemake.params.markers.get("excitatory_neuron", []),
    "inhibitory_neuron": snakemake.params.markers.get("inhibitory_neuron", []),
    "opc":               snakemake.params.markers.get("opc", []),
    "oligodendrocyte":   snakemake.params.markers.get("oligodendrocyte", []),
    "astrocyte":         snakemake.params.markers.get("astrocyte", []),
    "ependymal":         snakemake.params.markers.get("ependymal", []),
    # Additional TME cell types for annotation context
    "microglia":         ["AIF1", "CX3CR1", "TMEM119", "P2RY12", "PTPRC"],
    "t_cell":            ["CD3D", "CD3E", "CD4", "CD8A", "CD8B"],
    "endothelial":       ["PECAM1", "CDH5", "VWF", "ENG"],
    "tumor_gbm":         ["SOX2", "NES", "CD44", "PROM1", "MKI67"],
}

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

# --- Assign cell type per cluster by max mean score -------------------------
score_cols = [f"score_{s}" for s in scored_sets]

label_map: dict[str, str] = {}
if score_cols:
    cluster_means = (
        adata.obs[["leiden"] + score_cols]
        .groupby("leiden", observed=True)[score_cols]
        .mean()
    )
    for cluster in cluster_means.index:
        row = cluster_means.loc[cluster]
        label_map[cluster] = row.idxmax().replace("score_", "") if not row.empty else "unscored"
else:
    log_transformation(log, "scrna_annotate",
        "WARNING: no marker gene sets scored — all clusters labeled 'unscored'", status="WARNING")
    for cluster in adata.obs["leiden"].unique():
        label_map[cluster] = "unscored"

adata.obs["cell_type_predicted"] = adata.obs["leiden"].map(label_map).astype(str)

# Confidence: difference between top and second-best score per cell
if len(score_cols) >= 2:
    score_matrix = adata.obs[score_cols].values
    sorted_scores = np.sort(score_matrix, axis=1)
    adata.obs["cell_type_confidence"] = (sorted_scores[:, -1] - sorted_scores[:, -2]).astype(np.float32)
else:
    adata.obs["cell_type_confidence"] = np.float32(1.0)

# Record schema in uns for FAIR interoperability
adata.uns["annotation_schema"] = {
    "method":            "leiden_marker_scoring",
    "leiden_resolution": res,
    "marker_sets":       {k: v for k, v in MARKER_SETS.items() if k in scored_sets},
    "ontology":          "Allen Brain Atlas cell type taxonomy",
    "edam_operation":    "operation:3432",  # EDAM: Clustering
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
adata.write_h5ad(snakemake.output.h5ad)
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
    },
    description="Cell-type annotated AnnData with Leiden clusters, UMAP, and marker scores",
    ontology_operation="operation:3432",
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "scrna_annotate", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.summary])
