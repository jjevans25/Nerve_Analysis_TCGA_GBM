"""Assign per-cell canonical cell-type labels for scANVI anchoring.

Question answered: scANVI uses biological labels to push the latent space
away from patient-private structure. We need labels that are derived from
canonical markers (GFAP, MBP, OLIG2, PDGFRA, SYN1/SNAP25, …) — *not* from
the v1.0.0 X_scVI clusters — so using them to refine the embedding is not
circular.

Labels (4 categories + Unknown):
    neuron           = union of generic + excitatory + inhibitory markers
    opc              = PDGFRA / CSPG4 / SOX10 / OLIG1
    oligodendrocyte  = MBP / MOG / PLP1 / OLIG2
    astrocyte        = GFAP / S100B / AQP4 / VIM
    Unknown          = cells whose max marker score falls below the
                       configured percentile (default 20th)

scANVI handles "Unknown" as the unlabelled category: it still enters the
variational objective but is not used for the classification loss.
"""

import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, "workflow/scripts")
from fair_utils import (
    H5AD_COMPRESSION,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(int(snakemake.params.random_seed))
np.random.seed(int(snakemake.params.random_seed))

markers_cfg: dict = dict(snakemake.params.markers)
unknown_percentile = float(snakemake.params.unknown_percentile)

LABEL_GROUPS: dict[str, list[str]] = {
    "neuron": sorted(set(
        markers_cfg.get("neuron", [])
        + markers_cfg.get("excitatory_neuron", [])
        + markers_cfg.get("inhibitory_neuron", [])
    )),
    "opc": list(markers_cfg.get("opc", [])),
    "oligodendrocyte": list(markers_cfg.get("oligodendrocyte", [])),
    "astrocyte": list(markers_cfg.get("astrocyte", [])),
    # Motile-ciliated ventricular cells; markers cluster on DNAH/CFAP and
    # diverge from astrocyte despite the historical mislabel (cl19 in v1.0.0).
    "ependymal": list(markers_cfg.get("ependymal", [])),
}

log_transformation(log, "nerve_celltype_labels",
                   f"Loading counts h5ad: {snakemake.input.counts_h5ad}")
adata = ad.read_h5ad(snakemake.input.counts_h5ad)
log_transformation(log, "nerve_celltype_labels",
                   f"{adata.n_obs} cells x {adata.n_vars} genes")

# Map marker symbols -> ensembl IDs (var.index in this file is ensembl).
if "gene_symbol" not in adata.var.columns:
    raise RuntimeError("[FAIR-ALERT] var.gene_symbol missing — cannot resolve marker symbols")
symbol_to_ensembl = (
    adata.var.reset_index()
    .dropna(subset=["gene_symbol"])
    .drop_duplicates(subset=["gene_symbol"], keep="first")
    .set_index("gene_symbol")["index"]
    .to_dict()
)

resolved: dict[str, list[str]] = {}
missing: dict[str, list[str]] = {}
for label, symbols in LABEL_GROUPS.items():
    keep, drop = [], []
    for s in symbols:
        eid = symbol_to_ensembl.get(s)
        (keep if eid is not None else drop).append(eid if eid is not None else s)
    resolved[label] = keep
    missing[label] = drop
    log_transformation(log, "nerve_celltype_labels",
                       f"  {label}: {len(keep)} markers found, {len(drop)} missing "
                       f"{('(missing: ' + ','.join(drop) + ')') if drop else ''}")

for label, genes in resolved.items():
    if len(genes) < 2:
        raise RuntimeError(
            f"[FAIR-ALERT] cell-type '{label}' has < 2 resolvable markers; "
            f"cannot score reliably."
        )

# Normalize a working copy for scoring (do NOT mutate counts in .X).
adata_norm = adata.copy()
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)
log_transformation(log, "nerve_celltype_labels",
                   "Computed library-size normalisation + log1p on a working copy")

for label, genes in resolved.items():
    score_name = f"score_{label}"
    sc.tl.score_genes(
        adata_norm, gene_list=genes, score_name=score_name,
        random_state=int(snakemake.params.random_seed), use_raw=False,
    )
    # Copy back to the counts adata for persistence.
    adata.obs[score_name] = adata_norm.obs[score_name].values

# Argmax over the 4 score columns -> tentative label.
score_cols = [f"score_{lbl}" for lbl in LABEL_GROUPS]
score_matrix = adata.obs[score_cols].to_numpy()
argmax_idx = np.argmax(score_matrix, axis=1)
max_score = score_matrix[np.arange(len(argmax_idx)), argmax_idx]
tentative = np.asarray(list(LABEL_GROUPS.keys()))[argmax_idx]

threshold = float(np.percentile(max_score, unknown_percentile))
unknown_mask = max_score < threshold
final = tentative.copy()
final[unknown_mask] = "Unknown"

adata.obs["cell_type_score_max"] = max_score
adata.obs["cell_type"] = pd.Categorical(
    final, categories=list(LABEL_GROUPS.keys()) + ["Unknown"]
)

# ----------------------------------------------------------------------------
# Targets / verification gates
# ----------------------------------------------------------------------------
counts = adata.obs["cell_type"].value_counts()
total = int(counts.sum())
unknown_frac = float(counts.get("Unknown", 0)) / total
log_transformation(log, "nerve_celltype_labels",
                   f"Label distribution (n={total}):\n{counts.to_string()}")
log_transformation(log, "nerve_celltype_labels",
                   f"Unknown fraction = {unknown_frac:.3f} "
                   f"(threshold @ p{unknown_percentile:g} = {threshold:.4f})")

if unknown_frac >= 0.20 + 1e-6:
    log_transformation(log, "nerve_celltype_labels",
                       f"Unknown fraction {unknown_frac:.3f} exceeds soft cap 0.20 — "
                       f"consider lowering unknown_percentile.",
                       status="WARNING")

RARE_LABELS = {"ependymal"}  # anatomically rare; use a looser min-fraction floor
for lbl in LABEL_GROUPS:
    frac = float(counts.get(lbl, 0)) / total
    min_frac = 0.001 if lbl in RARE_LABELS else 0.01
    if frac < min_frac:
        raise RuntimeError(
            f"[FAIR-ALERT] cell-type '{lbl}' got only {frac:.3%} of cells "
            f"(floor {min_frac:.1%}) — label collapse; refine markers or threshold "
            f"before training scANVI."
        )

summary = pd.DataFrame(
    {
        "cell_type": list(counts.index),
        "n_cells": counts.values,
        "fraction": (counts.values / total).round(4),
        "n_markers_used": [len(resolved.get(lbl, [])) for lbl in counts.index],
        "n_markers_missing": [len(missing.get(lbl, [])) for lbl in counts.index],
    }
)
summary.to_csv(snakemake.output.summary, index=False)
verify_artifact(snakemake.output.summary, min_size_bytes=64)

adata.write_h5ad(snakemake.output.h5ad, compression=H5AD_COMPRESSION)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1_000_000)

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="nerve_celltype_labels",
    input_paths=[snakemake.input.counts_h5ad],
    tool_versions={
        "anndata": ad.__version__,
        "scanpy": sc.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
    },
    parameters={
        "n_cells": int(adata.n_obs),
        "unknown_percentile": unknown_percentile,
        "unknown_score_threshold": round(threshold, 4),
        "unknown_fraction": round(unknown_frac, 4),
        "label_groups": {k: v for k, v in LABEL_GROUPS.items()},
        "n_markers_resolved": {k: len(v) for k, v in resolved.items()},
        "n_markers_missing": {k: len(v) for k, v in missing.items()},
        "random_seed": int(snakemake.params.random_seed),
    },
    description=(
        "Per-cell canonical cell-type labels for scANVI anchoring. Labels "
        "derived from sc.tl.score_genes on library-size-normalised log1p "
        "expression of marker genes (independent of v1.0.0 X_scVI clusters)."
    ),
    ontology_operation="operation:2962",  # EDAM: Codon usage analysis (closest: classification)
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(
    log, "nerve_celltype_labels",
    "Complete",
    status="SUCCESS",
    artifact_paths=[
        snakemake.output.h5ad,
        snakemake.output.summary,
        snakemake.output.provenance,
    ],
)
