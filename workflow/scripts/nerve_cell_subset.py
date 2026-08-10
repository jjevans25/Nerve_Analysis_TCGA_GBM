"""Extract non-malignant neural/glial cells; re-cluster in nerve-cell subspace."""

import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from fair_utils import H5AD_COMPRESSION, log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

log_transformation(log, "nerve_cell_subset", f"Loading {snakemake.input.h5ad}")
adata = ad.read_h5ad(snakemake.input.h5ad)
n_total = adata.n_obs

# ---------------------------------------------------------------------------
# Step 1: Filter to non-malignant nerve cells
# ---------------------------------------------------------------------------
cell_types: list[str] = snakemake.params.cell_types
frozen_path: str | None = getattr(snakemake.params, "frozen_subset_file", None)
# Which malignancy flag excludes a cell from the nerve compartment. `is_malignant`
# is the balanced call used to BUILD the tumor compartment; a compartment that
# must be clean can instead key on a more sensitive flag, trading tumor-side
# recall for nerve-side purity. See scrna_malignancy.
malignant_flag: str = getattr(snakemake.params, "malignant_flag", None) or "is_malignant"
if malignant_flag not in adata.obs.columns:
    raise KeyError(
        f"[FAIR-ALERT] nerve_cells.malignant_flag='{malignant_flag}' is not a column "
        f"of the annotated input. Available: {sorted(c for c in adata.obs.columns if 'malig' in c)}"
    )


# Normalize type strings for comparison: lowercase, strip, collapse _↔space.
# scrna_annotate emits labels like "excitatory_neuron"; config lists them as
# "excitatory neuron" — both should match.
def _norm(s: str) -> str:
    return s.lower().strip().replace("_", " ")


if frozen_path:
    # v1.2.0 panel-tightening insulator: select cells by frozen barcode list
    # rather than by current cell_type_predicted / is_malignant. Keeps the
    # nerve subset's composition (and therefore nerve_leiden cluster IDs)
    # stable across ependymal-panel changes. Retire at v1.3.0 full rerun.
    frozen_barcodes = set(Path(frozen_path).read_text().split())
    log_transformation(log, "nerve_cell_subset",
        f"Frozen subset active: {len(frozen_barcodes)} barcodes from {frozen_path}")
    nerve_mask = adata.obs_names.isin(frozen_barcodes)
    missing = frozen_barcodes - set(adata.obs_names)
    if missing:
        log_transformation(log, "nerve_cell_subset",
            f"[FAIR-ALERT] {len(missing)} frozen barcodes missing from current "
            f"annotated input — upstream cell roster has changed since the "
            f"freeze. First 5: {list(missing)[:5]}", status="ERROR")
        raise RuntimeError(
            "Frozen subset incompatible with current annotated input. Either "
            "regenerate the freeze (scripts/freeze_nerve_subset_v1_1_0.py) "
            "against the current upstream artifact, or remove "
            "nerve_cells.frozen_subset_file from config to fall back to "
            "live cell_type_predicted-based subsetting."
        )
else:
    # Default path: derive nerve mask from current cell_type_predicted.
    #
    # Exact matching on a canonical label set, replacing substring matching.
    # The old `any(label in t or t in label ...)` test was silently bidirectional
    # and silently lenient: config said "oligodendrocyte precursor cell" while
    # annotate emitted "opc", neither is a substring of the other, and 70,881
    # cells were dropped without a word in the log. It was also over-permissive
    # in the other direction — "neuron" is a substring of "excitatory neuron",
    # so unrelated additions to the panel list could silently widen the mask.
    adata.obs["_ctype_norm"] = adata.obs["cell_type_predicted"].astype(str).map(_norm)
    type_targets = {_norm(t) for t in cell_types}
    observed_labels = set(adata.obs["_ctype_norm"].unique())

    # A configured cell type that matches no observed label is a typo or a stale
    # name, not an absent population. Fail loudly: silently selecting nothing is
    # precisely how this defect survived a full pipeline run.
    unmatched = sorted(type_targets - observed_labels)
    if unmatched:
        raise RuntimeError(
            f"[FAIR-ALERT] nerve_cells.cell_types entries match no observed "
            f"cell_type_predicted value: {unmatched}. Observed labels: "
            f"{sorted(observed_labels)}. Every configured type must match a "
            "canonical label emitted by scrna_annotate (see the panel names in "
            "config annotation_markers / nerve_cells.markers). Refusing to build "
            "a compartment from a config entry that selects nothing."
        )

    nerve_mask = (
        adata.obs["_ctype_norm"].isin(type_targets)
        & (~adata.obs[malignant_flag])
    )

adata_nerve = adata[nerve_mask].copy()
n_nerve = adata_nerve.n_obs

log_transformation(log, "nerve_cell_subset",
    f"Retained {n_nerve}/{n_total} non-malignant nerve cells "
    f"({100 * n_nerve / n_total:.1f}%)")
log_transformation(log, "nerve_cell_subset",
    f"Cell-type breakdown: {adata_nerve.obs['cell_type_predicted'].value_counts().to_dict()}")

# ---------------------------------------------------------------------------
# Step 1b: split the compartment into glia and neurons
#
# Neurons and glia are not interchangeable signalling partners, and lumping them
# into one "nerve" compartment lets a glial cluster's ligand profile be read as
# neuronal. They are separated here rather than at the LIANA step so the
# distinction is carried on the object and is auditable.
#
# Neurons stay a single group instead of being split per cluster: there are only
# ~3.4k of them across 170 donors, so per-cluster resolution would be noise.
# ---------------------------------------------------------------------------
neuron_labels = {_norm(t) for t in (getattr(snakemake.params, "neuron_labels", None) or [])}
_norm_pred = adata_nerve.obs["cell_type_predicted"].astype(str).map(_norm)
adata_nerve.obs["nerve_subcompartment"] = np.where(
    _norm_pred.isin(neuron_labels), "neuron", "glia"
)
_sub_counts = adata_nerve.obs["nerve_subcompartment"].value_counts().to_dict()
log_transformation(log, "nerve_cell_subset",
    f"Sub-compartments: {_sub_counts}")
if neuron_labels and _sub_counts.get("neuron", 0) == 0:
    log_transformation(log, "nerve_cell_subset",
        f"[FAIR-ALERT] no cells matched nerve_cells.neuron_labels {sorted(neuron_labels)} — "
        "the neuron group will be absent from the interaction tables entirely.",
        status="WARNING")

if n_nerve < 2:
    # Hard-fail instead of writing placeholders. A 0-cell nerve compartment
    # almost always means the integrated gene set lost its markers upstream
    # (every cell labeled 'unscored'), not that the biology is absent. The old
    # placeholder-and-exit-0 path laundered that upstream failure into a green
    # pipeline; CLAUDE.md requires each output be substantive.
    raise RuntimeError(
        f"[FAIR-ALERT] 0 of {n_total} cells selected as non-malignant nerve cells. "
        "An empty compartment almost always means the integrated gene set lost its "
        "markers upstream (all cells 'unscored'), not that the biology is absent. "
        "See markdowns/diagnosis_rerun_gene_intersection.md"
    )

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

# PCA on nerve-cell HVGs (kept as a reference embedding alongside the scVI
# integration; not used for the neighbourhood graph).
sc.tl.pca(adata_nerve, n_comps=min(30, n_nerve - 1, adata_nerve.n_vars - 1),
          random_state=snakemake.params.random_seed)
# Use the scVI batch-corrected latent space for the kNN graph; PCA on the
# nerve-cell subset re-introduces patient identity (verified by the
# nerve_batch_qc rule: PCA-based clustering produced 32/36 patient-pure
# clusters; X_scVI from the upstream `scrna_integration` rule preserves
# inter-patient mixing).
sc.pp.neighbors(adata_nerve, use_rep="X_scVI", random_state=snakemake.params.random_seed)
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
# Step 3b: Apply manual cluster split overrides (v1.3.0 cl15 surgery)
# ---------------------------------------------------------------------------
# The v1.2.0 rerun confirmed that ependymal-panel methodology alone cannot
# resolve cl15's MIXED status (cluster-level argmax is robust to per-cell score
# shifts). The fix is surgical: relabel cl15's four sub-Leiden subpopulations to
# separate nerve_leiden IDs. The split is frozen per-barcode in
# provenance/cl15_split_v1_3_0.csv (see scripts/freeze_cl15_split_v1_3_0.py) so
# no live re-clustering happens here — mirroring the frozen_subset_file pattern.
# Applied AFTER nerve_leiden so only the listed cl15 barcodes are relabeled;
# every other cluster's IDs (and their v1.0.0 verdicts) stay bit-exact.
split_path: str | None = getattr(snakemake.params, "split_assignments_file", None)
if split_path:
    split_df = pd.read_csv(split_path, dtype=str).set_index("barcode")
    log_transformation(log, "nerve_cell_subset",
        f"Cluster-split override active: {len(split_df)} barcodes from {split_path}")

    overlap = adata_nerve.obs_names.isin(split_df.index)
    n_overlap = int(overlap.sum())
    missing = set(split_df.index) - set(adata_nerve.obs_names)
    if missing:
        log_transformation(log, "nerve_cell_subset",
            f"[FAIR-ALERT] {len(missing)} split-assignment barcodes missing from "
            f"the nerve subset — the freeze that produced the split is stale. "
            f"First 5: {list(missing)[:5]}", status="ERROR")
        raise RuntimeError(
            "Cluster-split assignment incompatible with current nerve subset. "
            "Regenerate scripts/freeze_cl15_split_v1_3_0.py against the current "
            "nerve_cells.h5ad, or remove nerve_cells.cluster_overrides from config."
        )

    # Snapshot pre-split per-cluster sizes for the freeze-integrity check below.
    sizes_before = adata_nerve.obs["nerve_leiden"].astype(str).value_counts().to_dict()

    # Barcode-keyed remap: robust even if leiden numbering ever drifts, since we
    # relabel exactly the frozen barcodes rather than "whatever is in cl15".
    new_leiden = adata_nerve.obs["nerve_leiden"].astype(str).copy()
    target_by_barcode = split_df["target_cluster"]
    new_leiden.loc[adata_nerve.obs_names[overlap]] = target_by_barcode.reindex(
        adata_nerve.obs_names[overlap]
    ).to_numpy()
    adata_nerve.obs["nerve_leiden"] = new_leiden.astype("category")

    # Freeze-integrity check: every cluster NOT touched by the split must keep
    # its exact size (the v1.2.0 rerun guaranteed bit-exact cluster sizes; the
    # surgery must not perturb any cluster other than the split source).
    touched_sources = set(split_df["sub_id"].index)  # noqa: F841 (doc only)
    split_targets = set(split_df["target_cluster"].unique())
    sizes_after = adata_nerve.obs["nerve_leiden"].astype(str).value_counts().to_dict()
    # The source cluster(s) — clusters whose cells were reassigned — are those
    # present before but shrunk/absent after, excluding the new target IDs.
    drifted = []
    for cl, n in sizes_before.items():
        if cl in split_targets:
            continue
        if sizes_after.get(cl, 0) != n:
            drifted.append((cl, n, sizes_after.get(cl, 0)))
    # Any cluster fully consumed by the split (e.g. cl15 → 0) is expected; flag
    # only clusters that PARTIALLY changed, which would mean the split touched
    # cells outside its source cluster.
    partial = [(cl, b, a) for cl, b, a in drifted if a != 0]
    if partial:
        log_transformation(log, "nerve_cell_subset",
            f"[FAIR-ALERT] cluster-split perturbed non-source clusters "
            f"(cluster, before, after): {partial}. The split barcodes span more "
            f"than the intended source cluster.", status="ERROR")
        raise RuntimeError(
            "Cluster-split override altered sizes of clusters outside its source. "
            "Inspect provenance/cl15_split_v1_3_0.csv — every barcode must belong "
            "to the source cluster being split."
        )

    emptied = [cl for cl, _, a in drifted if a == 0]
    new_sizes = {t: sizes_after.get(t, 0) for t in sorted(split_targets, key=int)}
    log_transformation(log, "nerve_cell_subset",
        f"Split applied: source cluster(s) {emptied} → new IDs {new_sizes}; "
        f"{n_overlap} cells relabeled; all other clusters unchanged.")

    # Drop the now-empty source category so it doesn't linger in the dtype.
    adata_nerve.obs["nerve_leiden"] = (
        adata_nerve.obs["nerve_leiden"].cat.remove_unused_categories()
    )
    n_clusters = adata_nerve.obs["nerve_leiden"].nunique()

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
adata_nerve.write_h5ad(snakemake.output.h5ad, compression=H5AD_COMPRESSION)
verify_artifact(snakemake.output.h5ad, min_size_bytes=512)

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="nerve_cell_subset",
    input_paths=[snakemake.input.h5ad, snakemake.input.clinical],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__},
    parameters={
        "cell_types_selected":  cell_types,
        "frozen_subset_file":   frozen_path,
        "split_assignments_file": split_path,
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
