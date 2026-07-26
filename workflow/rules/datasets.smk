# workflow/rules/datasets.smk
# Cohort-namespaced replication track. Runs the reference cohort's analysis
# SCRIPTS unchanged, with all I/O redirected under a per-dataset namespace
# (<dir>/<dataset>/...). The reference cohort (top-level `samples:`) is never
# re-run. See markdowns/plan_second_gbm_replication_cohort.md.
#
# Scope note: the three-way nerve-tumor-immune interaction consumes the v1 nerve
# subset (nerve_cells.h5ad), NOT the scANVI-v2 latent, so the scANVI retrain
# branch is intentionally omitted here — it does not affect the concordance
# deliverable and would double the GPU cost.

import os
import re

DATASETS = config.get("datasets", {})


def _sheet(dataset: str) -> list[str]:
    """Committed per-dataset sample sheet (donor IDs), generated once by
    scripts/derive_dataset_sample_sheet.py. Must exist at parse time so the
    integration expand() can build the DAG."""
    sheet = os.path.join(config["dirs"]["data_raw"], dataset, "samples.txt")
    if not os.path.exists(sheet):
        return []
    return [ln.strip() for ln in open(sheet) if ln.strip()]


ALL_DS_SAMPLES = sorted({s for d in DATASETS for s in _sheet(d)})


# Constrain the dataset wildcard to known keys (reference rules have no
# {dataset} wildcard, so this is inert for them). Constrain sample to exclude
# '/' so no rule's {sample} can absorb the dataset directory boundary — this
# keeps the reference {sample} rules from matching dataset-scoped paths.
wildcard_constraints:
    dataset = "|".join(re.escape(d) for d in DATASETS) or "__no_datasets__",
    sample  = r"[^/]+",


def _dp(*parts) -> str:
    return os.path.join(config["dirs"]["data_processed"], "{dataset}", *parts)

def _tb(*parts) -> str:
    return os.path.join(config["dirs"]["tables"], "{dataset}", *parts)

def _fg(*parts) -> str:
    return os.path.join(config["dirs"]["figures"], "{dataset}", *parts)

def _pv(*parts) -> str:
    return os.path.join(config["dirs"]["provenance"], "{dataset}", *parts)

def _md(*parts) -> str:
    return os.path.join(config["dirs"]["models"], "{dataset}", *parts)

def _ext(*parts) -> str:
    return os.path.join(config["dirs"]["data_external"], "{dataset}", *parts)


def _entry(dataset: str) -> dict:
    return config["datasets"][dataset]


# =============================================================================
# Stage A — Generalized ingestion
# =============================================================================

rule ds_ingest_dataset:
    """Ingest one sample of a replication cohort to canonical h5ad (raw counts,
    ensembl_id+gene_symbol var, dataset/batch/sample_id obs)."""
    wildcard_constraints:
        # Exact donor alternation so {sample}.h5ad cannot absorb the _qc suffix.
        sample = "|".join(re.escape(s) for s in ALL_DS_SAMPLES) or "__no_samples__",
    input:
        raw = lambda wc: _entry(wc.dataset)["raw_file"],
    output:
        h5ad       = _dp("{sample}.h5ad"),
        provenance = _pv("{sample}_ingest_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_{sample}_ingest.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 2,
    params:
        dataset             = lambda wc: wc.dataset,
        sample_id           = lambda wc: wc.sample,
        source              = lambda wc: _entry(wc.dataset)["source"],
        raw_file            = lambda wc: _entry(wc.dataset)["raw_file"],
        filter              = lambda wc: _entry(wc.dataset).get("filter", {}),
        sample_key          = lambda wc: _entry(wc.dataset).get("sample_key", "donor_id"),
        subsample_per_donor = lambda wc: _entry(wc.dataset).get("subsample_per_donor"),
        batch_key           = lambda wc: _entry(wc.dataset).get("batch_key", "sample_id"),
        gene_id_type        = lambda wc: _entry(wc.dataset).get("gene_id_type", "ensembl"),
        random_seed         = config["scrna"]["random_seed"],
    script:
        "../scripts/ingest_dataset.py"


rule ds_gene_symbol_map:
    """Dataset-scoped Ensembl→symbol map from the cohort's own .var (covers the
    external cohort's IDs, which the reference MyGene cache does not)."""
    output:
        cache      = _dp("gene_symbol_map.tsv"),
        provenance = _pv("gene_symbol_map_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_gene_symbol_map.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 8000,
        threads = 1,
    params:
        dataset  = lambda wc: wc.dataset,
        raw_file = lambda wc: _entry(wc.dataset)["raw_file"],
    script:
        "../scripts/dataset_gene_symbol_map.py"


rule ds_clinical_stub:
    """Schema-compatible blank clinical TSV so nerve_cell_subset's clinical join
    works for a cohort without GDC clinical metadata."""
    output:
        clinical   = _ext("gdc_clinical.tsv"),
        provenance = _pv("clinical_stub_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_clinical_stub.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    params:
        dataset = lambda wc: wc.dataset,
        samples = lambda wc: _sheet(wc.dataset),
    script:
        "../scripts/dataset_clinical_stub.py"


# =============================================================================
# Stage B — QC + integration + annotation + malignancy
# =============================================================================

rule ds_scrna_qc:
    """Per-sample QC (reuses scrna_qc.py); is_nerve_marker + gene_presence."""
    input:
        h5ad = _dp("{sample}.h5ad"),
    output:
        h5ad          = _dp("{sample}_qc.h5ad"),
        qc_metrics    = _tb("{sample}_qc_metrics.csv"),
        gene_presence = _tb("{sample}_gene_presence.csv"),
        provenance    = _pv("{sample}_qc_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_{sample}_qc.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        sample_id    = lambda wc: wc.sample,
        min_genes    = config["scrna"]["min_genes"],
        max_genes    = config["scrna"]["max_genes"],
        max_pct_mito = config["scrna"]["max_pct_mito"],
        markers      = config["nerve_cells"]["markers"],
    script:
        "../scripts/scrna_qc.py"


def _ds_qc_h5ads(wc):
    return expand(_dp("{sample}_qc.h5ad"), dataset=wc.dataset, sample=_sheet(wc.dataset))


rule ds_scrna_integration:
    """Train this cohort's OWN scVI model (reuses scrna_integration.py); the
    per-cohort batch_key comes from the dataset entry."""
    input:
        h5ads = _ds_qc_h5ads,
    output:
        model_dir   = directory(_md("scvi_model")),
        latent_h5ad = _dp("integrated_latent.h5ad"),
        provenance  = _pv("scvi_integration_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_scvi_integration.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["gpu_threads"],
    params:
        device      = config["hardware"]["device"],
        precision   = config["hardware"]["precision"],
        n_latent    = config["scrna"]["n_latent"],
        n_layers    = config["scrna"]["n_layers"],
        batch_key   = lambda wc: _entry(wc.dataset).get("batch_key", "sample_id"),
        random_seed = config["scrna"]["random_seed"],
        min_cells   = config["scrna"]["min_cells"],
        # Census cohorts ship raw integer counts in .X → no recovery. A future
        # log1p-sourced cohort can opt in via `counts_from_log1p: true` in its
        # dataset entry.
        counts_from_log1p = lambda wc: _entry(wc.dataset).get("counts_from_log1p", False),
    script:
        "../scripts/scrna_integration.py"


rule ds_scrna_annotate:
    """Marker-scored cell-type annotation on this cohort's latent (reuses
    scrna_annotate.py); dataset-scoped symbol map guarantees marker mapping."""
    input:
        latent_h5ad = _dp("integrated_latent.h5ad"),
        model_dir   = _md("scvi_model"),
        symbol_map  = _dp("gene_symbol_map.tsv"),
    output:
        h5ad       = _dp("annotated.h5ad"),
        summary    = _tb("annotation_summary.csv"),
        provenance = _pv("annotation_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_annotate.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        leiden_resolution = config["nerve_cells"]["leiden_resolution"],
        markers           = config["nerve_cells"]["markers"],
        random_seed       = config["scrna"]["random_seed"],
        census_version    = config["databases"]["cellxgene_census_version"],
    script:
        "../scripts/scrna_annotate.py"


rule ds_scrna_malignancy:
    """CNV-based malignancy labeling (reuses scrna_malignancy.py)."""
    input:
        h5ad = _dp("annotated.h5ad"),
    output:
        h5ad       = _dp("malignancy_labeled.h5ad"),
        cnv_plot   = _fg("cnv_heatmap.png"),
        provenance = _pv("malignancy_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_malignancy.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        random_seed    = config["scrna"]["random_seed"],
        census_version = config["databases"]["cellxgene_census_version"],
        cnv_chunk_size = config["scrna"]["cnv_chunk_size"],
    script:
        "../scripts/scrna_malignancy.py"


# =============================================================================
# Stage C — Nerve + immune + three-way interaction
# =============================================================================

rule ds_nerve_cell_subset:
    """Fresh nerve-cell subset (reuses nerve_cell_subset.py); freezes OFF for a
    replication cohort (no frozen_subset_file / split_assignments_file params)."""
    input:
        h5ad     = _dp("malignancy_labeled.h5ad"),
        clinical = _ext("gdc_clinical.tsv"),
    output:
        h5ad        = _dp("nerve_cells.h5ad"),
        umap        = _fg("nerve_cells_umap.png"),
        composition = _tb("nerve_cluster_composition.csv"),
        provenance  = _pv("nerve_subset_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_cell_subset.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        leiden_resolution = config["nerve_cells"]["leiden_resolution"],
        cell_types        = config["nerve_cells"]["cell_types"],
        markers           = config["nerve_cells"]["markers"],
        random_seed       = config["scrna"]["random_seed"],
        n_top_genes       = config["scrna"]["n_top_genes"],
        # freezes intentionally omitted → getattr(..., None) → fresh cut
    script:
        "../scripts/nerve_cell_subset.py"


rule ds_nerve_batch_qc:
    """Per-cluster sample purity on the nerve subset (reuses nerve_batch_qc.py)."""
    input:
        h5ad = _dp("nerve_cells.h5ad"),
    output:
        purity          = _tb("nerve_cluster_sample_purity.csv"),
        umap_main       = _fg("nerve_cells_umap_by_sample.png"),
        umap_per_sample = _fg("nerve_cells_umap_per_sample_panel.png"),
        provenance      = _pv("nerve_batch_qc_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_batch_qc.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    params:
        random_seed              = config["scrna"]["random_seed"],
        dominant_fraction_max    = config["nerve_cells"]["batch_qc"]["dominant_fraction_max"],
        min_contributing_fraction= config["nerve_cells"]["batch_qc"]["min_contributing_fraction"],
        min_contributing_samples = config["nerve_cells"]["batch_qc"]["min_contributing_samples"],
    script:
        "../scripts/nerve_batch_qc.py"


rule ds_nerve_cell_heterogeneity:
    """Nerve DE markers + offline GSEA (reuses nerve_cell_heterogeneity.py);
    shares the reference MSigDB gmt download (gene sets are cohort-independent)."""
    input:
        h5ad       = _dp("nerve_cells.h5ad"),
        symbol_map = _dp("gene_symbol_map.tsv"),
        gmt_bp     = os.path.join(config["msigdb"]["download_dir"],
                                  config["msigdb"]["collections"]["gmt_bp"]["filename"]),
        gmt_mf     = os.path.join(config["msigdb"]["download_dir"],
                                  config["msigdb"]["collections"]["gmt_mf"]["filename"]),
        gmt_manifest = os.path.join(config["msigdb"]["download_dir"], "msigdb_manifest.json"),
    output:
        markers    = _tb("nerve_cluster_markers.csv"),
        enrichment = _tb("nerve_enrichment.csv"),
        dotplot    = _fg("nerve_dotplot.png"),
        abundance  = _fg("nerve_abundance_heatmap.png"),
        provenance = _pv("heterogeneity_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_cell_heterogeneity.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        markers        = config["nerve_cells"]["markers"],
        random_seed    = config["scrna"]["random_seed"],
        msigdb_release = config["msigdb"]["release"],
        bp_label       = config["msigdb"]["collections"]["gmt_bp"]["library_label"],
        mf_label       = config["msigdb"]["collections"]["gmt_mf"]["library_label"],
        prerank        = config["msigdb"]["prerank"],
    script:
        "../scripts/nerve_cell_heterogeneity.py"


rule ds_nerve_cluster_annotations:
    """Per-cluster biological labels for a replication cohort (reuses
    nerve_cluster_annotations.py).

    Maps each nerve Leiden id to its dominant cell type + top markers + canonical
    module scores. Without this a cohort's interaction results can only be read as
    opaque `nerve_c{N}` ids, which is what the TME explorer notebooks exist to
    avoid. The reference cohort's twin is pinned; this one is freshly runnable.
    """
    input:
        h5ad       = _dp("nerve_cells.h5ad"),
        markers    = _tb("nerve_cluster_markers.csv"),
        symbol_map = _dp("gene_symbol_map.tsv"),
    output:
        annotations = _tb("nerve_cluster_annotations.csv"),
        provenance  = _pv("nerve_cluster_annotations_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_cluster_annotations.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        markers     = config["nerve_cells"]["markers"],
        random_seed = config["scrna"]["random_seed"],
    script:
        "../scripts/nerve_cluster_annotations.py"


rule ds_nerve_tumor_interaction:
    """Two-way tumor-nerve LR (reuses nerve_tumor_interaction.py)."""
    input:
        malig = _dp("malignancy_labeled.h5ad"),
        nerve = _dp("nerve_cells.h5ad"),
    output:
        lr_table   = _tb("nerve_tumor_interactions.csv"),
        top_pairs  = _tb("nerve_tumor_top_pairs.csv"),
        heatmap    = _fg("nerve_tumor_sig_heatmap.png"),
        dotplot    = _fg("nerve_tumor_dotplot.png"),
        provenance = _pv("nerve_tumor_interaction_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_tumor_interaction.log"),
    conda:
        "../envs/scrna.yaml",
    threads: config["resources"]["default_threads"],
    resources:
        mem_mb = config["resources"]["default_mem_mb"],
    params:
        random_seed = config["scrna"]["random_seed"],
        # Replication cohorts carry raw UMIs in .X end to end (correct for scVI,
        # which wants counts), but LIANA assumes log1p. Default True; a cohort
        # whose .X is already log1p sets `normalize_counts: false` in its entry.
        normalize_counts = lambda wc: _entry(wc.dataset).get("normalize_counts", True),
    script:
        "../scripts/nerve_tumor_interaction.py"


rule ds_immune_cell_subset:
    """Immune compartment subset + reclustering (reuses immune_cell_subset.py)."""
    input:
        h5ad = _dp("malignancy_labeled.h5ad"),
    output:
        h5ad        = _dp("immune_cells.h5ad"),
        umap        = _fg("immune_cells_umap.png"),
        composition = _tb("immune_cluster_composition.csv"),
        purity      = _tb("immune_cluster_sample_purity.csv"),
        provenance  = _pv("immune_cell_subset_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_immune_cell_subset.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        source_label              = config["immune_cells"]["source_label"],
        leiden_resolution         = config["immune_cells"]["leiden_resolution"],
        n_top_genes               = config["scrna"]["n_top_genes"],
        random_seed               = config["scrna"]["random_seed"],
        dominant_fraction_max     = config["immune_cells"]["batch_qc"]["dominant_fraction_max"],
        min_contributing_fraction = config["immune_cells"]["batch_qc"]["min_contributing_fraction"],
        min_contributing_samples  = config["immune_cells"]["batch_qc"]["min_contributing_samples"],
    script:
        "../scripts/immune_cell_subset.py"


rule ds_immune_cluster_annotations:
    """Immune subtype labels per cluster (reuses immune_cluster_annotations.py)."""
    input:
        h5ad = _dp("immune_cells.h5ad"),
    output:
        h5ad          = _dp("immune_cells_labeled.h5ad"),
        annotations   = _tb("immune_cluster_annotations.csv"),
        subtype_purity= _tb("immune_subtype_sample_purity.csv"),
        provenance    = _pv("immune_cluster_annotations_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_immune_cluster_annotations.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        subtype_markers           = config["immune_cells"]["subtype_markers"],
        random_seed               = config["scrna"]["random_seed"],
        dominant_fraction_max     = config["immune_cells"]["batch_qc"]["dominant_fraction_max"],
        min_contributing_fraction = config["immune_cells"]["batch_qc"]["min_contributing_fraction"],
        min_contributing_samples  = config["immune_cells"]["batch_qc"]["min_contributing_samples"],
    script:
        "../scripts/immune_cluster_annotations.py"


rule ds_nerve_tumor_immune_interaction:
    """Three-way nerve-tumor-immune LR (reuses nerve_tumor_immune_interaction.py)."""
    input:
        malig  = _dp("malignancy_labeled.h5ad"),
        nerve  = _dp("nerve_cells.h5ad"),
        immune = _dp("immune_cells_labeled.h5ad"),
    output:
        lr_table   = _tb("nerve_tumor_immune_interactions.csv"),
        top_pairs  = _tb("nerve_tumor_immune_top_pairs.csv"),
        heatmap    = _fg("nerve_tumor_immune_sig_heatmap.png"),
        dotplot    = _fg("nerve_tumor_immune_dotplot.png"),
        provenance = _pv("nerve_tumor_immune_interaction_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_tumor_immune_interaction.log"),
    conda:
        "../envs/scrna.yaml",
    threads: config["resources"]["default_threads"],
    resources:
        mem_mb = 64000,
    params:
        random_seed = config["scrna"]["random_seed"],
        # See ds_nerve_tumor_interaction — same raw-UMI vs log1p asymmetry.
        normalize_counts = lambda wc: _entry(wc.dataset).get("normalize_counts", True),
    script:
        "../scripts/nerve_tumor_immune_interaction.py"


rule ds_annotate_cluster_qc:
    """Join batch-QC verdict onto the interaction tables (reuses
    annotate_cluster_qc.py). Replication cohorts start with empty exclude lists
    (no artifacts diagnosed yet)."""
    input:
        purity          = _tb("nerve_cluster_sample_purity.csv"),
        immune_purity   = _tb("immune_subtype_sample_purity.csv"),
        enrichment      = _tb("nerve_enrichment.csv"),
        markers         = _tb("nerve_cluster_markers.csv"),
        interactions    = _tb("nerve_tumor_interactions.csv"),
        top_pairs       = _tb("nerve_tumor_top_pairs.csv"),
        tw_interactions = _tb("nerve_tumor_immune_interactions.csv"),
        tw_top_pairs    = _tb("nerve_tumor_immune_top_pairs.csv"),
    output:
        enrichment      = _tb("nerve_enrichment_with_qc.csv"),
        markers         = _tb("nerve_cluster_markers_with_qc.csv"),
        interactions    = _tb("nerve_tumor_interactions_with_qc.csv"),
        top_pairs       = _tb("nerve_tumor_top_pairs_with_qc.csv"),
        tw_interactions = _tb("nerve_tumor_immune_interactions_with_qc.csv"),
        tw_top_pairs    = _tb("nerve_tumor_immune_top_pairs_with_qc.csv"),
        provenance      = _pv("annotate_cluster_qc_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_annotate_cluster_qc.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    params:
        exclude_clusters        = [],   # fresh cut: no diagnosed artifacts yet
        immune_exclude_subtypes = [],
    script:
        "../scripts/annotate_cluster_qc.py"


# =============================================================================
# Stage C-bis — scANVI-v2 nerve re-train (full parity with reference cohort).
# Side-branch: validates the nerve latent; does NOT feed the three-way
# interaction (which uses the v1 nerve_cells.h5ad). Own scVI+scANVI train.
# =============================================================================

rule ds_nerve_assemble_counts:
    """Reconstruct raw-count nerve AnnData from per-sample QC files (reuses
    nerve_assemble_counts.py)."""
    input:
        qc_h5ads  = _ds_qc_h5ads,
        nerve_ref = _dp("nerve_cells.h5ad"),
    output:
        h5ad       = _dp("nerve_cells_counts.h5ad"),
        provenance = _pv("nerve_assemble_counts_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_assemble_counts.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        random_seed = config["scrna"]["random_seed"],
        # Census cohorts ship raw integer counts in QC .X → passthrough (no expm1
        # recovery, no nCount_SCT). A log1p/SCT cohort opts in via the dataset
        # entry, mirroring ds_scrna_integration.
        counts_from_log1p = lambda wc: _entry(wc.dataset).get("counts_from_log1p", False),
    script:
        "../scripts/nerve_assemble_counts.py"


rule ds_nerve_celltype_labels:
    """Marker-scored cell-type labels for scANVI anchoring (reuses
    nerve_celltype_labels.py)."""
    input:
        counts_h5ad = _dp("nerve_cells_counts.h5ad"),
    output:
        h5ad       = _dp("nerve_cells_counts_labeled.h5ad"),
        summary    = _tb("nerve_celltype_label_summary.csv"),
        provenance = _pv("nerve_celltype_labels_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_celltype_labels.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = config["resources"]["default_threads"],
    params:
        markers            = config["nerve_cells"]["markers"],
        unknown_percentile = config["nerve_scanvi"]["unknown_percentile"],
        random_seed        = config["scrna"]["random_seed"],
    script:
        "../scripts/nerve_celltype_labels.py"


rule ds_nerve_scanvi_retrain:
    """scVI baseline + scANVI fine-tune on the nerve subset (reuses
    nerve_scanvi_retrain.py). Own model; MPS-accelerated."""
    input:
        counts_h5ad = _dp("nerve_cells_counts_labeled.h5ad"),
    output:
        model_dir   = directory(_md("nerve_scanvi_model")),
        latent_h5ad = _dp("nerve_cells_v2.h5ad"),
        curves      = _fg("nerve_scanvi_training_curves.png"),
        provenance  = _pv("nerve_scanvi_retrain_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_scanvi_retrain.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 64000,
        threads = config["resources"]["gpu_threads"],
    params:
        device              = config["hardware"]["device"],
        precision           = config["hardware"]["precision"],
        n_latent            = config["scrna"]["n_latent"],
        n_layers            = config["scrna"]["n_layers"],
        random_seed         = config["scrna"]["random_seed"],
        batch_key           = config["nerve_scanvi"]["batch_key"],
        labels_key          = config["nerve_scanvi"]["labels_key"],
        unlabeled_category  = config["nerve_scanvi"]["unlabeled_category"],
        scvi_max_epochs     = config["nerve_scanvi"]["scvi_max_epochs"],
        scanvi_max_epochs   = config["nerve_scanvi"]["scanvi_max_epochs"],
        n_samples_per_label = config["nerve_scanvi"]["n_samples_per_label"],
    script:
        "../scripts/nerve_scanvi_retrain.py"


rule ds_nerve_batch_qc_v2:
    """Per-cluster sample purity on the v2 scANVI latent (reuses
    nerve_batch_qc_v2.py)."""
    input:
        h5ad = _dp("nerve_cells_v2.h5ad"),
    output:
        purity          = _tb("nerve_cluster_sample_purity_v2.csv"),
        umap_main       = _fg("nerve_cells_umap_by_sample_v2.png"),
        umap_per_sample = _fg("nerve_cells_umap_per_sample_panel_v2.png"),
        provenance      = _pv("nerve_batch_qc_v2_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_nerve_batch_qc_v2.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = config["resources"]["default_mem_mb"],
        threads = 1,
    params:
        random_seed              = config["scrna"]["random_seed"],
        leiden_resolution        = config["nerve_cells"]["leiden_resolution"],
        dominant_fraction_max    = config["nerve_cells"]["batch_qc"]["dominant_fraction_max"],
        min_contributing_fraction= config["nerve_cells"]["batch_qc"]["min_contributing_fraction"],
        min_contributing_samples = config["nerve_cells"]["batch_qc"]["min_contributing_samples"],
    script:
        "../scripts/nerve_batch_qc_v2.py"


# =============================================================================
# Stage D — Cross-cohort concordance (the replication deliverable)
# =============================================================================

rule ds_cohort_concordance:
    """Compare this cohort's three-way _with_qc table against the reference
    cohort's; Jaccard overlap + Spearman on shared significant LR pairs."""
    input:
        dataset   = _tb("nerve_tumor_immune_interactions_with_qc.csv"),
        reference = os.path.join(config["dirs"]["tables"],
                                 "nerve_tumor_immune_interactions_with_qc.csv"),
    output:
        summary      = _tb("cohort_concordance_summary.json"),
        shared_pairs = _tb("cohort_concordance_shared_pairs.csv"),
        figure       = _fg("cohort_concordance.png"),
        provenance   = _pv("cohort_concordance_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "{dataset}_cohort_concordance.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 8000,
        threads = 1,
    params:
        dataset  = lambda wc: wc.dataset,
        pval_max = 0.05,
    script:
        "../scripts/cohort_concordance.py"
