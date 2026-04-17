# Plan — scRNA-seq Analysis of Nerve Cell Heterogeneity in the GBM TME

## Context

The researcher wants to characterize heterogeneity of **nerve cells within the Glioblastoma (GBM) tumor microenvironment**, using 18 TCGA GBM scRNA-seq loom files at `data/raw/gdc_extract/` (Seurat 1000×1000 format per sample). The project already has a FAIR-compliant, Snakemake-first infrastructure with scvi-tools/scanpy envs pre-configured for Apple MPS. However, **existing rules cannot yet execute end-to-end** against this data and there is **no annotation or nerve-cell-subset rule**. This plan extends the existing scaffold with the minimum set of new rules needed to answer the biological question while preserving the FAIR4RS and MPS constraints in `CLAUDE.md`.

### Key findings from exploration

1. **Infrastructure is partly scaffolded** (Snakefile:20-25): QC (qc.smk:7-28), scVI integration (integration.smk:6-32), RNA velocity, FAIR validation, and Marimo notebook rules already exist. Conda envs `scrna.yaml` (scvi-tools 1.4.2, scanpy 1.12.1, cellxgene-census 1.17.0, scvelo 0.3.4), `notebooks.yaml`, `proteomics.yaml`, `base.yaml` all exist.
2. **Critical input-path gap**: `qc.smk:10` expects `data/raw/{sample}.h5ad`, but data is at `data/raw/gdc_extract/{UUID}/{file}.seurat.1000x1000.loom`. A **loom→h5ad staging rule must be added** and QC input rewired.
3. **Sample count discrepancy**: MANIFEST.txt has 18 samples, `config/config.yaml:93-110` has 17, and one UUID in config (line 103: `8c0685b3-521e-45a6-9677-f8e3f246e09b`) does NOT match the MANIFEST (`8c0685b3-521e-45a6-9677-f8e3f186e09b`) — a typo (`f246` vs `f186`). Fix before execution.
4. **No annotation / no nerve-cell logic exists**. No CNV/malignancy classification exists. Both are essential for the question.
5. **Data scale caveat**: the `.seurat.1000x1000.loom` naming strongly suggests subsampled data (≈1000 cells × 1000 genes per sample) — GDC's reduced distribution format. Across 17 samples, this is ~17,000 cells and 1,000 genes max. Nerve markers (e.g., `SYN1`, `SNAP25`, `RBFOX3`, `MAP2`, `MBP`, `PDGFRA`, `OLIG2`) must be confirmed present before downstream work — not all may survive the gene selection. The loom-convert rule must emit a gene-presence report.
6. **No clinical metadata** at `data/raw/gdc_extract/` — only UUIDs. GDC clinical endpoint can be queried optionally for IDH/MGMT/primary-vs-recurrent annotation.

---

## Approach — Extend the existing Snakemake DAG with 7 new rules

**Scoping decisions (confirmed with researcher):**
- **Nerve-cell scope**: broad — neurons (excitatory + inhibitory) + OPCs + oligodendrocytes + non-malignant astrocytes. inferCNV step is load-bearing.
- **Data resolution**: proceed with the 1000×1000 subsampled GDC loom as an exploratory pass; code is written so the full-resolution matrix can be swapped in later without rule changes.
- **Clinical metadata**: pulled from the GDC `/cases` endpoint in a dedicated rule (Step 1b) and joined into `obs` downstream.

All new rules use existing envs (`scrna.yaml` for compute; `notebooks.yaml` for Marimo). All follow FAIR4RS: each emits a `*_provenance.json` via `workflow/scripts/fair_utils.py::log_transformation`, logs to `logs/`, and stores outputs in `data/processed/`, `results/tables/`, or `results/figures/`.

### Step 0 — Fix config (`config/config.yaml`)
- Correct the UUID typo at line 103: `f246` → `f186`.
- Add the missing 18th sample (`903b7af4-...` is present; confirm 18 entries match MANIFEST.txt line-for-line).
- Add a new `nerve_cells` section scoped to **neurons + full neural/glial lineage** (excitatory + inhibitory neurons, OPCs, oligodendrocytes, non-malignant astrocytes). Fields: `markers` (per cell type), `reference_dataset`, `leiden_resolution` (default 1.0), `include_astrocytes: true`.
- Add a `gdc_api` section: `base_url: "https://api.gdc.cancer.gov"`, `fields: ["case_id","primary_diagnosis","disease_type","diagnoses.prior_malignancy","diagnoses.tumor_grade","diagnoses.primary_diagnosis","samples.tissue_type","diagnoses.treatments.treatment_outcome"]`.

### Step 1 — NEW rule `loom_to_h5ad` (workflow/rules/ingest.smk)
- **Per-sample** rule: converts each GDC loom into h5ad with standardized `obs.batch = {sample_uuid}`, `obs.sample_id`, and Ensembl gene IDs mapped into `var`.
- Input: `data/raw/gdc_extract/{sample}/*.seurat.1000x1000.loom` (glob via helper).
- Output: `data/processed/{sample}.h5ad` + `results/tables/{sample}_gene_presence.csv` + `provenance/{sample}_ingest_provenance.json`.
- Uses `loompy` + `anndata` (both in `scrna.yaml`). Emits warning if any canonical nerve-cell marker is absent.
- **This also resolves the qc.smk input path**: after this rule, QC finds the expected `data/processed/{sample}.h5ad`. Edit `qc.smk:10` to reference `data_processed` instead of `data_raw`.

### Step 1b — NEW rule `gdc_clinical_fetch` (workflow/rules/ingest.smk)
- Single rule (not per-sample): queries GDC `/cases` endpoint for all 18 UUIDs in `config.samples` via `requests` (already available in `scrna.yaml` transitively).
- Pulls fields listed in `config.gdc_api.fields` (primary diagnosis, disease type, tissue type, prior malignancy, tumor grade, treatment outcome). IDH/MGMT are extracted from `diagnoses.primary_diagnosis` when available — note these are not consistently curated across TCGA-GBM; the script must tolerate missing values and log which samples lack each field.
- Input: none (pure API call).
- Output: `data/external/gdc_clinical.tsv` (one row per UUID), `provenance/gdc_clinical_provenance.json` (records endpoint URL, query timestamp, field map).
- Downstream: loaded in Steps 4, 6, 7 and joined into `adata.obs` on `sample_id` so clusters can be stratified by IDH/MGMT/recurrence.

### Step 2 — EXISTING rule `scrna_qc` (qc.smk:7-28) — rewire input only
- Change `qc.smk:10` input path to `data/processed/{sample}.h5ad`.
- No other change; thresholds from `config.scrna` (min_genes=200, max_genes=6000, max_pct_mito=20) apply.

### Step 3 — EXISTING rule `scrna_integration` (integration.smk:6-32) — no change
- Trains scVI VAE with `batch_key="batch"` (populated in Step 1), MPS backend, float32, seed=0.
- Output: `data/processed/integrated_latent.h5ad`, `results/models/scvi_model/`.

### Step 4 — NEW rule `scrna_annotate` (workflow/rules/annotation.smk)
- Reference-based annotation via scANVI + CELLxGENE Census brain/cortex reference (cellxgene-census 1.17.0 already pinned in `scrna.yaml`; `config.databases.cellxgene_census_version: "stable"` at config.yaml:83).
- Uses the scVI model from Step 3 as the base, then fits scANVI with reference labels.
- Input: `data/processed/integrated_latent.h5ad`, `results/models/scvi_model/`.
- Output: `data/processed/annotated.h5ad` (with `obs.cell_type_predicted`, `obs.cell_type_confidence`), `results/tables/annotation_summary.csv`, `provenance/annotation_provenance.json`.
- **Cell-type taxonomy** follows Allen Brain Atlas conventions and is written into `adata.uns['annotation_schema']` with EDAM operation terms per `CLAUDE.md` FAIR rules.

### Step 5 — NEW rule `scrna_malignancy` (workflow/rules/annotation.smk)
- CNV-based tumor/normal separation via **inferCNV or CopyKAT** (pure-Python inferCNVpy is lightweight and compatible with MPS-agnostic CPU; add to `scrna.yaml` as a new dep in the same rule). Uses non-glial cells (T-cells, endothelial) as reference.
- Input: `data/processed/annotated.h5ad`.
- Output: `data/processed/malignancy_labeled.h5ad` (with `obs.is_malignant: bool`), `results/figures/cnv_heatmap.png`, `provenance/malignancy_provenance.json`.
- **This is essential**: GBM cells often appear "astrocyte-like" or "OPC-like" by expression alone. Without CNV, nerve-cell counts are contaminated by tumor cells.

### Step 6 — NEW rule `nerve_cell_subset` (workflow/rules/nerve_cells.smk)
- Filter `malignancy_labeled.h5ad` to **non-malignant cells in the confirmed broad nerve-cell scope**: excitatory neurons, inhibitory neurons, OPCs, oligodendrocytes, and non-malignant astrocytes (where `is_malignant=False` AND scANVI label is in the nerve-cell taxonomy set).
- Re-run HVG selection, PCA, neighbors, UMAP, and Leiden clustering at a finer resolution (config: `nerve_cells.leiden_resolution`, default 1.0) in the nerve-cell subspace.
- Join clinical metadata from `data/external/gdc_clinical.tsv` onto `obs` so downstream heterogeneity can be stratified by IDH/MGMT/recurrence.
- Input: `data/processed/malignancy_labeled.h5ad`, `data/external/gdc_clinical.tsv`.
- Output: `data/processed/nerve_cells.h5ad`, `results/figures/nerve_cells_umap.png`, `results/tables/nerve_cluster_composition.csv`, `provenance/nerve_subset_provenance.json`.

### Step 7 — NEW rule `nerve_cell_heterogeneity` (workflow/rules/nerve_cells.smk)
- Core scientific output. For the nerve-cell subspace:
  - **Differential expression** (`scanpy.tl.rank_genes_groups`, Wilcoxon) per leiden cluster → top marker genes.
  - **Gene-set enrichment** via `gseapy` (add to `scrna.yaml`) against neuronal GO terms & MSigDB Hallmarks.
  - **Per-sample abundance**: cluster × sample heatmap (inter-patient heterogeneity).
  - **Marker dot plot** for canonical nerve-cell markers from `config.nerve_cells.markers`.
- Input: `data/processed/nerve_cells.h5ad`.
- Output: `results/tables/nerve_cluster_markers.csv`, `results/tables/nerve_enrichment.csv`, `results/figures/nerve_dotplot.png`, `results/figures/nerve_abundance_heatmap.png`, `provenance/heterogeneity_provenance.json`.

### Step 8 — Update `rule all` (Snakefile:43-49) and Marimo notebook (`notebooks/01_explore_gbm_data.py`)
- Add nerve-cell artifacts to `rule all` so Goal-backward verification catches missing outputs.
- Extend the existing Marimo notebook with a new section rendering the nerve-cell UMAP, marker dot plot, DE table (via DuckDB), and per-cluster GSEA summary.

---

## Critical files to modify / create

| Path | Action |
|---|---|
| `config/config.yaml` | Fix UUID typo line 103; add 18th sample; add `nerve_cells` + `gdc_api` sections |
| `workflow/rules/ingest.smk` | **NEW** — loom→h5ad staging + GDC clinical fetch |
| `workflow/rules/annotation.smk` | **NEW** — scANVI annotation + inferCNVpy malignancy |
| `workflow/rules/nerve_cells.smk` | **NEW** — subset + heterogeneity rules |
| `workflow/rules/qc.smk` | Edit line 10 input path → `data/processed/{sample}.h5ad` |
| `workflow/envs/scrna.yaml` | Add `infercnvpy`, `gseapy`, `requests` deps |
| `workflow/scripts/loom_to_h5ad.py` | **NEW** |
| `workflow/scripts/gdc_clinical_fetch.py` | **NEW** |
| `workflow/scripts/scrna_annotate.py` | **NEW** |
| `workflow/scripts/scrna_malignancy.py` | **NEW** |
| `workflow/scripts/nerve_cell_subset.py` | **NEW** |
| `workflow/scripts/nerve_cell_heterogeneity.py` | **NEW** |
| `Snakefile` | Add `include: "workflow/rules/ingest.smk"`, `annotation.smk`, `nerve_cells.smk`; extend `rule all` with nerve-cell + clinical artifacts |
| `notebooks/01_explore_gbm_data.py` | Append nerve-cell heterogeneity section |
| `CHANGELOG.md` | Log phase per `CLAUDE.md` requirement |

### Reused utilities (DO NOT duplicate)
- `workflow/scripts/fair_utils.py::log_transformation`, `stamp_artifact`, `sha256_file` — provenance.
- `workflow/rules/common.smk::p()` — path helper used by all rules.
- Existing `scrna_qc.py` threshold-logging pattern — reuse for new rules.

---

## Verification (end-to-end)

Run in order from project root, using `claude_science/bin/python3 -m snakemake` per memory note:

1. **Dry run (DAG validity)**:
   `claude_science/bin/python3 -m snakemake --use-conda --cores 8 -n`
   Expect no missing-input errors; 18 samples × (ingest, qc) + 1 × (integration, annotation, malignancy, subset, heterogeneity) + FAIR/notebook.

2. **Single-sample smoke test** (fast):
   `claude_science/bin/python3 -m snakemake --use-conda --cores 4 data/processed/06820e2c-9eb7-4e71-a1c3-976d561e659d_qc.h5ad`
   Verifies Step 1 + 2 on one sample before spending MPS time.

3. **Full pipeline**:
   `claude_science/bin/python3 -m snakemake --use-conda --cores all`
   Expect (per `rule all`): `fair_validation_report.json`, `qc_summary.csv`, `integrated_latent.h5ad`, plus new `nerve_cluster_markers.csv`, `nerve_dotplot.png`, `nerve_abundance_heatmap.png`.

4. **Goal-backward checks** (biology-facing):
   - Confirm `nerve_cluster_markers.csv` top markers per cluster include canonical neuronal genes (`SYN1`, `SNAP25`, `STMN2`, `RBFOX3`) for at least one cluster.
   - Confirm `annotation_summary.csv` reports non-zero neuron/OPC/oligodendrocyte counts.
   - Confirm `is_malignant=False` count in nerve-cell subset matches `nerve_cells.h5ad.n_obs`.
   - Run `claude_science/bin/python3 -m snakemake --report results/snakemake_report.html` for reproducibility trail.

5. **FAIR validation**:
   Confirm each new `*_provenance.json` is present and `fair_validation_report.json` shows no errors.

6. **Marimo interactive check**:
   `claude_science/bin/python3 -m marimo edit notebooks/01_explore_gbm_data.py` — verify new nerve-cell section renders UMAP + dot plot correctly.

---

## Caveats / residual risk to monitor during execution

1. **Subsampled 1000-gene matrix may drop key nerve markers.** Step 1's gene-presence CSV is the first checkpoint — if more than a couple of canonical markers (`SYN1`, `SNAP25`, `RBFOX3`, `MAP2`, `MBP`, `PDGFRA`, `OLIG2`) are missing across most samples, stop and reconsider fetching full-resolution loom files before running integration/annotation.
2. **scANVI transfer quality on 1k genes.** CELLxGENE brain references are trained on ~20k genes. Transfer accuracy may degrade — `obs.cell_type_confidence` should be inspected and low-confidence cells flagged in `annotation_summary.csv`.
3. **inferCNVpy requires contiguous genomic ordering.** Gene coordinates must be attached during Step 1 ingest (from Ensembl release 113 per `config.databases.ensembl_release`). If the loom lacks chromosome/position metadata, Step 5 will need an Ensembl coordinate lookup.
4. **TCGA-GBM IDH/MGMT curation is sparse.** Expect many NA values in `gdc_clinical.tsv`; stratification plots should include an "unknown" bucket rather than dropping samples.
