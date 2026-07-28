# GBM Nerve–Tumor–Immune Single-Cell Analysis

An autonomous, FAIR-compliant single-cell RNA-seq pipeline for characterizing nerve-cell heterogeneity in the Glioblastoma (GBM) tumor microenvironment (TME) — and the ligand–receptor crosstalk that connects the tumor and immune compartments to it. Built on Apple M4 Max hardware with Snakemake, scvi-tools, and scanpy, and tested for replication across two independent GBM cohorts.

---

## Scientific Questions

**1. Composition.** What subtypes of nerve cells (neurons, oligodendrocytes, OPCs, astrocytes) are present in the GBM TME, how heterogeneous are they across patients, and which molecular programs define each subtype?

**2. Crosstalk.** Which ligand–receptor interactions link malignant and immune cells to the nerve-cell compartment — the tumor→immune→nerve axis?

**3. Replication.** Do those interactions reproduce in an independent GBM cohort, or are they specific to the discovery data?

The analysis addresses these by:
1. Ingesting each cohort into a common AnnData format (GDC loom, or `.h5ad` from CELLxGENE Census)
2. Integrating samples with scVI to remove batch effects
3. Annotating cell types using canonical marker gene scoring
4. Separating malignant GBM cells from bona fide normal nerve cells via CNV inference
5. Subclustering the nerve-cell compartment and running differential expression + gene set enrichment
6. Subclustering the immune compartment into microglia / TAM / T / NK / dendritic subtypes
7. Scoring two-way (tumor→nerve) and three-way (tumor→immune→nerve) ligand–receptor interactions with LIANA
8. Re-running the whole chain on a replication cohort and scoring concordance against the reference

### Cohorts

| | Reference (discovery) | Replication |
|---|---|---|
| **Source** | TCGA-GBM via GDC | CELLxGENE Census (10x, primary data only) |
| **Scale** | 17 samples, subsampled loom | 1.29M cells × 61,497 genes, 174 donors |
| **Scope run** | full cohort | single largest study (`56c4912d`), capped at 5,000 cells/donor |
| **Matrix** | SCT log1p (counts recovered for scVI) | genuine raw UMIs |
| **Status** | frozen at **v1.3.0**, structurally pinned | active |

---

## Project Structure

```
Nerve_Analysis_TCGA_GBM/
│
├── CLAUDE.md                        # Project constitution — rules for all agentic work
├── README.md                        # This file
├── execution_instructions.md        # Step-by-step pipeline execution guide
├── CHANGELOG.md                     # Persistent lab notebook (all decisions logged here)
│
├── Snakefile                        # Pipeline entry point — defines rule all + imports
├── config/
│   └── config.yaml                  # All parameters, paths, hardware settings (no hardcoding)
│
├── workflow/
│   ├── rules/
│   │   ├── common.smk               # Path helper p(), BASELINE_PINNED, pinned_target()
│   │   ├── fair.smk                 # FAIR validation, report, provenance + reference freezes
│   │   ├── ingest.smk               # loom_to_h5ad (per-sample) + gdc_clinical_fetch
│   │   ├── qc.smk                   # scrna_qc (per-sample) + scrna_qc_report
│   │   ├── integration.smk          # scrna_integration (scVI VAE) + scrna_velocity
│   │   ├── annotation.smk           # scrna_annotate + scrna_malignancy (CNV scoring)
│   │   ├── nerve_cells.smk          # Nerve subset, heterogeneity, batch QC, scANVI-v2 retrain
│   │   ├── immune.smk               # Immune subset, subtype annotation, 3-way interaction
│   │   ├── datasets.smk             # Cohort-namespaced ds_* replication chain (Stage A→D)
│   │   ├── proteomics.smk           # AlphaPept MS rules (optional)
│   │   └── notebooks.smk            # Marimo notebook batch export
│   │
│   ├── scripts/                     # One module per biological operation (see below)
│   │   ├── fair_utils.py            # FAIR provenance utilities (UUID, SHA256, stamping)
│   │   ├── counts_utils.py          # Scale detection + log1p→counts recovery guards
│   │   ├── ingest_dataset.py        # Multi-source cohort loader (loom | h5ad | mtx | census)
│   │   ├── scrna_integration.py     # scVI VAE training (MPS-accelerated)
│   │   ├── scrna_malignancy.py      # CNV sliding-window scoring → is_malignant label
│   │   ├── nerve_cell_subset.py     # Filter + re-cluster nerve cells; join clinical metadata
│   │   ├── nerve_scanvi_retrain.py  # scANVI-v2 label-anchored latent space
│   │   ├── immune_cell_subset.py    # Immune compartment subset + re-clustering
│   │   ├── immune_cluster_annotations.py   # Microglia/TAM/T/NK/DC subtype calling
│   │   ├── nerve_tumor_interaction.py      # Two-way tumor→nerve LIANA scoring
│   │   ├── nerve_tumor_immune_interaction.py  # Three-way tumor→immune→nerve LIANA scoring
│   │   ├── cohort_concordance.py    # Reference vs replication agreement scoring
│   │   ├── freeze_pinned_reference.py / verify_pinned_reference.py   # v1.3.0 pin enforcement
│   │   └── …                        # QC, annotation, enrichment, velocity, proteomics
│   │
│   └── envs/
│       ├── scrna.yaml               # scRNA-seq conda env (scanpy, scvi-tools, infercnvpy, gseapy)
│       ├── notebooks.yaml           # Marimo/DuckDB env (isolated from scrna)
│       ├── proteomics.yaml          # AlphaPept env (isolated from bio)
│       └── base.yaml                # Minimal env for utility rules
│
├── scripts/                         # Out-of-band helpers (not pipeline rules)
│   ├── gbm_cellxgene_census_pull.py     # Census → cellxgene_data/gbm_10x_raw.h5ad
│   └── derive_dataset_sample_sheet.py   # Auto-derive per-cohort samples.txt
│
├── data/
│   ├── raw/gdc_extract/             # 17 GDC loom files ({UUID}/{file}.seurat.1000x1000.loom)
│   ├── processed/                   # Snakemake-produced AnnData files (h5ad)
│   └── external/                    # gdc_clinical.tsv (fetched from GDC API)
├── cellxgene_data/                  # Census replication cohort (raw UMI h5ad)
│
├── results/
│   ├── figures/                     # PNG + exported notebook HTML
│   │   └── {cohort}/                # Cohort-namespaced replication outputs
│   ├── tables/                      # CSV/JSON outputs (markers, enrichment, LR, concordance)
│   │   └── {cohort}/
│   └── models/scvi_model/           # Saved scVI / scANVI checkpoints
│
├── notebooks/                       # Marimo reactive explorers (01–05)
├── markdowns/                       # Plans, assessments, blocker write-ups
├── provenance/                      # FAIR provenance JSON per rule output
├── logs/                            # Snakemake rule logs
│
├── claude_science/                  # Python 3.12 virtual environment
└── requirements.txt                 # Top-level package list (pinned versions in envs/*.yaml)
```

### Marimo explorers

| Notebook | Scope |
|---|---|
| `01_explore_gbm_data.py` | QC explorer + nerve-cell viewer |
| `02_nerve_enrichment_explorer.py` | Nerve cluster DE / GSEA results |
| `03_nerve_tumor_immune_explorer.py` | Three-way interaction results |
| `04_tme_nerve_immune_explorer.py` | Full TME view across compartments |
| `05_census_nerve_immune_explorer.py` | Replication-cohort equivalent (per cohort) |

---

## Environment Setup

### Prerequisites

- macOS with Apple Silicon (M4 Max recommended; MPS GPU backend required for scVI training)
- [Conda/Mamba](https://github.com/conda-forge/miniforge) for environment management
- Python 3.12 virtual environment at `claude_science/`

### Activate the virtual environment

```bash
source claude_science/bin/activate
```

### Install top-level dependencies (first time only)

```bash
pip install -r requirements.txt
```

Conda environments for each Snakemake rule are built automatically on first run via `--use-conda`. This includes:

| Environment | Key packages |
|---|---|
| `scrna.yaml` | scanpy 1.12.1, scvi-tools 1.4.2, infercnvpy 0.4.3, gseapy 1.1.3, cellxgene-census 1.17.0 |
| `notebooks.yaml` | marimo 0.23.1, duckdb 1.5.2, scanpy 1.10.4 |
| `proteomics.yaml` | alphapept 0.5.3, numba |

---

## Running the Analysis

> Full step-by-step instructions with troubleshooting are in [`execution_instructions.md`](execution_instructions.md).

> ⚠️ **The v1.3.0 reference cohort is structurally pinned** (`baseline.pinned: true` in
> `config/config.yaml`). The 10 rules that would rewrite reference artifacts are *not defined*
> at parse time, so they cannot be scheduled by any invocation — see
> [Pinned reference](#pinned-reference-v130) before running anything that touches the
> reference cohort.

### 1. Validate the pipeline (dry run)

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores 8 -n
```

Job count depends on which cohorts are enabled in `config.yaml` (`samples:` and `datasets:`) and
on which pinned artifacts are already present on disk. Confirm the DAG resolves with no errors
and that no rule you did not intend to re-run appears in the plan.

### 2. Single-sample smoke test

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores 4 \
  data/processed/06820e2c-9eb7-4e71-a1c3-976d561e659d_qc.h5ad
```

Inspect the gene presence report before running the full pipeline:

```bash
cat results/tables/06820e2c-9eb7-4e71-a1c3-976d561e659d_gene_presence.csv
```

### 3. Full pipeline

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all
```

### 4. Replication cohort only

The replication chain is cohort-namespaced under `ds_*` rules; its terminal target pulls the
whole Stage A→D chain for one cohort without touching reference outputs:

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  results/tables/gbm_cellxgene_56c4912d/cohort_concordance_summary.json
```

### 5. Interactive exploration

```bash
claude_science/bin/python3 -m marimo edit notebooks/01_explore_gbm_data.py
```

---

## Pipeline Overview

```
GDC loom files (17 samples)                    Census h5ad (replication cohort)
        │                                                    │
        ▼                                                    ▼
  loom_to_h5ad ───────────── gdc_clinical_fetch       ds_ingest_dataset
  (per sample)                (GDC REST API)          (+ gene map, clinical stub)
        │                            │                       │
        ▼                            │                       ▼
   scrna_qc                          │                  ds_scrna_qc
  (per sample)                       │                       │
        │                            │                       ▼
        ▼                            │              ds_scrna_integration
 scrna_integration                   │                       │
 (scVI VAE, MPS)                     │                       ▼
        │                            │                ds_scrna_annotate
        ▼                            │                       │
 scrna_annotate                      │                       ▼
 (Leiden + marker scoring)           │               ds_scrna_malignancy
        │                            │                       │
        ▼                            │                       │
 scrna_malignancy ◄──────────────────┘                       │
 (CNV sliding-window → is_malignant)                         │
        │                                                    │
        ├──────────────────┬─────────────────┐               │
        ▼                  ▼                 │               │
 nerve_cell_subset   immune_cell_subset      │               │
 (non-malignant       (immune blob →         │               │
  neural/glial)        re-cluster)           │               │
        │                  │                 │               │
        ▼                  ▼                 │               │
 nerve_cell_          immune_cluster_        │               │
 heterogeneity        annotations            │               │
 (DE·GSEA·dotplot)    (microglia/TAM/T/…)    │               │
        │                  │                 │               │
        └────────┬─────────┘                 │               │
                 ▼                           ▼               │
   nerve_tumor_immune_interaction    nerve_tumor_interaction  │
   (three-way LIANA LR scoring)      (two-way tumor→nerve)    │
                 │                           │               │
                 └─────────┬─────────────────┘               │
                           ▼                                 │
                  annotate_cluster_qc                        │
                  (*_with_qc.csv tables)                     │
                           │                                 │
                           └──────────► ds_cohort_concordance ◄┘
                                        (reference vs replication)
```

A parallel **scANVI-v2 branch** (`nerve_celltype_labels` → `nerve_scanvi_retrain` →
`nerve_batch_qc_v2`) re-anchors the nerve latent space on marker-derived cell-type labels to
disentangle patient identity from cell-type biology. It runs for both cohorts and is not pulled
by the concordance target.

---

## Key Outputs

| File | Description |
|---|---|
| `data/external/gdc_clinical.tsv` | GDC clinical metadata per sample (primary diagnosis, tissue type, demographics) |
| `results/tables/annotation_summary.csv` | Cell-type counts and annotation confidence per cluster |
| `results/figures/cnv_heatmap.png` | CNV sliding-window heatmap — malignant/normal boundary |
| ~~`data/processed/nerve_cells.h5ad`~~ | Non-malignant nerve-cell AnnData — **lost 2026-07-21, not reproducible** for the reference cohort; exists per replication cohort under `data/processed/{cohort}/` |
| ~~`results/figures/nerve_cells_umap.png`~~ | Nerve-subspace UMAP — **lost with the above**; exists under `results/figures/{cohort}/` |
| `results/tables/nerve_cluster_markers.csv` | Wilcoxon DE markers per nerve-cell cluster (padj < 0.05) |
| `results/tables/nerve_enrichment.csv` | GO Biological Process / Molecular Function enrichment per cluster |
| `results/figures/nerve_dotplot.png` | Canonical marker expression dot plot across clusters |
| `results/figures/nerve_abundance_heatmap.png` | Per-sample cluster proportion heatmap |

**Immune compartment and crosstalk**

| File | Description |
|---|---|
| `results/figures/immune_cells_umap.png` | UMAP of the immune subspace by cluster and subtype |
| `results/tables/immune_cluster_annotations.csv` | Microglia / TAM / T / NK / dendritic subtype calls per cluster |
| `results/tables/nerve_tumor_interactions.csv` | Two-way tumor→nerve ligand–receptor pairs (LIANA) |
| `results/tables/nerve_tumor_immune_interactions.csv` | Three-way tumor→immune→nerve ligand–receptor pairs |
| `results/tables/*_with_qc.csv` | The above, joined to cluster-QC flags — **use these for interpretation** |
| `results/tables/nerve_tumor_immune_top_pairs_with_qc.csv` | Ranked top interaction pairs after QC filtering |

**Replication and reproducibility**

| File | Description |
|---|---|
| `results/tables/{cohort}/cohort_concordance_summary.json` | Reference vs replication agreement verdict |
| `results/tables/nerve_cluster_sample_purity_v2.csv` | Per-cluster donor purity after scANVI-v2 retrain |
| `results/figures/nerve_scanvi_training_curves.png` | scANVI train/val loss curves |
| `results/snakemake_report.html` | Full reproducibility report with DAG, rule stats, and provenance |
| `provenance/*.json` | FAIR provenance record for every rule output (artifact ID, SHA256, parameters) |

Replication-cohort outputs mirror the reference paths under a cohort namespace —
`results/tables/gbm_cellxgene_56c4912d/…` and `results/figures/gbm_cellxgene_56c4912d/…`.

---

## Configuration

All analysis parameters live in `config/config.yaml`. Key sections:

| Section | What it controls |
|---|---|
| `scrna` | QC thresholds, HVG count, scVI latent dims, batch key, random seed |
| `hardware` | MPS device, float32 precision, memory watermark ratio |
| `nerve_cells` | Leiden resolution, cell-type scope, canonical marker gene lists, batch-QC thresholds |
| `immune_cells` | Immune Leiden resolution, source label, subtype marker panels, QC exclusions |
| `nerve_scanvi` | scANVI-v2 batch/label keys, epoch budgets, per-label batch balance |
| `datasets` | Replication cohorts — source type, raw file, filters, sample key, per-donor cap |
| `baseline` | Semver of the citable result set + the v1.3.0 structural pin and its artifact list |
| `msigdb` | Gene-set libraries used for enrichment |
| `databases` | CELLxGENE Census version, Ensembl release |
| `gdc_api` | GDC REST API base URL and clinical fields to fetch |
| `dirs` | All input/output directory paths (no hardcoding elsewhere) |

To adjust the nerve-cell subclustering resolution:
```yaml
nerve_cells:
  leiden_resolution: 1.5   # increase for finer clusters
```

Then force re-run:
```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all --forcerun nerve_cell_subset
```

---

## Pinned reference (v1.3.0)

`data/processed/nerve_cells.h5ad` was destroyed by a failed `nerve_cell_subset` job on
2026-07-21 and **cannot be reproduced** — a retrained latent space re-clusters, so the frozen
`cl15` split no longer maps back to a single source cluster. The researcher decision was to pin
the surviving tables and never regenerate them.

That decision is enforced structurally, not by convention:

- `baseline.pinned: true` in `config/config.yaml`
- `workflow/rules/common.smk` exposes `BASELINE_PINNED`; the 10 affected rules sit behind
  `if not BASELINE_PINNED:` guards in `nerve_cells.smk` and `immune.smk`
- **A rule that is not defined cannot be scheduled** — so the pin holds for every invocation,
  including a bare `snakemake`, with no need to remember a command-line flag
- `pinned_target()` requests a frozen artifact only if it still exists on disk, so `rule all`
  never demands an output that was lost
- `freeze_pinned_reference` / `verify_pinned_reference` (in `fair.smk`) record and re-check the
  SHA256 of every pinned artifact

The scANVI-v2 chain is deliberately **not** pinned. To restore full reproducibility of the
reference cohort from scratch, set `baseline.pinned: false` — and accept that the resulting
clusters will differ from the published v1.3.0 labels.

---

## FAIR Compliance

This project implements FAIR (Findable, Accessible, Interoperable, Reusable) principles for all outputs:

- **Findable** — every artifact has a UUID5 persistent identifier and SHA256 hash in its `provenance/*.json`
- **Accessible** — all data accessed via standard APIs (GDC REST, CELLxGENE Census); no hardcoded local paths
- **Interoperable** — all matrices stored as AnnData (`.h5ad`) with Ensembl gene IDs and EDAM ontology operation codes
- **Reusable** — all transformations logged with tool versions and parameters; reproducible via `snakemake --report`

---

## Hardware Notes (Apple M4 Max)

- scVI training uses the **MPS backend** (`accelerator="mps", devices=1`) — do not change to CUDA or CPU unless MPS is unavailable
- Use **float32** (not float16) — MPS float16 is often slower on Apple Silicon
- Set `PYTORCH_MPS_HIGH_WATERMARK_RATIO=0.0` for long training runs (configured automatically in `scrna_integration.py`)
- If an MPS operation silently falls back to CPU, a `[MPS-ALERT]` warning is written to the rule log

---

## Caveats

1. **Subsampled reference data** — GDC distributes `seurat.1000x1000.loom` files (~1000 cells × 1000 genes per sample). Rare nerve cell populations may be under-represented. Check `*_gene_presence.csv` files after the smoke test.
2. **Malignancy classification** — uses sliding-window CNV scoring with T cells and endothelial cells as a reference. Confidence is proportional to the number of reference cells available per sample.
3. **GDC clinical metadata** — IDH and MGMT status are not consistently curated in TCGA-GBM. Expect `NA` values; all plots include an "unknown" category.
4. **GSEA** — requires internet access to query the Enrichr API. If offline, enrichment CSV will be empty but the pipeline continues without error.
5. **The reference cohort cannot be regenerated** — see [Pinned reference](#pinned-reference-v130). Reference nerve artifacts are frozen tables, not reproducible outputs.
6. **Matrix scale differs between cohorts** — the reference `.X` is SCT log1p (counts recovered for scVI via `counts_from_log1p`); the Census cohort is genuine raw UMIs. LIANA consumers normalize per-cohort via the `normalize_counts` parameter and hard-fail if the matrix is off-scale. Do not assume a shared scale when adding a new consumer — check `counts_utils.is_log1p_scale()`.
7. **Replication cohort is scoped, not complete** — only the single largest Census study (`56c4912d`) is run, capped at 5,000 cells per donor. Concordance is therefore a test against one independent study, not against all 174 donors.
8. **A substantial minority of clusters fail purity QC** — in the scANVI-v2 purity tables, the reference cohort passes 16 of 26 nerve clusters and the replication cohort 24 of 33. Failures are predominantly small, patient-dominated clusters. Always interpret from the `*_with_qc.csv` tables rather than the raw interaction tables.