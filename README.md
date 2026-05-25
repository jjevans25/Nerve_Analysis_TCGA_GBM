# TCGA GBM — Single-Cell Nerve Cell Heterogeneity Analysis

An autonomous, FAIR-compliant single-cell RNA-seq pipeline for characterizing the heterogeneity of nerve cells within the Glioblastoma (GBM) tumor microenvironment (TME), built on Apple M4 Max hardware with Snakemake, scvi-tools, and scanpy.

---

## Scientific Question

**What subtypes of nerve cells (neurons, oligodendrocytes, OPCs, astrocytes) are present in the GBM TME, how heterogeneous are they across patients, and which molecular programs define each subtype?**

The analysis addresses this by:
1. Converting 17 TCGA GBM loom files into a common AnnData format
2. Integrating samples with scVI to remove batch effects
3. Annotating cell types using canonical marker gene scoring
4. Separating malignant GBM cells from bona fide normal nerve cells via CNV inference
5. Subclustering the nerve-cell compartment and running differential expression + gene set enrichment

---

## Project Structure

```
Claude_TCGA_GBM/
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
│   │   ├── common.smk               # Shared path helper p()
│   │   ├── fair.smk                 # FAIR validation + Snakemake report rules
│   │   ├── ingest.smk               # loom_to_h5ad (per-sample) + gdc_clinical_fetch
│   │   ├── qc.smk                   # scrna_qc (per-sample) + scrna_qc_report
│   │   ├── integration.smk          # scrna_integration (scVI VAE) + scrna_velocity
│   │   ├── annotation.smk           # scrna_annotate + scrna_malignancy (CNV scoring)
│   │   ├── nerve_cells.smk          # nerve_cell_subset + nerve_cell_heterogeneity
│   │   ├── proteomics.smk           # AlphaPept MS rules (optional)
│   │   └── notebooks.smk            # Marimo notebook batch export
│   │
│   ├── scripts/
│   │   ├── fair_utils.py            # FAIR provenance utilities (UUID, SHA256, stamping)
│   │   ├── loom_to_h5ad.py          # GDC loom → AnnData conversion + gene presence check
│   │   ├── gdc_clinical_fetch.py    # GDC REST API clinical metadata fetch
│   │   ├── scrna_qc.py              # Per-sample QC filtering
│   │   ├── scrna_qc_report.py       # Aggregate QC summary
│   │   ├── scrna_integration.py     # scVI VAE training (MPS-accelerated)
│   │   ├── scrna_annotate.py        # Leiden clustering + marker-based annotation
│   │   ├── scrna_malignancy.py      # CNV sliding-window scoring → is_malignant label
│   │   ├── nerve_cell_subset.py     # Filter + re-cluster nerve cells; join clinical metadata
│   │   ├── nerve_cell_heterogeneity.py  # DE, GSEA, dot plot, abundance heatmap
│   │   └── scrna_velocity.py        # RNA velocity (scVelo)
│   │
│   └── envs/
│       ├── scrna.yaml               # scRNA-seq conda env (scanpy, scvi-tools, infercnvpy, gseapy)
│       ├── notebooks.yaml           # Marimo/DuckDB env (isolated from scrna)
│       ├── proteomics.yaml          # AlphaPept env (isolated from bio)
│       └── base.yaml                # Minimal env for utility rules
│
├── data/
│   ├── raw/
│   │   └── gdc_extract/             # 17 GDC loom files ({UUID}/{file}.seurat.1000x1000.loom)
│   ├── processed/                   # Snakemake-produced AnnData files (h5ad)
│   └── external/                    # gdc_clinical.tsv (fetched from GDC API)
│
├── results/
│   ├── figures/                     # PNG outputs (UMAP, CNV heatmap, dot plot, abundance)
│   ├── tables/                      # CSV outputs (QC summary, markers, enrichment)
│   └── models/
│       └── scvi_model/              # Saved scVI model checkpoint
│
├── notebooks/
│   └── 01_explore_gbm_data.py       # Marimo reactive notebook (QC explorer + nerve-cell viewer)
│
├── provenance/                      # FAIR provenance JSON per rule output
├── logs/                            # Snakemake rule logs
│
├── claude_science/                  # Python 3.12 virtual environment
└── requirements.txt                 # Top-level package list (pinned versions in envs/*.yaml)
```

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

### 1. Validate the pipeline (dry run)

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores 8 -n
```

Expect **44 jobs** across 12 rule types with no errors.

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

### 4. Interactive exploration

```bash
claude_science/bin/python3 -m marimo edit notebooks/01_explore_gbm_data.py
```

---

## Pipeline Overview

```
GDC loom files (17 samples)
        │
        ▼
  loom_to_h5ad ──────────────────────────────────── gdc_clinical_fetch
  (per sample)                                       (GDC REST API)
        │                                                    │
        ▼                                                    │
   scrna_qc                                                  │
  (per sample)                                               │
        │                                                    │
        ▼                                                    │
 scrna_integration                                           │
 (scVI VAE, MPS)                                            │
        │                                                    │
        ▼                                                    │
 scrna_annotate                                              │
 (Leiden + marker scoring)                                   │
        │                                                    │
        ▼                                                    │
 scrna_malignancy ◄──────────────────────────────────────────┘
 (CNV sliding-window → is_malignant)
        │
        ▼
 nerve_cell_subset
 (filter non-malignant neural/glial; re-cluster; join clinical)
        │
        ▼
 nerve_cell_heterogeneity
 (DE · GSEA · dot plot · abundance heatmap)
```

---

## Key Outputs

| File | Description |
|---|---|
| `data/external/gdc_clinical.tsv` | GDC clinical metadata per sample (primary diagnosis, tissue type, demographics) |
| `results/tables/annotation_summary.csv` | Cell-type counts and annotation confidence per cluster |
| `results/figures/cnv_heatmap.png` | CNV sliding-window heatmap — malignant/normal boundary |
| `data/processed/nerve_cells.h5ad` | Non-malignant nerve-cell AnnData with Leiden clusters and clinical obs |
| `results/figures/nerve_cells_umap.png` | UMAP of nerve-cell subspace colored by cluster and cell type |
| `results/tables/nerve_cluster_markers.csv` | Wilcoxon DE markers per nerve-cell cluster (padj < 0.05) |
| `results/tables/nerve_enrichment.csv` | GO Biological Process / Molecular Function enrichment per cluster |
| `results/figures/nerve_dotplot.png` | Canonical marker expression dot plot across clusters |
| `results/figures/nerve_abundance_heatmap.png` | Per-sample cluster proportion heatmap |
| `results/snakemake_report.html` | Full reproducibility report with DAG, rule stats, and provenance |
| `provenance/*.json` | FAIR provenance record for every rule output (artifact ID, SHA256, parameters) |

---

## Configuration

All analysis parameters live in `config/config.yaml`. Key sections:

| Section | What it controls |
|---|---|
| `scrna` | QC thresholds, HVG count, scVI latent dims, batch key, random seed |
| `hardware` | MPS device, float32 precision, memory watermark ratio |
| `nerve_cells` | Leiden resolution, cell-type scope, canonical marker gene lists |
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

1. **Subsampled data** — GDC distributes `seurat.1000x1000.loom` files (~1000 cells × 1000 genes per sample). Rare nerve cell populations may be under-represented. Check `*_gene_presence.csv` files after the smoke test.
2. **Malignancy classification** — uses sliding-window CNV scoring with T cells and endothelial cells as a reference. Confidence is proportional to the number of reference cells available per sample.
3. **GDC clinical metadata** — IDH and MGMT status are not consistently curated in TCGA-GBM. Expect `NA` values; all plots include an "unknown" category.
4. **GSEA** — requires internet access to query the Enrichr API. If offline, enrichment CSV will be empty but the pipeline continues without error.