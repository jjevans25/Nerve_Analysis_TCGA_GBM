# Execution Instructions — GBM Nerve Cell Heterogeneity Pipeline

## Prerequisites

All commands must be run from the project root:
```
/Users/jarrettevans/Documents/Biomedical Data Science/Projects/Nerve_Analysis_TCGA_GBM/
```

The virtual environment at `claude_science/` provides Python 3.12 and Snakemake. Use it in place of the system Python for all commands below.

Activate it by running: source claude_science/bin/activate 

---

## Step 1 — Validate the DAG (dry run)

Always run this first to confirm the pipeline is wired correctly before executing.

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores 8 -n
```

Expected output: **44 jobs** across 12 rule types, no errors.

---

## Step 2 — Single-sample smoke test

Runs only the ingest and QC rules on one sample to verify the loom→h5ad conversion and conda env build before committing to the full pipeline. The first run will also build the `scrna` conda environment (~5–10 min for new deps).

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores 4 \
  data/processed/06820e2c-9eb7-4e71-a1c3-976d561e659d_qc.h5ad
```

**After this step — check gene presence before proceeding:**

```bash
# Inspect which canonical nerve-cell markers survived the 1000-gene subsample
cat results/tables/06820e2c-9eb7-4e71-a1c3-976d561e659d_gene_presence.csv
```

If more than 3 canonical markers per cell type are absent across most samples, stop and consider re-downloading full-resolution loom files from GDC before running integration.

---

## Step 3 — Full pipeline

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all
```

This runs all 44 jobs in dependency order. On an M4 Max, expect:

| Stage | Approx. time |
|---|---|
| 17× loom_to_h5ad | 2–5 min |
| 17× scrna_qc | 5–10 min |
| scrna_integration (scVI VAE) | 20–60 min (MPS-accelerated) |
| gdc_clinical_fetch | 1–2 min |
| scrna_annotate | 5–15 min |
| scrna_malignancy (CNV scoring) | 10–20 min |
| nerve_cell_subset | 5–10 min |
| nerve_cell_heterogeneity + GSEA | 10–30 min |

---

## Step 4 — Goal-backward verification

After the full run, confirm all terminal artifacts are present and substantive:

```bash
# Check all rule-all targets exist and are non-empty
ls -lh results/tables/qc_summary.csv \
        results/tables/annotation_summary.csv \
        results/tables/nerve_cluster_markers.csv \
        results/tables/nerve_enrichment.csv \
        data/external/gdc_clinical.tsv

ls -lh results/figures/cnv_heatmap.png \
        results/figures/nerve_cells_umap.png \
        results/figures/nerve_dotplot.png \
        results/figures/nerve_abundance_heatmap.png

# Confirm nerve markers are biologically meaningful
head -20 results/tables/nerve_cluster_markers.csv

# Check malignancy classification stats
grep "pct_malignant" provenance/malignancy_*.json
```

---

## Step 5 — FAIR validation and Snakemake report

```bash
# Confirm all provenance JSONs are present
ls provenance/*.json | wc -l   # expect ≥ 20

# Generate reproducibility report
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  --report results/snakemake_report.html
```

---

## Step 6 — Interactive exploration (Marimo)

Launch the notebook to explore nerve-cell results interactively. The notebook renders the UMAP, dot plot, DE table, GSEA terms, and clinical metadata via DuckDB.

```bash
claude_science/bin/python3 -m marimo edit notebooks/01_explore_gbm_data.py
```

---

## Targeted re-runs

Re-run only nerve-cell analysis (e.g., after adjusting `leiden_resolution` in config.yaml):

```bash
# Force re-run of subset and heterogeneity rules only
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  --forcerun nerve_cell_subset
```

Re-run GDC clinical fetch only:

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores 2 \
  data/external/gdc_clinical.tsv
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `FileNotFoundError: No loom file found` | UUID mismatch or wrong `loom_dir` in config | Check `config/config.yaml` `loom_dir` and `samples` list against `data/raw/gdc_extract/MANIFEST.txt` |
| MPS silent CPU fallback | PyTorch MPS not available or op not supported | Check logs for `[MPS-ALERT]`; set `hardware.device: cpu` in config as fallback |
| GSEA timeout / no enrichment | gseapy Enrichr API unreachable | Check internet connectivity; enrichment CSV will be empty but pipeline continues |
| `scrna_annotate` produces only `tumor_gbm` cells | All 1000 genes from this sample are tumor-enriched | Normal — check `annotation_summary.csv` confidence scores; nerve-cell subset may be small |
| `nerve_cell_subset` produces < 20 cells | Nerve cells rare or absent in subsampled data | Review gene-presence CSVs; consider full-resolution data re-download |

---

## Key output files

| File | Description |
|---|---|
| `data/external/gdc_clinical.tsv` | GDC clinical metadata (IDH/MGMT/tissue type) per sample |
| `results/tables/annotation_summary.csv` | Cell-type counts and mean confidence scores |
| `results/figures/cnv_heatmap.png` | CNV sliding-window heatmap (malignant vs normal boundary) |
| `data/processed/nerve_cells.h5ad` | Non-malignant nerve-cell AnnData with Leiden clusters + clinical obs |
| `results/figures/nerve_cells_umap.png` | UMAP of nerve-cell subspace |
| `results/tables/nerve_cluster_markers.csv` | DE markers per nerve-cell Leiden cluster (Wilcoxon, padj < 0.05) |
| `results/tables/nerve_enrichment.csv` | GO enrichment per cluster (gseapy Enrichr) |
| `results/figures/nerve_dotplot.png` | Canonical marker expression dot plot |
| `results/figures/nerve_abundance_heatmap.png` | Per-sample cluster proportions heatmap |
| `provenance/*.json` | FAIR provenance records for every rule output |
| `results/snakemake_report.html` | Full reproducibility report with DAG and rule stats |
