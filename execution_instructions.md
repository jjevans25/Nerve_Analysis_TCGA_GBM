# Execution Instructions — GBM Nerve–Tumor–Immune Pipeline

> ## ⚠️ Read first: the v1.3.0 reference cohort is pinned
>
> `config/config.yaml` sets `baseline.pinned: true`. The 10 rules that would rewrite reference
> nerve/immune artifacts are **not defined at parse time**, so no invocation can schedule them —
> including a bare `snakemake`. This is deliberate: `data/processed/nerve_cells.h5ad` was
> destroyed on 2026-07-21 and **cannot be reproduced** (a retrained latent re-clusters, so the
> frozen `cl15` split no longer maps to one source cluster).
>
> Consequences for everything below:
> - Any rule that *reads* the reference nerve tables must be run with
>   **`--allowed-rules <rule>`** plus `--rerun-triggers mtime`. Without it, Snakemake plans a
>   9-job nerve-cascade rebuild — verified by dry run.
> - `--forcerun nerve_cell_subset` and friends will fail with "no such rule". That is the pin
>   working, not a bug.
> - Set `baseline.pinned: false` only if you accept that new clusters will not match published
>   v1.3.0 labels.

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

Job count varies with which cohorts are enabled (`samples:` and `datasets:` in `config.yaml`)
and which pinned artifacts still exist on disk, so there is no fixed expected number. What to
check instead:

- the DAG resolves with **no errors**
- **no reference-cohort rule you did not intend appears in the plan** — if you see a nerve
  cascade you did not ask for, add `--allowed-rules` before running for real
- replication work appears as `ds_*` rules only

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

This runs every enabled target in dependency order. On an M4 Max, expect roughly:

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
| immune_cell_subset + immune_cluster_annotations | 5–15 min |
| nerve_tumor_interaction (2-way LIANA) | 5–15 min |
| nerve_tumor_immune_interaction (3-way LIANA) | 5–15 min |
| nerve_scanvi_retrain (scANVI-v2) | 30–90 min (MPS) |

With `baseline.pinned: true`, the nerve/immune stages above are skipped for the reference
cohort — their outputs are frozen tables.

---

## Step 3b — Replication cohort

The replication chain is cohort-namespaced (`ds_*`). Its terminal concordance target pulls the
whole Stage A→D chain for one cohort and leaves reference outputs untouched:

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  results/tables/gbm_cellxgene_56c4912d/cohort_concordance_summary.json
```

Scale note: the Census source file is 1.29M cells × 61,497 genes. The configured run takes the
single largest study (`56c4912d`) capped at `subsample_per_donor: 5000`. Removing that cap
changes the compute profile by orders of magnitude — do not do it casually on a laptop.

The scANVI-v2 side branch is not pulled by concordance; request it explicitly:

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  results/tables/gbm_cellxgene_56c4912d/nerve_cluster_sample_purity_v2.csv
```

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
        results/figures/nerve_dotplot.png \
        results/figures/nerve_abundance_heatmap.png

# Immune compartment + crosstalk
ls -lh results/figures/immune_cells_umap.png \
        results/tables/immune_cluster_annotations.csv \
        results/tables/nerve_tumor_interactions.csv \
        results/tables/nerve_tumor_immune_interactions.csv

# Replication cohort (per enabled dataset)
ls -lh results/tables/gbm_cellxgene_56c4912d/cohort_concordance_summary.json

# Confirm nerve markers are biologically meaningful
head -20 results/tables/nerve_cluster_markers.csv

# Check malignancy classification stats
grep "pct_malignant" provenance/malignancy_provenance.json

# Verify the pinned reference still matches its recorded SHA256s.
# This is the ONLY drift detector for those files: while baseline.pinned is true
# their producing rules do not exist, so Snakemake cannot notice a modified or
# deleted pinned table on its own.
claude_science/bin/python3 -m snakemake --use-conda --cores 1 verify_pinned_reference
# -> results/pinned_reference_verification.json
```

> `results/figures/nerve_cells_umap.png` and `data/processed/nerve_cells.h5ad` are **absent by
> design** for the reference cohort — lost 2026-07-21, not reproducible. Do not treat their
> absence as a failed run. The replication cohort has both under its namespace.

**LIANA sanity check.** Both interaction scripts hard-fail if `.X` is off-scale, but confirm the
tables are populated rather than empty-ranked:

```bash
# specificity_rank must not be empty across all rows
head -3 results/tables/nerve_tumor_immune_interactions.csv
```

---

## Step 5 — FAIR validation and Snakemake report

```bash
# Confirm all provenance JSONs are present
ls provenance/*.json | wc -l   # 65 as of v1.3.0 + replication cohort

# Generate reproducibility report
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  --report results/snakemake_report.html
```

---

## Step 6 — Interactive exploration (Marimo)

Launch a notebook to explore results interactively. Each renders figures, DE/LR tables, and
clinical metadata via DuckDB.

```bash
claude_science/bin/python3 -m marimo edit notebooks/01_explore_gbm_data.py
```

| Notebook | Scope |
|---|---|
| `01_explore_gbm_data.py` | QC explorer + nerve-cell viewer |
| `02_nerve_enrichment_explorer.py` | Nerve cluster DE / GSEA results |
| `03_nerve_tumor_immune_explorer.py` | Three-way interaction results |
| `04_tme_nerve_immune_explorer.py` | Full TME view across compartments |
| `05_census_nerve_immune_explorer.py` | Replication-cohort equivalent (per cohort) |

Notebooks 03 and 04 read the **reference** tables; 05 reads a replication cohort. All are also
exported to HTML by the pipeline into `results/figures/` — editing is only needed for
interactive work.

---

## Targeted re-runs

> ⚠️ **`--forcerun nerve_cell_subset` no longer works and should not be made to work.** While
> `baseline.pinned: true`, that rule is not defined; Snakemake will report no such rule. Adjusting
> `nerve_cells.leiden_resolution` cannot re-cut the reference cohort. Tune resolution on a
> **replication cohort** instead, where the equivalent `ds_*` rules are live.

Re-run a replication cohort's nerve analysis (e.g. after adjusting `leiden_resolution`):

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  --forcerun ds_nerve_cell_subset \
  results/tables/gbm_cellxgene_56c4912d/cohort_concordance_summary.json
```

Re-run only the LIANA interaction rules for a cohort (the pattern used for the 2026-07-26
normalisation fix — 5 jobs, all `ds_*`, reference untouched, ~7 min):

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  --rerun-triggers mtime \
  --forcerun ds_nerve_tumor_interaction ds_nerve_tumor_immune_interaction
```

Re-run a rule that **reads** pinned reference tables — the `--allowed-rules` insulator is
mandatory, not optional. Without it Snakemake plans a 9-job nerve-cascade rebuild:

```bash
claude_science/bin/python3 -m snakemake --use-conda --cores all \
  --rerun-triggers mtime \
  --allowed-rules tme_nerve_immune_notebook \
  results/figures/04_tme_nerve_immune_explorer.html
```

Always dry-run (`-n`) these first and read the job list before dropping `-n`.

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
| `There is no rule named nerve_cell_subset` | `baseline.pinned: true` — rule not defined | Expected. Work on a replication cohort's `ds_*` rule instead |
| Snakemake plans a 9-job nerve cascade you didn't ask for | Missing `--allowed-rules` on a rule that reads pinned tables | Add `--allowed-rules <rule> --rerun-triggers mtime`; verify with `-n` first |
| LIANA rows all have empty `specificity_rank` | `.X` on the wrong scale at the LIANA call | Check the `normalize_counts` param for that cohort; the scripts now hard-fail rather than proceed silently |
| Interaction script raises on matrix scale | `.X` not log1p and `normalize_counts` false | Reference `.X` is SCT log1p (`false`); Census is raw UMIs (`true`). See `counts_utils.is_log1p_scale()` |
| Pinned artifact hash mismatch | A frozen v1.3.0 table was modified or deleted | `verify_pinned_reference` reports which; restore from git — do **not** re-freeze to absorb drift |
| Census ingest exhausts memory | Full 1.29M-cell file loaded | Confirm `filter.dataset_id` and `subsample_per_donor` are set in the `datasets:` entry |

---

## Key output files

| File | Description |
|---|---|
| `data/external/gdc_clinical.tsv` | GDC clinical metadata (IDH/MGMT/tissue type) per sample |
| `results/tables/annotation_summary.csv` | Cell-type counts and mean confidence scores |
| `results/figures/cnv_heatmap.png` | CNV sliding-window heatmap (malignant vs normal boundary) |
| ~~`data/processed/nerve_cells.h5ad`~~ | **Lost 2026-07-21, not reproducible** for the reference; exists per cohort under `data/processed/{cohort}/` |
| ~~`results/figures/nerve_cells_umap.png`~~ | **Lost with the above**; exists under `results/figures/{cohort}/` |
| `results/tables/nerve_cluster_markers.csv` | DE markers per nerve-cell Leiden cluster (Wilcoxon, padj < 0.05) |
| `results/tables/nerve_enrichment.csv` | GO enrichment per cluster (gseapy Enrichr) |
| `results/figures/nerve_dotplot.png` | Canonical marker expression dot plot |
| `results/figures/nerve_abundance_heatmap.png` | Per-sample cluster proportions heatmap |
| `results/figures/immune_cells_umap.png` | UMAP of the immune subspace by cluster and subtype |
| `results/tables/immune_cluster_annotations.csv` | Microglia / TAM / T / NK / dendritic subtype calls |
| `results/tables/nerve_tumor_interactions.csv` | Two-way tumor→nerve ligand–receptor pairs (LIANA) |
| `results/tables/nerve_tumor_immune_interactions.csv` | Three-way tumor→immune→nerve ligand–receptor pairs |
| `results/tables/*_with_qc.csv` | The above joined to cluster-QC flags — **interpret from these, not the raw tables** |
| `results/tables/{cohort}/cohort_concordance_summary.json` | Reference vs replication agreement verdict |
| `results/tables/nerve_cluster_sample_purity_v2.csv` | scANVI-v2 donor purity (reference: 16/26 pass; replication: 24/33) |
| `provenance/pinned_reference_v1.3.0.json` | SHA256 manifest of every pinned v1.3.0 artifact |
| `provenance/*.json` | FAIR provenance records for every rule output |
| `results/snakemake_report.html` | Full reproducibility report with DAG and rule stats |


## To Pickup Where it Left Off
claude --continue  # or use shorthand: claude -c
