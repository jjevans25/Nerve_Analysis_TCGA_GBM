# TCGA-GBM Nerve-Cell Pipeline — Project Overview

*Landing page for the project. Current cut: `v1.1.0` (2026-05-24). Original frozen baseline: `v1.0.0` (2026-05-09).*

---

## At a glance

This repository implements an autonomous, FAIR-compliant single-cell RNA-seq pipeline for characterising the heterogeneity of **non-malignant nerve cells inside the glioblastoma (GBM) tumor microenvironment**. The stack is Snakemake (workflow), scvi-tools (batch correction), scanpy (clustering / DE), gseapy + MSigDB (offline GSEA), LIANA (cell–cell signaling), and marimo (interactive notebooks), targeting Apple M4 Max via the MPS PyTorch backend. The headline deliverables are three interactive marimo notebooks — exported to self-contained HTML — plus the tabular and image artifacts they read.

## Scientific question

> *"What subtypes of nerve cells (neurons, oligodendrocytes, OPCs, astrocytes) are present in the GBM TME, how heterogeneous are they across patients, and which molecular programs define each subtype?"*
> — `README.md`

The pipeline approaches this in three passes: (1) characterise the **distribution** of neural / glial cells across patients; (2) identify the **molecular programs** active in each cluster via differential expression and GO BP / GO MF enrichment; (3) surface candidate **paracrine signaling axes** linking malignant GBM cells to the surrounding nerve-cell compartment via ligand–receptor inference.

**Cohort:** 17 TCGA-GBM scRNA-seq samples (GDC loom files) → 184,494 cells post-QC → **106,603 non-malignant nerve cells across 24 Leiden clusters** 

## Pipeline architecture

Wired in `Snakefile` and `workflow/rules/*.smk`. Forty-four jobs in dependency order:

```
loom_to_h5ad (×17)        ─┐
                           ├─→ scrna_qc (×17)
gdc_clinical_fetch ────────┤
                           ├─→ scrna_integration   (scVI VAE on M4 Max MPS)
                           ├─→ scrna_annotate      (Leiden + marker scoring)
                           ├─→ scrna_malignancy    (CNV sliding-window)
                           ├─→ nerve_cell_subset   (re-cluster on X_scVI, 24 clusters)
                           ├─→ nerve_cell_heterogeneity   (DE + offline GSEA)
                           ├─→ nerve_cluster_annotations
                           ├─→ nerve_clinical_association
                           ├─→ nerve_tumor_interaction    (LIANA consensus L–R)
                           ├─→ nerve_batch_qc      (sample-purity verdict)
                           └─→ explore_gbm_notebook + nerve_enrichment_notebook
```

Every rule writes a JSON provenance record to `provenance/<rule>_provenance.json` (UUID5 + SHA-256 + tool versions + parameters), and the `freeze_baseline_provenance` rule bundles all of those into `provenance/baseline_<version>.json` for citable snapshots. GSEA runs **offline** against pinned MSigDB C5 GO BP+MF `.gmt` files (release `2024.1.Hs`, SHA-256-verified) to keep results deterministic and reproducible.

## Current state — `v1.1.0` cut (v1.0.0 frozen baseline retained)

**v1.1.0 deltas (2026-05-24, post failing-cluster diagnosis):** cl21 dropped from downstream `_with_qc.csv` tables as a dissociation artifact (stress-marker dominant). Ependymal added as the 5th nerve cell type (FOXJ1 / RFX3 / DNAH / CFAP / PIFO / RSPH1 markers + added to `nerve_cells.cell_types` so cells survive the nerve_cell_subset filter). After re-clustering with the +5,713 ependymal cells, **three** ependymal-dominant clusters emerged: **cl11** (72% ependymal, 2,585 cells; ependymal regulators PARAIL/GLIS3/DTNA/YAP1), **cl15** (89%, 1,854 cells), and **cl19** (83%, 1,020 cells; clean motile-cilia signature CFAP54/DNAH7/DNAH9). scANVI v2 retrained on 5 labels (was 4) with classifier accuracy 0.85. The `02_nerve_enrichment_explorer.py` theme regex now includes `ependymal` (CILIUM/CILIARY/AXONEMAL/DYNEIN/EPENDYM patterns). v1.0.0 baseline (`provenance/baseline_v1.0.0.json`) remains frozen and citable.

### v1.0.0 baseline (unchanged)

- **24 nerve-cell Leiden clusters** computed on the scVI latent space (`X_scVI`).
- **Batch-correction QC: 18 / 24 clusters PASS.** Median dominant-sample fraction `0.313`; median normalised entropy `0.722` of the uniform-distribution maximum `4.087` bits.
- **Per-cluster GSEA** (480 rows; 312 with FDR < 0.05) — sample top-line hits from the new clustering: cluster 7 → `GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS` (FDR 0.025); cluster 0 → DNA replication / cell-cycle (proliferating state); cluster 2 → axon-extension / axon-guidance neurons.
- **LIANA tumor↔nerve interactions:** 1,207 significant L–R pairs across the 24 clusters (`results/tables/nerve_tumor_interactions.csv`).
- **Provenance frozen:** `provenance/baseline_v1.0.0.json` (49 rule records + git commit + SHA-256s).
- **Known caveat — 6 batch-QC-"failing" clusters re-diagnosed (2026-05-23):** clusters **13, 15, 19, 21, 22, 23** fail the dominance test, but post-hoc inspection of markers + per-cell QC + clinical metadata showed that **5 of 6 are real biology, not batch-correction failures.** See `markdowns/failing_cluster_diagnosis.md` for the full evidence; headline:
   - **cl13, cl15, cl22** (synaptic-adhesion markers CNTNAP5 / NBEA / FGF14 / RIMS2) and **cl23** (perivascular GPC5 / LAMA2 / SLC4A4) are **patient-anatomy-specific neural subtypes** preserved primarily in donors whose `tissue_type` is `Normal` (tumor-adjacent intact brain). Their patient-dominance is *because* one donor's intact tissue captured the subtype, not because the subtype is patient-private. **Retained** in `_with_qc.csv`, flagged by `batch_qc_pass = False`.
   - **cl19** is **ependymal cells** (motile-cilia markers DNAH7 / 9 / 11, CFAP54, RFX3) mislabeled as astrocyte by upstream marker scoring that lacked `score_ependymal`. The next pipeline iteration adds `nerve_cells.markers.ependymal` in `config/config.yaml` and both per-cell scorers (`scrna_annotate.py`, `nerve_celltype_labels.py`) pick it up; cl19 will be re-labelled on the next full run.
   - **cl21 is a true technical artifact** — dissociation-stress signature (5 / 10 top markers are MT-RNR2 / MT-RNR1 / MT-CO1 / RPL41 / RPS12; `n_genes_median = 1099`, lowest of all 24 clusters). **Dropped** from `_with_qc.csv` variants via the new `nerve_cells.batch_qc.exclude_clusters = ["21"]` config key; the source `nerve_cluster_markers.csv`, `nerve_enrichment.csv`, and LIANA tables retain cl21 for audit. The marimo explorers no longer surface cl21 because they read the filtered `_with_qc.csv` tables.
   - **Why coarsening Leiden / adding scANVI didn't fix the "failures":** the sweep (`results/tables/nerve_leiden_resolution_sweep.csv`) showed no resolution clears the ≥ 90 % pass-rate gate; the scANVI v2 retrain (`results/tables/nerve_cluster_sample_purity_v2.csv`) made it worse (16 / 25 = 64 % pass) because a 4-label classifier head cannot distinguish neuron *subtypes*. The dominance test cannot distinguish "real biology preserved by one donor's intact tissue" from "patient bias in many-sample shared biology" — that's why a marker + QC + clinical cross-reference is needed before drop / retain decisions on dominance-failing clusters.

Detailed history is in `CHANGELOG.md` under the `[2026-05-09]` entries.

---

## The three notebooks

Each is a marimo reactive notebook (`.py` source) with a published HTML export. Edit interactively with `marimo edit notebooks/<name>.py`; or just open the HTML in any browser.

### 6.1 `notebooks/01_explore_gbm_data.py` — TCGA-GBM Sample Explorer

Sample-level QC plus a gallery of downstream nerve-cell artifacts. Use this to triage the cohort and confirm the pipeline produced everything.

**What you'll see (in scroll order):**
- 17-row sample manifest, DuckDB-queried from `MANIFEST.txt` (file, size, state).
- Per-sample QC histograms with config thresholds overlaid (`min_genes`, `max_genes`, `max_pct_mito`).
- Sample-selector dropdown, per-sample filter pass/fail summary table.
- The full nerve-cell artifact gallery: UMAP image, top-200 markers table (FDR < 0.05), top-100 enrichment table (FDR < 0.05), canonical-marker dot plot, per-sample abundance heatmap.
- Clinical-metadata table joined from `data/external/gdc_clinical.tsv`.
- FAIR provenance callout (run UUID, timestamp, tool versions).

**Reads:** `data/raw/gdc_extract/*.loom`, `data/processed/nerve_cells.h5ad`, `results/tables/nerve_cluster_markers.csv`, `results/tables/nerve_enrichment.csv`, `data/external/gdc_clinical.tsv`, plus pre-rendered figures in `results/figures/`.

**Exported HTML:** `results/figures/01_explore_gbm_data.html` (via the `explore_gbm_notebook` rule).

### 6.2 `notebooks/02_nerve_enrichment_explorer.py` — Nerve-Cell GSEA Enrichment Explorer

Four interactive views over `results/tables/nerve_enrichment.csv`. Use this to compare clusters by their pathway signatures and to drill into any one cluster.

**What you'll see:**
1. **Term-by-cluster heatmap** of `-log10(FDR)`, with library picker (BP / MF), top-N slider, FDR threshold slider, and an optional row-clustering switch (correlation linkage).
2. **Top-K terms per cluster** (interactive) — bar charts for GO BP and GO MF side-by-side, cluster-picker dropdown, companion table showing the leading-edge / set-size overlap.
3. **Cluster similarity from enrichment profiles** — pairwise cluster–cluster heatmap (cosine on `-log10(FDR)` weights or Jaccard on the binarised matrix) plus a hierarchical-linkage dendrogram.
4. **Term-category roll-up** across nine curated themes (axon/dendrite, synapse, myelin/glia, neuron development, immune, metabolic, cell-cycle, signaling, transcription, vasculature). Per-cluster theme score = mean `-log10(FDR)` over significant matching terms; companion table lists each cluster's dominant theme.

**Reads:** `results/tables/nerve_enrichment.csv` (and `nerve_cluster_markers.csv` as a companion reference).

**Exported HTML:** `results/figures/02_nerve_enrichment_explorer.html` (via the `nerve_enrichment_notebook` rule).

### 6.3 `notebooks/nerve_tumor_exploration.py` — Tumor–Nerve Ligand–Receptor Explorer

Interactive exploration of LIANA-computed L–R pairs between malignant cells and nerve-cell clusters. Use this to surface candidate paracrine axes and quantify directional asymmetry.

**What you'll see:**
- **Glossary accordion** explaining the LIANA columns (`cellphone_pvals`, `lrscore`, `magnitude_rank`, `specificity_rank`, `direction`).
- **Overview-statistics table** (counts, % significant, per-cluster / per-ligand / per-receptor breakdowns).
- **Filter panel** — toggle full vs. top_pairs, direction radio (both / malignant→nerve / nerve→malignant), nerve-cluster multiselect, magnitude_rank + p-value sliders, ligand / receptor text search.
- **Filtered results table** (paginated, 25 rows / page), sorted by magnitude_rank.
- **Cluster × direction significance heatmap** (cell value = number of L–R pairs passing both thresholds).
- **Top-K L–R dotplot** — clusters on the x-axis, "ligand → receptor" labels on the y-axis; dot size = `-log10(magnitude_rank)`, colour = `lrscore`.
- **Per-cluster bidirectional scatter** — `lrscore(malignant→nerve)` vs. `lrscore(nerve→malignant)`; on-diagonal points are bidirectional, off-axis points are unidirectional. Companion bar chart highlights the top-K most asymmetric pairs.
- **Export button** — saves the filtered view as CSV with a `.provenance.txt` sidecar (timestamp, source SHA-256[:12], filter values).

**Reads:** `results/tables/nerve_tumor_interactions.csv`, `results/tables/nerve_tumor_top_pairs.csv`.

**Exported HTML:** *not currently a Snakemake rule* — interactive only via `marimo edit notebooks/nerve_tumor_exploration.py`. (Wiring an export rule mirroring the other two is a clean follow-up if a static deliverable is needed.)

---

## Where to go next

| What you want | Where to look |
|---|---|
| Run / re-run the pipeline | `execution_instructions.md` (project root) |
| Pending work + remediation status | `markdowns/next_steps_interpretation.md` (§6 checklist) |
| Stratification follow-ups | `markdowns/stratifying_open_issues.md` |
| Lab-notebook session-by-session audit | `CHANGELOG.md` |
| Verify provenance of v1.1.0 outputs | `provenance/baseline_v1.1.0.json` |
| Verify provenance of v1.0.0 outputs | `provenance/baseline_v1.0.0.json` |
| Project rules / conventions | `CLAUDE.md` (project constitution) |

## Reproduce

```bash
# v1.1.0 (current cut)
git checkout v1.1.0
snakemake --use-conda --cores all

# v1.0.0 (original frozen baseline)
git checkout v1.0.0
snakemake --use-conda --cores all
```
