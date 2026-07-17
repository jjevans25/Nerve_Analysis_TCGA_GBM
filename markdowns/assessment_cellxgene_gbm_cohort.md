# Assessment: Suitability of `cellxgene_data/gbm_10x_raw.h5ad` as a Second GBM Replication Cohort

> **Status:** Assessment only — no pipeline execution, no config/loader changes.
> All findings came from **metadata-only** reads (h5py structure + categorical value-counts +
> a 200k-nnz sample of `.X`). Nothing compute-intensive was run.
> Companion to `plan_second_gbm_replication_cohort.md`.

## Context

The v1.3.0 baseline produced a tumor→immune→nerve ligand–receptor result on a single 17-sample
TCGA-GBM cohort. We want an **independent GBM scRNA-seq cohort** to test whether those interactions
reproduce. A prior candidate (SCP393 / Neftel, in `broad_data/`) had a **thin nerve arm (~500 oligo)**
and is **Smart-seq2 logTPM** — a poor match for the scVI/scANVI (raw-count) pipeline. This new file is
a CELLxGENE Census pull (`scripts/gbm_cellxgene_census_pull.py`: `disease == 'glioblastoma' and
is_primary_data == True`, 10x-only).

---

## Verdict: **HIGH suitability — the strongest replication candidate available.**

Materially better than SCP393 on every axis that matters to the pipeline.

### What the file is
- **1,288,507 cells × 61,497 genes**, CSR sparse, 3.28B non-zeros, 7.6 GB (gzip `.h5ad`).
- 174 donors, 4 Census `dataset_id`s (studies), 4 10x chemistries. `is_primary_data == True` throughout
  (de-duplicated — no meta-analysis double-counting).

### Why it fits (evidence)

| Requirement (from CLAUDE.md / v1.3.0 pipeline) | This file | Verdict |
|---|---|---|
| **Raw counts** for scVI (`.X`) | `.X` sample = integers, min 1, values 1/2/4/… | ✅ genuine raw UMIs |
| **Ensembl + symbol** gene IDs (`.var`) | `feature_id` = ENSG…, `feature_name` = symbols | ✅ canonical, no remap needed |
| **Tumor arm** | malignant 340,756 + neoplastic 58,327 (~399k) | ✅ well-powered |
| **Immune arm** | macrophage 336k, microglia 240k, T 99k, monocyte 42k, DC/NK/B/neutrophil | ✅ very well-powered |
| **Nerve arm** (the SCP393 weakness) | oligodendrocyte 67,664 + OPC 23,881 + neuron 3,413 (~95k) | ✅ **~190× the SCP393 nerve arm** |
| **Marker panels transfer** (nerve + immune) | 9/10 panels 100% present; ependymal 5/6 (only PIFO missing — already absent in the reference cohort, not a regression) | ✅ subtype calling will work |
| **Per-donor batch key** (`batch_key: sample_id`) | 174 donors; pull already sets `sample_id = donor_id` | ✅ transfers cleanly |
| **Ingestion path exists in the plan** | `source: h5ad` branch of the proposed `ingest_dataset.py` | ✅ no new loader logic |

The nerve arm is the decisive improvement: the SCP393 memo flagged only ~500 oligodendrocytes and no
usable neuron population; here the tumor→immune→**nerve** third leg is backed by ~95k nerve-lineage cells
including 3,413 neurons — enough to actually re-derive nerve clusters rather than annotate a handful.

---

## Caveats / things to resolve before a run (not blockers)

1. **Compute scale is the real constraint.** 1.29M cells is orders of magnitude above the 17-sample
   cohort (scVI ~5h there). Chosen scope = **single largest study** — but study `56c4912d` is itself
   **1,020,902 cells across 170 donors and 3 chemistries**, so it does **not** remove batch structure.
   At ~1M cells scVI/scANVI on the M4 Max is still a day-scale run. **Recommendation:** even within the
   chosen study, add a per-donor subsample cap (e.g. ≤3–5k cells/donor → ~150–300k cells) for the first
   replication pass; the nerve arm stays well-powered (~15–20k oligo/OPC/neuron) and scVI fits an
   overnight run. This is a config/ingest knob, not a science change.

   Composition of the chosen study `56c4912d` (1,020,902 cells, 170 donors, 10x 3' v2/v3 + 5' v1):
   malignant 340,756 · macrophage 253,683 · microglia 191,490 · mature T 99,840 · monocyte 39,963 ·
   oligodendrocyte 34,723 · OPC 21,549 · DC 8,901 · neuron 3,413 · (+ endothelial/NK/B/plasma/mast/…).
   All three interaction arms are intact in this single study.

2. **`batch_key` must still be `donor_id`, not `dataset_id`.** The chosen study spans 170 donors and 3
   chemistries — keep the per-donor batch covariate (the pull already maps `sample_id = donor_id`).
   Consider chemistry as a second covariate if integration looks driven by 3'/5' differences.

3. **Census `cell_type` labels are provided but the pipeline re-annotates.** The pipeline does its own
   marker-scored annotation + scANVI. The Census ontology labels are a strong **cross-check / potential
   scANVI seed**, not a substitute — expect to still run `scrna_annotate` + `scrna_malignancy`. The
   provided `malignant cell` label is useful for validating the malignancy call.

4. **Gene space is the full 61,497-feature Census panel** (61,378 measured). HVG selection in
   `scrna_integration` handles this; no pre-filter required, but ingest should carry the `nnz`/
   `n_measured_obs` `.var` columns for provenance.

5. **Freezes off.** A replication cohort gets a fresh cut — `freeze_nerve_subset: false`,
   `freeze_cl15_split: false` (those SHA-pinned freezes belong to the TCGA reference cohort only).

---

## How this maps onto the existing plan

This file exercises the **`source: h5ad`** branch of `plan_second_gbm_replication_cohort.md` — the
cohort-namespaced parallel track (own scVI/scANVI, own `results/tables/<dataset>/…`, then
`cohort_concordance` vs the reference cohort's `nerve_tumor_immune_interactions_with_qc.csv`). No new
design is needed beyond that plan; the only addition this data argues for is a **per-donor subsample
knob** in the ingest step to keep the ~1M-cell study tractable.

## Recommended next step (deferred — nothing heavy now)

First light action: a **Stage-A ingest smoke test** on a single donor from study `56c4912d` (subset the
`.h5ad`, confirm `.X` raw counts + symbol mapping) — no scVI. The full overnight scVI/scANVI run is the
subsequent, explicitly-gated step.
