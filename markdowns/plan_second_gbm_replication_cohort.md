# Plan: Add a Second GBM scRNA-seq Cohort as a Replication Track

> **Status:** DRAFT for researcher review. Nothing has been executed or changed in the pipeline.
> This document only describes the proposed approach.

## Context

**Why:** The project has produced a tumor→immune→nerve ligand–receptor interaction result on a single
17-sample TCGA-GBM cohort (the v1.3.0 baseline). The researcher wants to run the *same* analysis on
**other GBM scRNA-seq data** and check whether the tumor→immune→nerve pathways **align/reproduce** —
i.e. use a second, independent GBM cohort as a **replication cohort** to add robustness to the published
interactions.

**Decisions locked in (from clarifying questions):**
- **Data source:** "Both / not sure yet" → the ingestion layer must be generalized to accept either more
  GDC looms *or* an external GBM cohort (GEO / CELLxGENE / raw `.h5ad`).
- **Integration:** "Separate parallel cohort" → the new data gets its **own** scVI/scANVI model and its own
  end-to-end run. It is **not** co-trained with the frozen 17-sample atlas. Results are then **compared**.

### Is the project conducive to this today?

| Scenario | Verdict |
|---|---|
| **More GDC TCGA looms** (same portal/format, versioned-Ensembl IDs) | **Yes — nearly config-only.** Drop looms under `data/raw/gdc_extract/<uuid>/`, add UUIDs to `samples:`, retire the nerve freezes, rerun. |
| **A genuinely external GBM cohort** (GEO/CELLxGENE, `.h5ad`/`.mtx`, HGNC-symbol genes, non-UUID IDs) | **Not yet.** Three hard blocks in the code prevent it (below). |

The three blocks for external data:
1. **Loader is loom-only.** `workflow/rules/ingest.smk:_find_loom` globs `*.loom`; `loom_to_h5ad.py` uses
   `loompy.connect`.
2. **Sample wildcard is UUID-locked.** `ingest.smk:51-52` constrains `sample` to a UUID regex — a GEO
   `GSM…`/sample name won't match.
3. **Single flat namespace, no dataset dimension.** `config.yaml` has one `samples:` list + scalar
   `loom_dir`; every output is a flat singleton (`data/processed/integrated_latent.h5ad`, `annotated.h5ad`,
   `results/tables/nerve_tumor_immune_interactions.csv`). A second cohort would collide with the first.

This plan removes all three **without touching the existing cohort's outputs** — critical because the
v1.0.0–v1.3.0 semver baselines and the nerve freezes (`nerve_subset_v1_1_0.txt`, `cl15_split_v1_3_0.csv`)
are pinned to exact paths and SHA-256s.

---

## Design principle: cohort-namespaced parallel track

Keep the existing pipeline **byte-for-byte unchanged** (frozen baselines stay valid). Add a new
**dataset-scoped** replication pipeline that reuses the existing analysis *scripts* but redirects all I/O
under a per-dataset namespace:

```
data/raw/<dataset>/            # new cohort's raw inputs (+ MANIFEST.txt)
data/processed/<dataset>/      # <sample>.h5ad, <sample>_qc.h5ad, integrated_latent.h5ad, annotated.h5ad, …
results/models/<dataset>/      # scvi_model, nerve_scanvi_model
results/tables/<dataset>/      # nerve_tumor_immune_interactions.csv, …
results/figures/<dataset>/
provenance/<dataset>/          # <dataset>_<sample>_ingest_provenance.json, …
```

The current TCGA-GBM cohort is treated as the **reference** and is *not re-run* — the concordance step
just reads its existing `results/tables/nerve_tumor_immune_interactions_with_qc.csv`.

---

## Config changes (`config/config.yaml`)

Add a `datasets:` block. Each entry is self-describing so the loader can dispatch:

```yaml
datasets:
  gbm_validation_01:                 # <dataset> wildcard value
    description: "External GBM scRNA-seq replication cohort"
    disease: "glioblastoma"
    source: "h5ad"                   # gdc_loom | h5ad | mtx | cellxgene_census
    raw_dir: "data/raw/gbm_validation_01"
    gene_id_type: "symbol"           # symbol | ensembl | ensembl_versioned
    samples:                         # arbitrary IDs (no UUID requirement)
      - GSMxxxxxx1
      - GSMxxxxxx2
    batch_key: "sample_id"           # per-cohort scVI batch covariate
    # Replication cohorts get a FRESH cut — disable the reference cohort's freezes:
    freeze_nerve_subset: false
    freeze_cl15_split: false
```

- Reuse the existing `scrna:`, `nerve_cells.markers`, `nerve_cells.cell_types`, `immune_cells.subtype_markers`,
  `nerve_scanvi:` blocks **as shared defaults** (same disease → same panels). Only the freeze/override files
  are per-dataset (off for replication).
- Legacy top-level `samples:`/`loom_dir:` stay for the reference cohort — untouched.

---

## Selected replication cohort: `gbm_cellxgene_56c4912d` (CHOSEN 2026-07-16)

The abstract "which second cohort?" is now resolved. We use the **CELLxGENE Census GBM pull** already on
disk at `cellxgene_data/gbm_10x_raw.h5ad` (1,288,507 cells × 61,497 genes; raw 10x UMIs; `.var` carries
Ensembl `feature_id` + symbol `feature_name`). Full suitability write-up:
`markdowns/assessment_cellxgene_gbm_cohort.md` (verdict: HIGH).

**Scoping decisions:**
- **Restrict to the single largest Census study** `dataset_id == 56c4912d-2bae-4b64-98f2-af8a84389208`
  (~1,020,902 cells; 170 donors; 10x 3′ v2/v3 + 5′ v1). Drops cross-study batch effects but **not**
  within-study structure — 170 donors / 3 chemistries remain, so **`batch_key` = `donor_id`** (chemistry
  optionally added as a second covariate).
- **Per-donor subsample cap = 5,000 cells/donor** (seed = `scrna.random_seed` = 0). Takes the study from
  ~1.02M cells to roughly ~150–300k while keeping the nerve arm well-powered (~15–20k oligo/OPC/neuron),
  so scVI/scANVI fit an overnight run instead of a day-scale one. This is a config/ingest knob, not a
  science change.
- **`sample` == `donor`.** The pull sets `sample_id = donor_id`, so each donor is a "sample": the
  Stage-A per-sample ingest emits `data/processed/gbm_cellxgene_56c4912d/<donor>.h5ad` and the 5k cap
  applies per that file. The `samples:` list is **auto-derived** from `obs["donor_id"]` (filtered to the
  target study) rather than hand-enumerated (170 IDs).

**Concrete config entry (replaces the `gbm_validation_01` placeholder above):**

```yaml
datasets:
  gbm_cellxgene_56c4912d:
    description: "CELLxGENE Census GBM 10x replication cohort — single largest study"
    disease: "glioblastoma"
    source: "h5ad"
    raw_file: "cellxgene_data/gbm_10x_raw.h5ad"       # single file, not a dir
    filter: { dataset_id: "56c4912d-2bae-4b64-98f2-af8a84389208" }
    gene_id_type: "ensembl"          # .var feature_id = Ensembl; feature_name = symbol
    sample_key: "donor_id"           # donor == sample; auto-derive `samples` from obs
    subsample_per_donor: 5000        # per-donor cap; random_state = scrna.random_seed
    batch_key: "sample_id"           # = donor_id
    freeze_nerve_subset: false       # replication cohort → fresh cut
    freeze_cl15_split: false
```

**`ingest_dataset.py` `h5ad` branch (behavior confirmed by smoke test):** read `raw_file` backed → subset
to `filter` + one `sample_key` value → `.to_memory()` → if `> subsample_per_donor`, `sc.pp.subsample`
(seeded) → map `feature_id`→`ensembl_id` (canonical `var_names`) and `feature_name`→`gene_symbol` → flag
MT via `MT-` symbol prefix → set `obs` `dataset`/`batch`/`sample_id` from `donor_id` → write
`data/processed/<dataset>/<donor>.h5ad` + provenance. **Stage-A smoke test (2026-07-16, donor `BT389`,
5,028→5,000 cells): raw counts OK (int, min 1, max 8725), 100% symbol mapping, 37 MT genes, all marker
panels present, round-trips in scanpy.**

---

## New / changed files

**New rules file — `workflow/rules/datasets.smk`** (modeled on the `immune.smk` additive-stage precedent),
`include:`-d in `Snakefile` after `notebooks.smk`. Defines the dataset-scoped chain with a `{dataset}`
wildcard driving every path. Helper `_find_dataset_input(wildcards)` resolves a sample's raw file by the
entry's `source`/`raw_dir` (replacing the loom-only `_find_loom`), and `wildcard_constraints` derives the
`sample` pattern per-dataset (relaxed from the UUID regex).

**New generalized loader — `workflow/scripts/ingest_dataset.py`.** Dispatches on `source`:
- `gdc_loom` → reuse existing `loom_to_h5ad.py` logic (extract into a shared function).
- `h5ad` → `ad.read_h5ad`, standardize `.X` to raw counts, set `obs["batch"|"sample_id"|"dataset"]`.
- `mtx` → `sc.read_10x_mtx`.
- `cellxgene_census` → query by SOMA `value_filter` (disease/assay) using the pattern already in
  `scrna_malignancy.py:26-55`.
Then **gene-space harmonization**: map every dataset to the project canonical (`ensembl` + `gene_symbol`
columns in `.var`) so downstream marker scoring works. For `gene_id_type: symbol`, add the reverse
(symbol→Ensembl) direction to `build_gene_symbol_map.py`. Stamps provenance via `fair_utils` exactly like
`loom_to_h5ad.py:94-109`.

**New concordance step — `workflow/scripts/cohort_concordance.py`** + rule `cohort_concordance`. The
robustness deliverable: read the new cohort's `nerve_tumor_immune_interactions_with_qc.csv` and the
reference cohort's, then compute:
- overlap of significant LR pairs (Jaccard, shared/unique pair table),
- rank correlation of interaction magnitudes on shared pairs (Spearman),
- shared vs. cohort-specific pathway summary,
- a comparison figure + a small `cohort_concordance_summary.json`.
Optionally wire a marimo explorer (`notebooks/04_cohort_concordance_explorer.py`) following the existing
notebook-export pattern in `notebooks.smk`.

**`Snakefile`:** add `include: "workflow/rules/datasets.smk"`; add `DATASETS = config.get("datasets", {})`;
append dataset-scoped targets to `rule all` gated on `if DATASETS` (mirroring the `if SAMPLES` pattern),
including the final `cohort_concordance_summary.json`.

**Reused as-is (I/O redirected via dataset-scoped rules, no script logic changes beyond param wiring):**
`scrna_qc.py`, `scrna_integration.py`, `scrna_annotate.py`, `scrna_malignancy.py`, `nerve_cell_subset.py`
(with freezes disabled), `nerve_scanvi_retrain.py`, `immune_cell_subset.py`,
`nerve_tumor_immune_interaction.py`, and all of `fair_utils.py`.

---

## Staged implementation

Deliver in stages so each is verifiable before the next (Snakemake-first; each rule's artifact verified
via `fair_utils.verify_artifact` before proceeding):

- **Stage A — Generalized ingestion:** `datasets:` config block, `ingest_dataset.py`, `datasets.smk`
  ingest rules, gene-space harmonization. Output: `data/processed/<dataset>/<sample>.h5ad` with
  `ensembl_id`+`gene_symbol` var and `dataset`/`batch`/`sample_id` obs.
- **Stage B — Cohort QC + integration + annotation:** dataset-scoped `scrna_qc` → `scrna_integration`
  (own scVI model, `batch_key` per entry) → `scrna_annotate` → `scrna_malignancy`. Output:
  `data/processed/<dataset>/malignancy_labeled.h5ad`.
- **Stage C — Nerve + immune + three-way interaction:** dataset-scoped `nerve_cell_subset` (freezes off)
  → `nerve_scanvi_retrain` → `immune_cell_subset` → `nerve_tumor_immune_interaction`. Output:
  `results/tables/<dataset>/nerve_tumor_immune_interactions_with_qc.csv`.
- **Stage D — Cross-cohort concordance:** `cohort_concordance` rule → robustness summary + figure
  (+ optional explorer). This is the deliverable that answers "do the pathways align?".

**Compute note:** each cohort re-runs scVI (~hours on M4 Max/MPS per the v1.2.0 changelog: ~5h for the
full 54-step chain) plus scANVI. Budget a full overnight run per cohort.

---

## FAIR / provenance conventions (kept consistent)

- Every new rule emits a provenance JSON via `fair_utils.stamp_artifact`/`write_provenance` under
  `provenance/<dataset>/`, naming `{dataset}_{sample}_ingest_provenance.json` and
  `{dataset}_{rule}_provenance.json` (extends the existing two-pattern scheme).
- New cohort's raw inputs get a `data/raw/<dataset>/MANIFEST.txt` (id/filename/md5/size) like
  `gdc_extract/MANIFEST.txt`.
- No hardcoded paths — all via `config["dirs"]` + `{dataset}` wildcard.
- When a cohort's result set is publication-ready, freeze its own `provenance/baseline_<dataset>_vX.Y.Z.json`
  via the existing `freeze_baseline_provenance` mechanism. Reference baselines stay untouched.
- CHANGELOG.md `[STATUS]` + a session entry updated; plan committed under `.claude/plans/` (per CLAUDE.md).

---

## Verification (end-to-end)

1. **Dry run:** `snakemake -n --use-conda` with the new `datasets:` block present → DAG includes the
   `<dataset>`-scoped rules and the reference cohort's targets are unchanged (no forced re-run).
2. **Reference integrity:** confirm `snakemake -n` reports **0** changed jobs for existing outputs
   (`integrated_latent.h5ad`, nerve/immune tables) and that the freeze SHA-256s in `config.yaml` still match.
3. **Stage A smoke test:** run only `data/processed/<dataset>/<sample>.h5ad` for one sample; assert it
   loads in scanpy, `.X` is raw counts, and `gene_symbol` mapping rate is reasonable (log the %, as
   `loom_to_h5ad.py:77` does).
4. **Full cohort run:** `snakemake --use-conda --cores all` to the dataset targets; verify each artifact
   exists and is non-empty (goal-backward, per CLAUDE.md).
5. **Concordance check:** open `cohort_concordance_summary.json` — confirm the shared-LR-pair overlap and
   Spearman correlation are computed and the comparison figure renders.

---

## Implementation status (2026-07-16) — CODE COMPLETE, not yet run

All code for the `gbm_cellxgene_56c4912d` track is written and the DAG is validated
(`snakemake -n`). Nothing heavy has been executed.

**New files**
- `config/config.yaml` — `datasets:` block (concrete entry).
- `scripts/derive_dataset_sample_sheet.py` — setup util; generated
  `data/raw/gbm_cellxgene_56c4912d/{samples.txt (170 donors), MANIFEST.txt}`.
- `workflow/scripts/ingest_dataset.py` — Stage-A loader (`h5ad` branch; other sources stubbed).
- `workflow/scripts/dataset_gene_symbol_map.py` — Ensembl→symbol map from the cohort's own `.var`
  (the reference MyGene cache does not cover Census unversioned IDs).
- `workflow/scripts/dataset_clinical_stub.py` — blank clinical TSV so `nerve_cell_subset`'s clinical
  join works without GDC metadata.
- `workflow/scripts/cohort_concordance.py` — Stage-D Jaccard + Spearman deliverable.
- `workflow/rules/datasets.smk` — full cohort-namespaced chain (16 `ds_*` rules), reusing every
  reference analysis script unchanged.
- `Snakefile` — `include` + `DATASETS` + terminal `cohort_concordance_summary.json` target.

**Scope decisions baked in**
- **scANVI-v2 retrain INCLUDED** (full parity, per researcher 2026-07-16): `ds_nerve_assemble_counts →
  ds_nerve_celltype_labels → ds_nerve_scanvi_retrain → ds_nerve_batch_qc_v2`. This is a nerve-latent
  validation side-branch — it does NOT feed the three-way interaction (which uses the v1
  `nerve_cells.h5ad`), but is run for apples-to-apples parity with the reference cohort. Adds a second
  MPS scVI+scANVI train (extra GPU time).
- **Freezes off** (fresh cut); **`exclude_clusters`/`exclude_subtypes` empty** (no artifacts diagnosed
  on this cohort yet — revisit after inspecting its purity tables).
- **`batch_key = donor_id`** (within-study batch, both scVI and scANVI).

**Reference protection — MANDATORY run flag.** The default `snakemake` would try to *re-run the 17-sample
reference cohort* (provenance/code/env rerun-triggers fire on the committed baseline that concordance
reads). Running with **`--rerun-triggers mtime`** isolates the replication track: verified that only the
16 `ds_*` rules run (354 jobs) and the reference deliverable reports "Nothing to be done."

**Run commands** (researcher-initiated; nothing auto-run). `caffeinate -i` keeps the Mac awake for the
duration; `--rerun-triggers mtime` is MANDATORY (protects the reference baseline):
```bash
# (optional) light Stage-A smoke via Snakemake — one donor, no scVI:
snakemake --use-conda --cores 4 --rerun-triggers mtime -- \
  data/processed/gbm_cellxgene_56c4912d/BT389.h5ad

# full replication run — concordance deliverable + scANVI-v2 parity branch
# (heavy; overnight — two MPS scVI/scANVI trains on ~150–300k cells):
caffeinate -i snakemake --use-conda --cores all --rerun-triggers mtime -- \
  results/tables/gbm_cellxgene_56c4912d/cohort_concordance_summary.json \
  results/tables/gbm_cellxgene_56c4912d/nerve_cluster_sample_purity_v2.csv
```

---

## Open items for execution (non-blocking)

- ~~**Which concrete second GBM dataset?**~~ **RESOLVED 2026-07-16** → `gbm_cellxgene_56c4912d`
  (`cellxgene_data/gbm_10x_raw.h5ad`, study `56c4912d`, `source: h5ad`, 5k/donor cap). See the
  "Selected replication cohort" section above. Stage-A ingest smoke-tested and passing.
- **Marker availability:** on a non-TCGA cohort, re-verify the nerve + immune marker symbols are present in
  `var_names` (the immune block comments already track presence rates) before trusting subtype calls.
