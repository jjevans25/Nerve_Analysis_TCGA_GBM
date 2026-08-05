# Plan — Compartment Integrity Fix (`scrna_annotate` → downstream)

**Status:** reviewed against code, approved for execution · **Arms:** `gbm_cellxgene_56c4912d_full` + `gbm_cellxgene_56c4912d`
**Trigger:** S1PR1 localization audit, `markdowns/s1pr1_localization_report.md` (2026-08-05)

> **Approved 2026-08-05.** A second pass verified every structural claim below against the code
> and the on-disk artifacts. All seven defects (D1–D7) are real and correctly located. That pass
> also found **four things that would have broken the fix as originally written** — see
> [Corrections after code review](#corrections-after-code-review-2026-08-05), which supersede the
> conflicting text in the Design section. Execution is gated: researcher review after Phase 1,
> researcher approval before Phase 5.

---

## Context

A question about whether S1PR1 sits on microglia or on malignant cells returned "neither — it is
endothelial," and in doing so exposed a compartment-labelling failure upstream of most of the
project's nerve-side conclusions.

**I reproduced every claim in that report independently against the on-disk artifacts**, and the
numbers match exactly. I also found the problem is **wider than the report states**.

### Verified findings

Joining `nerve_cells.h5ad` back to the Census author annotation (`cell_type`, preserved intact in
`malignancy_labeled.h5ad`):

| "nerve" compartment, full arm | n | share |
|---|---|---|
| malignant cell | 222,198 | **58.9%** |
| myeloid (macrophage/mono/microglia/DC/neutrophil/mast) | 103,924 | **27.5%** |
| genuinely neural | 41,809 | **11.1%** |
| vascular (endothelial + mural) | 7,058 | 1.9% |

Only **5 of 40** nerve clusters are dominantly neural (`c2, c25, c26, c28, c37` — 33,354 cells).
`nerve_c24` is **92.9% endothelial cell**, and both S1PR1 rows in
`nerve_tumor_immune_top_pairs_with_qc.csv` are `nerve_c24` (SPP1→S1PR1, immune→nerve and
tumor→nerve). The report's mechanism is exactly right.

### One finding beyond the report

The report says immune–tumor axes are unaffected. **Half true.** They do not use the nerve mask, but
the `tumor` compartment in *every* LIANA rule is `is_malignant == True`
(`nerve_tumor_immune_interaction.py:88`), and I measured that flag against the Census annotation:

- **recall 18.0%** (60,027 of 333,681 census-malignant cells flagged)
- **precision 47.9%** — the "tumor" compartment is **~41% myeloid** (25.5% macrophage, 10.6%
  microglia, 5.4% oligodendrocyte, 4.3% monocyte)

So immune–tumor axes are substantially immune–myeloid, i.e. within-compartment. The immune *mask*
itself is clean — I measured it at **99.5% pure** — but the compartment it is paired against is not.

### Root causes (four stacked defects, plus three carried)

| # | Defect | Evidence |
|---|---|---|
| **D1** | `sc.tl.score_genes` runs on **raw UMI counts** in the Census arms (`scrna_annotate.py:89`) | `mean_confidence` spans 0.13–24.2 across panels vs 0.11–0.73 in the reference arm. Already documented, still open: `markdowns/blocker_census_annotation_scoring.md` |
| **D2** | Labels assigned by **per-Leiden-cluster argmax** across uncalibrated panel scores, at res 1.0 → **47 clusters for 1,006,344 cells** (`scrna_annotate.py:108-115`) | 10 of 47 global clusters are <90% pure; cluster 4 (70,881 cells) is 82.5% malignant yet labelled `opc` |
| **D3** | Astrocyte panel = `GFAP, S100B, AQP4, VIM` — cannot separate normal astrocyte from AC-/MES-like malignant. No panel exists at all for macrophage, mural, mast, B, NK, neutrophil, DC | `cell_type_predicted` calls **322,152** cells astrocyte; the Census authors call **304** of them astrocyte. Nerve clusters dominated by mast (c33), mural (c21), B cells (c22) all forced to "astrocyte" |
| **D4** | CNV caller broken twice: gene order from an **Ensembl-ID numeric-suffix proxy, not chromosomal coordinates** (`scrna_malignancy.py:42-48`); and the "normal reference" set is `cell_type_predicted ∈ {t_cell, endothelial}` — drawn from the same broken annotation (`:71-72`) | 18% recall / 48% precision, measured above |
| D5 | OPC string-match bug — config says `"oligodendrocyte precursor cell"`, annotate emits `"opc"`, neither substring matches → 70,881 cells silently dropped (`nerve_cell_subset.py:59-75`) | No `opc` category exists in either arm's `nerve_cells.h5ad` |
| D6 | `nerve_celltype_labels.py:126-128` **overwrites `cell_type`**, destroying the Census CL-ontology column — the only external ground truth | `nerve_cells_v2.h5ad` `cell_type` holds pipeline strings; audit required joining back to the parent object |
| D7 | Placeholder-on-empty still live in `nerve_cell_heterogeneity.py:81` and `nerve_cluster_annotations.py:88` | Same laundering pattern already hard-failed in `nerve_cell_subset.py:86-97` |

> **⚠ D5 must never be fixed alone.** Global Leiden cluster 4 — the 70,881 cells labelled `opc` — is
> **82.5% malignant**. Repairing the string match without D4 would *add* ~58,000 malignant cells to
> the nerve compartment and make the problem worse.

### Decisions taken (researcher, 2026-08-05)

1. **Full fix from annotate down** — repair nerve + tumor + immune together.
2. **Neuroglial compartment with neurons broken out** as a separate LIANA group.
3. **Both Census arms** (full + capped), so the depth comparison stays valid.
4. Gene coordinates for D4 sourced from a **new Ensembl 113 GTF rule** + `infercnvpy` (C1 below).
5. **Bundle** `markdowns/task_conda_env_enforcement.md` (Defect 1) into this work, in Phase 4.
6. Approval gates **after Phase 1** and **before Phase 5**; Phases 2–4 run through.

---

## Corrections after code review (2026-08-05)

Verified against `workflow/rules/datasets.smk`, the scripts, and the on-disk artifacts. **These
supersede the conflicting text in Design and Execution below.**

Confirmations worth recording, since they are load-bearing:

- `annotation_summary.csv` reproduces this plan's numbers exactly — full arm 1,006,344 cells,
  `opc` = 70,881, `mean_confidence` spanning 0.13 → 24.2 across panels.
- The **capped arm assigned only 5 labels** (astrocyte, microglia, t_cell, oligodendrocyte,
  neuron). No `endothelial`, no `opc`, no `ependymal` — so its CNV normal reference
  (`scrna_malignancy.py:71`, `{t_cell, endothelial}`) silently degraded to T cells alone.
- `results/pinned_reference_verification.json` currently reads `pass: true, 37/37`.

### C1 — There is no gene-coordinate source in this project (blocks D4)

The D4 section below says the symbol map "already carries a `chromosome` column … extend
`ds_gene_symbol_map` to emit `start` too." **The column exists but is a literal placeholder:**
`dataset_gene_symbol_map.py:35` writes `chromosome = "unknown"` for every row, and the on-disk map
confirms it — 61,497/61,497 rows read `unknown`. Census `.var` carries `feature_id` /
`feature_name` / `feature_length` only, no coordinates. D4's hard gate is unreachable without a
real source.

**Fix (approved):** new rule **`ds_gene_positions`**, modelled on the existing
`download_msigdb_gmt` pattern in `fair.smk` (download → SHA-256 manifest → committed digest):

- Download the Ensembl **release 113** GTF (`config.databases.ensembl_release` already pins 113).
- Emit `data/external/ensembl113_gene_positions.tsv` — `ensembl_id, chromosome, start, end`.
- New `gene_positions:` config block (`url`, `filename`, `sha256`, `request_timeout`,
  `max_retries`), mirroring `msigdb:`.
- `scrna_malignancy.py` joins on unversioned Ensembl ID, sorts by `(chromosome, start)`, and hands
  off to **`infercnvpy`** — already declared and unused at `workflow/envs/scrna.yaml`
  (`infercnvpy==0.4.3`), and already claimed by that script's docstring.
- Genes that fail to map are dropped from the CNV pass and the count logged, not silently
  `-1`-ordered as today (`scrna_malignancy.py:56`).

### C2 — New marker panels must NOT go into `nerve_cells.markers`

The D1/D2/D3 section says to move the four hardcoded panels and the seven new ones into
`config.yaml` "alongside" the existing ones — i.e. into `nerve_cells.markers`. **That key is not
private to annotation.** `datasets.smk:165` passes it as `params.markers` to **`ds_scrna_qc`**,
which runs per-sample across 170 donors × 2 arms; changing it invalidates every QC artifact, which
invalidates `ds_scrna_integration` — **the 11 h 13 m scVI train this whole plan is built around
preserving.** It also feeds `ds_nerve_cell_heterogeneity`, `ds_nerve_cluster_annotations` and
`ds_nerve_celltype_labels`. It is semantically wrong besides: `scrna_qc` uses these markers to set
`is_nerve_marker`, and `CD163`/`MRC1` are not nerve markers.

**Fix:** a **separate top-level `annotation_markers:` block**, consumed *only* by
`ds_scrna_annotate` (and its reference twin in `annotation.smk`), holding the four
currently-hardcoded panels plus the seven new ones. `scrna_annotate.py` merges `params.markers`
(neural, unchanged) with `params.annotation_markers`. `nerve_cells.markers` stays
**byte-identical**, so no QC job is invalidated under any rerun-trigger setting.

### C3 — A `macrophage` panel silently guts the immune compartment

`immune_cell_subset.py:55` selects `cell_type_predicted == source_label` — an **exact match on a
single string**, `immune_cells.source_label: "microglia"`. Today all 361,112 myeloid cells land in
that one bucket, which is *why* the immune mask measures 99.5% pure. The moment `macrophage`,
`dc`, `neutrophil`, `mast`, `b_cell` and `nk_cell` panels exist, those cells get their own labels
and **fall out of the immune subset entirely**. Verification item 5 below ("immune purity ≥95%")
would still pass — purity is not a size check — while the compartment quietly loses most of its
cells.

**Fix:** `immune_cells.source_label` (str) → `source_labels` (list) covering every myeloid/lymphoid
label the new panels can emit, with a str fallback for backward compatibility. Add a **size** gate
to `ds_compartment_audit` alongside the purity gate: no more than ~15% loss vs today.

### C4 — Per-arm scoping and headroom

- **`nerve_scanvi.labels_key` is global.** Repointing it to `cell_type_marker_label` (D6) also hits
  the reference cohort's `nerve_scanvi_retrain`, which is *not* pinned and whose on-disk
  `nerve_cells_counts_labeled.h5ad` carries `cell_type`. Make it a per-arm override
  (`datasets.<arm>.scanvi_labels_key`, defaulting to the global value).
- **`nerve_cells.leiden_resolution` is read by four rules** — `ds_scrna_annotate`,
  `ds_nerve_cell_subset`, `ds_nerve_batch_qc_v2`, plus reference twins. The split must introduce
  `annotate_leiden_resolution` and leave the other three on the existing key.
- **Per-arm expected sizes.** The "~55–62k" figure below is a full-arm number. The capped arm has
  no `opc`/`ependymal` today and 615k cells; its gate needs its own range, derived from its own
  Census `cell_type` counts in Phase 1 rather than assumed now.
- **Disk headroom is thin: 76 GB free**, against 117 GB across the two arms. Rewrites are in place
  so the run is net-neutral, but preservation must stay SHA-256-only for the `.h5ad`s.

---

## Answering "does this force a re-cluster?"

**Yes — but not a full pipeline re-run.**

`nerve_leiden` is computed *inside* the nerve subset (`nerve_cell_subset.py:153-156`, on `X_scVI`),
so any change to the mask changes the kNN graph and therefore every cluster ID. Everything keyed on
cluster ID invalidates.

**What is preserved:** ingest, QC, and — critically — `ds_scrna_integration`, the **11 h 13 m scVI
train**. `X_scVI` is unsupervised and label-free, so it is unaffected by any of D1–D7. Re-entry is at
`ds_scrna_annotate`, which recomputes neighbors/UMAP/Leiden from the existing latent in ~15 min.

Estimated **6–10 h full arm + 4–7 h capped arm**, against ~31 h for the original run.

---

## Design

### Corrected compartment definition

`config/config.yaml` → `nerve_cells:`

```yaml
cell_types:                      # neuroglial compartment
  - neuron
  - excitatory_neuron
  - inhibitory_neuron
  - opc                          # canonical label, matches what annotate emits
  - oligodendrocyte
  - astrocyte
  - ependymal
  - radial_glial
neuron_labels: [neuron, excitatory_neuron, inhibitory_neuron]   # broken out as its own LIANA group
```

`nerve_cell_subset.py` gains an `obs["nerve_subcompartment"]` column (`glia` / `neuron`).
In `nerve_tumor_interaction.py:83` and the three-way equivalent, the label mint becomes:

- glial cells → `glia_c{nerve_leiden}` (per-cluster, as today)
- neurons → a single `neuron` group (too few for per-cluster resolution)

**Expected post-fix size, checkable against Census counts (full arm):** oligodendrocyte 34,160 +
OPC 21,460 + neuron 3,401 + radial glial 3,264 + astrocyte 347 ≈ **62,600 before malignancy
exclusion**; expect **~55–62k** after. Down from 377,343.

**Keep the `compartment` string as `"nerve"`.** The concordance pair key is
`(source_compartment, target_compartment, ligand_complex, receptor_complex)`; renaming it to
`neuroglia` would break every comparison against the pinned v1.3.0 reference. Record the definition
change in provenance and in the audit table instead.

### Fixes, per defect

**D1/D2/D3 — `workflow/scripts/scrna_annotate.py`**
- Score on a normalized working copy (`normalize_total(1e4)` + `log1p`), guarded by the **existing**
  `counts_utils.is_log1p_scale()` helper built for the 2026-07-26 LIANA fix — reuse it, do not
  rewrite. Leave `.X` raw for downstream consumers.
- Z-score each panel score across cells before argmax so panels with high-expression genes stop
  winning by scale.
- Replace VIM in the astrocyte panel with `SLC1A2, SLC1A3, ALDH1L1, GJA1`; keep AQP4.
- Add the missing panels so non-neural cells have somewhere to go: macrophage
  (`CD163, MRC1, MSR1, F13A1`), mural (`PDGFRB, RGS5, ACTA2, DCN`), mast (`TPSAB1, CPA3, KIT`),
  B/plasma (`CD79A, MS4A1, MZB1`), NK (`NKG7, GNLY, KLRD1`), neutrophil (`FCGR3B, CSF3R, S100A8`),
  DC (`FLT3, CLEC9A, LAMP3`). Move the four currently hardcoded panels
  (`scrna_annotate.py:45-49`) into `config.yaml` alongside them — **into the new
  `annotation_markers:` block, NOT `nerve_cells.markers`; see C2.**
- **Widen `immune_cells.source_label` to `source_labels` (list) in the same change — see C3.**
  Adding a macrophage panel without this drops most of the immune compartment on the floor.
- Raise annotate `leiden_resolution` for the Census arms (1.0 → 2.0) so mixed clusters can separate;
  make it a per-arm key rather than sharing `nerve_cells.leiden_resolution` across two different
  clustering passes.
- Add an `ambiguous` label when top-minus-second margin falls below a configured floor, instead of
  forcing a call.

**D4 — `workflow/scripts/scrna_malignancy.py`** *(highest-risk item)*
- Order genes by **real chromosome + start coordinate**. ~~The symbol map already carries a
  `chromosome` column (`scrna_annotate.py:65`); extend `ds_gene_symbol_map` to emit `start` too.~~
  **Superseded by C1** — that column is the literal string `"unknown"` for all 61,497 genes in both
  Census arms. Coordinates come from the new `ds_gene_positions` (Ensembl 113 GTF) rule instead.
- Draw the normal reference set from post-fix high-confidence immune cells, and assert it is
  non-empty (the capped arm never assigned `endothelial` at all, silently degrading the reference to
  T cells only).
- **Preferred:** replace the hand-rolled smoother with `infercnvpy`, which the docstring already
  claims is in use and which is declared but unused in `workflow/envs/scrna.yaml`.
- **Hard gate:** must reach **≥0.80 recall and ≥0.85 precision** against Census `cell_type` before
  any downstream rule runs.

**D5/D6/D7**
- Replace substring matching with an explicit canonical label set, plus an assertion that every
  configured `cell_types` entry matches at least one observed label — fail loudly on a typo rather
  than silently selecting nothing.
- `nerve_celltype_labels.py` writes to `cell_type_marker_label`; set
  `config.nerve_scanvi.labels_key` to match. **Never write to `cell_type`.**
- Hard-fail the two remaining placeholder paths, consistent with `nerve_cell_subset.py:86-97`. This
  matters more now: with a ~62k compartment and a ~3.4k neuron group, empty/tiny groups are far more
  likely than before.

### New rule — `ds_compartment_audit` (the Test Oracle)

CLAUDE.md mandates comparing against a known reference before marking complete. This rule makes the
S1PR1 finding **reproducible** — its three companion CSVs were never written to disk, so the audit
currently exists only as prose.

Cross-tabs every mask against Census `cell_type`; runs on existing artifacts, no re-run needed.
Outputs: `compartment_audit.csv`, `nerve_compartment_cluster_audit.csv`, `malignancy_confusion.csv`.

Config-driven gates (`compartment_audit:` block), failing the rule when unmet:

| gate | now | required |
|---|---|---|
| nerve compartment neural fraction | 11.1% | ≥85% |
| tumor compartment malignant fraction | 47.9% | ≥85% |
| malignancy recall | 18.0% | ≥80% |
| immune compartment purity | 99.5% | ≥95% (regression guard) |

### Bundle these two open items now

- **`markdowns/task_conda_env_enforcement.md`** — blocked solely because fixing it "may shift Leiden
  cluster IDs → needs sign-off." **This re-run discards every cluster ID anyway**, so its cost has
  dropped to roughly zero. Closing it here also clears the standing `[FAIR-ALERT]` that every
  artifact in this arm was produced by `claude_science/bin/python` rather than the declared
  `scrna.yaml` env.
- **Pin `numba` in `workflow/envs/scrna.yaml`** — currently unpinned; it caused both SIGSEGVs.

---

## Execution phases

| # | Phase | Compute | Gate |
|---|---|---|---|
| 0 | Sync this document with corrections C1–C4 | — | Atomic commit |
| 1 | Snapshot current results; add `ds_compartment_audit` (`enforce: false`); run it on existing artifacts to record the failing baseline numerically | ~10 min | Audit reproduces 11.1% / 58.9% / 27.5%. **⟵ STOP FOR RESEARCHER REVIEW** |
| 2 | D1/D2/D3 annotate fixes (incl. `annotation_markers:` per C2, `source_labels` per C3); flake8; dry-run on one sample | ~30 min | `flake8` exit 0 |
| 3 | `ds_gene_positions` (C1) + D4 CNV rebuild on `infercnvpy`; validate against oracle | ~1 h | **recall ≥0.80, precision ≥0.85 — stop if unmet** |
| 4 | D5/D6/D7 + compartment redefinition + per-arm scoping (C4) + conda-env enforcement + numba pin; flip `compartment_audit.enforce: true` | — | `flake8` exit 0; conda smoke test passes |
| 5 | Re-run **full** arm from `ds_scrna_annotate`; then **capped** arm | 6–10 h / 4–7 h | **⟵ STOP FOR RESEARCHER APPROVAL BEFORE STARTING.** All audit gates pass |
| 6 | Re-derive shortlist; re-check S1PR1, CXCR4, LRP1; refresh concordance | ~1 h | — |

### Preservation before Phase 5

- Snapshot `results/tables/<arm>/` for both arms (small) and record SHA256 for the large `.h5ad`s —
  do not duplicate 77 GB.
- **Pinned v1.3.0 reference must remain untouched.** It is structurally protected (producing rules
  are undefined at parse time when `baseline.pinned: true`), but verify
  `results/pinned_reference_verification.json` reads `pass: true, 37/37` both before and after.
- Commit atomically per CLAUDE.md; the branch is kept, never deleted.

### Operational notes (from the 2026-08-02 run — do not rediscover)

- `NUMBA_THREADING_LAYER=workqueue` must be set **before** any import that pulls in numba, in every
  script that calls `sc.pp.neighbors` after torch/MPS init.
- `--allowed-rules` takes `nargs='+'` and swallows targets placed after it — **pass targets first**.
- `--rerun-triggers=mtime` needs the `=`.
- tmux is not installed; use `nohup` + `caffeinate -ims`.
- `resources: mem_mb` is a scheduler gate only; it cannot cap a single process. Never declare above
  36000.
- Long-training rules should write to an **undeclared sidecar** — Snakemake deletes declared outputs
  of failed jobs.

---

## Verification (Goal-backward)

1. `compartment_audit.csv` — all gates pass, **both arms**.
2. Nerve compartment ≥85% neural by Census annotation; n ≈ 55–62k for the **full** arm, capped-arm
   range fixed from the Phase 1 numbers (C4).
3. **No nerve cluster exceeds 20% endothelial** — the `nerve_c24` failure mode is gone.
4. **Zero nerve-side S1PR1 rows** in `nerve_tumor_immune_top_pairs_with_qc.csv`.
5. Immune compartment purity still ≥95% (regression guard against today's 99.5%) **and ≥85% of its
   current size** — purity alone does not catch the C3 failure mode.
6. Malignancy recall ≥0.80, precision ≥0.85 against Census `cell_type`.
7. Neuron group present, n ≈ 2.5–3.4k (full arm), carried as its own LIANA group.
8. `cell_type` intact in every artifact including `nerve_cells_v2.h5ad`.
9. `results/pinned_reference_verification.json` still reads `pass: true, 37/37`.
10. `flake8` exit 0 on every edited script; provenance JSON per rule; `CHANGELOG.md` entry.

## Expected outcomes to interpret, not treat as regressions

- **Cohort concordance against the v1.3.0 reference will move, probably down.** The reference's own
  nerve compartment was built by this same uncorrected logic on a different data source. A corrected
  compartment scored against an uncorrected reference is not a like-for-like comparison. Report it,
  but demote it from headline until the reference's status is resolved — it carries no author
  annotation, so it is currently *unknown*, not cleared.
- **The 40-axis curated shortlist (`results/tables/nerve_crosstalk_lead_targets.csv`) will need
  re-derivation.** Note S1PR1, CXCR4 and LRP1 are **not on it** — they appear only in the raw LIANA
  tables, and S1PR1 in just 2 rows of the full arm. Whatever shortlist named them was produced
  outside this repo.
- Batch-purity QC will flag the `neuron` group. That is expected — it is a single group of ~3k cells
  across 170 samples, not a cluster. Exempt it explicitly with a written justification rather than
  letting it read as a failure.

## Known limits of this fix

- Census `cell_type` is an **external author annotation**, harmonized to CL ontology from
  per-study labels of varying provenance — not a re-derivation from these counts. It is the best
  available oracle, not truth.
- It likely **under-calls normal astrocytes** (347 in a 1M-cell cohort), because most GBM study
  authors label AC-like cells "malignant cell". A corrected pipeline may legitimately find more
  astrocytes than the oracle does; that specific disagreement should not be auto-failed.
- Separating malignant AC-/OPC-like cells from normal glia is **intrinsically hard** — they share
  GFAP, PTPRZ1, SLC1A3. This is why D4 (a genuine CNV caller) carries the weight, and why marker
  panels alone cannot close it.
- The cohort contains only **3,401 neurons (0.34%, ~20/sample)**. Any strictly neuronal
  neuro-oncology hypothesis is depth-limited in this data regardless of how well the compartment is
  fixed.
