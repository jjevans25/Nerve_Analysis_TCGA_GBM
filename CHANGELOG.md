# CHANGELOG.md — Lab Notebook & Persistent Memory
# Autonomous Agentic Biomedical Research Environment

This file is the **persistent lab notebook** for all agentic sessions. It serves as long-term memory across context windows. Every significant computational decision, tool failure, parameter choice, and artifact verification must be recorded here.

Format each entry with: date, phase, action taken, outcome, and any open issues.

---

## How to Use This File

- **Agents:** Append an entry at the start and end of every session. Record failures, not just successes.
- **Researcher:** Use this as an audit trail. Entries here are the ground truth for what the system has done.
- **GSD Orchestrator:** Read the most recent `[STATUS]` block to resume after context resets.

---

## Entry Format

```
### [YYYY-MM-DD] | Phase: <phase name> | Status: IN-PROGRESS | COMPLETE | FAILED
**Action:** <what was attempted>
**Outcome:** <what actually happened>
**Artifacts:** <file paths produced or modified>
**Tool Versions:** <key package versions used>
**Open Issues:** <unresolved problems or follow-ups>
**FAIR Notes:** <any FAIR compliance actions taken or deferred>
```

---

## Current Project State

```
[STATUS]
Phase:          Branch merge to main + documentation resync — COMPLETE
Last Updated:   2026-07-28
Active Agent:   lead-researcher
Current Task:   feat/nerve-tumor-immune-interaction merged to main and pushed (main @ 2cb6776; reconciled with the stale PR #1 merge commit — no content change). README.md and execution_instructions.md rewritten to current scope: immune arm, LIANA crosstalk, replication cohort, structural v1.3.0 pin. Removed guidance that `--forcerun nerve_cell_subset` can re-cut the reference cohort (that rule is undefined while pinned). No analysis code ran this session.
Blocked On:     Nothing. Prior science state unchanged: both explorers (04 reference, 05 Census) runnable, inputs valid.
Next Action:    Researcher review of the 30/40 curated-axis replication and the refreshed concordance (Jaccard 0.4691 / ρ 0.6249 — the old 0.4248/0.5357 figures are void). Then decide on markdowns/blocker_census_annotation_scoring.md (NEW, open — expensive: re-clusters the Census cohort and invalidates the above). Still open from 2026-07-25: review of the 9 flagged micro-clusters. Note `marimo` on PATH has a broken matplotlib; run notebooks via the snakemake notebooks env or `claude_science/bin/python -m marimo`.
```

Prior status (v1.3.0 baseline, retained): cl15 surgical sub-cluster split COMPLETE; freeze insulator retained; cluster set {0-14, 16-27}.
```

---

## Session Log

---

### [2026-07-28] | Phase: Branch merge to main + documentation resync | Status: COMPLETE
**Action:** Merge `feat/nerve-tumor-immune-interaction` into `main`, push to origin, and bring
`README.md` and `execution_instructions.md` back in line with the current pipeline.
**Outcome:** Merge into `main` was a fast-forward locally, but the push was rejected: GitHub PR #1
had merged the same branch on 2026-07-11 from an **older tip** (`076bd9f`), so the 5 later commits
(`d8b384e`…`3ad59c5`) had never reached the remote while the PR's merge commit was absent locally.
Verified `37db19d` was already an ancestor of `076bd9f` and `git diff 076bd9f 9d3b314` was empty —
i.e. the merge commit carried no unique content — then reconciled with a merge commit (`2cb6776`)
rather than a rebase, since the 5 commits were already published on the remote feature branch.
`git diff 3ad59c5 HEAD` empty, confirming no file content changed. Pushed.

Both docs were stale at the *scope* level: they described only the nerve-heterogeneity arm of a
17-sample TCGA cohort, with no mention of the immune compartment, LIANA crosstalk, the replication
cohort, or the structural v1.3.0 pin. Rewrote scope, structure, DAG, outputs, config, and caveats.
Three factual corrections made against live data rather than prior notes: (1) the Census LIANA
raw-counts blocker is RESOLVED, not open; (2) the 24/33 purity figure is the *replication* cohort —
the reference is 16/26; (3) `nerve_cells.h5ad` and `nerve_cells_umap.png` were listed as available
outputs but do not exist, so they are now struck through and annotated as unreproducible.

Most consequential doc fix: `execution_instructions.md` instructed `--forcerun nerve_cell_subset`
for resolution tuning. Under `baseline.pinned: true` that rule is not defined, so the command fails —
and the guidance implied the reference cohort could be re-cut, which it cannot. Replaced with the
`ds_*` replication-cohort equivalent plus the mandatory `--allowed-rules` insulator pattern.
**Artifacts:** `README.md` (+192/−61), `execution_instructions.md`, `CHANGELOG.md`; `main` at
`2cb6776` on origin.
**Tool Versions:** git 2.x; no analysis code executed — documentation and VCS only.
**Open Issues:** PR #1 remains merged-but-stale on GitHub (shows only the pre-`076bd9f` work);
`scrna_annotate.py` marker scoring on the Census cohort is still unaudited
(`markdowns/blocker_census_annotation_scoring.md`). The README does not enumerate the 21 `ds_*`
rules individually.
**FAIR Notes:** No artifacts produced or modified; provenance untouched (65 JSONs). Documentation
now states the pinned-reference drift-detection procedure (`verify_pinned_reference`) and the
per-cohort `.X` scale contract, both of which were previously undocumented and had already caused
one silent defect.

---

### [2026-07-26] | Phase: Census cohort LIANA normalization fix + replication explorer | Status: COMPLETE

**Action:** Unblock exploration of the CELLxGENE Census replication cohort. Its ligand–receptor
tables were void — computed on raw UMI counts against LIANA's log1p assumption (see the
`[2026-07-26] TME × Nerve × Immune` entry below for the discovery).

**Root cause, corrected framing.** Raw `.X` is *correct* for most of this pipeline: scVI wants
counts, which is why `scrna_integration` sets `counts_from_log1p=True` for the reference
(recovering counts from SCT log1p) and passes Census raw straight through. Integration and
clustering were never wrong. The defect is that two *consumers* requiring log1p —
`nerve_tumor_interaction` and `nerve_tumor_immune_interaction` — never normalized for themselves
and only **logged** the assumption. Fixing it per-consumer is therefore the right design, not
just the cheap one: normalizing upstream would break scVI's input contract.

**Outcome:**
- **`counts_utils.py`** — new `is_log1p_scale()` + `LOG1P_MAX_PLAUSIBLE = 50.0`.
- **Both LIANA scripts** — `normalize_counts` param runs
  `normalize_total(target_sum=1e4)` + `log1p` on the **combined** AnnData (post-concatenation, so
  all compartments share one scale) just before the LIANA call. The passive "assuming log1p" log
  line is replaced by a **hard guard that raises** on an off-scale matrix, for either cohort.
  This is the durable part: the silent-wrong-answer mode is now impossible.
- **Wiring** — `normalize_counts = _entry(wc.dataset).get("normalize_counts", True)` on the two
  `ds_*` rules; pinned `False` on the reference twins (already SCT log1p).
- **New `ds_nerve_cluster_annotations`** — the cohort had no cluster → cell-type map, so its
  results could only be read as opaque `nerve_c{N}` ids. 35 clusters annotated.
- **New `notebooks/05_census_nerve_immune_explorer.py`** + `ds_census_nerve_immune_notebook`.
  A separate notebook rather than a cohort switch on 04, per researcher preference and because
  the cohorts differ in what exists (no clinical metadata, no per-cohort curated list, 169 donors
  vs 17). Panels: compartment census, per-cluster patient purity (bar charts, not the 27×169
  heatmap that would be unreadable at this cohort size), biology-labelled LR browser,
  cell-type × immune-subtype matrix, **curated-axis replication**, and whole-table concordance.
  It self-checks for the raw-counts signature at load and refuses to vouch for defective data.
- **`run_notebook_export.py`** — optional `params.env` passthrough so a cohort-namespaced rule can
  tell the notebook which dataset to read (marimo has no argv passthrough). Additive; the five
  existing notebook rules are unaffected.

**Verified (re-run: 5 jobs, all `ds_*`, reference untouched, ~7 min):**

| Marker | Before | After |
|---|---:|---:|
| `X.max()` at the LIANA call | 53,027.000 | **8.773** |
| three-way rows with empty `specificity_rank` | 71,189 / 71,189 | **0** |
| two-way rows with empty `specificity_rank` | 25,672 / 25,672 | **0** |
| `lr_logfc = inf` (three-way) | 13,492 | **0** |
| Jaccard vs reference | 0.4248 | **0.4691** |
| Spearman ρ (shared pairs) | 0.5357 | **0.6249** |
| Census significant pairs | 3,850 | **4,583** |

**Scientific result.** Of the reference cohort's 40 curated lead axes, **30 replicate** in this
independent 169-donor cohort at `magnitude_rank ≤ 0.05`, 9 are present but not significant, and 1
(`EGF | RHBDL2`) is absent. Strongest replicating axes: `CD44|SPP1` (133 significant rows),
`APP|CD74` (133), `ABCA1|APOE` (95), `PTN|PTPRZ1` (75), `NCAM1|PTPRZ1` (50), `NLGN1|NRXN1` (44) —
the synaptic-adhesion and PTPRZ1 programs carry across cohorts.

**Artifacts:** refreshed `results/tables/gbm_cellxgene_56c4912d/{nerve_tumor_interactions,
nerve_tumor_immune_interactions,*_top_pairs,*_with_qc}.csv`, `cohort_concordance_summary.json`,
`cohort_concordance_shared_pairs.csv`, new `nerve_cluster_annotations.csv`,
`results/figures/gbm_cellxgene_56c4912d/05_census_nerve_immune_explorer.html` (881 KB, 4 figures,
0 errors), `notebooks/05_census_nerve_immune_explorer.py`,
`markdowns/blocker_census_annotation_scoring.md`.

**Open Issues:**
- **NEW blocker split out:** `markdowns/blocker_census_annotation_scoring.md`. The same raw-counts
  family affects `scrna_annotate`'s `score_genes` — Census `mean_confidence` 1.41–15.29 vs the
  reference's 0.11–0.73, and only 5 cell types resolved vs 10 (no `endothelial`, `opc`,
  `ependymal`, `tumor_gbm`). **Not fixed:** the re-run re-clusters the cohort and would invalidate
  everything above, including the replication finding. Filed so exploration is not blocked behind
  that decision. Consequence already visible: 55,477 marker-labelled T cells (9.0% of the cohort)
  sit outside the immune compartment, because `immune_cell_subset` takes only `microglia`-labelled
  cells. Surfaced in notebook 05 Panel A.
- Cluster-level cell-type interpretation for this cohort is provisional until that is resolved.
  The LR results do not depend on the labels being correct, only on the compartment split being
  stable, so the replication finding stands on its own.
- Pre-existing `E302` in `nerve_tumor_interaction.py` fixed in passing (file was already being
  edited); no other reformatting.

**FAIR Notes:** every re-run artifact re-stamped with provenance JSON recording
`normalize_counts`. The baseline pin held throughout — no pinned rule was scheduled at any point,
confirmed by dry run before the re-run and by unchanged reference mtimes after.

---

### [2026-07-26] | Phase: Baseline pin — make the v1.3.0 freeze structural | Status: COMPLETE

**Action:** The decision to never regenerate the v1.3.0 reference cohort was enforced only by
remembering to pass `--allowed-rules` on every invocation. `Snakefile:57` requested
`results/figures/nerve_cells_umap.png` unconditionally — one of the four outputs of the failed
`nerve_cell_subset` job, and also missing — so a bare `snakemake` scheduled the producer and
cascaded through ~60 artifacts, overwriting the pinned tables. Four of the five notebook targets
did the same (7–9 jobs each). Requested fix was `ancient()`; that turned out not to work.

**Finding — `ancient()` cannot pin a missing file.** Verified against the installed Snakemake
9.20.0 rather than assumed. The flag has exactly two semantic consumers and both concern
timestamps: `io/__init__.py:763-769` (`is_newer` → `False`) and `dag.py:1630-1640`. The
missing-file decision runs down a separate, ancient-blind path — `dag.py:1620-1629` queues a
producer purely from `job_.missing_output(files)`, and `jobs.py:725-737` decides that on
`not await f.exists()` alone. Even the one propagation-skip that reads the flag is gated on
existence: `if all([f.is_ancient and await f.exists() for f in files])` (`dag.py:1633`). Net:
ancient input + missing file + a producing rule ⇒ the producer is scheduled normally.

**Outcome — parse-time conditional rule definition instead.** A rule that is never defined cannot
be scheduled by any invocation, and Snakemake treats an existing file with no producer as a source
file (`dag.py:1300-1310`). The project already depended on that: `nerve_crosstalk_lead_targets.csv`
has no producing rule anywhere and resolves fine.

- **`config.yaml`** — new `baseline.pinned: true` plus `baseline.pinned_artifacts` enumerating the
  **37** surviving outputs of the pinned rules.
- **`common.smk`** — `BASELINE_PINNED` flag and a `pinned_target()` helper that omits a frozen
  artifact from `rule all` when it no longer exists (otherwise Snakemake raises "No rule to
  produce" for the three files lost with the h5ad).
- **10 rules wrapped in `if not BASELINE_PINNED:`** — nine in `nerve_cells.smk`
  (`nerve_cell_subset`, `nerve_cell_heterogeneity`, `nerve_cluster_annotations`,
  `nerve_clinical_association`, `nerve_batch_qc`, `nerve_leiden_resolution_sweep`,
  `nerve_tumor_interaction`, `nerve_assemble_counts`, `annotate_cluster_qc`) and
  `nerve_tumor_immune_interaction` in `immune.smk`.
- **`Snakefile`** — 19 `rule all` targets converted to `pinned_target()`.
- **Integrity manifest** — `freeze_pinned_reference` / `verify_pinned_reference` rules in
  `fair.smk` + two new scripts, writing `provenance/pinned_reference_v1.3.0.json` and
  `results/pinned_reference_verification.json`. A fresh manifest was required because
  `baseline_v1.3.0.json` is dated 2026-05-25, several tables were legitimately regenerated in
  July, and it stores `artifact_sha256: null` for `annotate_cluster_qc`.

**Measured result (dry runs, no `--allowed-rules`):**

| Target | Before | After |
|---|---:|---:|
| `04_tme_nerve_immune_explorer.html` | 9 jobs | 0 (up to date) |
| `03_nerve_tumor_immune_explorer.html` | 7 | 0 |
| `02_nerve_enrichment_explorer.html` | 7 | 1 (own export only) |
| `nerve_tumor_exploration.html` | 7 | 1 (own export only) |
| bare `snakemake -n --rerun-triggers mtime` | schedules `nerve_cell_subset` | 5 jobs, zero pinned |
| bare `snakemake -n` (default triggers) | 420 | 410, zero pinned |

**Two honest qualifications.**
1. **A bare `snakemake` is still not cheap.** With default rerun-triggers it plans ~410 jobs —
   full re-ingest of both cohorts, scVI/scANVI retrains, the entire `ds_*` chain. That is stale
   Snakemake metadata, pre-existing and unrelated to this change. The pin removes exactly the 10
   pinned rules; it protects the v1.3.0 tables, not the compute budget. `--rerun-triggers mtime`
   still reduces the same command to 5 jobs.
2. **One verification step was wrong on the first pass.** The `snakemake --list` check used
   `grep -x <rulename>`, but `--list` appends each rule's docstring, so it matched nothing in
   *either* pin state and appeared to pass vacuously. Re-run with an anchored pattern it is a
   real test: 10 rules listed with `pinned: false`, 0 with `pinned: true`.

**Verification:** workflow parses; anchored `--list` check passes in both states; reversibility
confirmed by flipping `pinned: false` (rules return, 04 target back to 9 jobs) and back;
`freeze_pinned_reference` → 37 artifacts hashed; `verify_pinned_reference` → PASS 37/37 unchanged;
**negative test** — injected a bad sha and deleted a manifest entry, verifier correctly hard-failed
on both (`1 modified, 0 missing, 1 configured but not in manifest`) and the manifest was restored
to a clean pass. `flake8` clean on both new scripts and notebook 04. All 202 tracked mtimes under
`results/`, `data/processed/` and `provenance/` unchanged except the intended new artifacts.

**Artifacts:** `workflow/scripts/{freeze,verify}_pinned_reference.py`,
`provenance/pinned_reference_v1.3.0.json` (37 entries),
`results/pinned_reference_verification.json`, `markdowns/plan_pin_v1_3_0_reference.md`.

**Open Issues / residual risk:**
- `annotate_cluster_qc` is pinned, so edits to `nerve_cells.batch_qc.exclude_clusters` or
  `immune_cells.batch_qc.exclude_subtypes` no longer regenerate the six `_with_qc` tables for the
  reference cohort. Unpinning is the escape hatch; the `ds_*` track is unaffected.
- The scANVI-v2 chain (`nerve_celltype_labels`, `nerve_scanvi_retrain`, `nerve_batch_qc_v2`) stays
  runnable by choice and *is* scheduled under default rerun-triggers — it would retrain scANVI and
  rewrite `nerve_cluster_sample_purity_v2.csv`. Reproducible, but not free.
- `results/pinned_reference_verification.json` is deleted by Snakemake when the verifier fails, so
  the report is also embedded in the FAILURE entry of `logs/verify_pinned_reference.log`.
- `datasets.<name>.freeze_nerve_subset` and `freeze_cl15_split` (`config.yaml:145-146`) are read by
  nothing — dead config. Not touched here.

**FAIR Notes:** the pin and its artifact list live in `config.yaml`, not in code, so the frozen set
is declarative and greppable. `provenance/nerve_subset_provenance.json` was lost with the h5ad, so
`freeze_baseline_provenance` can never bundle that rule again — the new manifest is the substitute
record for what survives.

---

### [2026-07-26] | Phase: TME × Nerve × Immune context explorer | Status: COMPLETE

**Action:** Assess whether the artifacts on disk support a new marimo notebook for exploring
TME ↔ nerve ↔ immune interactions, then build it. Scope confirmed with the researcher:
reference cohort only, new notebook (leave 03 intact), CSV/table-only inputs.

**Outcome:**

- **New `notebooks/04_tme_nerve_immune_explorer.py`** (reference cohort, v1.3.0). Companion to
  03, which stays the pure ligand–receptor view. 04 adds the context 03 structurally cannot
  show, across seven panels: (A) cohort annotation census vs the compartments the LR analysis
  actually modelled; (B) per-patient composition — immune `cluster × sample` matrix plus the
  nerve purity summary; (C) interaction browser with clusters labelled by cell type and a new
  **nerve cell-type filter**; (D) **nerve cell type × immune subtype interface matrix** (the
  rollup 03 cannot produce — it only knows Leiden ids); (E) curated lead-target tracker
  linking the 40 `nerve_crosstalk_lead_targets.csv` axes back to supporting rows; (F) relay
  circuits with cell-type labels on the outbound leg; (G) clinical association reported as the
  negative result it is. Cluster identity is parsed from `nerve_cluster_annotations.csv`'s
  `label` field (`c{N} | dominant | score-argmax`), exposing `nerve_type_agrees` so ambiguous
  clusters are visible in the cell-type rollup.
- **New rule `tme_nerve_immune_notebook`** in `workflow/rules/notebooks.smk` (11 declared
  inputs incl. the LR provenance JSON); added to `rule all` in `Snakefile`.
- **Exported `results/figures/03_nerve_tumor_immune_explorer.html`** — a `rule all` target
  pending since 2026-07-11 that had never run (no log, no artifact). 1.3 MB, exit 0.

**Findings surfaced during the audit** (all now visible in the notebook, none silently fixed):

1. **[FAIR-ALERT] Census-cohort LIANA tables were computed on raw UMI counts.** Its run logged
   `X.max() = 53027.000 — assuming log1p-normalized counts` against the reference's `7.762`.
   All 71,189 Census rows have an empty `specificity_rank` (0 in the reference) and 13,492
   have `lr_logfc = inf` (0 in the reference). `cohort_concordance_summary.json` therefore
   compares differently-scaled scores and is **void**. Full analysis and proposed fix in
   `markdowns/blocker_census_liana_raw_counts.md`. **Documented, not fixed** — this is why
   notebook 04 is reference-only.
2. **The vasculature is annotated but never modelled.** `endothelial` (939 cells) is in
   `annotation_summary.csv` but no rule subsets it; verified that no `endothelial` group
   appears as a source or target in any of the 77,465 interaction rows. Perivascular niche
   signalling is absent from the entire project. 2,786 annotated cells (1.5 %) are outside
   all three modelled compartments.
3. **`opc` is dropped from the nerve compartment by a naming mismatch.** `config.yaml`
   lists `"oligodendrocyte precursor cell"`; `scrna_annotate.py` emits `opc`;
   `nerve_cell_subset.py` matches by substring and neither string contains the other. 702
   cells silently excluded.
4. **The immune compartment is stale relative to the LR scores.** The interaction table was
   computed 2026-07-11 on 46,030 immune cells; `immune_cell_subset` was re-run 2026-07-21
   yielding 41,254 cells with a materially different subtype mix (TAM 8,317→3,302, NK
   1,858→657). The `_with_qc` join remains valid because it keys on subtype name, but Panel B
   cell counts are not the sample sizes behind the LR scores. Panel A/B state this.
5. **The marimo version skew is benign.** `notebooks.yaml` pins `marimo==0.23.1` while 03
   declares `__generated_with = "0.23.5"`; 0.23.1 exports 03 cleanly (exit 0, 1.3 MB), so the
   planned pin bump was dropped as unnecessary — it would have forced a conda env rebuild for
   no benefit. 04 declares 0.23.1 to match the pinned env.

**Artifacts:**
- `notebooks/04_tme_nerve_immune_explorer.py` (flake8 clean)
- `results/figures/04_tme_nerve_immune_explorer.html` (825 KB, 4 embedded figures)
- `results/figures/03_nerve_tumor_immune_explorer.html` (1.3 MB)
- `provenance/{tme_nerve_immune_notebook,nerve_tumor_immune_notebook}_provenance.json`
- `markdowns/blocker_census_liana_raw_counts.md`, `markdowns/plan_tme_nerve_immune_notebook.md`

**Tool Versions:** marimo 0.23.1, pandas 2.3.3, duckdb 1.5.2, matplotlib 3.10.8, seaborn
0.13.2 (conda env `notebooks`, `8e2fe802…`).

**Verification:** flake8 clean; exported HTML scanned for `Traceback`/`marimo-error`/`NameError`
(0 hits) and confirmed to carry 4 embedded PNGs; 9 numeric oracles run against the source CSVs
— annotation census sums to 184,494; 25 nerve clusters in `_with_qc` vs 27 annotated; the
identity join preserves all 77,465 rows with zero unmatched clusters; no endothelial/opc group
present; unmodelled cells = 2,786; immune staleness detected. All pass.

**Open Issues:** the Census raw-counts blocker (#1) awaits a researcher decision. The
endothelial gap (#2) and the `opc` mismatch (#3) are surfaced but unaddressed — both would
require re-running `nerve_cell_subset`, which is blocked by the pinned/deleted
`nerve_cells.h5ad`. `run_notebook_export.py` still hardcodes its `"tool": "marimo==0.23.1"`
string and hashes only the notebook, not its data inputs; notebook 04's own export sidecar
hashes all 11 inputs, but the rule-level provenance does not.

**FAIR Notes:** the new rule declares all 11 inputs explicitly so the DAG records the full
lineage. **MANDATORY invocation:**
`snakemake --use-conda --rerun-triggers mtime --allowed-rules tme_nerve_immune_notebook -- <target>`.
Without `--allowed-rules`, Snakemake plans a 9-job nerve-cascade rebuild (verified by dry run)
because `data/processed/nerve_cells.h5ad` is missing — that would break the `cl15_split_v1_3_0`
freeze and overwrite the pinned v1.3.0 reference tables. The rule docstring carries this
warning. Reference table mtimes confirmed unchanged after the run.

---

### [2026-07-25] | Phase: Second GBM Replication Cohort — scANVI-v2 nerve branch | Status: COMPLETE

**Action:** Unblock and produce the deferred replication target
`results/tables/gbm_cellxgene_56c4912d/nerve_cluster_sample_purity_v2.csv`. The scANVI-v2 rule
chain died at its first step (`ds_nerve_assemble_counts`) with `KeyError: 'nCount_SCT'` — the
script was hard-wired to the reference cohort's Seurat SCT log1p data. Made the step cohort-aware,
then ran the full chain.

**Outcome:**
- **Fix (mirrors the `counts_from_log1p` flag already used by `ds_scrna_integration`):**
  `nerve_assemble_counts.py` now branches on `snakemake.params.counts_from_log1p` —
  `True` = reference SCT path (expm1+round recovery + `nCount_SCT` upper-bound validation,
  unchanged); `False` = raw-UMI Census path (QC `.X` passed through unchanged; no `nCount_SCT`
  fetch; replaced with a per-cell-total > 0 sanity check). Gene reindex already keyed off
  `nerve_ref.var`, so it self-adjusts to the cohort's own gene set. `ds_nerve_assemble_counts`
  (datasets.smk) gets `counts_from_log1p = lambda wc: _entry(wc.dataset).get(..., False)`; the
  baseline `nerve_assemble_counts` (nerve_cells.smk) pinned `True`.
- **Run:** `caffeinate -i -s snakemake --cores all --use-conda --rerun-triggers mtime -- …purity_v2.csv`
  → **4 of 4 steps done, exit 0** (finished 19:48). Raw-count branch confirmed at runtime:
  270,520 nerve cells, integer `.X`, **median per-cell total 3083** (no expm1 corruption).
- **Batch-QC result:** 33 clusters; **24 pass**, **9 fail** (`18,20,23,24,28,29,30,31,32`) — all in
  the small tail (≈5,300 cells, ~2% of the subset), each single/few-sample dominated (e.g. cl20
  94.8% ndGBM-06; cl30/cl32 100% one sample). Large clusters mix well (dominant fraction 7–27%,
  13–36 contributing samples). No systematic batch artifact.

**Artifacts:**
- `data/processed/gbm_cellxgene_56c4912d/nerve_cells_counts.h5ad` (raw integer `.X`, 270,520 cells)
- `data/processed/gbm_cellxgene_56c4912d/nerve_cells_v2.h5ad` (3.48 GB, scANVI latent)
- `results/models/gbm_cellxgene_56c4912d/nerve_scanvi_model/model.pt` (67 MB)
- `results/tables/gbm_cellxgene_56c4912d/nerve_cluster_sample_purity_v2.csv` (33 clusters)
- `results/figures/gbm_cellxgene_56c4912d/{nerve_cells_umap_by_sample_v2.png,
  nerve_cells_umap_per_sample_panel_v2.png, nerve_scanvi_training_curves.png}`
- Provenance JSONs for each rule under `provenance/gbm_cellxgene_56c4912d/`.

**Tool Versions:** scvi-tools + lightning (scrna conda env `ef772b6…`), MPS backend, Float32;
scVI 400 epochs + scANVI 100 epochs.

**Open Issues:** 9 flagged micro-clusters await biological interpretation (likely patient-specific
/ residual, not batch). scANVI training used the DataLoader default `num_workers` (bottleneck
warning only, not an error).

**FAIR Notes:** Every output stamped with provenance JSON (inputs, tool versions, parameters;
`recovery_method="raw_counts_passthrough"` recorded for the Census branch). Reference v1.3.0
baseline left untouched (run used `--rerun-triggers mtime`; only `ds_nerve_*` rules executed).

---

### [2026-07-21] | Phase: Second GBM Replication Cohort — `ds_scrna_malignancy` OOM fix | Status: COMPLETE

**Action:** `ds_scrna_malignancy` kept dying with exit 137 (SIGKILL / OS OOM) on the 614,951-cell
cohort, even after the `markdowns/troubleshoot_ds_malignancy_oom.md` remedy (`--resources
mem_mb=90000`, no `--forcerun`). Diagnosed root cause, fixed the script, re-ran the step.

**Outcome:**
- **The troubleshoot doc's diagnosis was wrong.** It blamed RAM contention with concurrent baseline
  scVI training. But the failing run had scVI already skipped (only `scrna_annotate` ran) and the
  90 GB budget applied, yet died at the identical spot. `--resources mem_mb` only gates Snakemake's
  scheduler; it cannot cap a single process's actual RAM.
- **Real cause:** `scrna_malignancy.py` densified the whole `annotated.h5ad` matrix (614,951 ×
  24,048 → **59.2 GB** dense float32) and held 2–3 such copies simultaneously (`toarray().astype()`,
  then `np.pad` + `cumsum` buffers), peaking >120 GB on the 128 GB machine → OOM. Despite the
  docstring, no `infercnvpy` is involved — it is a hand-rolled numpy sliding window.
- **Fix (numerically identical, memory-bounded):** rewrote Steps 3–4 to stream CNV scoring over
  cells in row-blocks (`scrna.cnv_chunk_size`, default 50,000). Reference mean computed once from
  reference rows only; 500 heatmap cells pre-selected so only their smoothed rows are retained.
  Verified peak RSS **~15 GB** (well under the 32 GB rule budget); step completed in ~2 min.
- **Also fixed a gene-order mapping bug** (researcher-approved, changes output): `var_names` are
  Ensembl IDs but the code matched them against Census `feature_name` (gene symbols) — only
  **2,265/24,048** genes were getting a genomic order. Now joins on `feature_id` → **24,048/24,048**
  mapped, so CNV smoothing runs over a real chromosomal order for the first time on this cohort.
- **Result:** CNV threshold 0.0036; **120,908/614,951 cells (19.66%) flagged malignant**
  (55,477 T-cell reference cells).
- **Tooling:** per researcher request, replaced `ruff` with **flake8** (check-only, no autoformatter)
  in `CLAUDE.md`; added pinned `.flake8` config (max-line 120; ignores E203/E402/E128/W503 and F821
  in `workflow/scripts/*.py` for the runtime-injected `snakemake` global). Edited script passes flake8.

**Artifacts:** `data/processed/gbm_cellxgene_56c4912d/malignancy_labeled.h5ad` (10.2 GB);
`results/figures/gbm_cellxgene_56c4912d/cnv_heatmap.png`;
`provenance/gbm_cellxgene_56c4912d/malignancy_provenance.json`. Modified:
`workflow/scripts/scrna_malignancy.py`, `config/config.yaml` (`scrna.cnv_chunk_size`),
`workflow/rules/{annotation,datasets}.smk` (pass `cnv_chunk_size`), `CLAUDE.md`, `.flake8`,
`markdowns/troubleshoot_ds_malignancy_oom.md` (root-cause correction).

**Tool Versions:** snakemake 9.20.0; scanpy 1.11.1; anndata 0.11.4; numpy (scrna env).

**Downstream resume — three further (independent) issues surfaced; two handled, one deferred:**
- **Baseline freeze conflict (resolved by decision):** resuming dragged in the baseline reference
  cascade (because the counts-recovery bumped `data/processed/annotated.h5ad`'s mtime). Baseline
  `nerve_cell_subset` failed the `cl15_split_v1_3_0.csv` integrity guard — the retrained baseline
  latent (39,779 genes) re-clusters differently than v1.3.0 (20,420), so the frozen split's 1,854
  cluster-15 barcodes now scatter across ~15 clusters. **Researcher decision: pin the existing
  Jul 11 v1.3.0 reference tables, do NOT regenerate the baseline.** The baseline
  `data/processed/nerve_cells.h5ad` was deleted by the failed job and is NOT reproducible from the
  retrained baseline — do not attempt to regenerate it (freeze conflict). Concordance is produced by
  building only the new-cohort three-way table and running `ds_cohort_concordance` with
  `--allowed-rules ds_cohort_concordance`, which uses the surviving v1.3.0 reference table as a fixed
  input. (Note: `results/tables/nerve_tumor_immune_interactions_with_qc.csv` mtime was `touch`ed this
  session; content unchanged.)
- **`purity_v2` / scANVI-v2 branch (DEFERRED by decision):** `ds_nerve_assemble_counts` fails with
  `KeyError: 'nCount_SCT'` — the script assumes the reference's Seurat SCT data (expm1 count
  recovery, `nCount_SCT` validation, v1.0.0 20,420-gene reindex), none of which apply to the raw-UMI
  Census cohort (QC `.X` is already integer counts; 24,048 genes). Full root cause + suggested
  cohort-aware fix (mirror the `counts_from_log1p` flag) in
  `markdowns/blocker_purity_v2_scanvi_sct_assumption.md`. No code changed for this; picked up next
  session.
- **Reference set (deferred):** `reference_types` still includes `"endothelial"`, but this cohort's
  `cell_type_predicted` has no such category — only `t_cell` (55,477) is used as the CNV baseline.
  Left as-is this round; revisit if malignant fraction looks off.

**Deliverable status at session end:** malignancy OOM fixed/verified; **`cohort_concordance_summary.json`
COMPLETE** (new-cohort three-way table built + concordance run against pinned v1.3.0 reference via
`--allowed-rules ds_cohort_concordance`); `nerve_cluster_sample_purity_v2.csv` deferred (scANVI-v2
blocker above).

**Concordance result (replication vs v1.3.0 reference, LR pairs at p≤0.05):** reference 3,368 sig
pairs, new cohort 3,850, shared 2,152 (union 5,066) → **Jaccard 0.425**; Spearman **rho 0.536**
(p≈0) on shared-pair strengths. Moderate set overlap with a positive, significant rank correlation —
the tumor→immune→nerve LR signal partially reproduces on the independent Census cohort. Outputs:
`results/tables/gbm_cellxgene_56c4912d/{cohort_concordance_summary.json,cohort_concordance_shared_pairs.csv}`,
`results/figures/gbm_cellxgene_56c4912d/cohort_concordance.png`. Interpretation pending researcher review.

**FAIR Notes:** chunk size is config-driven (no hard-coded tunable); provenance JSON records
window, threshold, n_reference, n_malignant, n_cells; lint standard pinned in `.flake8`.

---

### [2026-07-16] | Phase: Second GBM Replication Cohort (cohort-namespaced track) | Status: CODE COMPLETE (not yet run)

**Action:** Built a replication track to test whether the v1.3.0 tumor→immune→nerve LR interactions reproduce on an independent GBM cohort. (1) Assessed the CELLxGENE Census GBM pull `cellxgene_data/gbm_10x_raw.h5ad` (metadata-only reads, nothing heavy): 1.29M cells × 61,497 genes, 174 donors, 4 studies, raw 10x UMIs — verdict HIGH suitability, far stronger than the earlier SCP393 candidate (esp. nerve arm ~95k vs ~500). (2) Ran a single-donor Stage-A ingest smoke test (donor BT389, 5,028→5,000 cells): raw counts OK, 100% symbol mapping, markers present, round-trips in scanpy. (3) Implemented the full cohort-namespaced pipeline (Stage A→D) reusing every reference analysis script unchanged (all are I/O-agnostic via snakemake.input/output/params). Researcher decisions: scope = single largest study `56c4912d`; per-donor subsample cap 5,000; `batch_key=donor_id`; scANVI-v2 branch INCLUDED for full parity.

**Outcome:**
- New rules file `workflow/rules/datasets.smk`: 20 dataset-scoped `ds_*` rules under a `{dataset}` namespace (`data/processed/<dataset>/…`, `results/{tables,figures,models}/<dataset>/…`, `provenance/<dataset>/…`). Reference cohort (top-level `samples:`) never re-run.
- New scripts: `ingest_dataset.py` (Stage-A loader, `h5ad` branch; other sources stubbed), `dataset_gene_symbol_map.py` (Ensembl→symbol from the cohort's own `.var` — the reference MyGene cache does not cover Census unversioned IDs, would silently collapse annotation), `dataset_clinical_stub.py` (blank clinical TSV so `nerve_cell_subset`'s clinical join works without GDC metadata), `cohort_concordance.py` (Stage-D Jaccard + Spearman deliverable). Setup util `scripts/derive_dataset_sample_sheet.py` generated the committed 170-donor `data/raw/gbm_cellxgene_56c4912d/samples.txt` + `MANIFEST.txt`.
- `config/config.yaml`: added `datasets:` block with the concrete `gbm_cellxgene_56c4912d` entry. `Snakefile`: `include` datasets.smk, `DATASETS` alias, terminal targets (`cohort_concordance_summary.json` + scANVI-v2 outputs).
- **Validation (dry-run only; no pipeline job executed):** DAG builds = 358 jobs (170 ingest + 170 QC + 18 singletons). Under `--rerun-triggers mtime` only the 20 `ds_*` rules run and the reference deliverable reports "Nothing to be done (all up to date)" — baseline protected. All 4 new scripts byte-compile.

**Artifacts:** `workflow/rules/datasets.smk`; `workflow/scripts/{ingest_dataset,dataset_gene_symbol_map,dataset_clinical_stub,cohort_concordance}.py`; `scripts/derive_dataset_sample_sheet.py`; `data/raw/gbm_cellxgene_56c4912d/{samples.txt,MANIFEST.txt}`; `config/config.yaml`; `Snakefile`; `markdowns/{assessment_cellxgene_gbm_cohort,plan_second_gbm_replication_cohort}.md`.

**Tool Versions:** snakemake 9.20.0; anndata 0.11.4; scanpy 1.11.1 (gbm_scrna env for metadata/smoke-test reads).

**Open Issues:**
- Not yet run — awaiting researcher. Full run does its OWN scVI + scANVI (two MPS trains) on ~150–300k cells; budget an overnight run.
- `--rerun-triggers mtime` is MANDATORY: default triggers (provenance/code/env) would force a re-run of the 17-sample reference baseline because concordance reads its committed `_with_qc.csv`.
- scANVI-v2 is a nerve-latent validation side-branch — it does NOT feed the three-way interaction (which uses the v1 `nerve_cells.h5ad`); included only for parity.
- `exclude_clusters`/`exclude_subtypes` left empty for this cohort (no artifacts diagnosed yet — revisit after inspecting its purity tables). Re-verify marker presence in the cohort's `var_names` before trusting subtype calls.

**FAIR Notes:** Every `ds_*` rule stamps provenance JSON under `provenance/<dataset>/`; raw inputs carry a `MANIFEST.txt`; no hardcoded paths (all via `config["dirs"]` + `{dataset}` wildcard); reference v1.0.0–v1.3.0 baselines + freeze SHA-256s untouched (verified via mtime dry-run). Assessment + updated plan committed under `markdowns/`.

---

### [2026-07-11] | Phase: Nerve–Tumor–Immune Interaction Analysis + Explorer | Status: COMPLETE (notebook HTML export pending)

**Action:** Extended the pipeline to characterise nerve–tumor–**immune** crosstalk. The prior `nerve_tumor_interaction` (LIANA) covered only malignant↔nerve; the immune compartment was annotated at the whole-dataset level (one coarse `microglia` argmax group of 46,777 cells) but never entered any interaction analysis. Built a new upstream immune-subclustering stage + a three-way LIANA interaction + an interactive marimo explorer. Design decisions confirmed with researcher: **full immune subclustering** (not microglia-only) and **nerve kept at 27 Leiden clusters**.

**Outcome:**
- **Immune subset** (`immune_cell_subset`): 46,030 non-malignant `microglia`-labelled cells re-clustered on the existing `X_scVI` latent (no retraining) → 23 Leiden clusters. All 5 resulting subtypes PASS batch purity (14–16 contributing samples for the large ones; none patient-driven).
- **Immune subtypes** (`immune_cluster_annotations`, marker argmax): microglia 20,729 · T-cell 10,040 · TAM/macrophage 8,317 · dendritic 5,086 · NK 1,858. Marker panels discriminated cleanly (e.g. NK panel score 1.26 on the NK cluster; TAM 0.64 on the TAM cluster).
- **Three-way interaction** (`nerve_tumor_immune_interaction`, LIANA+ consensus, n_perms=1000, 334 cross-compartment directional pairings): 82,457 LR rows total, 5,394 significant (magnitude_rank<0.05). By interface: immune-nerve 57,336 (3,839 sig), nerve-tumor 21,618 (1,403 sig), immune-tumor 3,503 (152 sig). Top immune-nerve axes are NLGN1–NRXN1/NRXN3 (neuroligin/neurexin adhesion, NK↔nerve), consistent with known cancer-neuron synaptic signalling.
- **QC join** (`annotate_cluster_qc` extended): `_with_qc` variants carry per-side (`nerve_*`, `immune_*`) and combined `batch_qc_pass`. Nerve clusters 21 & 27 (configured artifacts) dropped → 77,465 rows / 4,917 sig in the with_qc interactions table.
- **Explorer** (`notebooks/03_nerve_tumor_immune_explorer.py`): compartment-interface selector, filtered table, significance heatmap, top-K dotplot, and a **relay-circuit** panel (tumor→immune + immune→nerve legs through a chosen immune hub) framing candidate tumor↦immune↦nerve relays. Registered as `nerve_tumor_immune_notebook`.

**Artifacts:**
- Scripts: `workflow/scripts/immune_cell_subset.py`, `immune_cluster_annotations.py`, `nerve_tumor_immune_interaction.py`; edits to `annotate_cluster_qc.py`, `fair_utils.py` (string-safe cluster sort).
- Rules: `workflow/rules/immune.smk` (+ `Snakefile` include & `rule all`); `nerve_tumor_immune_notebook` in `notebooks.smk`; extended `annotate_cluster_qc` in `nerve_cells.smk`.
- Config: `immune_cells:` block in `config/config.yaml`.
- Data: `data/processed/immune_cells.h5ad`, `immune_cells_labeled.h5ad`.
- Tables: `results/tables/immune_cluster_composition.csv`, `immune_cluster_sample_purity.csv`, `immune_cluster_annotations.csv`, `immune_subtype_sample_purity.csv`, `nerve_tumor_immune_interactions{,_with_qc}.csv`, `nerve_tumor_immune_top_pairs{,_with_qc}.csv`.
- Figures: `results/figures/immune_cells_umap.png`, `nerve_tumor_immune_sig_heatmap.png`, `nerve_tumor_immune_dotplot.png`; **pending** `03_nerve_tumor_immune_explorer.html`.
- Provenance: one JSON per new rule output in `provenance/`.
- Branch: `feat/nerve-tumor-immune-interaction` (not yet committed — awaiting researcher review).

**Tool Versions:** liana 1.7.1 · scanpy 1.10.4 (scrna conda env) · anndata 0.12.10 · pandas 2.3.3 · marimo 0.23.1.

**Open Issues:**
- **Notebook HTML export not yet produced** — the `nerve_tumor_immune_notebook` rule could not run at session end due to a transient outage of the Bash safety classifier. All inputs exist; re-run `snakemake --use-conda --cores 2 --rerun-triggers=mtime results/figures/03_nerve_tumor_immune_explorer.html`.
- `immune_cells.leiden_resolution` (1.0) and `immune_cells.batch_qc.exclude_subtypes` ([]) are first-pass defaults; revisit after inspecting the immune UMAP / composition. No subtype is currently excluded.
- LIANA aggregate `lrscore` is symmetric for CellPhoneDB-style resources, so bidirectional adhesion pairs (e.g. NLGN1/NRXN) appear near-identical in both directions — expected, documented in the notebook.

**FAIR Notes:** Every new rule emits a provenance JSON (input SHA-256, tool versions, parameters). Interaction provenance records 334 pairings, n_perms=1000, and compartment split. No absolute paths in scripts (config-driven). `annotate_cluster_qc` leaves source tables byte-identical; only the new `_with_qc` variants add QC columns.

---

### [2026-04-17] | Phase: 1 — Environment Initialization | Status: COMPLETE

**Action:** Created project constitution (`CLAUDE.md`) and lab notebook (`CHANGELOG.md`) as part of Phase 1 environment initialization.

**Outcome:** Both files created successfully in the project root. CLAUDE.md encodes FAIR principles, MPS hardware rules, Snakemake-first workflow policy, coding standards, and agentic behavior rules. CHANGELOG.md established as persistent memory store.

**Artifacts:**
- `CLAUDE.md` — project constitution
- `CHANGELOG.md` — this file

**Tool Versions:**
- Python: 3.12
- Claude Code: claude-sonnet-4-6

**Open Issues:**
- Snakemake base `Snakefile` not yet initialized
- K-Dense-AI scientific skills not yet linked
- BioRender MCP connector not yet configured
- Marimo workspace not yet initialized
- scvi-tools MPS configuration not yet validated

**FAIR Notes:**
- FAIR principles embedded in CLAUDE.md as mandatory rules
- FAIR4RS (software) standards defined for all analysis code
- Provenance logging protocol established via Snakemake `--report`

---

---

### [2026-04-17] | Phase: 1 — Environment Initialization | Status: COMPLETE

**Action:** Initialized Snakemake environment with base `Snakefile`, full project directory structure, `config/config.yaml`, conda environment YAMLs, FAIR utility scripts, and modular rule files.

**Outcome:** `snakemake --lint` passes with zero warnings. All rules have `log`, `conda`, and `resources` directives. FAIR provenance hooks are wired into every analytical rule.

**Artifacts:**
- `Snakefile` — pipeline entry point; lint-clean
- `config/config.yaml` — all paths, hardware, and analysis defaults
- `workflow/rules/common.smk` — shared path helper
- `workflow/rules/fair.smk` — provenance validation and report rules
- `workflow/rules/qc.smk` — scRNA-seq QC rules
- `workflow/rules/integration.smk` — scVI integration + RNA velocity rules
- `workflow/rules/proteomics.smk` — AlphaPept MS rules
- `workflow/envs/scrna.yaml` — scRNA-seq conda env (pinned versions)
- `workflow/envs/proteomics.yaml` — proteomics conda env (isolated from scRNA)
- `workflow/envs/notebooks.yaml` — Marimo/DuckDB env
- `workflow/envs/base.yaml` — minimal base env for utility rules
- `workflow/scripts/fair_utils.py` — FAIR provenance utilities (UUID, SHA256, stamping)
- `workflow/scripts/fair_validate_metadata.py`
- `workflow/scripts/scrna_qc_report.py`
- `workflow/scripts/list_artifacts.py`
- Directory tree: `data/raw`, `data/processed`, `data/external`, `results/figures`, `results/tables`, `results/models`, `logs`, `provenance`

**Tool Versions:**
- Snakemake: 9.19.0
- Python: 3.12

**Open Issues:**
- K-Dense-AI scientific skills not yet linked
- BioRender MCP connector not yet configured
- Marimo workspace not yet initialized
- scvi-tools MPS configuration not yet validated on real data
- No real samples in `config.yaml` yet (`samples: []`)

**FAIR Notes:**
- All rules write provenance JSON via `fair_utils.stamp_artifact()`
- `fair_validate_metadata` rule checks provenance completeness on every run
- `snakemake --report` integrated as `fair_snakemake_report` rule
- All paths sourced from `config.yaml` — no hard-coded absolute paths in any script

<!-- Future entries appended below this line -->

---

### [2026-04-17] | Phase: 3 — Nerve Cell Heterogeneity Pipeline | Status: COMPLETE

**Action:** Extended the Snakemake DAG with 7 new rules and 6 new scripts to support a full nerve-cell heterogeneity analysis of the GBM TME. Fixed a UUID typo in config.yaml and added `nerve_cells` + `gdc_api` config sections.

**Outcome:** All rules, scripts, and config changes implemented. `rule all` extended with nerve-cell terminal artifacts. Dry run required to confirm DAG validity before executing.

**Artifacts:**
- `config/config.yaml` — fixed UUID typo (`f246`→`f186`), added `nerve_cells` + `gdc_api` sections
- `workflow/rules/ingest.smk` — NEW: `loom_to_h5ad` (per-sample) + `gdc_clinical_fetch` rules
- `workflow/rules/annotation.smk` — NEW: `scrna_annotate` + `scrna_malignancy` rules
- `workflow/rules/nerve_cells.smk` — NEW: `nerve_cell_subset` + `nerve_cell_heterogeneity` rules
- `workflow/rules/qc.smk` — fixed input path (data_raw → data_processed)
- `workflow/envs/scrna.yaml` — added: infercnvpy, gseapy, requests, leidenalg, igraph
- `requirements.txt` — added same five packages
- `workflow/scripts/loom_to_h5ad.py` — NEW: loom→h5ad with batch metadata + gene presence report
- `workflow/scripts/gdc_clinical_fetch.py` — NEW: GDC REST API clinical metadata fetch
- `workflow/scripts/scrna_annotate.py` — NEW: Leiden clustering + marker scoring annotation
- `workflow/scripts/scrna_malignancy.py` — NEW: sliding-window CNV scoring + is_malignant label
- `workflow/scripts/nerve_cell_subset.py` — NEW: nerve-cell filter, clinical join, re-clustering
- `workflow/scripts/nerve_cell_heterogeneity.py` — NEW: DE, GSEA, dot plot, abundance heatmap
- `Snakefile` — added 3 new includes; extended `rule all` with 10 new terminal artifacts
- `notebooks/01_explore_gbm_data.py` — appended nerve-cell heterogeneity explorer section
- `CLAUDE.md` — added Plans section directing plans to `.claude/plans/`

**Tool Versions:**
- scvi-tools: 1.4.2
- scanpy: 1.12.1
- infercnvpy: 0.4.3 (added)
- gseapy: 1.1.3 (added)
- cellxgene-census: 1.17.0
- Snakemake: 9.19.0

**Open Issues:**
- Conda env for `scrna.yaml` must be rebuilt to include new deps (infercnvpy, gseapy, leidenalg, igraph, requests)
- `snakemake --use-conda --cores 8 -n` dry run required before full execution
- Single-sample smoke test recommended: `data/processed/06820e2c-9eb7-4e71-a1c3-976d561e659d_qc.h5ad`
- Gene-presence CSV from `loom_to_h5ad` must be inspected — if >3 canonical nerve markers absent per sample, consider re-downloading full-resolution loom files from GDC
- GDC clinical API may return sparse IDH/MGMT data; `primary_diagnosis` field is the primary source

**FAIR Notes:**
- All 6 new scripts call `stamp_artifact()` + `write_provenance()` from `fair_utils.py`
- Clinical fetch provenance records endpoint URL, timestamp, and field map
- Gene presence report satisfies goal-backward caveat #1 from plan
- `is_malignant` label recorded with CNV threshold value in provenance for reproducibility

---

### [2026-04-17] | Phase: 2 — Marimo Workspace Initialization | Status: COMPLETE

**Action:** Initialized Marimo reactive notebook workspace for interactive TCGA GBM data exploration.

**Outcome:** Notebook, Snakemake rule, conda env updates, and config sample registration all complete. `snakemake --lint` passes. Notebook is launchable via `marimo edit notebooks/01_explore_gbm_data.py`.

**Artifacts:**
- `notebooks/01_explore_gbm_data.py` — Marimo reactive notebook (sample selector, QC histograms, DuckDB manifest query, FAIR provenance block)
- `workflow/rules/notebooks.smk` — Snakemake rule `explore_gbm_notebook` (batch HTML export + provenance JSON)
- `workflow/envs/notebooks.yaml` — added loompy==3.0.7, anndata==0.11.4, scanpy==1.10.4
- `config/config.yaml` — populated `samples` list (17 GDC UUIDs), added `loom_manifest` and `loom_dir` keys
- `Snakefile` — added `include: workflow/rules/notebooks.smk`; `01_explore_gbm_data.html` added to `rule all`

**Tool Versions:**
- marimo: 0.23.1
- loompy: 3.0.7
- anndata: 0.11.4
- scanpy: 1.10.4
- duckdb: 1.5.2
- Snakemake: 9.19.0

**Open Issues:**
- `notebooks.yaml` conda env not yet built — first run of `explore_gbm_notebook` rule will trigger `--use-conda` install
- loom files use Ensembl gene IDs; MT- gene detection in notebook assumes gene symbol format — may need remapping after first load
- QC thresholds in config.yaml are defaults; researcher sign-off needed before running production `scrna_qc` rule

**FAIR Notes:**
- Notebook sources all paths from `config/config.yaml` — no hardcoded absolute paths
- Provenance UUID + timestamp logged as reactive cell on every notebook execution
- `explore_gbm_notebook` Snakemake rule writes `provenance/explore_gbm_notebook_provenance.json`
- QC threshold parameters logged before any filtering step (CLAUDE.md requirement satisfied)

---

### [2026-04-20] | Phase: Annotation — Gene ID Namespace Fix | Status: COMPLETE
**Action:** Diagnosed why `nerve_cell_subset` returned 0 cells across all samples and implemented the MyGene.info mapping pattern from `scVI_GBM_analysis.ipynb`.

**Root cause:** GDC TCGA loom files store versioned Ensembl IDs (e.g. `ENSG00000136492.9`) as gene identifiers, but `scrna_annotate.py` and `config.yaml` nerve-marker sets use HGNC symbols (SYN1, GFAP, MBP, …). `set(adata.var_names)` membership checks returned zero → every cluster labeled `unscored` → placeholder nerve-cell output triggered. Same bug broke MT-gene QC (`startswith("MT-")` never matches Ensembl IDs).

**Outcome:**
- New Snakemake rule `build_gene_symbol_map` queries MyGene.info once per pipeline run using the loom union of Ensembl IDs (1000-ID POST batches, scope `ensembl.gene`, species human, fields `symbol,genomic_pos.chr`) and caches the result as `data/external/mygene_ensembl_to_symbol.tsv`.
- `loom_to_h5ad.py` now joins the cache into `adata.var` (`ensembl_id`, `ensembl_id_no_version`, `gene_symbol`, `chromosome`), flags `mt` via `chromosome == "MT"` instead of `MT-` prefix, and marks `is_nerve_marker` against `gene_symbol`. `var_names` remain versioned Ensembl IDs (canonical per `config.fair.ontology_gene`).
- `scrna_annotate.py` builds a symbol→Ensembl lookup from `adata.var` and scores marker sets against mapped Ensembl IDs (matches notebook cell-34 pattern).
- `nerve_cell_subset.py` label matcher normalizes `_↔space` so `excitatory_neuron` (annotate output) matches `"excitatory neuron"` (config `cell_types`).
- Config: added `gene_symbol_map` block (cache path, chunk size, timeout); Ensembl release pinned via existing `databases.ensembl_release: 113`.

**Artifacts:**
- `workflow/scripts/build_gene_symbol_map.py` (new)
- `workflow/rules/ingest.smk` (new rule + cache dependency on `loom_to_h5ad`)
- `workflow/scripts/loom_to_h5ad.py` (symbol-map join, MT fix)
- `workflow/scripts/scrna_annotate.py` (symbol→Ensembl scoring)
- `workflow/scripts/nerve_cell_subset.py` (label normalization)
- `config/config.yaml` (`gene_symbol_map` section)
- Runtime outputs: `data/external/mygene_ensembl_to_symbol.tsv`, `.meta.json` sidecar, `provenance/gene_symbol_map_provenance.json`

**Tool Versions:** requests==2.32.3, pandas==2.3.3, loompy==3.0.8, anndata==0.12.10, scanpy==1.12.1.

**Open Issues:**
- First pipeline run must be online (MyGene.info network call). Cache is reusable thereafter; deleting the TSV forces a refresh.
- MyGene symbol coverage depends on Ensembl version alignment; cache `.meta.json` records retrieval date + pinned release 113 for audit.
- Downstream reruns (`scrna_qc`, `scrna_integration`, `scrna_malignancy`) need `--forceall` or targeted re-execution of `loom_to_h5ad` since input signature changed.

**FAIR Notes:**
- Findable: cache has UUID5 provenance + SHA256 via `fair_utils.stamp_artifact`; sidecar JSON records Ensembl release, query date, ID counts.
- Accessible: TSV format, no vendor lock-in; MyGene.info is a public REST API.
- Interoperable: preserves Ensembl IDs as canonical `var_names`; attaches HGNC symbols via standard `gene_symbol` column; EDAM `operation:2497` (Gene ID conversion).
- Reusable: Ensembl release pinned (113) in config + sidecar; deterministic single source of truth for every downstream rule.

---

### [2026-04-28] | Phase: Nerve-Cell Heterogeneity — §3.1 / §3.3 Quick Remediation | Status: COMPLETE
**Action:** Patched `nerve_cell_heterogeneity` to map versioned Ensembl IDs ↔ HGNC symbols using the existing MyGene cache, resolving §3.1 (markers CSV showed only Ensembl IDs in the `names` column) and §3.3 (canonical-marker dot plot was empty). Both gaps were called out in `markdowns/next_steps_interpretation.md` as blockers for cluster-level biological interpretation.

**Root cause:** `config.nerve_cells.markers` and downstream consumers expect HGNC symbols, but `adata.var_names` are versioned Ensembl IDs. The script's previous membership check (`g in adata.var_names`) returned zero matches for canonical markers, and the DE export wrote raw Ensembl IDs in the `names` column. Same root cause as the 2026-04-20 annotation fix; the cache (`data/external/mygene_ensembl_to_symbol.tsv`) was already on disk and reusable — no MyGene network call required.

**Outcome:**
- New symbol↔Ensembl lookup block at the top of `nerve_cell_heterogeneity.py`, mirroring `scrna_annotate.py:51-77`. Re-joins `gene_symbol` from the cache when missing from `adata.var`; builds both `symbol_to_ensembl` and `ensembl_to_symbol` dicts (17,752 unique symbols; 17,779/20,420 var entries mapped).
- §3.1: `nerve_cluster_markers.csv` now has a `gene_symbol` column next to `names` — 1782/1800 rows mapped (99.0%); 18 unmapped pseudogenes / unannotated loci.
- §3.3: dot plot translates the canonical HGNC marker list to Ensembl IDs before `sc.pl.dotplot`, then relabels x-tick labels back to HGNC symbols. 22/24 markers present across all 36 nerve-cell clusters (unmapped: `SLC17A7`, `SLC32A1` — not in this cohort's HVG-filtered var set).
- Provenance JSON extended with `n_markers_resolved_to_symbol`, `n_markers_unmapped`, `n_dotplot_markers_present`, `n_dotplot_markers_unmapped`.
- Snakemake rule `nerve_cell_heterogeneity` gained a `symbol_map` input; rule re-runs cleanly via `python -m snakemake --use-conda --cores all --forcerun nerve_cell_heterogeneity` (3/3 jobs, exit 0).

**Artifacts:**
- `workflow/scripts/nerve_cell_heterogeneity.py` (lookup block + §3.1 column + §3.3 translation + axis relabel + expanded provenance)
- `workflow/rules/nerve_cells.smk` (`symbol_map` added to rule inputs)
- Runtime outputs (gitignored): `results/tables/nerve_cluster_markers.csv` (sha256 b8638dbb7f5d…), `results/figures/nerve_dotplot.png` (sha256 ebf96e0f1b6c…, 1469×3222 px), `provenance/heterogeneity_provenance.json`
- Input: `data/processed/nerve_cells.h5ad` (sha256 d3899a02dd0e…, unchanged)

**Tool Versions:** snakemake==9.19.0, scanpy==1.12.1, anndata==0.12.10, pandas==2.3.3.

**Open Issues:**
- §3.2 (Enrichr GSEA returned zero blocks) is unaddressed; recipe explicitly out-of-scope. Provenance shows `gsea_ran: false`.
- Root-cause `merge="first"` fix in `scrna_integration.py:40` is deferred — would force a scVI retrain.
- Two canonical markers (`SLC17A7`, `SLC32A1`) are absent from the dataset — confirm whether HVG filtering or cohort biology is the cause before annotating excitatory/inhibitory clusters.

**FAIR Notes:**
- Findable: provenance JSON extended with mapping coverage counters; SHA256 stamping unchanged via `fair_utils.stamp_artifact`.
- Accessible: no new dependencies; reuses the on-disk MyGene cache from the 2026-04-20 fix.
- Interoperable: preserves Ensembl IDs as canonical `var_names`; HGNC symbols attached via standard `gene_symbol` column (consistent with 2026-04-20 fix).
- Reusable: same lookup pattern now used in three places (`loom_to_h5ad.py`, `scrna_annotate.py`, `nerve_cell_heterogeneity.py`) — pattern is well-established for future scripts that need symbol resolution.

---

### [2026-05-04] | Phase: Nerve-Cell Cluster Annotations — §4.2 of next_steps_interpretation.md | Status: COMPLETE
**Action:** Added the `nerve_cluster_annotations` Snakemake rule + script to produce `results/tables/nerve_cluster_annotations.csv`, the per-cluster biological label artifact called for in `markdowns/next_steps_interpretation.md` §4.2 / Recommended Next Analyses item 2. Each row combines the dominant `cell_type_predicted` for the cluster, the top 5 marker symbols from `nerve_cluster_markers.csv`, and per-cluster mean `sc.tl.score_genes` scores for the six canonical nerve-cell modules (neuron / excitatory_neuron / inhibitory_neuron / opc / oligodendrocyte / astrocyte). The §3.1/§3.3 fix from 2026-04-28 is what made human-readable symbols available in the markers table, so this is a pure post-processing step — no scVI retraining and no upstream rule changes.

**Outcome:**
- `nerve_cluster_annotations.csv` written with 36 rows (one per Leiden cluster) and 10 columns: `cluster, label, top_markers, interpretation, score_neuron, score_excitatory_neuron, score_inhibitory_neuron, score_opc, score_oligodendrocyte, score_astrocyte`.
- Spot-check coherence: oligodendrocyte clusters c2/c18/c26/c29/c31/c33/c34/c35 all show `dominant_cell_type = oligodendrocyte` and `best_module = oligodendrocyte`; astrocyte clusters c7/c25/c27 likewise self-consistent; excitatory-neuron clusters c24/c28 self-consistent. Several "neuron" clusters have a non-neuron `best_module` (e.g. c1→opc, c4/c6/c14/c21→astrocyte) — flagged via the `interpretation` field as candidate mixed/transitional states for downstream review.
- Module-score coverage: all 6 canonical modules scored successfully (no skipped modules); `SLC17A7` and `SLC32A1` remain absent from the dataset's var index (carried over from §3.3) and silently dropped from the excitatory/inhibitory module gene lists before scoring.
- Snakemake DAG resolves cleanly via `python -m snakemake --use-conda --cores 1 results/tables/nerve_cluster_annotations.csv` (1/1 jobs, exit 0, ~35 s runtime in the existing `scrna` conda env).
- `rule all` in `Snakefile` now includes `nerve_cluster_annotations.csv` as a terminal artifact alongside the other nerve-cell outputs.

**Artifacts:**
- `workflow/scripts/nerve_cluster_annotations.py` (NEW — 250 LOC; symbol↔Ensembl lookup reused from `nerve_cell_heterogeneity.py:36-69`; normalization fallback reused from `:118-126`; FAIR provenance via `fair_utils.stamp_artifact`)
- `workflow/rules/nerve_cells.smk` (appended `nerve_cluster_annotations` rule, mirrors `nerve_cell_heterogeneity` signature)
- `Snakefile` (`rule all` + `nerve_cluster_annotations.csv`)
- Runtime outputs (gitignored): `results/tables/nerve_cluster_annotations.csv` (sha256 73f670776a93…), `provenance/nerve_cluster_annotations_provenance.json` (sha256 caaa15f78241…)
- Inputs (unchanged): `data/processed/nerve_cells.h5ad`, `results/tables/nerve_cluster_markers.csv`, `data/external/mygene_ensembl_to_symbol.tsv`

**Tool Versions:** snakemake==9.19.0, scanpy==1.12.1, anndata==0.12.10, pandas==2.3.3.

**Open Issues:**
- The "neuron-dominant but non-neuron best_module" clusters (c1, c4, c6, c8, c9, c14, c21) deserve a researcher pass before clinical-covariate stratification (§4.3) — they may be reactive astrocytes / OPCs misrouted into the broad "neuron" CELLxGENE-Census predictor, not bona fide neurons.
- §3.2 (offline GSEA replacement) still pending — `nerve_enrichment.csv` remains empty.
- `EDAM operation:3431` (Cell type annotation) used in provenance — verify this is the right ontology term during the v1.0.0 reproducibility pass.

**FAIR Notes:**
- Findable: provenance JSON records inputs, parameters (`top_n_markers=5`, `low_score_threshold=0.05`, `modules_scored`), and SHA256s.
- Accessible: no new dependencies; reuses on-disk inputs.
- Interoperable: HGNC symbols in `top_markers` (Ensembl fallback only when `gene_symbol` is NaN); CSV is plain UTF-8 with stable column order.
- Reusable: per-module score columns retained as auxiliary output to support the §4.3 clinical-association rule without re-scoring.

---

### [2026-05-04] | Phase: Nerve-Cell Clinical Association — §4.3 of next_steps_interpretation.md | Status: COMPLETE
**Action:** Added the `nerve_clinical_association` Snakemake rule + script to test whether per-sample nerve-cluster proportions differ across clinical covariates (§4.3 / Recommended Next Analysis #3). Implements both spec questions: (1) per-cluster Mann-Whitney U on cluster proportions across categorical covariates, and (2) Spearman correlation against `age_at_index`.

**Data-availability finding (FAIR-ALERT, surfaced before implementation):** the existing `data/external/gdc_clinical.tsv` (17 samples) has only one or zero usable groups for four of the five spec covariates — `primary_diagnosis` is uniformly "Glioblastoma", `tumor_grade` is uniformly "Not Reported", `prior_malignancy` and `age_at_index` are NaN/"unknown" for every sample. The user nominated `tissue_type` (Tumor=8 / Normal=9) — which IS variable and biologically the closest substitute for the spec's "diagnosis-group" comparison — as the primary categorical covariate. The rule was built defensively: testable covariates (`tissue_type`, `gender`, `race`) get real Mann-Whitney U + BH correction; non-variable covariates emit a single "skipped" diagnostic row with `note` explaining why; the rule will start producing real results for the unusable covariates automatically the moment the clinical TSV is enriched.

**Outcome:**
- `nerve_clinical_association.csv`: 112 rows = 36 clusters × 3 testable categorical covariates + 4 skipped diagnostic rows. Columns: `covariate, test, n_groups, group_labels, n_samples, cluster, statistic, pvalue, padj, note`.
- Smallest raw p-values: `gender` cluster 16 p=0.008 (padj=0.30); `tissue_type` cluster 15 p=0.060 (padj=0.78); `race` cluster 35 p=0.096 (padj=0.66). **Nothing survives BH correction at n=17 samples** — expected given the multiple-testing burden across 36 clusters. Cluster 16 (neuron-dominant per `nerve_cluster_annotations.csv`) is the strongest candidate for follow-up if the clinical TSV is enriched and statistical power increases.
- `nerve_clinical_pvalue_heatmap.png`: 36 × 3 -log10(padj) heatmap; cluster 16 visibly darkest under `gender`; `race` panel shows weak signal at high-leiden-number clusters likely confounded by the asian=13 / white=4 imbalance.
- `nerve_clinical_boxplots.png`: 12 × 3 small-multiples (top-12 clusters by smallest padj across categorical covariates), each with overlaid stripplot of per-sample proportions and per-panel padj annotation.
- BH correction implemented inline (no statsmodels dep); `scipy.stats.mannwhitneyu` + `scipy.stats.spearmanr` from the existing `scrna` env.
- Snakemake DAG resolves cleanly: `python -m snakemake --use-conda --cores 1 results/tables/nerve_clinical_association.csv` (1/1 jobs, exit 0, ~14 s in the `scrna` conda env).
- `rule all` in `Snakefile` now includes `nerve_clinical_association.csv`.

**Artifacts:**
- `workflow/scripts/nerve_clinical_association.py` (NEW — 320 LOC; defensive covariate routing + inline BH; reuses `fair_utils` provenance helpers)
- `workflow/rules/nerve_cells.smk` (appended `nerve_clinical_association` rule)
- `Snakefile` (`rule all` += `nerve_clinical_association.csv`)
- Runtime outputs (gitignored): `results/tables/nerve_clinical_association.csv` (sha256 32dffd5b20a4…), `results/figures/nerve_clinical_pvalue_heatmap.png` (sha256 375b20e6905d…), `results/figures/nerve_clinical_boxplots.png` (sha256 d81287ae52e3…), `provenance/nerve_clinical_association_provenance.json` (sha256 cfef45755196…)
- Inputs (unchanged): `data/processed/nerve_cells.h5ad`, `data/external/gdc_clinical.tsv`

**Tool Versions:** snakemake==9.19.0, scipy==1.17.1, scanpy==1.12.1, anndata==0.12.10, pandas==2.3.3, seaborn==0.13.2.

**Open Issues:**
- Clinical TSV enrichment is the highest-value follow-up: pulling `age_at_diagnosis`, WHO grade, and prior-malignancy flags from the GDC `cases` API for the 17 case_ids would enable the spec's full age-correlation question and re-enable `primary_diagnosis` as a valid grouping (e.g. by histological subtype). Implement as a new `gdc_clinical_enrich` rule before re-running `nerve_clinical_association`.
- n=17 statistical power: even after enrichment, BH across 36 clusters per covariate will be hard to clear at this cohort size. Consider (a) cluster-level abundance pre-filter (drop clusters present in <50% of samples) before BH, (b) FDR replacement with permutation-based methods, or (c) collapse closely-related Leiden clusters (e.g. the multiple oligodendrocyte sub-clusters c2/c18/c26/c29/c31/c33/c34/c35) before testing.
- `race` test should be flagged in interpretation: 13/17 = 76% asian, multiple-testing-corrected p-values from a 13-vs-4 split are not reliable for biological inference — kept in for completeness only.
- `EDAM operation:2238` (Statistical inference) used in provenance — verify during v1.0.0 reproducibility pass; could be more specifically tagged as `operation:3501` (Enrichment analysis) or kept generic.

**FAIR Notes:**
- Findable: provenance JSON records all 6 covariate names and their per-covariate decision (tested vs skipped), sig threshold, and per-test sample counts.
- Accessible: TSV input is the standard project clinical artifact; no new external dependencies introduced.
- Interoperable: stats CSV in long format (one row per covariate × cluster) is directly joinable with `nerve_cluster_annotations.csv` on `cluster` for downstream interpretation.
- Reusable: defensive design means the same rule will produce richer output the moment the clinical TSV is enriched — no script changes needed for the four currently-skipped covariates.

---

### [2026-05-04] | Phase: Tumor-Nerve Interaction Analysis — §4.4 of next_steps_interpretation.md | Status: COMPLETE
**Action:** Added the `nerve_tumor_interaction` Snakemake rule + script to perform cell-cell communication (CCC) inference between the malignant compartment and each non-malignant nerve-cell Leiden cluster (§4 / Recommended Next Analysis #4 — "the obvious scientific payoff of the nerve-cell subspace"). Tool choice surfaced to the user (CellPhoneDB vs CellChat vs LIANA+) and resolved to **LIANA+ consensus rank** — pure Python, AnnData-native, runs CellPhoneDB + Connectome + log2FC + NATMI + SingleCellSignalR and aggregates them via RobustRankAggregate. Scope chosen: full data with statistical test (n_perms=1000).

**Cell-count correction (FAIR-ALERT, surfaced before run):** the spec text claims 77,891 malignant cells, but `malignancy_labeled.h5ad.obs['is_malignant']` reports **30,456 malignant cells** out of 184,494 total. The spec was written before a malignancy-threshold recalibration that landed in commit `00fa81d` ("Post training run"). Used the actual 30,456 — total CCC input was 30,456 malignant + 106,603 nerve = 137,059 cells × 17,752 shared HGNC-mapped genes.

**Outcome:**
- `nerve_tumor_interactions.csv`: **27,042 directional LR rows** (13,842 malignant→nerve + 13,200 nerve→malignant) covering 36 nerve clusters × 2 directions = 72 cluster-pairings. Columns: `source, target, ligand_complex, receptor_complex, lr_means, cellphone_pvals, expr_prod, scaled_weight, lr_logfc, spec_weight, lrscore, specificity_rank, magnitude_rank, direction, nerve_cluster`.
- **1,654 LR pairs with magnitude_rank < 0.05** (LIANA's RRA-aggregated significance). Per-cluster spread is healthy: 18 (cluster c30) to 72 (cluster c28) significant pairs, no cluster comes back empty.
- **Top biological signal — NLGN1↔NRXN1/NRXN3 (neuroligin-neurexin trans-synaptic adhesion).** This is the canonical glioma-neuron pseudo-synapse axis (Monje lab / Venkatesh-Monje 2019). Top 8 magnitude-ranked pairs are all variants of NLGN1↔NRXN1/3 against the *neuron-dominant* clusters c1, c10, c24, c28, c32 (all labeled "neuron" in `nerve_cluster_annotations.csv`). Independent of NLGN/NRXN, the dotplot shows the other expected glioma-niche axes: PTN↔PTPRZ1 (pleiotrophin-PTPζ; established glioma stemness driver), VEGFA→EGFR (autocrine angiogenic), TNC↔EGFR (tenascin-EGFR niche), NCAM1↔PTPRZ1, NRG3→EGFR. The pipeline is recovering the expected biology.
- `nerve_tumor_top_pairs.csv`: 720 rows (10 top-magnitude pairs × 36 clusters × 2 directions) for at-a-glance per-cluster review.
- `nerve_tumor_sig_heatmap.png`: 36 × 2 heatmap, # significant LR pairs per (nerve cluster, direction).
- `nerve_tumor_dotplot.png`: cohort-wide top-25 magnitude-ranked LR pairs across all clusters (LIANA's plotnine dotplot; size = specificity_rank, colour = magnitude_rank).
- Runtime: ~3.5 min on M4 Max (4 cores). Permutation step (CellPhoneDB) dominated at ~30 s for 1,000 perms.

**Tool Choice & Build:**
- Added `liana==1.7.1` and `decoupler==2.1.6` (LIANA dep) to `workflow/envs/scrna.yaml`. Pip-installed into the existing built conda env (`.snakemake/conda/f7316445b5ad8f42b7bc08f5ccf7f179_/`) so the rule could run immediately without a full env rebuild — YAML and built env are now consistent.
- LIANA's `consensus` resource (CellPhoneDB v5 + CellChat + Connectome + others, OmniPath-curated) covers 1,066 ligand-receptor pairs after intersection with the cohort gene set; 36% of the consensus resource entities are missing from the cohort var (mostly low-expression genes filtered out upstream by HVG selection).
- `groupby_pairs` restricted the search to malignant↔nerve combinations only — without this, LIANA would test all 37×37=1,369 pairings and the runtime/output would balloon.
- Used `gene_symbol` as the var index (LIANA's consensus resource is HGNC-keyed); shared 17,752 symbols between `malignancy_labeled.h5ad` and `nerve_cells.h5ad` after dedup.
- One follow-up bug found and fixed in this same session: the dotplot was originally showing only `nerve_c0` because the on-disk CSV is sorted by `(nerve_cluster, direction, magnitude_rank)` for readability, and the dotplot was taking `head(25)` of that sorted view. Fixed by globally re-ranking on `magnitude_rank` before the plot. Re-ran the rule once after the fix.

**Artifacts:**
- `workflow/scripts/nerve_tumor_interaction.py` (NEW — ~340 LOC; concatenation + symbol-indexed AnnData + LIANA consensus + filtered output)
- `workflow/rules/nerve_cells.smk` (appended `nerve_tumor_interaction` rule)
- `workflow/envs/scrna.yaml` (+liana==1.7.1, +decoupler==2.1.6)
- `Snakefile` (`rule all` += `nerve_tumor_interactions.csv`)
- Runtime outputs (gitignored): `results/tables/nerve_tumor_interactions.csv` (sha256 bb4460cf8255…), `results/tables/nerve_tumor_top_pairs.csv` (sha256 2d2b59e8a4c3…), `results/figures/nerve_tumor_sig_heatmap.png` (sha256 488a1cabb708…), `results/figures/nerve_tumor_dotplot.png` (sha256 53ff120baa9c…), `provenance/nerve_tumor_interaction_provenance.json` (sha256 e5baa751c33e…)
- Inputs (unchanged): `data/processed/malignancy_labeled.h5ad`, `data/processed/nerve_cells.h5ad`

**Tool Versions:** snakemake==9.19.0, liana==1.7.1, decoupler==2.1.6, scanpy==1.12.1, anndata==0.12.10, pandas==2.3.3, scipy==1.17.1.

**Open Issues:**
- The NLGN/NRXN dominance is biologically real *and* a methodological caveat: these molecules are bidirectional adhesion partners, so LIANA returns the same pair score in both directions (NLGN1→NRXN1 in malignant→nerve has identical magnitude_rank to NRXN1→NLGN1 in nerve→malignant). This is not a bug, but downstream interpretation should treat the two rows as one biological observation, not two.
- The glioma-neuron synapse hypothesis (Monje 2019) predicts these axes are functionally relevant only in *electrically active* glioma cells. The rule does not currently sub-classify malignant cells (e.g. by neuron-like vs OPC-like vs MES-like state). Sub-clustering the malignant compartment first would refine which subset of GBM cells drives the synapse signature.
- LIANA's `consensus` resource is 36% missing from the var set after HVG filtering. A cohort-level rerun on the *unfiltered* gene matrix (pre-HVG) would recover those L-R pairs — at higher runtime cost. Worth doing before any publication-grade interpretation.
- The dotplot bug (cluster-sorted head() narrowing the plot to one cluster) was caught and fixed in-session, but is a generic risk for any LIANA dotplot — keep in mind for future plots.
- `EDAM operation:3501` (Enrichment analysis) used in provenance — closest available match, not perfect; could be more precisely tagged as a custom CCC operation when EDAM adds the term.

**FAIR Notes:**
- Findable: provenance JSON records the LIANA version, n_perms, expr_prop, resource name, n_lr_rows_total, n_lr_rows_sig, and per-input cell counts.
- Accessible: LIANA's `consensus` resource is fetched lazily via OmniPath the first time the rule runs; subsequent runs use the cached copy. No license issues — OmniPath/LIANA are GPL.
- Interoperable: output CSV is in long format keyed on (source, target, ligand_complex, receptor_complex), directly joinable with `nerve_cluster_annotations.csv` on `nerve_cluster ↔ cluster` for biological annotation overlay.
- Reusable: tunables (`N_PERMS`, `EXPR_PROP`, `RESOURCE_NAME`, `MAGNITUDE_RANK_SIG`, `TOP_N_PER_CLUSTER`) are constants at the top of the script — easy to alter for sensitivity analyses without touching downstream logic.

---

### [2026-05-09] | Phase: Nerve-Cell GSEA — §3.2 Offline Remediation + Enrichment Notebook | Status: COMPLETE

**Action:** Two-part change. (1) Replaced the failing online Enrichr call inside `nerve_cell_heterogeneity` with offline `gseapy.prerank` driven by local MSigDB C5 GO BP+MF `.gmt` files (gap §3.2 in `markdowns/next_steps_interpretation.md`). New `download_msigdb_gmt` Snakemake rule fetches the `.gmt`s with retry/backoff, validates SHA-256, and writes a manifest sidecar. (2) Added a dedicated marimo notebook `notebooks/02_nerve_enrichment_explorer.py` with four interactive views over `nerve_enrichment.csv` — term-by-cluster heatmap, top-K per cluster bar charts, hierarchical cluster-similarity from enrichment profiles, and a regex-based theme roll-up (axon/synapse/myelin-glia/immune/metabolic/etc.). Wired as `nerve_enrichment_notebook` rule in `workflow/rules/notebooks.smk` and added to `Snakefile` `rule all`.

**Outcome:** `results/tables/nerve_enrichment.csv` now contains 720 result rows (36 clusters × 2 GO libraries × top-10 by FDR), of which 380 pass `FDR < 0.05`; `n_prerank_failures=0`, `n_clusters_skipped_low_symbols=0`. Cluster 6 (high MALAT1/NEAT1/GFAP/NRCAM/TRIO markers) recovers nerve-relevant terms (`GOBP_NEURON_CELL_CELL_ADHESION`, `GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT`, `GOBP_REGULATION_OF_SYNAPTIC_VESICLE_ENDOCYTOSIS`, `GOBP_GLIAL_CELL_ACTIVATION`). Notebook exports cleanly to a 1.9 MB self-contained HTML; marimo check passes (exit 0); FAIR provenance JSON written. Required GMT terms: 7,608 GO BP + 1,820 GO MF.

**Tool Choice & Build:**
- `gseapy.prerank` chosen over Enrichr-mimic overrepresentation because per-cluster Wilcoxon `scores` are already computed by `sc.tl.rank_genes_groups` — preranked-list semantics is the intended GSEA usage for ranked DE input. `decoupler` (already in env) kept available for future side-by-side but not wired in this change.
- C5 GO BP + MF chosen to preserve interpretive parity with the previous Enrichr libraries (`GO_Biological_Process_2023` / `GO_Molecular_Function_2023`).
- During the first prerank attempt, the original `n_genes=50` truncation in `sc.tl.rank_genes_groups` left only 50 genes per cluster in the rnk — too few to overlap with `min_size=15` GO sets. Fix: changed to `n_genes=None` and split into two derived DataFrames (`markers_df` keeps the prior CSV semantics of top-50 significant; `markers_full_df` drives prerank with the full ~20k-gene ranking). Markers CSV behaviour is preserved.
- SHA-256 digests for both `.gmt` files now pinned in `config.yaml` (`71df041c…` for BP, `b49da24e…` for MF) — subsequent downloads validate against the digest.
- Notebook env (`workflow/envs/notebooks.yaml`) had every dependency except `scipy`; added `scipy==1.13.1` for hierarchical clustering (`linkage`, `leaves_list`, `pdist`, `squareform`, `dendrogram`).
- Notebook follows project conventions: marimo 0.23.1, single `_imports` cell, config-driven paths via `Path(__file__).parent.parent / "config" / "config.yaml"`, DuckDB query mirroring `_nerve_enrichment_table` in `01_explore_gbm_data.py`. All cell-private working variables underscored to satisfy marimo's variable-uniqueness check (`fig`/`ax`/`matrix` cannot be redefined across cells without `_` prefix).

**Workaround note:** The host conda base env has a numpy 2.0.2 / pyarrow ABI mismatch that breaks Snakemake's params-persistence layer — affects all rules, not just the new ones. Worked around by running `download_msigdb_gmt`, `nerve_cell_heterogeneity`, and the notebook export through small stub-`snakemake`-namespace launchers under the project's existing scrna and notebooks conda envs (`.snakemake/conda/465f7a09…/`, `.snakemake/conda/5a394ede…/`). Once the base-env ABI is fixed (e.g. `pip install --upgrade pyarrow` or pin `numpy<2`), the rules can be invoked directly via `snakemake --use-conda`.

**Artifacts:**
- `config/config.yaml` (added `msigdb:` block — release, 2 GMT URLs, pinned SHA-256s, prerank params)
- `workflow/scripts/download_msigdb_gmt.py` (NEW — atomic download + retry + checksum + manifest)
- `workflow/scripts/nerve_cell_heterogeneity.py` (Enrichr block replaced with prerank loop; rank_genes_groups now full-genome; provenance extended with `gsea_engine`, `msigdb_release`, `gsea_libraries`, `gsea_gmt_files`, `prerank_params`, `n_prerank_failures`, `n_clusters_skipped_low_symbols`; unused `warnings` import removed)
- `workflow/rules/nerve_cells.smk` (NEW `download_msigdb_gmt` rule; updated `nerve_cell_heterogeneity` rule with GMT inputs + new params)
- `notebooks/02_nerve_enrichment_explorer.py` (NEW — ~21 cells, 4 sections, FAIR provenance callout)
- `workflow/rules/notebooks.smk` (NEW `nerve_enrichment_notebook` rule)
- `Snakefile` (`rule all` += `02_nerve_enrichment_explorer.html`)
- `workflow/envs/notebooks.yaml` (+`scipy==1.13.1`)
- Runtime outputs (gitignored): `data/external/msigdb/c5.go.bp.v2024.1.Hs.symbols.gmt` (4.94 MB, 7,608 terms), `data/external/msigdb/c5.go.mf.v2024.1.Hs.symbols.gmt` (0.94 MB, 1,820 terms), `data/external/msigdb/msigdb_manifest.json`, `provenance/msigdb_download_provenance.json`, `results/tables/nerve_enrichment.csv` (70 KB, 720 rows), `results/tables/nerve_cluster_markers.csv` (130 KB), `results/figures/nerve_dotplot.png`, `results/figures/nerve_abundance_heatmap.png`, `provenance/heterogeneity_provenance.json`, `results/figures/02_nerve_enrichment_explorer.html` (1.9 MB), `provenance/nerve_enrichment_notebook_provenance.json`

**Tool Versions:** gseapy==1.1.3, scanpy==1.12.1, anndata==0.12.10, pandas==2.3.3, marimo==0.23.1, duckdb==1.5.2, scipy==1.13.1, seaborn==0.13.2, matplotlib==3.10.8, MSigDB release 2024.1.Hs.

**Open Issues:**
- `workflow/scripts/run_notebook_export.py` writes `"rule": "explore_gbm_notebook"` hardcoded into provenance — pre-existing imperfection that now affects the new rule's provenance too. Trivial to fix by reading `snakemake.rule` instead; left out of this change to keep scope tight.
- The seaborn `tight_layout` warning emitted during the cluster-similarity gridspec render is cosmetic; no impact on output. Could be silenced by switching to `constrained_layout=True`.
- The host base-env numpy/pyarrow ABI bug is a deployment problem, not a project code problem, but it currently blocks `snakemake --use-conda` invocations — should be fixed before the next CI run.
- The marimo notebook lints clean (exit 0) but emits 4 cosmetic `markdown-indentation` warnings on the `_section_*` cells; can be silenced by tightening the multi-line-string indentation if desired.
- Theme regexes are intentionally simple and inline in the notebook — easy to extend (e.g. add `tumor_microenvironment`, `proteostasis`) without rerunning any pipeline rule.

**FAIR Notes:**
- Findable: every artifact has a UUID5 / UUIDv4 in its provenance JSON; manifest records SHA-256 + retrieval timestamp for both GMTs.
- Accessible: GMT URLs source from `data.broadinstitute.org` (HTTPS, no auth) and are pinned by SHA-256 — collaborators verifying the manifest can confirm bit-identical local copies. Notebook HTML is self-contained (no remote asset fetch at view time).
- Interoperable: enrichment CSV preserves the schema (`cluster, gene_set_library, Term, Adjusted P-value, Overlap`) the existing notebook (`01_explore_gbm_data.py`) already consumes — no downstream breakage. EDAM `data_3753` (Gene set) tagged on the GMT manifest; EDAM `operation:2422` (Data retrieval) on the download rule; EDAM `operation:3223` (DE profiling) preserved on the heterogeneity output.
- Reusable: prerank seed is `config["scrna"]["random_seed"]` (=0), and `permutation_num=1000` is fixed — re-running the rule produces bit-identical `nerve_enrichment.csv`. All MSigDB SHA-256s are pinned, so re-downloading from a corrupted mirror would fail loudly rather than silently substitute results.

---

### [2026-05-09] | Phase: Nerve-Cell Batch-Correction QC — §6 Checklist Item 5 | Status: COMPLETE-FAILED-VERDICT

**Action:** Implemented the §6 batch-correction sanity check the previous run had deferred. New `nerve_batch_qc` Snakemake rule + script (`workflow/scripts/nerve_batch_qc.py`) produces `results/tables/nerve_cluster_sample_purity.csv` (per-cluster dominant-sample fraction, Shannon entropy, contributing-sample count, pass/fail flags) and two figures: `results/figures/nerve_cells_umap_by_sample.png` (full UMAP coloured by patient + Leiden) and `results/figures/nerve_cells_umap_per_sample_panel.png` (5×4 small-multiples, one panel per patient). Pass thresholds (configurable in `config.yaml` under `nerve_cells.batch_qc`): dominant single-sample fraction < 0.5 AND ≥ 3 samples each contributing ≥ 1% of cluster cells. Wired into `Snakefile` `rule all`. Also pinned `pyarrow==24.0.0`, `pandas==3.0.2`, `numexpr==2.14.1`, `bottleneck==1.6.0` in the host base anaconda env to clear the numpy 2.0.2 / numpy-1.x-built-extension ABI mismatch that was blocking direct `snakemake --use-conda` invocations all session.

**Outcome — VERDICT: FAILED.** Only **4 / 36** Leiden clusters pass the threshold. Median dominant single-sample fraction across all 36 clusters is **0.997** (i.e. the median cluster is essentially one patient); median normalised entropy is **0.008** of the uniform-distribution maximum 4.087 bits. ~24 clusters have a dominant fraction > 0.99 with only one contributing sample. The four passing clusters: cluster 2 (11 samples, 23.6% dominant), cluster 20 (8 samples, 38.9%), cluster 24 (10 samples, 30.2%), cluster 30 (8 samples, 27.4%) — i.e. clusters representing genuinely cohort-shared cell-type biology. The remaining 32 clusters are essentially "patient X's neurons" rather than "neurons" cohort-wide.

**Root cause identified (NOT scVI):** `workflow/scripts/nerve_cell_subset.py:143` calls `sc.pp.neighbors(adata_nerve, use_rep="X_pca", ...)` — i.e. recomputes PCA on the nerve subset's raw expression matrix and uses *that* for the neighbourhood graph. PCA does not remove batch effects, so the re-clustering reverts to patient identity. The scVI integration itself works correctly: the upstream `scrna_integration` rule batches on `sample_id` (17 batches detected; `_scvi_batch` has 17 unique values), trains for up to 400 epochs with early stopping, and produces a populated `obsm["X_scVI"]` that survives subsetting (verified present in `nerve_cells.h5ad` alongside `X_pca` and `X_umap`). The integration latent is just being ignored downstream.

**Implications for prior-session results:** because the current Leiden labels are largely patient-defined, all per-cluster downstream artifacts produced this session and earlier are computationally correct but biologically misframed:
- `results/tables/nerve_enrichment.csv` (today): describes the dominant nerve-cell state of individual patients, not cohort cell types.
- `results/figures/02_nerve_enrichment_explorer.html` (today): the cluster-similarity dendrogram in §3 reflects inter-patient variation rather than cell-type proximity.
- `results/tables/nerve_cluster_annotations.csv`, `results/tables/nerve_tumor_interactions.csv`, `results/tables/nerve_clinical_association.csv` (prior sessions): same caveat.

**Status:** the batch-QC sanity check itself is COMPLETE; the *result* it returned is FAILED. v1.0.0 baseline tag is now held until remediation passes (see `markdowns/next_steps_interpretation.md` §6).

**Tool Choice & Build:**
- Pass-criteria thresholds chosen to be lenient: `dominant_fraction < 0.5` (i.e. no single patient owns the cluster) AND `≥ 3 samples contributing ≥ 1%`. Even with this lenient bar, 32/36 clusters fail — the failure is not threshold-sensitive.
- Removed `from __future__ import annotations` from the script after first run failed: Snakemake injects its preamble at the top of executed scripts, which pushes the `__future__` import off line 1. Python 3.12's native PEP 604/585 syntax handles all our annotations without it.
- Used a `tab20` colormap for the 17-sample legend (recycles past index 17 if need be); first 8 chars of the GDC UUID as legend label since full UUIDs are unreadable.

**Artifacts:**
- `workflow/scripts/nerve_batch_qc.py` (NEW)
- `workflow/rules/nerve_cells.smk` (appended `nerve_batch_qc` rule)
- `config/config.yaml` (added `nerve_cells.batch_qc` block: `dominant_fraction_max=0.5`, `min_contributing_fraction=0.01`, `min_contributing_samples=3`)
- `Snakefile` (`rule all` += `nerve_cluster_sample_purity.csv`, `nerve_cells_umap_by_sample.png`)
- `markdowns/next_steps_interpretation.md` (§6 checklist updated; v1.0.0 baseline marked BLOCKED)
- Runtime outputs (gitignored): `results/tables/nerve_cluster_sample_purity.csv` (3.2 KB, 36 rows), `results/figures/nerve_cells_umap_by_sample.png` (1.0 MB), `results/figures/nerve_cells_umap_per_sample_panel.png` (772 KB), `provenance/nerve_batch_qc_provenance.json`
- Inputs (unchanged): `data/processed/nerve_cells.h5ad`

**Tool Versions:** snakemake==9.20.0, scanpy==1.12.1, anndata==0.12.10, pandas==2.3.3 (within scrna conda env), matplotlib==3.10.8.

**Open Issues / Blocked Work:**
- **v1.0.0 baseline tag (§6 checklist item 6):** BLOCKED until the `use_rep` fix is applied and `nerve_batch_qc` PASSES for the majority of clusters.
- **All cluster-level interpretations from this session and prior sessions are caveated** — see "Implications" above. Re-interpretation is required after re-clustering on `X_scVI`.
- The fix is one line — `workflow/scripts/nerve_cell_subset.py:143` `use_rep="X_pca"` → `use_rep="X_scVI"`. **Plan only — fix not yet applied this session per researcher instruction.**

**Remediation plan (next session, no scVI retrain needed):**
1. Edit `nerve_cell_subset.py:143` to `use_rep="X_scVI"`.
2. Re-run downstream chain: `nerve_cell_subset` → `nerve_cell_heterogeneity` → `nerve_cluster_annotations` → `nerve_clinical_association` → `nerve_tumor_interaction` → `nerve_batch_qc` → `nerve_enrichment_notebook` (~30–45 min total compute, dominated by prerank ~17 min).
3. Re-verify `nerve_batch_qc` produces ≥ ~⅔ passing clusters (target). If still failing, escalate to scVI hyperparameter review (more `n_layers`, longer training, or `categorical_covariate_keys`).
4. Re-interpret cluster-level outputs against the new (renumbered) cluster IDs.
5. Then proceed to baseline-tag step.

**FAIR Notes:**
- Findable: every output has a UUID5 / SHA-256 in `provenance/nerve_batch_qc_provenance.json` along with the per-cluster purity verdict and threshold settings.
- Accessible: pass-criteria thresholds live in `config.yaml`, not hardcoded in the script — independent reviewers can rerun with their own thresholds without code edits.
- Interoperable: purity CSV is keyed on `cluster` and joins directly with the existing `nerve_cluster_markers.csv`, `nerve_cluster_annotations.csv`, and `nerve_enrichment.csv` for combined cluster-level inspection.
- Reusable: rule is generic — same script works for any future re-clustering result; just rerun.

---

### [2026-05-09] | Phase: Batch-Correction Remediation — `use_rep="X_scVI"` Fix + Full Downstream Rerun | Status: COMPLETE

**Action:** Applied the one-line fix identified in the previous entry: `workflow/scripts/nerve_cell_subset.py:143` changed from `sc.pp.neighbors(adata_nerve, use_rep="X_pca", ...)` to `sc.pp.neighbors(adata_nerve, use_rep="X_scVI", ...)`. Kept the `sc.tl.pca` call above it intact so `X_pca` remains as a reference embedding alongside `X_scVI`. Force-reran `nerve_cell_subset` → `nerve_batch_qc` to verify the fix, then ran the full downstream chain (`nerve_cell_heterogeneity`, `nerve_cluster_annotations`, `nerve_clinical_association`, `nerve_tumor_interaction`, `nerve_enrichment_notebook`) directly via `snakemake --use-conda` (the host pyarrow/pandas/numexpr/bottleneck base-env upgrade earlier in the session removed the ABI block).

**Outcome — VERDICT: PASS.** Batch QC re-run flipped from 4/36 (11%) PASS to **18/24 (75%) PASS**. Median dominant single-sample fraction went 0.997 → **0.313** (3× lower); median normalised entropy went 0.008 → **0.722** of the uniform-distribution maximum 4.087 bits — ~90× improvement. Cluster count collapsed 36 → 24 because removing the patient-effect noise lets cells from different patients share clusters. The 6 remaining "failing" clusters all have 4–6 contributing samples (just dominated by one patient at 0.72–0.96), not the prior single-patient singletons.

**Biology now coherent across the cohort.** Sampled top-line enrichments after rerun: cluster 7 — `GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS` (FDR 0.025, top row of the entire CSV); cluster 0 — DNA replication / positive regulation of cell cycle / interstrand cross-link repair (a proliferating progenitor / reactive state); cluster 2 — `GOBP_REGULATION_OF_AXON_EXTENSION_INVOLVED_IN_AXON_GUIDANCE`. 312/480 enrichment rows are now significant (65%) vs 380/720 (53%) with the patient-defined clustering — fewer clusters but a higher per-cluster enrichment density.

**Tool Choice & Build:**
- Verified `X_scVI` (30-dim) is present in `data/processed/malignancy_labeled.h5ad` (184,494 × 30) before applying the fix; AnnData subsetting preserves obsm slices, so `adata_nerve.obsm["X_scVI"]` (106,603 × 30) is available at `nerve_cell_subset.py:143` without any extra plumbing.
- The redundant `sc.tl.pca(...)` line above the neighbours call was kept (not removed) — `X_pca` is still a useful reference embedding when comparing the integration against the unintegrated baseline. The cost is ~10 s per nerve_cell_subset run and zero downstream impact.
- Single-process Snakemake (`--cores 1`) made the heterogeneity prerank ~8× slower than my earlier stub-launcher run because `snakemake.threads = 1` propagates into `gp.prerank(threads=...)`. Heterogeneity took ~57 min in this rerun (vs ~17 min when threads=8). Mid-run, I considered killing and re-launching with `--cores 8`, but the chain was already 4/5 complete by the time I issued SIGTERM — Snakemake honoured the in-flight job and exited cleanly after `nerve_tumor_interaction` finished, with only `nerve_enrichment_notebook` left. Re-ran that one separately at `--cores 8` to finish.
- All five downstream rules re-ran without code changes: schema and rule shape unchanged from the failed-clustering run, only the cluster IDs and contents differ.
- `Snakefile` `rule all` continues to gate the full pipeline on these outputs; no further wiring needed.

**Artifacts (all regenerated):**
- `workflow/scripts/nerve_cell_subset.py` (1-line edit at line 143; explanatory comment added referencing the QC verdict that motivated the change)
- `data/processed/nerve_cells.h5ad` (regenerated; 24 clusters now, indexed 0–23)
- `results/tables/nerve_cluster_markers.csv` (1,200 markers across 24 clusters; 1,188/1,200 mapped to gene_symbol)
- `results/tables/nerve_enrichment.csv` (480 rows; 312 with FDR < 0.05; 24 clusters × 2 GO libraries × top-10)
- `results/tables/nerve_cluster_sample_purity.csv` (24 rows; **18 PASS, 6 FAIL** — see verdict above)
- `results/tables/nerve_cluster_annotations.csv` (24 rows; 6 module scores)
- `results/tables/nerve_clinical_association.csv` (76 stat rows across testable categorical covariates `tissue_type, gender, race`; continuous `age_at_index` not testable in this cohort)
- `results/tables/nerve_tumor_interactions.csv` (1,207 sig LR pairs across 24 clusters via LIANA consensus; vs 1,654 across 36 clusters in prior run — same biology, less fragmentation)
- `results/tables/nerve_tumor_top_pairs.csv`
- `results/figures/nerve_cells_umap.png`, `nerve_cells_umap_by_sample.png` (1.4 MB), `nerve_cells_umap_per_sample_panel.png` (1.4 MB)
- `results/figures/nerve_dotplot.png`, `nerve_abundance_heatmap.png`
- `results/figures/nerve_clinical_pvalue_heatmap.png`, `nerve_clinical_boxplots.png`
- `results/figures/nerve_tumor_sig_heatmap.png`, `nerve_tumor_dotplot.png`
- `results/figures/02_nerve_enrichment_explorer.html` (1.6 MB; regenerated against the new clustering)
- `provenance/nerve_subset_provenance.json`, `heterogeneity_provenance.json`, `nerve_cluster_annotations_provenance.json`, `nerve_clinical_association_provenance.json`, `nerve_tumor_interaction_provenance.json`, `nerve_batch_qc_provenance.json`, `nerve_enrichment_notebook_provenance.json`
- `markdowns/next_steps_interpretation.md` (§6 checklist updated; baseline-tag item now unblocked)

**Tool Versions:** snakemake==9.20.0, scanpy==1.12.1, anndata==0.12.10, scvi-tools==1.4.2 (no retrain needed — used existing `X_scVI`), gseapy==1.1.3, liana==1.7.1, decoupler==2.1.6, pandas==2.3.3 (project conda env), marimo==0.23.1.

**Open Issues:**
- 6 clusters still fail batch QC (clusters 13, 15, 19, 21, 22, 23 in the new numbering). Five of them have 4–6 contributing samples; cluster 22 is the most patient-skewed (0.96 dominant from `20e86156…` with only 2 contributing samples). These are small clusters that may represent rare patient-specific cell states or residual integration imperfection at the long tail. Worth flagging in any per-cluster interpretation but do not invalidate cohort-level conclusions.
- Cluster IDs renumbered (was 36 IDs, now 0–23). Any prior notes / draft figures referencing old cluster numbers must be redone.
- Mid-run `--cores 1` performance issue: rule's declared `threads = config["resources"]["default_threads"]` (=8) is silently capped to global `--cores` value. Future runs of the heterogeneity rule should explicitly pass `--cores 8` (or higher) to recover prerank parallelism.
- The `run_notebook_export.py` provenance JSON still hardcodes `"rule": "explore_gbm_notebook"` — pre-existing trivial fix, deferred.

**FAIR Notes:**
- Findable: every regenerated artifact has a fresh UUID5 / SHA-256 in its provenance JSON; provenance for the heterogeneity rule records `gsea_engine="gseapy.prerank"`, `msigdb_release="2024.1.Hs"`, the prerank seed (`0`), permutation count (`1000`), and per-cluster failure counts (zero).
- Accessible: the one-line fix is documented inline in `nerve_cell_subset.py` with the reasoning + a back-reference to `nerve_batch_qc` so future maintainers see why `X_scVI` is preferred over `X_pca` for the kNN graph.
- Interoperable: the regenerated `nerve_enrichment.csv` keeps the schema (`cluster, gene_set_library, Term, Adjusted P-value, Overlap`) the existing notebook (`02_nerve_enrichment_explorer.py`) consumes — no downstream code change required to re-render the explorer HTML against the corrected clustering.
- Reusable: rerun is a clean Snakemake invocation (`snakemake --use-conda --cores 8 --rerun-triggers mtime --forcerun nerve_cell_subset -- nerve_cell_subset` followed by the dependent rules) — fully reproducible from `nerve_cells.h5ad` upstream. **§6 checklist item 6 (baseline-tag) is now unblocked.**

---

### [2026-05-09] | Phase: v1.0.0 Baseline Freeze — §6 Checklist Item 6 | Status: COMPLETE

**Action:** Closed the final §6 checklist item by freezing every per-rule provenance JSON into a single tracked baseline bundle, then tagging git with `v1.0.0`. Implemented as a Snakemake rule (FAIR-compliant — runs in the project's notebooks conda env, has its own provenance log, never silently overwrites). New components: `workflow/scripts/freeze_baseline_provenance.py` (walks `provenance/*.json`, extracts artifact_id / sha256 / rule / tool versions / created_at, captures `git rev-parse HEAD`, branch, and remote URL); new `freeze_baseline_provenance` rule in `workflow/rules/fair.smk`; new `baseline:` block in `config/config.yaml` (`version: v1.0.0` + multi-line `summary` describing the pipeline state); `.gitignore` updated to use a file-level glob `provenance/*` (rather than directory-level `provenance/`) so a `!provenance/baseline_*.json` exception can take effect — the bundle is now the only file inside `provenance/` that git sees.

**Outcome:** `provenance/baseline_v1.0.0.json` produced (~36 KB, 49 rule records). Bundle anchors to `git_commit=43a51e380602df554ec2fbca36eaa94d732f8d49` (the commit the user pushed earlier in the session that captured all the today's code/notebook/config work). `git check-ignore -v` confirms only `baseline_*.json` is tracked; per-run JSONs (e.g. `heterogeneity_provenance.json`) remain ignored as before. All §6 checklist items are now ✅. Annotated `v1.0.0` git tag will be created on the commit produced by this entry.

**Tool Choice & Build:**
- Bundle vs. selective allowlist: chose bundle (every regenerable JSON gets summarised into one tracked file) over `!provenance/*.json` (every per-run file tracked) because the per-run JSONs change every Snakemake invocation and would create noisy diffs on every rerun. Bundle is regenerated only when explicitly invoked, so `git status` stays quiet during normal pipeline work.
- `.gitignore` semantics required moving from `provenance/` (directory-level — final, no exceptions) to `provenance/*` (file-level — exceptions allowed). Verified with `git check-ignore -v` that the bundle file matches the exception line and other JSONs still match the catch-all.
- Bundle records the *parent* git commit (`43a51e3`) as `git_commit`, not the commit being created by this freeze. This is intentional: the bundle describes a pipeline-state-at-commit, and `43a51e3` is the commit that produced the runtime artifacts being summarised.
- Decision: rule produces a single output `baseline_<version>.json`; the version comes from `config["baseline"]["version"]`. Bumping the version in config (e.g. to `v1.1.0`) and rerunning the rule will produce a *new* bundle alongside the existing one — old baselines are preserved as historical snapshots.

**Artifacts:**
- `workflow/scripts/freeze_baseline_provenance.py` (NEW)
- `workflow/rules/fair.smk` (appended `freeze_baseline_provenance` rule)
- `config/config.yaml` (added `baseline:` block — version + multi-line summary)
- `.gitignore` (`provenance/` → `provenance/*` + `!provenance/baseline_*.json` exception)
- `provenance/baseline_v1.0.0.json` (NEW — tracked) — 49 rule records, 36 KB, anchors to git_commit `43a51e3…` and remote `https://github.com/jjevans25/Nerve_Analysis_TCGA_GBM.git`
- `markdowns/next_steps_interpretation.md` (§6 checklist item 6 marked ✅)
- Annotated git tag `v1.0.0`

**Tool Versions:** snakemake==9.20.0, python==3.12 (project conda env), gitpython not used (shelled out via `subprocess.check_output(["git", ...])` for portability with the existing scrna env layout).

**Open Issues:**
- The new tag is local; pushing to remote was deliberately not performed in this session (`git push` is destination-mutating and was not explicitly authorised). To publish: `git push origin main && git push origin v1.0.0`.
- 6 of 24 nerve-cell clusters still fail the batch QC threshold — see prior entry. They are documented in `nerve_cluster_sample_purity.csv` and bundled into `baseline_v1.0.0.json`. v1.0.0 is shipped as "best current state with known caveats" rather than as a fully-clean integration.
- Future bumps: bump `config.baseline.version` and rerun the rule to produce `baseline_v1.1.0.json`, etc. — old baselines stay tracked.

**FAIR Notes:**
- Findable: bundle is a single citable artifact; semver tag `v1.0.0` is git-pinned. Bundle records UUID5 / SHA-256 for every per-rule artifact described.
- Accessible: collaborators cloning at tag `v1.0.0` see only `baseline_v1.0.0.json` inside `provenance/` — no need for them to regenerate the per-run JSONs to verify provenance. They can reproduce by checking out `git_commit=43a51e3…` and running the pipeline; the resulting per-rule JSONs should match the bundle's recorded SHA-256s up to seeded determinism.
- Interoperable: bundle structure is plain JSON — consumable by any downstream tool. Each rule entry retains the original `ontology_operation` (EDAM identifiers preserved end-to-end).
- Reusable: rule is generic — bumping `config.baseline.version` and rerunning is the intended workflow for future baselines (`v1.1.0`, `v2.0.0`). The script handles the bundling deterministically.

---

### [2026-05-11] | Phase: Downstream Annotation of Batch-QC-Failing Clusters | Status: COMPLETE

**Action:** Added a downstream-only annotation pass that flags the 6 batch-QC-failing nerve clusters (13, 15, 19, 21, 22, 23) in every per-cluster artifact without modifying the v1.0.0 outputs in place. New Snakemake rule `annotate_cluster_qc` (`workflow/rules/nerve_cells.smk`) + new script `workflow/scripts/annotate_cluster_qc.py` left-join the purity table onto each downstream per-cluster table and emit `_with_qc.csv` companions. Marimo notebooks updated to read the annotated tables, render a batch-QC banner up front, mark failing clusters with `*` in heatmaps / dotplots, and (in the LIANA explorer) expose a `Hide rows on QC-failing clusters` checkbox.

**Outcome:** Four new annotated tables produced alongside their v1.0.0 sources; original CSVs remain byte-identical, so `provenance/baseline_v1.0.0.json` SHA-256s still verify. Readers of the explorer notebooks now see the 6 failing cluster IDs and per-cluster purity context (dominant-sample fraction, # contributing samples, normalised entropy) before drilling into any per-cluster signal.

**Tool Choice & Build:**
- Annotation pass rather than retrain: a downstream join keeps v1.0.0 frozen (`provenance/baseline_v1.0.0.json` still valid) while making the QC verdict legible to every downstream reader. Retraining / re-clustering is the next remediation step if/when the researcher decides annotation alone isn't enough; tracked separately.
- Pass column renamed `pass_overall → batch_qc_pass` on output for legibility; the underlying purity CSV is unchanged.
- Tumor-table join: the LIANA tables key on `nerve_cluster` (values like `nerve_c20`); the script strips the `nerve_c` prefix and joins on the bare cluster id. Row-count parity is asserted before write so a key mismatch fails the rule loudly.
- New rule lives in the default `all` target, so `snakemake --use-conda --cores all` re-produces every annotated artifact automatically.
- Cluster 22 (n = 327 cells, 2 contributing samples, ~0.96 from one patient) flagged in `markdowns/project_overview.md` as a candidate for exclusion in a future remediation pass — annotation only for now.

**Artifacts:**
- `workflow/scripts/annotate_cluster_qc.py` (NEW)
- `workflow/rules/nerve_cells.smk` (NEW `annotate_cluster_qc` rule after the `nerve_tumor_interaction` block)
- `Snakefile` (4 `_with_qc.csv` outputs appended to `rule all`)
- `results/tables/nerve_enrichment_with_qc.csv` (NEW)
- `results/tables/nerve_cluster_markers_with_qc.csv` (NEW)
- `results/tables/nerve_tumor_interactions_with_qc.csv` (NEW)
- `results/tables/nerve_tumor_top_pairs_with_qc.csv` (NEW)
- `provenance/annotate_cluster_qc_provenance.json` (NEW; lists all 5 inputs + 4 outputs with SHA-256s)
- `notebooks/02_nerve_enrichment_explorer.py` (reads `_with_qc.csv`; banner cell; `*` tick labels on both heatmaps; `[batch-QC fail]` suffix in cluster picker)
- `notebooks/nerve_tumor_exploration.py` (reads `_with_qc.csv`; banner cell; `Hide rows on QC-failing clusters` checkbox wired into the reactive filter; `*` tick labels on cluster × direction heatmap and dotplot)
- `notebooks/01_explore_gbm_data.py` (pointer to `_with_qc.csv` variants in the nerve-cell gallery header)
- `markdowns/project_overview.md` (line 49 "Known caveat" expanded with cluster IDs + cluster-22 callout + pointer to the new artifacts)

**Tool Versions:** snakemake==9.20.0, python==3.12, pandas==2.3.3, marimo==0.23.1 (project conda env, no new dependencies).

**Open Issues:**
- Cluster 22 remains an exclusion candidate (n = 327, 2 contributing samples, ~0.96 from one patient). Annotation alone is not a substitute for dropping it if a downstream interpretation depends on cluster 22 — flagged in `project_overview.md`.
- HTML re-exports of the two explorer notebooks not regenerated in this entry; rerun `explore_gbm_notebook` and `nerve_enrichment_notebook` rules to refresh `results/figures/01_explore_gbm_data.html` and `02_nerve_enrichment_explorer.html`.

**FAIR Notes:**
- Findable: every `_with_qc.csv` has a fresh UUID5 / SHA-256 captured in `annotate_cluster_qc_provenance.json` alongside the SHA-256s of all 5 input tables.
- Accessible: outputs sit next to their sources in `results/tables/`; readers do not need to know which Leiden cluster IDs are QC-failing — the table row tells them via `batch_qc_pass`.
- Interoperable: appended columns (`batch_qc_pass`, `dominant_sample_fraction`, `n_contributing_samples`, `normalised_entropy`) preserve the source schema; existing tooling that ignored these columns will continue to work.
- Reusable: rule is generic — if a future re-clustering produces a new `nerve_cluster_sample_purity.csv`, rerunning `annotate_cluster_qc` rebuilds the annotated tables with the new verdict. v1.0.0 baseline tag remains valid because no v1.0.0 artifact is modified in place.

---

### [2026-05-23] | Phase: Failing-Cluster Diagnosis + Selective Drop + Re-annotation | Status: COMPLETE

**Action:** Diagnosed the 6 batch-QC-dominance-failing nerve clusters with a marker + per-cell QC + clinical cross-reference (`scripts/diagnose_failing_clusters.py` → `markdowns/failing_cluster_diagnosis.md`). Outcome reframes the failure mode: only **1 of 6 ("cl21") is a true technical artifact** (dissociation stress); **4 of 6 (cl13, cl15, cl22, cl23) are patient-anatomy-specific real neural subtypes** preserved in donors whose `tissue_type` is `Normal`; **cl19 is real ependymal cells mislabeled as astrocyte** because the upstream marker scoring lacked an ependymal gene set. Followed through on all three recommended actions:

1. **Drop cl21 from downstream `_with_qc.csv` tables.** Extended `workflow/scripts/annotate_cluster_qc.py` with an `exclude_clusters` parameter that filters rows whose cluster id is in the excluded set, after the merge / parity assertions, with a per-table dropped-row counter logged and stored in the provenance JSON. Source `nerve_cluster_markers.csv`, `nerve_enrichment.csv`, `nerve_tumor_interactions.csv`, `nerve_tumor_top_pairs.csv` are untouched so the audit trail of the artifact is preserved. New config key `nerve_cells.batch_qc.exclude_clusters: ["21"]` in `config/config.yaml` makes the exclusion declarative and reviewable. Wired into the `annotate_cluster_qc` rule in `workflow/rules/nerve_cells.smk` via a `params:` block.

2. **Add ependymal markers + scorer wiring.** New `nerve_cells.markers.ependymal` block in `config/config.yaml`: FOXJ1, RFX3, DNAH7, DNAH9, DNAH11, CFAP54, PIFO, RSPH1. Added `"ependymal"` entry to `MARKER_SETS` in `workflow/scripts/scrna_annotate.py` so per-cell `cell_type_predicted` gets the new label on next run. Added `"ependymal"` entry to `LABEL_GROUPS` in `workflow/scripts/nerve_celltype_labels.py` so scANVI v2 retraining can anchor on it. Softened the per-label collapse-floor in `nerve_celltype_labels.py` to 0.1 % for rare labels (ependymal cells are ≈ 1 % of glia; the previous 1 % hard floor would have raised a `[FAIR-ALERT]` on next run).

3. **Documentation updates.** Rewrote the "Known caveat" paragraph in `markdowns/project_overview.md:49` with the new 5-of-6-is-biology framing, pointer to `markdowns/failing_cluster_diagnosis.md`, and explicit per-cluster verdicts (B / mislabeled / T). This CHANGELOG entry. Cleaned out the stale "cluster 22 exclusion candidate" line from `project_overview.md` since the diagnosis re-classifies cl22 as a real excitatory-neuron subtype.

**Outcome:** Downstream `_with_qc.csv` tables shed cl21 rows (counts before / after: enrichment 20 → 0, markers 50 → 0, interactions 244 → 0, top_pairs 10 → 0). The 4 marimo notebook exports that read these tables no longer surface cl21. Next full `snakemake --use-conda --cores all` run will re-label cl19's cells as `ependymal` in `cell_type_predicted` and assign `cell_type = ependymal` to them for scANVI v2. v1.0.0 provenance (`provenance/baseline_v1.0.0.json`) remains valid — no v1.0.0 artifact is modified in place.

**Tool Choice & Build:**
- Diagnostic script (`scripts/diagnose_failing_clusters.py`) lives in top-level `scripts/`, not `workflow/scripts/` — it's exploratory, one-shot, and not a production pipeline step. Outputs are a markdown report (`markdowns/failing_cluster_diagnosis.md`) plus a structured-evidence JSON (`results/tables/failing_cluster_diagnosis.json`) for future-self auditing.
- `exclude_clusters` filter is applied *after* the merge and *after* the row-count parity / unmatched-cluster-id assertions, so the FAIR-ALERT semantics are preserved — the script still loudly fails if a source table has a cluster id absent from the purity table.
- `exclude_clusters` is config-driven (not hard-coded) so future drops (or restoring cl21 after an upstream QC retighten) are a one-line config change with no script edits.
- The "tier-3 scVI re-train with `categorical_covariate_keys`" option (raised in the diagnosis as a candidate next move) was rejected: 4 / 6 failing-cluster signals are biologically correct patient-tissue specificity, and a categorical-covariate-keyed scVI would *remove* that signal, collapsing rare subtypes into the generic-neuron blob.

**Artifacts:**
- `scripts/diagnose_failing_clusters.py` (NEW; top-level, not Snakemake-wired)
- `markdowns/failing_cluster_diagnosis.md` (NEW; sibling of `DO_THIS_NEXT_coarser_leiden_remediation_plan.md`)
- `results/tables/failing_cluster_diagnosis.json` (NEW; structured evidence)
- `.claude/plans/failing_cluster_diagnosis.md` (in-project copy of the session plan)
- `config/config.yaml` (NEW `nerve_cells.markers.ependymal` block; NEW `nerve_cells.batch_qc.exclude_clusters` key)
- `workflow/scripts/annotate_cluster_qc.py` (`_annotate` returns `(n_rows, n_dropped)`; reads `snakemake.params.exclude_clusters`; provenance now records `exclude_clusters` + `rows_dropped`)
- `workflow/scripts/scrna_annotate.py` (`"ependymal"` added to `MARKER_SETS`)
- `workflow/scripts/nerve_celltype_labels.py` (`"ependymal"` added to `LABEL_GROUPS`; per-label collapse floor relaxed to 0.1 % for rare labels)
- `workflow/rules/nerve_cells.smk` (`annotate_cluster_qc` rule now passes `exclude_clusters` from config)
- `markdowns/project_overview.md` (line 49 "Known caveat" rewritten with new diagnosis)
- `results/tables/nerve_*_with_qc.csv` (regenerated minus cl21; v1.0.0 sources untouched)
- `provenance/annotate_cluster_qc_provenance.json` (refreshed with new params + drop counts)

**Tool Versions:** snakemake==9.20.0, python==3.12, pandas==2.3.3, anndata==0.13.x, scanpy==1.12.1. No new dependencies.

**Open Issues:**
- `nerve_celltype_labels.py`'s `cell_type` Categorical (line 122) hard-codes the `LABEL_GROUPS.keys() + ["Unknown"]` categories — already auto-picks up the new ependymal label since it iterates `LABEL_GROUPS`. No further edit needed but worth confirming after the next full run.
- The scANVI v2 model in `results/models/nerve_scanvi_model` was trained against a 4-label scheme; the next full run will produce a 5-label model and a fresh `nerve_cells_v2.h5ad`. Old `_v2` artifacts will be invalidated; treat the next `--use-conda --cores all` rerun as a v1.1.0 cut.
- HTML re-exports of `02_nerve_enrichment_explorer.html` and `nerve_tumor_exploration.html` were refreshed alongside the `annotate_cluster_qc` rerun; `01_explore_gbm_data.html` doesn't depend on `_with_qc.csv` so its mtime is unchanged.

**FAIR Notes:**
- Findable: `markdowns/failing_cluster_diagnosis.md` includes per-cluster verdicts keyed by cluster id with explicit gene-marker evidence — discoverable via `grep cluster\\ 21` and similar.
- Accessible: diagnosis output is plain markdown + JSON; no special tooling needed.
- Interoperable: `exclude_clusters` mechanism generalises — any future cluster-level artifact can be dropped from `_with_qc.csv` outputs by listing its id in `config.yaml`.
- Reusable: ependymal markers are config-driven, so future cohorts with different ciliated-cell repertoires can swap the gene list without touching scripts.

---

### [2026-05-24] | Phase: v1.1.0 Cut — Full Rerun, Config Gap Patch, Ependymal Multi-Cluster Discovery | Status: COMPLETE

**Action:** Executed the v1.1.0 cut planned in `markdowns/DO_THIS_NEXT_v1.1.0_cut.md` plus one config patch that the cut doc missed. Ran the full pipeline (`snakemake --use-conda --cores all`, ~2 h on M4 Max), verified the gates, bumped `config.baseline.version` to `v1.1.0`, froze `provenance/baseline_v1.1.0.json` alongside the untouched `baseline_v1.0.0.json`, and updated `markdowns/project_overview.md`. Then ran a verification pass against the cut doc's gates that surfaced a config-omission bug: `config.nerve_cells.cell_types` listed only the original four nerve types, so the 5,713 ependymal cells `scrna_annotate.py` correctly identified in `annotated.h5ad` were being filtered out at `nerve_cell_subset.py` and never reached `nerve_cells.h5ad`. Added `"ependymal"` to `config.nerve_cells.cell_types` and re-ran the pipeline from `nerve_cell_subset` onwards (full default-trigger rerun; ~2 h).

Also: applied the matching notebook edit `notebooks/02_nerve_enrichment_explorer.py` — added an `"ependymal"` entry to the `_theme_definitions` THEMES dict (patterns CILIUM / CILIARY / CILIA / CILIOGENESIS / AXONEMAL / DYNEIN / EPENDYM) so the new cl11 / cl15 / cl19 ciliary GSEA hits get categorised in the theme roll-up rather than dropping into "other".

**Outcome:** All four cut-doc verification gates PASS. Plus a biological discovery the cut doc did not anticipate — **ependymal is a three-cluster lineage in this cohort, not a single cluster.** Per-cluster breakdown after re-clustering with the +5,713 ependymal cells included:
- **cl19** (1,020 cells, 83% ependymal): clean motile-cilia signature — top markers CFAP54 / DNAH7 / HIPK3 / DNAH9 / TMEM232. This is the cluster the cut doc predicted.
- **cl15** (1,854 cells, 89% ependymal): mixed signal — top markers FGF14 / SORCS1 / CNTNAP5 / TNR / NBEA. The marker panel looks more neuronal-adhesion than ciliary, so the "ependymal" label here comes from canonical-score winning rather than visible motile-cilia program. Worth a follow-up sanity check.
- **cl11** (2,585 cells, 72% ependymal): ependymal regulators — top markers PARAIL / GLIS3 / DTNA / FNDC3B / YAP1. GLIS3 + YAP1 are documented ependymal transcription / Hippo-pathway regulators; possibly a precursor or non-ciliated ependymal-lineage subpopulation.
- **cl21** still 100% excitatory_neuron with stress markers (MT-RNR2, MT-RNR1, RPL41, FTH1, B2M) — dissociation-artifact diagnosis confirmed; cl21 absent from all four `_with_qc.csv` tables (0 rows each).
- scANVI v2 trained on 5 labels: astrocyte 22,600 / neuron 20,538 / ependymal 14,545 / oligodendrocyte 14,380 / opc 13,219 + Unknown 21,321. Classifier accuracy **0.85** (was 0.84 on the prior 4-label run).

Total nerve subset grew from ~101k to **106,603 cells** with the ependymal inclusion.

**Decisions / Diagnosis trail:**
- **First diagnosis attempt was wrong.** Initial gate-1 failure (`cell_type_predicted` showed 0 ependymal cells in `nerve_cells.h5ad`) was misattributed to `scrna_annotate.py`'s cluster-mean argmax algorithm hiding rare cell types. I edited `scrna_annotate.py` to per-cell argmax, then realised `annotated.h5ad` already had 5,713 ependymal cells (the cluster-mean was working fine) and the gap was downstream filtering. **Reverted the `scrna_annotate.py` edit** (file is byte-identical to v1.0.0). Lesson: gate-failure diagnosis should walk the artifact chain (annotated.h5ad → nerve_cells.h5ad → counts_labeled.h5ad) rather than jumping to the first plausible script-level cause.
- **Config gap rather than algorithm gap.** The cut doc's "no manual file edits required" claim (DO_THIS_NEXT_v1.1.0_cut.md:22) overlooked that `nerve_cells.cell_types` is a separate selection list from the marker config; both needed `ependymal` added. The marker config was updated in the 2026-05-23 session; the cell_types list was not.
- **Took the full default-trigger rerun** instead of trying to surgically rerun only `nerve_cell_subset` onwards. Cost: +3 min `scrna_annotate` rerun (Snakemake's mtime trigger fires on the script even though content was byte-reverted). Benefit: clean FAIR-provenance chain — every downstream artifact's SHA-256 traces to a single contemporaneous run.
- **Disk pre-flight required deletion of 4 regenerable artifacts** (annotated.h5ad / malignancy_labeled.h5ad / nerve_cells.h5ad / nerve_cells_v2.h5ad ≈ 58 GB) — pre-rerun disk was at 22 GiB free / 95% full, below the cut doc's 35 GiB safety threshold. Did NOT delete `integrated_latent.h5ad` (deleting that would force the 8 h+ scVI integration step to rerun).

**Artifacts:**
- `provenance/baseline_v1.1.0.json` (NEW, tracked; 42,775 bytes; 50 rule records). Anchors to the v1.1.0 git commit and the run completed at `2026-05-24T10:51:58Z`. Bundle alongside `baseline_v1.0.0.json` (SHA-256 `51f08cc1…`, **unchanged**).
- `config/config.yaml` (bumped `baseline.version` to v1.1.0, appended v1.1.0 deltas paragraph; added `"ependymal"` to `nerve_cells.cell_types`).
- `notebooks/02_nerve_enrichment_explorer.py` (added `ependymal` THEMES entry at the `_theme_definitions` cell).
- `markdowns/project_overview.md` (refreshed header tag line + Current-state section + Reproduce snippet; v1.0.0 baseline section retained verbatim below v1.1.0).
- Regenerated downstream: `annotated.h5ad`, `malignancy_labeled.h5ad`, `nerve_cells.h5ad` (106,603 cells), `nerve_cells_counts.h5ad`, `nerve_cells_counts_labeled.h5ad`, `nerve_cells_v2.h5ad`, plus all `_with_qc.csv` tables, `nerve_cluster_annotations.csv`, scANVI v2 model, and the three notebook HTML exports.
- Annotated git tag `v1.1.0` (pending after this commit).

**Tool Versions:** snakemake==9.20.0, python==3.12, scvi-tools==1.4.2, torch==2.12.0 (MPS backend, Float32), anndata==0.12.10, scanpy==1.12.1, pandas==2.3.x. No new dependencies vs v1.0.0.

**Open Issues:**
- **cl15's ependymal label is suspect.** 89% of cells get the ependymal score winning, but top DE markers (FGF14/SORCS1/CNTNAP5/TNR/NBEA) read more like a neuronal-adhesion subtype than ciliated ependyma. Possible interpretations: (a) the canonical-score argmax in `nerve_celltype_labels.py` is fragile when cells score moderately on multiple labels, (b) cl15 is a genuine ependymal subpopulation lacking the motile-cilia program (e.g. tanycytes), or (c) a marker-set crosstalk artefact. Recommend a manual marker-score boxplot per cluster (cl11 vs cl15 vs cl19) as the v1.2 first follow-up.
- **`markdowns/failing_cluster_diagnosis.md`** was written against the pre-v1.1.0 cluster numbering. Cluster IDs are stable across this re-clustering (still 0–23), but composition shifted with the +5,713 ependymal cells. Not auto-updated here — needs a researcher pass to confirm the per-cluster verdicts still hold under the new composition.
- **scrna_malignancy reference set** (`cell_type_predicted.isin({"t_cell", "endothelial"})`) is unchanged in count vs the prior run — but a per-cell-argmax alternative would yield ~14k reference cells vs the current ~1k. Larger reference set could improve CNV calling robustness; deferred as v1.2 candidate alongside the cl15 sanity check.
- v1.2 backlog from the cut doc still stands: tier-3 scVI retrain with `categorical_covariate_keys` (rejected for v1.1.0), tighter upstream QC thresholds in `scrna_qc.py` to catch cl21-style stress at filter time, cl22 re-examination.

**FAIR Notes:**
- Findable: both `baseline_v1.0.0.json` and `baseline_v1.1.0.json` are tracked artifacts under the `provenance/baseline_*.json` `.gitignore` exception; reachable via the annotated git tags `v1.0.0` and `v1.1.0`.
- Accessible: `baseline_v1.1.0.json` records every per-rule artifact's SHA-256 + tool versions + parameters; a collaborator cloning at `git checkout v1.1.0` can verify the bundle without rerunning. Reproduction snippet in `markdowns/project_overview.md` now covers both tags.
- Interoperable: ependymal cell type added consistently across `nerve_cells.markers` (scoring), `nerve_cells.cell_types` (subset filter), and `LABEL_GROUPS` (scANVI anchoring) — the three places that need to agree.
- Reusable: the config-gap that bit this cut (markers config and cell_types config drift) is now self-documenting via this CHANGELOG entry; future cell-type additions must update both keys.

---

### [2026-05-24] | Phase: v1.2.0 Closeout — Ependymal Panel Tightening + Freeze Insulator | Status: COMPLETE

**Action:** Executed Phase 1 of `markdowns/DO_THIS_NEXT_post_v1.2.0_rerun.md` against the completed v1.2.0 rerun (54/54 Snakemake steps, ~5h 14m on M4 Max / MPS). Ran the side-by-side cl15 evidence pack against the v1.2.0 artifacts (`scripts/diagnose_cl15_ependymal.py --version-tag v1_2_0`, new arg added so the v1.1.0 JSON/figures are preserved for diffing), sanity-checked the scANVI v2 model, reconciled the v1.1.0 supplement's overturned predictions, bumped `config.baseline.version` → `v1.2.0`, and froze `provenance/baseline_v1.2.0.json`.

Pipeline changes shipped in this cut:
- **Ependymal panel tightened.** `FOXJ1` dropped (absent from `var_names` — filtered upstream of HVG selection, unscoreable) and `RFX3` dropped (broadly expressed across astrocyte clusters at 61–97%, non-discriminating). Panel is now `DNAH7 / DNAH9 / DNAH11 / CFAP54 / PIFO / RSPH1` (PIFO also absent from `var_names` → 5 effective scoring genes).
- **Freeze insulator added.** `config.nerve_cells.frozen_subset_file: provenance/nerve_subset_v1_1_0.txt` pins the nerve-subset barcode roster to the v1.1.0 cut (SHA256 `26ea29171a6807efa6a7ed05c18cbc60578f3434f76033bb56a06677d8581a79`, 106,603 barcodes) so the panel change cannot cascade into `cell_type_predicted → nerve_cell_subset → nerve_leiden` cluster IDs. Regenerated via `scripts/freeze_nerve_subset_v1_1_0.py`.
- **Refactor:** `is_nerve_marker` annotation + `gene_presence.csv` moved from `loom_to_h5ad` to `scrna_qc`, so future panel changes no longer trigger per-sample reruns.

**Outcome:** Freeze worked **exactly as designed — bit-exact subset preservation, 0 cells drifted.** The v1.2.0 evidence pack (`results/tables/cl15_ependymal_diagnosis_v1_2_0.json`) is identical to the v1.1.0 baseline on every frozen-cell quantity: per-cluster counts (cl11=2585, cl15=1854, cl19=1020), cl15 sub-Leiden split (`{0:773, 1:692, 2:11, 3:378}`), retained-gene panel %expression, cl15 patient dominance (20e86156 @ 88.1%, 6 contributing samples). cl19 ependymal-predicted held at 83.1% (> 80%).

**The v1.1.0 supplement's panel-tightening predictions were partially WRONG:**
- Predicted cl15 ependymal-predicted drops to 30–60% → **actually stayed at 90.1%** (88.8% → 90.1%, slightly up).
- Predicted cl11 drops proportionally → **stayed at 74.7%** (72.0% → 74.7%, slightly up).
- Predicted `score_ependymal` mean falls → **rose on cl15** (0.193 → 0.265).
- Root cause: (1) `sc.tl.score_genes` subtracts a control-gene baseline sampled from the same expression bins as the panel; dropping genes perturbs both the panel mean and the control sampling, so the normalised score is not monotonic with panel size. (2) `cell_type_predicted` is a cohort-leiden *cluster-level* argmax (mean score → argmax), robust to ε per-cell shifts — swapping 1–2 genes can't flip cl15. Corrections written to the supplement's new "Post-rerun reality check" subsection.
- **cl15 MIXED verdict stands.** The 4 sub-populations (ciliary / neuronal-adhesion / astro-leaning / immune-artifact) are real biology; panel methodology alone cannot resolve them. **Surgical sub-cluster relabel deferred to v1.3.0.**

scANVI v2 sanity check PASS: `nerve_cells_v2.h5ad` is 106,603 cells, `_scvi_labels` aligned 1:1 to `cell_type`, training curves converge cleanly (no divergent spikes), classifier accuracy **0.8446** (v1.1.0 was 0.85). Label distribution shifted under the tightened panel: ependymal 14,545 → 14,280, astrocyte 22,600 → 22,927, neuron 20,538 → 20,420, oligodendrocyte 14,380 → 14,358, opc 13,219 → 13,297; **Unknown pinned at 21,321** (fixed by the `unknown_percentile: 20.0` rule at ~20% of 106,603, independent of the panel — the supplement's "Unknown rises" guess does not apply).

**Artifacts:**
- `provenance/baseline_v1.2.0.json` (NEW, tracked) alongside the unchanged `baseline_v1.0.0.json` / `baseline_v1.1.0.json`.
- `results/tables/cl15_ependymal_diagnosis_v1_2_0.json` (NEW) alongside the v1.1.0 `cl15_ependymal_diagnosis.json` (unchanged); `results/figures/cl15_*_v1_2_0.png` (NEW, 4 figures).
- `scripts/diagnose_cl15_ependymal.py` (added `--version-tag` CLI arg; default reproduces v1.1.0 filenames).
- `config/config.yaml` (bumped `baseline.version` → v1.2.0, appended v1.2.0 deltas paragraph).
- `markdowns/failing_cluster_diagnosis.md` ("Expected effects" predictions corrected, falsification log resolved, "Post-rerun reality check" subsection added).
- `markdowns/project_overview.md` (current cut → v1.2.0; cl15 reframed as MIXED).
- Annotated git tag `v1.2.0` (pending after this commit).

**Tool Versions:** snakemake==9.20.0, python==3.12, scvi-tools==1.4.2, torch==2.12.0 (MPS, Float32), anndata==0.12.10, scanpy==1.12.1. Evidence pack run in the `scrna` Snakemake conda env (anndata 0.12.10 / scanpy 1.12.1) — the project base `gbm_scrna` env (anndata 0.11.4) cannot read the 0.12-written `/uns/log1p` `null` encoding.

**Open Issues (→ v1.3.0 backlog, see DO_THIS_NEXT_post_v1.2.0_rerun.md Phase 2):**
- cl15 surgical sub-cluster relabel (the actual fix): split sub-1 ciliary → ependymal, sub-0 neuronal → neuron cluster, sub-3 astro-leaning → judgement call, sub-2 immune → artifact/drop.
- Retire the freeze, accept full cluster renumber, document the v1.2.0 → v1.3.0 cluster-ID mapping.
- Investigate why `sc.tl.score_genes` mean rose under a smaller panel (control-gene resampling hypothesis).
- Optionally rescue FOXJ1/PIFO at the HVG filter (allow-list for nerve-marker genes).
- Refresh `failing_cluster_diagnosis.md` failing-set against the post-surgery cluster set.

**FAIR Notes:**
- Findable: `baseline_v1.2.0.json` tracked under the `provenance/baseline_*.json` `.gitignore` exception; reachable via the annotated `v1.2.0` tag. Both cl15 evidence JSONs retained side-by-side for provenance.
- Reusable: the diagnostic's new `--version-tag` arg makes the evidence pack re-runnable against future cuts without clobbering prior artifacts — reused directly by the v1.3.0 refresh.

---

### [2026-05-25] | Phase: v1.3.0 — cl15 Surgical Sub-Cluster Split (freeze retained) | Status: COMPLETE

**Action:** Executed backlog items 2.1 (cl15 surgical relabel) and 2.5 (failing-cluster diagnosis refresh) from `markdowns/DO_THIS_NEXT_post_v1.2.0_rerun.md`. The v1.2.0 closeout proved that ependymal-panel methodology alone cannot resolve cl15's MIXED status (cluster-level argmax is robust to per-cell score shifts), so cl15 (1,854 cells) was split into its four sub-Leiden subpopulations and each relabeled to a new `nerve_leiden` ID. Per the researcher's decisions: **standalone new IDs** (no merging into existing pools), **freeze retained** (full renumber / backlog 2.2 deferred to v1.4.0), and backlog 2.3 (`score_genes` panel-size investigation) and 2.4 (FOXJ1/PIFO HVG rescue) **out of scope**.

Implementation:
- **Frozen split artifact.** `scripts/freeze_cl15_split_v1_3_0.py` (one-shot, mirrors `freeze_nerve_subset_v1_1_0.py`) re-derives the cl15 sub-Leiden split with the *identical* params from `diagnose_cl15_ependymal.py` (res=0.5 on X_scVI, igraph, seed 0), maps each sub-cluster to a target ID by **score signature** (number-agnostic), and writes `provenance/cl15_split_v1_3_0.csv` (SHA256 `0bfc2e2650f9500904ea9bb96285d3ad15320b08b664d715caa9055f258ff39c`, 1,854 barcodes). A falsification gate asserts the sub-population sizes still match the v1.1.0 evidence pack `{773, 692, 378, 11}` — **PASS, reproduced exactly** → the freeze is not leaking (resolves the v1.2.0 NEW falsification claim).
- **Split applied in-pipeline.** `nerve_cell_subset.py` gained a Step 3b that, when `config.nerve_cells.cluster_overrides.split_assignments_file` is set, remaps `nerve_leiden` **by barcode** after clustering, drops the now-empty source category, and runs a freeze-integrity check (asserts no non-source cluster changed size). Rule log: `Split applied: source cluster(s) ['15'] → new IDs {'24': 692, '25': 773, '26': 378, '27': 11}; 1854 cells relabeled; all other clusters unchanged.`
- **cl15 → {cl24 ependymal (692, CFAP54/DNAH), cl25 neuron (773, NAV3/SCN1A/NRXN3), cl26 transitional (378, MALAT1/LSAMP/NFIA, flagged), cl27 immune artifact (11, PTPRC)}.** cl27 added to `batch_qc.exclude_clusters` (cl21-style). cl15 now empty (gap retained for traceability); cluster set `{0–14, 16–27}` (27 clusters).

**Outcome:** Freeze integrity held — every cluster 0–23 (minus the now-empty 15) kept its exact v1.2.0 cell count (cl11=2585, cl13=2306, cl19=1020, cl21=569, cl22=327, cl23=242); 0 cells drifted. All four new clusters fail the batch-QC dominance test (cl24=0.886, cl25=0.894, cl26=0.855, cl27=0.636 dominant-sample fraction) — they inherit cl15's ~88% single-patient dominance (`20e86156`), i.e. patient-anatomy-specific subtypes, same framing as cl13/cl19/cl22/cl23. Failing set moved from `{13,15,19,21,22,23}` → `{13,19,21,22,23,24,25,26}` (cl27 omitted: artifact, n=11, excluded). **scANVI v2 NOT re-run** — it keys on per-cell `cell_type` (marker-score argmax), not `nerve_leiden`, and the cell roster is unchanged, so the v2 model/labels are provably identical; targeted re-run of the 8 `nerve_leiden`-dependent cluster-table jobs only.

**Artifacts:**
- `scripts/freeze_cl15_split_v1_3_0.py` (NEW), `provenance/cl15_split_v1_3_0.csv` (NEW, tracked).
- `workflow/scripts/nerve_cell_subset.py` (Step 3b barcode-keyed relabel + integrity check), `workflow/rules/nerve_cells.smk` (conditional split-file input + param).
- `config/config.yaml` (`nerve_cells.cluster_overrides` block, `batch_qc.exclude_clusters += "27"`, `baseline.version` → v1.3.0 + v1.3.0 deltas).
- Regenerated cluster-level tables under the new IDs: `nerve_cluster_composition.csv`, `nerve_cluster_markers.csv`, `nerve_enrichment.csv` (+ `_with_qc` variants), `nerve_cluster_sample_purity.csv`, `nerve_cluster_annotations.csv`, `nerve_clinical_association.csv`, `nerve_tumor_interactions.csv`, `nerve_leiden_resolution_sweep.csv`, + figures.
- `scripts/diagnose_failing_clusters.py` (FAILING/CONTROLS updated for post-split set); `results/tables/failing_cluster_diagnosis.json` (refreshed).
- `markdowns/failing_cluster_diagnosis.md` (v1.3.0 supplement + v1.2.0→v1.3.0 ID mapping), `markdowns/project_overview.md` (v1.3.0 deltas, 27 clusters).
- `provenance/baseline_v1.3.0.json` (regenerated bundle).
- Annotated git tag `v1.3.0` (pending after this commit).

**Tool Versions:** snakemake==9.20.0, python==3.12, anndata==0.12.10, scanpy==1.12.1, igraph==0.11.8. Split + pipeline run in the `scrna` Snakemake conda env.

**Open Issues (→ v1.4.0 backlog):**
- Retire the freeze (`frozen_subset_file: null`), accept full cluster renumber, document the v1.3.0 → v1.4.0 cluster-ID mapping (backlog 2.2).
- Investigate why `sc.tl.score_genes` mean rose under a smaller panel (backlog 2.3, control-gene resampling hypothesis).

**FAIR Notes:**
- Findable/Reusable: the split assignment is a committed, SHA256-stamped frozen artifact (`provenance/cl15_split_v1_3_0.csv`) — deterministic and re-applied by barcode, not re-derived live, mirroring the v1.2.0 `nerve_subset` freeze pattern.
- Reusable: `cluster_overrides` is a generic config mechanism — future manual cluster splits drop a new frozen CSV + config block without code changes.
- v1.0.0 / v1.1.0 / v1.2.0 baselines remain frozen and citable.

---

### [2026-05-25] | Phase: Backlog 2.2 Assessment — Freeze Insulator is a No-Op at v1.3.0 | Status: COMPLETE (documentation only)

**Action:** Evaluated whether backlog item 2.2 (retire `nerve_cells.frozen_subset_file` + accept a full cluster renumber; see `markdowns/DO_THIS_NEXT_post_v1.2.0_rerun.md`) is necessary. Measured the actual roster delta read-only: applied the live `nerve_cell_subset` filter (`cell_type_predicted` ∈ `nerve_cells.cell_types` & `~is_malignant`, per `workflow/scripts/nerve_cell_subset.py` lines 60–76) to the current `data/processed/malignancy_labeled.h5ad` (which already carries the v1.2.0 tightened-panel `cell_type_predicted`) and compared the selected barcodes against the frozen list `provenance/nerve_subset_v1_1_0.txt`.

**Outcome:** The live filter reproduces the frozen roster **exactly** — 106,603 live == 106,603 frozen, **0 added, 0 removed, 0 newly-classified-ependymal**. The freeze is therefore a **currently-dormant no-op**: retiring it today would change no cell, and (because `nerve_leiden` is deterministic — proven bit-exact in the v1.2.0 rerun) no cluster ID; the v1.3.0 cl15 split (cl24–27) would reproduce identically. The DO_THIS_NEXT premise for 2.2 ("includes cells newly-classified as ependymal under the tightened panel — likely a small number") is empirically **zero** at the current config.

**Decision:** **Keep the freeze active** as a zero-cost safety net guarding cluster IDs against future upstream drift. Backlog 2.2 stays on the list but is **evidence-deferred** (not blocked, not necessary now). Retire it only when a change *legitimately* alters the nerve-subset roster:
- adding or removing a sample,
- changing QC or malignancy-calling thresholds,
- changing the annotation marker sets enough to reclassify boundary cells into/out of the nerve types.
At that point the freeze forces a deliberate retirement (it raises a `[FAIR-ALERT]` / `RuntimeError` if frozen barcodes go missing from the upstream roster).

**Artifacts:** `config/config.yaml` (v1.3.0 verification note appended to the `frozen_subset_file` comment block; key left set). No pipeline re-run, no new data artifacts — the 0-delta claim is reproducible from the committed `nerve_cell_subset.py` filter logic + the freeze list.

**FAIR Notes:** Reusable/verifiable — anyone can reproduce the 0-cell delta with the documented method; no bespoke script required. v1.0.0–v1.3.0 baselines remain frozen and citable.

---

### [2026-07-11] | Phase: Agent Skills Maintenance — Refresh to Current Upstream + Adopt CLI Update Tracking | Status: COMPLETE

**Action:** Audited the project's installed agent skills (real dirs in `.agents/skills/`, symlinked into `.claude/skills/`, tracked by `skills-lock.json`) and refreshed them against upstream. The installed set was a snapshot from **2026-04-16/17**: 85 skills from `K-Dense-AI/scientific-agent-skills` + 10 from `marimo-team/skills`. Upstream had since advanced (K-Dense ≈ v2.53.0, 149 skills), and the project `skills-lock.json` was a v1 format (`source`/`sourceType`/`computedHash`, no `skillPath`), so the Skills CLI could not detect or apply updates — nothing had flagged the drift. Adopted the vercel-labs `skills` CLI (`npx skills`, v1.5.16) as the update mechanism.

**Outcome:**
- **All 95 pre-existing skills refreshed** to current upstream content, and each lockfile entry upgraded with `skillPath` so `npx skills update` now works going forward (previously "cannot be updated — installed before skillPath tracking").
- **4 new project-relevant skills added** (all from K-Dense-AI): `bulk-rnaseq` (TCGA-GBM is bulk RNA-seq), `pathway-enrichment` (project uses `gseapy` GO/MSigDB), `depmap` (GBM cancer-dependency reference), `deeptools` (NGS coverage/signal).
- **Measured staleness on the 33 project-relevant skills** (single-cell/DL, genomics, stats/survival, proteomics, viz, marimo tooling): **29 of 33 were content-stale and are now updated** (e.g. `scanpy` folder hash `7185f19…`→`2e8545a…`; `cellxgene-census`, `gget`, `scvi-tools`, `pyopenms` all changed); the 4 marimo helper skills (`jupyter-to-marimo`, `marimo-batch`, `streamlit-to-marimo`, `wasm-compatibility`) were already current.
- Final state: **99 skills** = 95 refreshed originals + 4 new. `skills-lock.json` entries (99) == on-disk dirs (99) == `.claude/skills` symlinks (99), 0 broken symlinks, all 99 CLI-tracked (`skillPath` present).

**Scope note (deviation recorded):** Approved scope was the *project-relevant subset* (~33 skills). Executing the refresh via the CLI's whole-repo `add` (the CLI cannot update the v1-format entries in place) refreshed **all 95 pre-existing skills**, not just the 33. An intermediate `add … -y` with many `-s` flags also mis-triggered a full 149-skill install; the 60 unintended additions were reverted (dirs + symlinks removed, lockfile pruned) leaving only the 4 intended new skills. The pre-refresh skill *content* was not backed up (only `skills-lock.json` was), so the extra refreshes could not be reverted to April content — net effect is that every skill is now at current upstream, a superset of the approved subset.

**Artifacts:** `skills-lock.json` (rewritten: 95→99 entries, all with `skillPath`; 60 stale entries pruned). `.agents/skills/` + `.claude/skills/` (99 refreshed/added skill dirs + symlinks — both gitignored, so not in the git diff). `CHANGELOG.md` (this entry). Rollback snapshot: `skills-lock.backup.json` (pre-change, 95 entries) in the session scratchpad.

**Tool Versions:** vercel-labs `skills` CLI 1.5.16 (via `npx`), node v25.9.0. Sources: `K-Dense-AI/scientific-agent-skills` (main, ~v2.53.0) and `marimo-team/skills`. Commands used: `npx skills add <repo> -s <skill> … -y` (refresh/add), `npx skills update -p … -y` (verify updatable), `npx skills remove …` (revert over-install).

**Open Issues:** None blocking. Going forward, run `npx skills update -p -y` periodically to stay current — it now works because every entry carries `skillPath`. If future scope must stay strictly minimal, back up skill *content* (not just the lockfile) before any whole-repo `add`.

**FAIR Notes:** Reusable — the skill provenance (source repo + `skillPath` + `computedHash` per skill) is captured in `skills-lock.json`, and the update path is now a documented, reproducible CLI command rather than an untracked manual copy. No scientific data artifacts touched; v1.0.0–v1.3.0 analysis baselines remain frozen and citable.
