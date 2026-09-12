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
Phase:          Pre-publication (blog). Analysis COMPLETE on both Census arms; housekeeping + write-up.
Last Updated:   2026-09-11
Repo:           github.com/jjevans25/GBM_Nerve_Tumor_Immune_Single_Cell_Analysis (renamed 2026-07-28; local dir intentionally still Nerve_Analysis_TCGA_GBM — see CLAUDE.md)
Active Agent:   lead-researcher
Current Task:   NONE. Both Census arms are COMPLETE and pass 10/10 enforcing compartment gates (nerve 95.39%/94.76% neural, tumor 92.96%/93.47% malignant, immune purity 99.63%, malignancy F1 0.911). Cross-arm lead-axes shortlist delivered: 183 rows / 144 unique axes, tiers 34/18/57/57/17, now carrying per-row `selection_arm` + `qc_filter_applied`. Four Census notebooks render on both arms. Project is in **pre-publication write-up**: engineering post drafted at `markdowns/GBM_TME_Crosstalk_Analysis.md`, science post at `markdowns/blog_01_science_nerve_immune_crosstalk.md`. Dashboard hosting is the researcher's, out of scope here.
Blocked On:     Nothing. **Defect 1 CLOSED 2026-09-11, fully.** `--use-conda` is enforced (`scripts/run_snakemake.sh`, smoke test 14/14) and `fair_validate_metadata` parses all 797 provenance records against the declared pins, failing the build on any unreasoned conflict. The former blanket [FAIR-ALERT] is now ONE named, reasoned exception: the full arm's scvi_integration ran under torch 2.11.0 vs the pinned 2.12.0 (2026-07-30). **Researcher decided 2026-09-11 to RETAIN it — no re-derivation** — and the pin deliberately stays 2.12.0, the env all 15 post-fix artifacts were built under. Capped arm is pin-clean, so cross-arm replication does not rest on it.
Next Action:    (1) ~~Defect 1 / Option B~~ — **DONE 2026-09-11** (closed; latent retained by decision). (2) ~~Make the repo linkable~~ — **DONE**: 92-file / 49 MB results slice committed and pushed on `docs/blog-prep-2026-09`. Remaining: `results/` is 100% gitignored (0 tracked files); a reader following either post can obtain no table. Zenodo deposit (already earmarked, and the only off-machine copy of the 26 MB v1.3.0 archive + the 1.4 GB unreproducible `nerve_cells_counts.h5ad`) or a committed slice. (3) Push `docs/blog-prep-2026-09` and `chore/disk-reclamation-2026-08` — both unpushed, on the single disk flagged as un-backed-up. (4) Review purity_v2 failing clusters (16/43 full, 9/33 capped). (5) The five recommended analyses in `markdowns/post_compartment_fix_next_steps.md` §3.1–3.5, none executed. Note `marimo` on PATH has a broken matplotlib; run notebooks via the snakemake notebooks env. **Invocation: always `scripts/run_snakemake.sh`, targets first, and pin `--allowed-rules` — `--forcerun` alone pulled the network-only `refresh_drug_annotation` into the DAG on 2026-09-11.**
```

Prior status (v1.3.0 baseline, retained): cl15 surgical sub-cluster split COMPLETE; freeze insulator retained; cluster set {0-14, 16-27}.
```

---

## Session Log

---

### [2026-09-11] | Phase: Researcher decision — retain the torch 2.11.0 latent | Status: COMPLETE

**Action:** Researcher decision on the last open item from the Defect 1 closure: whether to
re-derive `data/processed/gbm_cellxgene_56c4912d_full/integrated_latent.h5ad`, the one artifact in
797 provenance records that ran under an unpinned dependency (torch 2.11.0 vs pinned 2.12.0).

**Outcome: RETAIN. No re-derivation.** Recorded as a decision in three places — the
`ACCEPTED_VERSION_EXCEPTIONS` entry in `fair_validate_metadata.py`, §0 of
`markdowns/task_conda_env_enforcement.md`, and the engineering blog draft.

**The pin stays at `torch==2.12.0`.** This was the non-obvious half of the decision and is
deliberate: 2.12.0 is the environment all fifteen post-fix artifacts were built under, so
downgrading the pin to match this one historical file would make those fifteen newly
non-conformant — inverting the defect rather than closing it. The accepted exception covers **one
historical artifact, not a go-forward environment choice.** A future re-run of
`ds_scrna_integration` will record 2.12.0 and stop matching the exception on its own; there is no
hand-maintained suppression to un-stick.

**What stays unknown, and should be stated rather than glossed in any write-up:** whether scVI
under torch 2.11.0 produced a materially different embedding than 2.12.0 would have. The only way
to establish that is the retrain that was just declined. Mitigating evidence: the capped arm's own
integration ran under 2.12.0 and is pin-clean, and this project's replication claim is cross-arm
agreement — so the headline result does not rest on the drifted file alone.

**Cost avoided:** a multi-hour scVI retrain on 1,006,344 cells which would renumber every Leiden
cluster and thereby invalidate every per-cluster table keyed on `nerve_c{N}` — annotation, purity,
compartment audit, both interaction tables and the lead-axes shortlist. (Correction to an earlier
note in this session: the "~11 h" figure quoted for this job was carried over from a `rule all`
comment describing the **17-sample v1.3.0 reference** retrain, not the 170-donor full arm. The full
arm's scVI time has never been measured in isolation; the 2026-08-02 full-cohort run was ~31 h
total with scVI named as the long pole.)

**Artifacts:** `workflow/scripts/fair_validate_metadata.py`,
`markdowns/task_conda_env_enforcement.md` (untracked by design),
`markdowns/GBM_TME_Crosstalk_Analysis.md`, `results/fair_validation_report.json`.

**Outcome verification:** `fair_validate_metadata` re-run, exit 0, `pass: true`, 797 records, 0
unaccepted conflicts, 3 accepted exceptions. flake8 clean.

**Open Issues:** None on Defect 1 — it is fully closed. Remaining project items are unchanged:
Zenodo deposit (still the only off-machine plan for 740 MB of models and the 1.4 GB unreproducible
`nerve_cells_counts.h5ad`); `rule all` schedules 365 jobs rather than a no-op since the Tier-3
reclamation; 4 malformed `tool_versions` strings remain as non-blocking warnings.

**FAIR Notes:** The project's environment story is now precisely stateable: one artifact of 797 is
reproducible from a *recorded* version rather than a *declared* one, the gate names it, and
retaining it was an on-the-record decision rather than an unexamined fact.

---

### [2026-09-11] | Phase: Defect 1 closure + results publication | Status: COMPLETE

**Action:** Researcher directed: fix Defect 1 (conda-env shadowing), and publish `results/` if
the size is manageable.

**Outcome — Defect 1 is CLOSED.**

The *remediation* shipped 2026-08-05 in `a375811` (`scripts/run_snakemake.sh`); the task doc's
`OPEN — not started` status line was simply never updated, and `conda_env_smoke_test` has passed
14/14 against its own acceptance criteria ever since. Neither candidate approach in the doc was
needed — Snakemake was not reinstalled, the venv was not rebuilt.

What was genuinely missing was any *check that the pins stay enforced*, and any measurement of how
far the original defect reached. `fair_validate_metadata` globbed `provenance/*.json`, counted the
files and wrote the count. It never opened one, and never saw the ~650 records under
`provenance/<dataset>/` — so it reported on the retired v1.3.0 reference while the live Census arms
went unaudited. It now recurses (**70 -> 797 records**), compares every recorded `tool_versions`
entry against every `pkg==version` pin in `workflow/envs/*.yaml`, and **fails the build with
`[FAIR-ALERT]`** on any conflict not in a reasoned `ACCEPTED_VERSION_EXCEPTIONS` entry.

**Measured blast radius.** The task doc claimed "every `script:`-rule artifact this repo has ever
produced". Across 797 records there is **exactly one**: the full arm's `scvi_integration`,
2026-07-30, **torch 2.11.0** against the pinned 2.12.0. The capped arm's own integration recorded
**2.12.0** and is clean — so, since the replication evidence here is cross-arm agreement, the
result does not rest on the drifted artifact. Every load-bearing rule downstream (annotate,
malignancy, both subsets, both interaction rules, compartment audit, scANVI, concordance) was
rebuilt **2026-08-06** with versions matching the pins. Re-deriving that one file is an ~11 h scVI
retrain that renumbers Leiden clusters, so per the doc's own Risk section it stays a **researcher
decision** and is recorded as an explicit accepted exception, not a silent pass.

**The new gate caught two lying provenance records on its first run**, both fixed at source:
`nerve_clinical_association.py` recorded `"scanpy": ad.__version__` — anndata's version filed under
scanpy's key, with a comment calling it a "stand-in for the env" — and `nerve_assemble_counts.py`
recorded the literal string `"scipy.sparse@n/a"`, a version field naming no version. Their
reference-cohort artifacts cannot be regenerated (reference pinned, `nerve_cells.h5ad` deleted), so
those are allowlisted as known-bad historical records with the reason stated.

**Outcome — `results/` is published and pushed.** 92 files / 49 MB (~22 MB packed), force-added
past the wholesale `results/` ignore; repo 13 MB -> 37 MB. Both arms' analysis tables, figures and
the 8 notebook HTMLs, the lead-axes shortlist, the audit gates/snapshots, and the FAIR +
smoke-test reports. `results/README.md` documents the slice, the four ways these tables are easy to
misread, and the invocation rules.

**Excluded by researcher decision:** the ~40 root-level v1.3.0 reference tables. That cohort's
nerve compartment was 11% neural / 59% malignant; its nerve claims are void and uncorrectable, and
publishing them beside the corrected tables with nothing in the filename to warn a reader would
hand out void numbers. Also excluded: `results/models/` (740 MB) and the pre-fix snapshot (153 MB),
both Zenodo material; 680 per-sample QC/gene-presence CSVs (51 MB); the pre-QC-join duplicates
(verified strictly redundant — same 44,140 rows, 19 cols -> 28); and four 170-panel UMAP contact
sheets (26 MB).

**Pushed:** `docs/blog-prep-2026-09` -> origin, carrying 10 commits, including the four
disk-reclamation commits that had been local-only on a single un-backed-up disk.

**Artifacts:** `workflow/scripts/fair_validate_metadata.py` (rewritten),
`workflow/scripts/nerve_clinical_association.py`, `workflow/scripts/nerve_assemble_counts.py`,
`results/README.md`, `.gitignore`, `markdowns/task_conda_env_enforcement.md` (closed, untracked by
design), `markdowns/GBM_TME_Crosstalk_Analysis.md`, `results/fair_validation_report.json`.

**Tool Versions:** validated against 31 declared pins across `workflow/envs/*.yaml`. flake8 clean
on all three scripts (also fixed a pre-existing E501 in `nerve_clinical_association.py`).

**Open Issues:** (1) The one accepted exception — re-deriving the full arm's scVI latent under
torch 2.12.0 — remains a researcher call. (2) **`rule all` is no longer a no-op:** the Tier-3 disk
reclamation removed the per-sample intermediates under `data/processed/<arm>/`, so a dry run now
schedules **365 jobs** (169 ingest + 170 QC + the downstream chain), not the 7 recorded on
2026-08-27. Terminal artifacts are all present; this only matters if someone runs a bare `rule all`
expecting a cheap check. (3) Zenodo deposit still the only off-machine plan for the 740 MB of
models and the 1.4 GB unreproducible `nerve_cells_counts.h5ad`. (4) 4 malformed `tool_versions`
strings remain as warnings in the FAIR report.

**FAIR Notes:** The standing `[FAIR-ALERT]` is retired as a blanket claim and replaced by a named,
reasoned, machine-checked exception. Provenance records are now parsed rather than counted — the
project's Reusable claim is enforced by a rule that can fail, which is what it always asserted.

---

### [2026-09-11] | Phase: Pre-publication fact-lock (blog) | Status: COMPLETE

**Action:** Researcher asked whether the project is ready for a blog post ending in the
interactive marimo dashboard. Assessed state, then locked down the factual items that a
reader checking the repo would find. Scope agreed: full arm, two posts (science +
engineering), dashboard hosting handled by the researcher.

**Outcome:**

1. **`rule nerve_immune_lead_axes` now emits its own selection rule.** The shortlist is two
   halves under genuinely different rules (oligodendrocyte = full arm, batch QC required;
   neuron = capped arm, no QC filter). The rule/module docstrings said so, the emitted row
   did not — `compartment_side` was the only hint. Added per-row `selection_arm` and
   `qc_filter_applied`, surfaced in notebook 05 Panel E, recorded in provenance as
   `selection_rule_counts`. **183 rows before and after, every pre-existing column
   byte-identical**, only the two columns added (full/qc=True 128, capped/qc=False 55).
   Panel E, which recomputes the shortlist independently, still asserts 183/183 on both arms.

2. **`markdowns/blocker_census_annotation_scoring.md` was stale by five weeks and is now
   closed.** It was fixed 2026-08-05 by `c9874d4` as defects D1/D2/D3 of the compartment
   work, which absorbed it under a different filename, so its status line was never updated.
   All three evidence strands verified dead: the cohort resolves **19** cell types, not 5
   (opc 53,805 / ependymal 37,760 / endothelial 4,375); the `mean_confidence` comparison is
   invalid post-fix because the column is now a **z-score margin**, not a raw `score_genes`
   margin, so it is not comparable to the reference's 0.11–0.73; and the T-cell exclusion is
   gone (`source_label` str -> `source_labels` list of 6, immune compartment 1.80x).
   What survives: the two arms disagree on the annotation — capped has **zero** endothelial
   and 61,717 `ambiguous` vs the full arm's 4,375 / 7,172, and `tumor_gbm` differs 10x on a
   1.6x cell-count difference. Documented as a real limit.

3. **Two corrections + one new caveat in the engineering draft.** "9 of 24 nerve groups" fail
   batch QC matched neither arm (actual: 8 of 20 full, 7 of 21 capped); the interface split
   30,622/11,280/2,237 read as though it described the 2,673 significant rows when it
   describes all 44,139 (significant split: 2,118/459/96). New section-9 caveat: five of 23
   nerve clusters sit below the 0.80 neural bar (1,677 cells, 4.4%), four reach the
   interaction tables carrying 12.8% of rows, and **cl12 — 83.5% malignant against the
   oracle — passes `batch_qc_pass`**, because 7 donors at 0.4368 dominant is a healthy donor
   spread. Purity and donor diversity are independent failure modes and only one has a flag.
   Exposure bounded: 27 of 144 axes draw some support from a purity-failing cluster, 3 rest
   on one entirely (`APOE|SCARB1`, `BCAN|EGFR`, `IGSF11|VSIR`, all cl12); `BCAN|EGFR` reaches
   tier 1 on a single LR row.

4. **Science post drafted** — `markdowns/blog_01_science_nerve_immune_crosstalk.md`.

**Incident (recovered, no data loss):** the first re-run used `--forcerun nerve_immune_lead_axes`
*without* `--allowed-rules`, which pulled the opt-in, network-dependent `refresh_drug_annotation`
into the DAG — it would have overwritten the adopted 2026-08-19 snapshot that `config.lead_axes.
drug_annotation` pins. Killed mid-run; Snakemake then deleted that rule's outputs per its usual
behaviour. Both files were tracked and restored via `git checkout`; the restored
`nerve_immune_axis_drug_annotation_refreshed_2026-08-19.csv` hashes to
`0280e1f82a9bd7a4d02774003de7c2830f193ebb88243ea0ad5747da1981c726`, matching `MANIFEST.json`.
`nerve_immune_lead_axes_postfix.csv` was verified byte-identical to a pre-run copy. **The
documented invocation pins `--allowed-rules` for exactly this reason; `--forcerun` alone is not
a safe substitute.**

**Artifacts:** `results/tables/nerve_immune_lead_axes_postfix.csv` (+2 cols),
`provenance/nerve_immune_lead_axes_provenance.json`,
`results/figures/{gbm_cellxgene_56c4912d,gbm_cellxgene_56c4912d_full}/05_census_nerve_immune_explorer.html`,
`workflow/scripts/nerve_immune_lead_axes.py`, `notebooks/05_census_nerve_immune_explorer.py`,
`markdowns/GBM_TME_Crosstalk_Analysis.md`, `markdowns/blog_01_science_nerve_immune_crosstalk.md`,
`markdowns/blocker_census_annotation_scoring.md` (untracked by design).

**Tool Versions:** marimo 0.23.1, pandas 2.3.3, liana 1.7.1 (artifacts unchanged), flake8 clean.

**Open Issues:** Defect 1 (conda-env shadowing) still OPEN and still the blocker for claiming
reproducibility in print — the engineering post's thesis is reproducibility discipline, so it
must be fixed or stated. `results/` is 100% gitignored (0 tracked files), so neither post is
linkable yet; Zenodo deposit or a committed slice is a researcher decision. Branch
`docs/blog-prep-2026-09` is unpushed, as is `chore/disk-reclamation-2026-08` before it.

**FAIR Notes:** The shortlist is now self-describing at row level — a consumer no longer needs
the generator source to know how a row was selected. The standing **[FAIR-ALERT]** on conda-env
shadowing is unchanged and unaddressed.

---

### [2026-08-01 → 2026-08-02] | Phase: Full-cohort rerun unblock (Option A) + Option B tracked | Status: COMPLETE

**Action:** Resume the failed `gbm_cellxgene_56c4912d_full` run (died 2026-07-30 at 343/360 in
`ds_scrna_annotate` on `ImportError: Please install the igraph package`) by taking **Option A** from
`markdowns/blocker_full_cohort_annotate_venv_conda_shadowing.md` — patch the `claude_science` venv —
and open **Option B** (fix the conda-env shadowing) as a separately tracked task.

**[FAIR-ALERT] — READ BEFORE REUSING THESE ARTIFACTS.**
Every `script:`-rule artifact in `data/processed/gbm_cellxgene_56c4912d_full/`,
`results/*/gbm_cellxgene_56c4912d_full/` and `provenance/gbm_cellxgene_56c4912d_full/` was produced
by **`claude_science/bin/python`**, *not* by the `scrna.yaml` conda environment the rules declare.
`--use-conda` was passed and the log prints `Activating conda environment:
.snakemake/conda/ef772b6a…`, but that environment never reaches `sys.path` (Defect 1 below). These
artifacts are **not reproducible from the declared environment spec.** Do not treat them as such.
Re-deriving them under a correctly-enforced `scrna.yaml` may shift Leiden cluster IDs and cascade
through annotation, purity, concordance, and every downstream table. Known unclosed drift: venv
`torch 2.11.0` vs the `scrna.yaml` pin `torch 2.12.0`.

**Outcome:**

*Option A — venv patch (dry-run gated).* `pip install --dry-run` was run first and the resolver plan
inspected; it proposed only new packages and would not move `numpy`/`pandas`/`scanpy`/`scvi-tools`/
`torch`/`anndata`/`scipy`/`matplotlib`, so the install proceeded. Installed at the `scrna.yaml`
pins: `igraph==0.11.8`, `leidenalg==0.10.2`, `gseapy==1.1.3`, `liana==1.7.1`, `decoupler==2.1.6`
(plus transitive `adjustText`, `legendkit`, `marsilea`, `mizani`, `plotnine`, `texttable`). The
`pip freeze` delta is **11 pure additions, zero modified pins** — see
`logs/venv_freeze_pre_optionA_20260801.txt` vs `logs/venv_freeze_post_optionA_20260801.txt`. That
diff, not this prose, is the authoritative record of what Option A changed.

*Scope correction — it was four missing packages, not two.* The blocker doc proposed installing only
`igraph`+`leidenalg`. An audit of every module-level import across all 17 remaining jobs found
`gseapy` (`workflow/scripts/nerve_cell_heterogeneity.py:192`) and `liana`
(`nerve_tumor_interaction.py:40`, `nerve_tumor_immune_interaction.py:42`) also absent and unguarded.
The two-package fix would have failed twice more. `infercnvpy` is absent but genuinely unused —
`scrna_malignancy.py` uses a custom CNV path — and was deliberately **not** installed.

*Root-cause correction — the blocker doc's regression window was wrong.* It guessed commit `c2a3d31`
reinstalled from `requirements.txt` and dropped `igraph`. `CHANGELOG.md:1345` records that the repair
only rewrote 76 hardcoded `Claude_Setup/` path prefixes in `claude_science/bin/` and installed
nothing. Before it, `claude_science/bin/snakemake` was broken and `/opt/anaconda3/bin/snakemake` was
the real driver — and anaconda base has `igraph 1.0.0` / `leidenalg 0.11.0`, which is why the capped
arm's Leiden passed on 2026-07-20. `c2a3d31` made the venv's snakemake runnable; the 2026-07-29
runbook opened with `source claude_science/bin/activate`; the shadowing interpreter silently changed
from anaconda base to the venv. **Base is not a fallback** — base `scanpy` now fails to import on a
numpy/h5py ABI mismatch (`numpy.dtype size changed, Expected 96 from C header, got 88`).

*`ds_scrna_annotate` — re-run in isolation as a smoke gate, exit 0.* 47 Leiden clusters at
resolution 1.0 over **1,006,344 cells** (matches the post-QC count exactly); 24,139/24,139 genes
mapped to symbols. **Peak RSS 27.22 GB** (27,218,329,600 B), wall clock **892.84 s = 14 m 53 s**.
This settles the attribution question the blocker doc raised: annotate alone accounts for 27.22 of
the failed run's 29.40 GB whole-run peak, so `ds_scrna_annotate` — **not** `ds_scrna_integration` —
is the memory ceiling, and the 2026-07-28 integration remediation is not implicated. It is still
well above the runbook's "expect < 20 GB" and leaves ~3 GB against the ~30 GB usable budget.

*Remaining 16 jobs — running.* Dry-run confirmed 16 jobs, all `ds_*` (the pinned v1.3.0 insulator
holds). Deviation from the documented resume command, taken on the 27.22 GB measurement and
confirmed with the researcher: **`--resources mem_mb=28000` added.** `ds_nerve_cell_subset` and
`ds_immune_cell_subset` both depend only on `malignancy_labeled.h5ad` and are mutually independent,
so they could otherwise co-schedule. Every remaining rule declares `mem_mb=28000`, so the flag
forces strict one-at-a-time execution. The 2026-07-29 plan omitted it to avoid serializing the 340
ingest/QC jobs; those are complete, so the cost is now zero. Per `CLAUDE.md` this is a **scheduler
gate only** — it cannot cap a single process's RAM.

**Tooling defect found (`--allowed-rules` argument order).** `--allowed-rules` takes `nargs='+'` and
silently swallows target paths placed after it, leaving Snakemake with no target; it then falls back
to rule `all` and raises a misleading `MissingInputException` naming the five final artifacts. Pass
targets **first**. Same class of footgun as `--rerun-triggers=mtime` requiring the `=`. Confirmed on
snakemake 9.19.0. Documented in the blocker doc, `markdowns/task_conda_env_enforcement.md`, and
memory — this pattern is used as the pinned-reference insulator throughout the project, so getting
the order wrong is a live risk to that guarantee.

*Resume batch — 13/16 jobs completed, then `ds_nerve_scanvi_retrain` **SEGFAULTED**.* Wall clock
36,160.88 s = **10 h 02 m 41 s**; **peak RSS 32.78 GB** (32,782,401,536 B) on a 36 GB machine.

The crash is **not** an ImportError — Option A held through all 13 completed rules, including the
first-ever executions of `liana` (both LIANA rules) and `gseapy` in this environment. It is:

```
EXC_BAD_ACCESS (SIGSEGV), KERN_INVALID_ADDRESS at 0x580
  libomp.dylib  __kmp_suspend_initialize_thread
  libomp.dylib  __kmp_fork_barrier / __kmp_launch_worker
```

— an OpenMP worker thread failing to initialise inside `sc.pp.neighbors()`
(`nerve_scanvi_retrain.py:147`), 21 s after the "Computing neighbors + UMAP" line. Thread-stack
allocation failing under memory exhaustion presents exactly this way, as a null deref rather than a
clean OOM. The capped arm runs this same code path fine at 270,520 cells (~1 m 39 s); this arm has
**377,343**. Crash report: `~/Library/Logs/DiagnosticReports/python3.12-2026-08-02-001047.ips`.

**`--resources mem_mb=28000` did not and cannot prevent this** — per `CLAUDE.md` it is a scheduler
gate only and cannot cap a single process's RAM. It did do its actual job (no co-scheduling).

*Script-ordering defect — the real cost.* scANVI had **finished training** (classifier accuracy
**0.5979**; capped arm 0.6369) before the crash, but `scanvi_model.save()` sat at line 153, *after*
the `neighbors`/`umap` calls at 147–148. A crash in a purely visual step therefore destroyed
**2 h 22 m** of completed training. Line 99 shows the author had already reasoned about this exact
failure mode for stage 1 ("Save the baseline so a scANVI-stage crash doesn't lose the long
pretrain") but did not extend the guard to stage 2. Compounding it, the stage-1 checkpoint is
**written but never read** — line 84 unconditionally constructs a fresh `SCVI` and line 91
unconditionally trains — so the protective intent in that comment was never actually implemented and
a re-run pays the full ~7 h again.

**CODE CHANGE (researcher-approved, option "reorder save + free memory"):**
`workflow/scripts/nerve_scanvi_retrain.py` —
1. `scanvi_model.save()` moved to immediately after `predict()`, before any visualization.
2. `del scanvi_model` + `gc.collect()` + `torch.mps.empty_cache()` inserted before
   `sc.pp.neighbors`, releasing model weights and torch's cached MPS blocks so the kNN graph has
   headroom.
3. Removed a pre-existing unused `import pandas as pd` (the file's only flake8 finding; confirmed
   pre-existing at HEAD, and `pd.` has zero occurrences).

`flake8` exit 0. **No scientific impact** — identical model, identical seed, merely persisted
earlier with memory released; `scanvi_history` is captured at line 124 into a plain dict, so the
training-curves figure does not need the model after the `del`. Dry-run confirmed the edit does
**not** schedule the pinned baseline `nerve_scanvi_retrain` rule — `--rerun-triggers=mtime`
correctly suppresses the code trigger, so the v1.3.0 freeze is intact. Caveat recorded honestly:
the memory release is a well-founded mitigation, **not a guarantee**. If the segfault was a genuine
OpenMP runtime conflict (torch's `libomp` vs numba/pynndescent's) rather than memory pressure, it
can recur — but it would then cost only the UMAP step, since the model is already on disk.

*Attempt 2 also SEGFAULTED — and disproved the memory hypothesis.* The mitigation cut peak RSS to
**16.05 GB** (from 32.78 GB) and `sc.pp.neighbors` crashed anyway with the **identical** signature.
Memory was therefore never the cause. The crash report's loaded images identified it:
`libomp.dylib` (torch/sklearn) + `libtbb.12.8.dylib` + **`omppool.cpython-312-darwin.so`** — numba's
OpenMP threading pool — coexisting in one process. numba-jitted pynndescent inside `sc.pp.neighbors`
spins up an OpenMP worker that collides with the `libomp` torch/MPS already initialised.

**The save-before-UMAP fix worked at the script level but was defeated by Snakemake.** The 66 MB
model was verified on disk mid-run, then removed: `Removing output files of failed job
ds_nerve_scanvi_retrain since they might be corrupted`. It is a declared `output:`, so Snakemake
cleans it on failure. **The stage-1 `nerve_scvi_baseline` survived precisely because it is NOT a
declared output** — the undeclared checkpoint lived, the declared one died. This is a general lesson
for long-training rules in this project.

*Root cause proven by controlled experiment*, not inference — 200,000 synthetic cells × 30 dims, no
project data, torch/MPS initialised then `sc.pp.neighbors`:

| Arm | `NUMBA_THREADING_LAYER` | Result |
|---|---|---|
| Control | default (OpenMP) | **exit 139 — SIGSEGV** |
| Treatment | `workqueue` | **exit 0**, `NEIGHBORS OK (200000, 200000)` |

**SECOND CODE CHANGE (researcher-approved):** `workflow/scripts/nerve_scanvi_retrain.py` —
1. `os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")` set **before** any import that
   pulls in numba (scanpy → pynndescent/umap). `E402` is project-ignored per `.flake8`.
2. Stage-1 baseline is now **loaded when present** instead of always retraining — implementing what
   the line-99 comment had promised since it was written but never did (line 84 unconditionally
   constructed a fresh `SCVI`). Cut the retry from ~7 h to 2 h 35 m.
3. Stage-2 model also written to an **undeclared sidecar** (`nerve_scanvi_stage2`) so Snakemake's
   failure cleanup cannot destroy it again.

`flake8` exit 0. `history_` is persisted inside `model.pt` (verified: 400 epochs of loss series), so
the training-curves figure works from a loaded checkpoint.

*Attempt 3 — SUCCESS, exit 0.* Baseline reuse confirmed in the log (`Reusing baseline scVI
checkpoint … (skipping stage-1 pretrain)`), stage-1 skipped in <1 s vs 4 h 36 m 55 s of retraining.
`sc.pp.neighbors` + UMAP completed in **3 m 18 s** (22:02:50 → 22:06:08) — the step that segfaulted
twice. Wall clock **9,345.54 s = 2 h 35 m 46 s**, **peak RSS 15.02 GB**. Classifier accuracy
**0.5979** reproduced bit-identically across all three attempts, confirming seed control.

**GOAL-BACKWARD VERIFICATION — all five targets present and non-empty:**

| Target | Size |
|---|---|
| `cohort_concordance_summary.json` | 4.0K |
| `nerve_cluster_sample_purity_v2.csv` | 4.0K |
| `nerve_celltype_label_summary.csv` | 4.0K |
| `nerve_scanvi_training_curves.png` | 104K |
| `05_census_nerve_immune_explorer.html` | 896K |

- **Capped arm untouched** — zero files under `data/processed/gbm_cellxgene_56c4912d/` or
  `results/tables/gbm_cellxgene_56c4912d/` modified since 2026-07-31.
- **Pinned v1.3.0 reference intact — freshly re-verified 2026-08-02T22:19 UTC: `pass: true`,
  37 checked / 37 unchanged / 0 modified / 0 missing.** (The on-disk verification JSON had been
  stale from Jul 28; it was force-re-run via `--allowed-rules verify_pinned_reference --forcerun`,
  dry-run confirmed 1 job.)
- `nerve_cluster_sample_purity_v2`: full arm **43 clusters, 27 pass / 16 fail**; capped arm 33
  clusters, 24 pass / 9 fail. More clusters at the same resolution and a higher absolute fail count
  — consistent with the finer partition at 1.66× depth; the small-tail failures are the same class
  flagged on 2026-07-25 and remain a researcher-review item.

**SCIENTIFIC RESULT — the arm comparison (`cohort_concordance_summary.json`), the payoff:**

| Metric | Capped (614,951 cells) | Full (1,020,902 cells) | Δ |
|---|---|---|---|
| Reference significant pairs | 3,368 | 3,368 | — (same reference) |
| Dataset significant pairs | 4,583 | 5,034 | +451 (+9.8%) |
| Shared pairs | 2,539 | 2,619 | +80 (+3.2%) |
| Union pairs | 5,412 | 5,783 | +371 |
| **Jaccard overlap** | **0.4691** | **0.4529** | **−0.0162** |
| **Spearman ρ (shared)** | **0.6249** | **0.6083** | **−0.0166** |

**1.66× the cells did not improve concordance — both summary metrics moved slightly down.** But
Jaccard is symmetric and penalises the larger set; decomposing it: reference **recall improved**
(2,539/3,368 = 75.4% → 2,619/3,368 = **77.8%**, +2.4 pp) while **precision fell** (55.4% → **52.0%**,
−3.4 pp). The added depth surfaced 451 more significant pairs, only 80 of which were reference
pairs. Spearman ρ on shared pairs is essentially flat (−2.7% relative), i.e. the depth changed
*which* pairs clear significance, not how the shared ones rank. Whether the 371 extra non-reference
pairs are depth-limited discoveries or higher-power false positives is a **scientific call for the
researcher**; note the reference is the 17-sample v1.3.0 baseline, so "absent from the reference"
carries limited evidential weight against a 170-sample cohort.

**Other results from the completed 13:** `ds_scrna_malignancy` CNV threshold 0.0056, **125,305/
1,006,344 malignant (12.5%)**, 117,134 reference cells. `ds_immune_cell_subset` **329,608 microglia
(32.8%)** → 28 clusters, batch purity 18/28 PASS. `ds_nerve_cell_subset` **377,343 non-malignant
nerve cells (37.5%)** — the silent-placeholder path (`[[nerve-subset-silent-placeholder]]`) did
**not** trigger. Both LIANA rules passed the raw-counts scale guard identically
(`X.max() 53027.000 → 8.740`), confirming the 2026-07-26 per-consumer normalization fix holds:
35,255 tumor–nerve LR rows / 2,090 sig pairs; 84,422 three-way rows (immune–nerve 46,777,
nerve–tumor 35,255, immune–tumor 2,390) / 5,223 sig pairs. GSEA: **80 result blocks, 0 clusters
skipped**. `ds_scrna_annotate` 24,139/24,139 genes mapped.

**Note on `Unknown` = 20%:** `nerve_celltype_labels.py:120` sets `Unknown` as the bottom
`unknown_percentile` (default 20th) of max marker score, so it is **20.0% by construction in every
cohort** and carries no biological signal — do not compare the Unknown *fraction* between arms. Two
label groups do score on reduced panels in **both** arms: `oligodendrocyte` 3/4 markers,
`ependymal` 5/6. That caveat is real and is not depth-related.

**Artifacts:** `data/processed/gbm_cellxgene_56c4912d_full/annotated.h5ad` (5.19 GB, sha256
`5c32c067…`, artifact_id `1b2f56e9-578b-53f4-aa26-53d65f5a6fd8`);
`results/tables/gbm_cellxgene_56c4912d_full/annotation_summary.csv`;
`provenance/gbm_cellxgene_56c4912d_full/annotation_provenance.json`;
`data/processed/gbm_cellxgene_56c4912d_full/malignancy_labeled.h5ad`;
`logs/full_cohort_annotate_20260801_135036.log`; `logs/venv_freeze_{pre,post}_optionA_20260801.txt`;
`markdowns/task_conda_env_enforcement.md` (new).

**Tool Versions:** snakemake 9.19.0 (`claude_science/bin/snakemake`); scanpy 1.12.1; anndata 0.12.10;
scvi-tools 1.4.2; torch 2.11.0 (**pin says 2.12.0**); igraph 0.11.8; leidenalg 0.10.2; gseapy 1.1.3;
liana 1.7.1; decoupler 2.1.6; numpy 2.3.5; pandas 2.3.3.

**Open Issues:**
- **Defect 1 (conda-env shadowing) remains OPEN** → `markdowns/task_conda_env_enforcement.md`.
  Blast radius is every `script:`-rule artifact in the repo's history, both arms — not just this run.
- **`numba` is unpinned in `scrna.yaml`** (0 matches) — the threading conflict is downstream of
  Defect 1: the env that actually executes is the venv, which drifted. The capped arm passed this
  same code on 2026-07-25 under anaconda base (a different numba build); base is now broken outright
  (`numpy.core.multiarray failed to import`). Pin `numba` when Option B is actioned.
- **Scientific call for the researcher:** are the 371 extra non-reference LR pairs in the full arm
  genuine depth-limited discoveries or higher-power false positives? See the concordance table.
- **Snakemake deletes a failed rule's declared outputs.** Any mid-rule checkpoint meant to survive a
  crash must be written to an undeclared sidecar path. Applied here to stage 2; worth auditing other
  long-training rules for the same exposure.
- No cluster was assigned `tumor_gbm` or `inhibitory_neuron` at resolution 1.0 despite both being
  scored, and `astrocyte` took 391,369 cells at mean confidence 6.09 — while malignancy flagged only
  125,305 cells (12.5%) overall, so the large astrocyte compartment is mostly **not** malignant by
  CNV. Plausible for a multi-donor Census cohort carrying non-tumour brain tissue, but worth a look
  rather than an assumption. Researcher review.
- 47 clusters here vs 33 in the capped arm — cluster IDs are **not** comparable one-to-one between
  arms. Use the concordance summary, not cluster-number matching.
- **`--allowed-rules` argument order** (see above) is a live risk to the pinned-reference insulator,
  since that flag is exactly how the v1.3.0 freeze is protected across this project.

**FAIR Notes:** The `[FAIR-ALERT]` above is the substantive one. Provenance JSON with input/output
sha256 and `artifact_id` was written for every completed rule. `provenance/baseline_v1.*.json` was
not touched. No `workflow/envs/*.yaml`, `Snakefile`, or analysis-script changes were made — the
`scrna.yaml` spec is correct; it simply never loads.

---

### [2026-07-28] | Phase: Repository rename + constitution path fix | Status: COMPLETE
**Action:** Rename the GitHub repo `Nerve_Analysis_TCGA_GBM` →
`GBM_Nerve_Tumor_Immune_Single_Cell_Analysis` to match project scope, and correct the stale
project-root declaration in `CLAUDE.md`.
**Outcome:** Renamed via GitHub API (HTTP 200; `gh` is not installed — used the osxkeychain git
credential). `origin` updated with `git remote set-url`; fetch verified. GitHub redirects the old
URL, so existing clones keep working — **the old name must never be reused** or redirects break.

**Deliberate decision: the local directory was NOT renamed.** It stays
`Projects/Nerve_Analysis_TCGA_GBM/`. Renaming it would (a) invalidate ~6.5 GB of
`.snakemake/conda/` environments, whose hashes include the prefix path, forcing a full rebuild,
and (b) orphan the Claude Code project memory and session history, which are keyed on the
filesystem path. Neither cost buys anything: `config/config.yaml` contains zero absolute paths and
all 37 `baseline.pinned_artifacts` are relative, so no rule, target, or drift check depends on the
directory name. `CLAUDE.md` now documents the name mismatch as intentional.

`CLAUDE.md:7` had declared the root as `…/Projects/CLAUDE_TCGA_GBM/`, which never existed, and its
path-correction rule mapped stale names onto `NERVE_ANALYSIS_TCGA_GBM` (wrong case, also not the
real directory). Both corrected to the actual `Nerve_Analysis_TCGA_GBM`, with an explicit
exception added: `provenance/*.json` absolute paths and `git_remote_url` values are historical
FAIR records and must be left stale rather than rewritten.
**Artifacts:** `CLAUDE.md`, `CHANGELOG.md`; GitHub repo renamed; `.git/config` remote URL updated.
**Tool Versions:** GitHub REST API 2022-11-28; git 2.x. No analysis code executed.
**Open Issues:** `provenance/baseline_v1.0.0`–`v1.3.0.json` record the old remote URL and old
absolute artifact paths — intentionally left as-is per the rule above. PR #1 still merged-but-stale.
**FAIR Notes:** No artifacts produced or modified; provenance deliberately untouched (65 JSONs).
Findability is preserved by GitHub's permanent redirect from the old URL.

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

---

### [2026-07-28] | Phase: Local directory rename — reverted | Status: COMPLETE

**Action:** The local project directory had been renamed `Nerve_Analysis_TCGA_GBM` → `GBM_Nerve_Tumor_Immune_Single_Cell_Analysis` to match the repo. Audited what actually depends on that path, then reverted the rename.

**Outcome:**
- **Reverted** the local dir to `Nerve_Analysis_TCGA_GBM`. The local/remote name mismatch is restored and remains deliberate per `CLAUDE.md:9-15`.
- **Why it mattered — conda envs.** Snakemake 9.20 mixes `realpath(envs_dir)` into the env hash on purpose (`snakemake/deployment/conda.py:236-240`: *"moving the working directory around automatically invalidates all environments… hardcoded absolute RPATHs"*). Recomputed both ways: `scrna.yaml`→`ef772b6a…` (1.6 GB) and `notebooks.yaml`→`8e2fe802…` (878 MB) resolve to existing env dirs **only** under the old name. The rename stranded ~2.5 GB of live envs; the revert recovered them without a rebuild.
- **Why it mattered — FAIR PIDs.** `workflow/scripts/fair_utils.py:19` seeds `uuid5` on the *resolved absolute path*, so every artifact_id shifted (`nerve_cells.h5ad`: `6119a339-…` → `c8198df0-…`). Post-revert it returns `6119a339-a1c1-501d-a543-c5b5eb7a25f1` again, matching existing provenance.
- **Not affected** (verified, no action needed): `Snakefile`/`*.smk` (relative `configfile`, no `workdir:`), all `config.yaml` paths, all five marimo notebooks (`Path(__file__).parent.parent`), git, and Snakemake's output metadata (records `conda_env` as base64 of env **YAML content**, not path — so no spurious re-runs).
- **Corrected the "6.5 GB" figure in `CLAUDE.md`:** only ~2.5 GB of `.snakemake/conda/` is live. The other three env dirs (~4 GB: `465f7a09…`, `5a394ede…`, `f7316445…`) were already orphaned by earlier edits to `base.yaml`/`proteomics.yaml`, which currently resolve to no env on disk and will rebuild on their next run regardless.
- **Repaired the `claude_science/` venv** (pre-existing breakage, unrelated to the rename): 76 files in `claude_science/bin/` hardcoded `…/Projects/Claude_Setup/claude_science/`, a directory that has never existed here. `./claude_science/bin/snakemake` failed outright, so `execution_instructions.md` was pointing at a venv that could not run. Rewritten to the correct prefix; now returns 9.19.0. Note this is **older than the 9.20.0 on PATH** at `/opt/anaconda3/bin/snakemake`, which is what has actually been running the pipeline.

**Artifacts:** `CLAUDE.md` (plans path corrected `Claude_TCGA_GBM` → `Nerve_Analysis_TCGA_GBM`; committed `81ff3e0`). `.claude/plans/project_writeup_skeleton.md` (placeholder repo URL resolved; gitignored). `claude_science/bin/*` (76 files; gitignored). `../TCGA_GBM_scRNA-seq_2.code-workspace` (folder path fixed; outside repo). `results/pinned_reference_verification.json` + `provenance/pinned_reference_v1.3.0.json` (see FAIR Notes). `CHANGELOG.md` (this entry).

**Tool Versions:** snakemake 9.20.0 (`/opt/anaconda3/bin`), python 3.12.

**Open Issues:**
- **Root cause unfixed.** `fair_utils.py:19` seeding PIDs on absolute paths means any future move silently re-issues every identifier. Making `artifact_id`/`stamp_artifact` project-relative would fix it permanently, but re-issues every PID once and requires re-freezing v1.3.0 — deferred to its own change.
- `README.md:40` and `execution_instructions.md:24` were left untouched: they already read `Nerve_Analysis_TCGA_GBM` and became correct again on revert.

**FAIR Notes:** All ~60 files under `provenance/` were left stale per `CLAUDE.md:24-27` — the 9 pinned provenance JSONs still carry the original absolute path, which is precisely why their sha256 values did not drift. `verify_pinned_reference` passes: **37/37 unchanged, 0 modified, 0 missing**, and it activated the recovered env `8e2fe802…` rather than building one — direct confirmation the envs were reused. Caveat on that check: the run first re-executed `freeze_pinned_reference` (pre-existing "code has changed" trigger), rewriting the manifest's `frozen_at_utc`, so the pass is partly self-referential. The independent evidence is stronger: the newest recorded artifact mtime is **2026-07-22**, six days before this session, so none of the 37 pinned artifacts were touched. The 5 tracked `baseline_v1.*.json` files are unmodified (`git status` clean).

---

### [2026-07-28] | Phase: Memory-efficiency remediation (items 1-7) + uncapped cohort namespace | Status: COMPLETE (code only — no pipeline run)
**Action:** Investigated whether polars could relieve the memory pressure that forced
`subsample_per_donor: 5000` (which discards 47.7% of the Census cohort). Answer: no — the footprint
is a sparse count matrix, not a DataFrame. Measured `annotated.h5ad`: `X` = 9.82 GB vs `obs` = 121 MB,
so polars addresses ~1.3% of it. Instead applied the seven memory-hygiene fixes from
`markdowns/assessment_pipeline_memory_efficiency.md` and added an uncapped dataset namespace.

**Headline finding — the documented hardware budget was wrong by 3.5x.** `system_profiler` reports
**36 GB** unified memory (`sysctl hw.memsize` = 38654705664). `CLAUDE.md:47`, `config/config.yaml:36`
and `project_plan_orchestration.md:20` all claimed 128 GB. Consequently `CHANGELOG.md:442`'s account
of the 2026-07-21 OOM ("peaking >120 GB on the 128 GB machine") describes a peak this machine cannot
reach — the process was SIGKILLed far earlier. The chunking fix applied then was still correct.

**Outcome (all seven items applied):**
1. Hardware facts corrected in the three docs; `max_memory_gb: 120 -> 30`; `default_mem_mb: 32000 -> 28000`.
   The two unsatisfiable `mem_mb = 64000` declarations (`datasets.smk:480`, `:595`) now use the default.
2. `scrna_integration.py` — `adatas.clear(); del adatas` after `ad.concat`. The per-sample list
   (~12 GB Census / ~21 GB baseline) was previously held alongside the concatenated copy.
3. `scrna_integration.py` — the identity `layers["counts"]` copy is no longer materialized when `.X`
   already holds raw counts; scVI now gets `layer=None` for those cohorts. **Numerically identical**
   (the layer was a bit-identical copy of `.X`), and it matches how the existing Census
   `integrated_latent.h5ad` was actually built per `markdowns/plan_recover_baseline_raw_counts.md`.
   The baseline arm is untouched — `integration.smk:34` hard-codes `counts_from_log1p = True`.
   Added a CSR-index dtype guard that warns if nnz ever promotes indices to int64.
4. `scrna_malignancy.py` — the reference block is now streamed in row-blocks like the main loop
   (was 55,477 cells x 24,048 genes = 5.3 GB in a single `.toarray()`, the last unchunked
   densification). `cnv_chunk_size: 50000 -> 10000` (~9.6 GB -> ~2 GB working set).
5. Both LIANA scripts release `malig_full` / `nerve` / `immune` after subsetting/concat; cell counts
   captured beforehand for the provenance blocks.
6. `backed="r"` for `nerve_batch_qc`, `nerve_clinical_association`, `nerve_leiden_resolution_sweep`
   — all three read only `obs`/`obsm`. **`nerve_batch_qc_v2` was deliberately NOT converted**: it
   rewrites its own input file in place, and writing to a file held open in backed mode is unsafe.
7. All 14 `write_h5ad` call sites now pass `compression=H5AD_COMPRESSION` (new constant in
   `fair_utils.py`). `data/processed` was 143 GB of uncompressed HDF5 on a disk that is 74% full.

**Uncapped cohort:** new `datasets:` entry `gbm_cellxgene_56c4912d_full` — same study filter, same 170
donors, `subsample_per_donor: null` -> **1,020,902 cells** (vs 614,951 capped). Deliberately a separate
namespace so the capped cohort's 77 GB of artifacts survive as a comparison arm.
`data/raw/gbm_cellxgene_56c4912d_full/{samples.txt,MANIFEST.txt}` generated via
`scripts/derive_dataset_sample_sheet.py` (backed metadata read).

**Feasibility (measured, not estimated):** sampled 3,000 cells from `annotated.h5ad` — 1,995 nnz/cell,
51.6% of nonzeros in the top-3,000 genes. Full study at 24,048 genes ~= 2.04e9 nnz = **16.3 GB** in
memory, and 2.04e9 < 2**31 so int32 CSR indices remain legal. Fits in 36 GB *only because* of items
2+3, which freed ~22 GB. Disk: uncompressed the run would need ~128 GB against 116 GiB free — item 7
is load-bearing for this run, not cosmetic.

**Artifacts:** `markdowns/assessment_pipeline_memory_efficiency.md` (new); `CLAUDE.md`,
`project_plan_orchestration.md`, `config/config.yaml`, `workflow/rules/datasets.smk`,
`workflow/scripts/{fair_utils,scrna_integration,scrna_malignancy,nerve_tumor_interaction,
nerve_tumor_immune_interaction,nerve_batch_qc,nerve_clinical_association,
nerve_leiden_resolution_sweep}.py` + 14 write-site edits across the AnnData-emitting scripts;
`data/raw/gbm_cellxgene_56c4912d_full/`.

**Verification:** `flake8` output diffed against `HEAD` for every changed file — **identical**, no new
errors (pre-existing ones left alone per `CLAUDE.md` "do not mass-reformat"). All changed `.py` files
byte-compile. `snakemake --dry-run` builds a valid 772-job DAG. Scoped dry-run on the five `_full`
targets with `--rerun-triggers=mtime` yields **360 jobs, every one `ds_*` under the `_full`
namespace**; the only non-namespaced paths in the DAG are read-only *inputs* (the pinned v1.3.0
concordance comparator and `nerve_crosstalk_lead_targets.csv`). No existing artifact is an output.

**Tool Versions:** anndata 0.12.10, scanpy 1.12.1, scvi-tools 1.4.2, scipy 1.17.1, numpy 2.3.5,
h5py 3.16.0, torch 2.12.0.

**Open Issues:**
- **The pipeline has NOT been run.** Researcher asked for the edits only. See the invocation below.
- **A bare `snakemake` is now unsafe.** Script edits tripped the `code` rerun-trigger, so the default
  invocation schedules 772 jobs including the 17-sample baseline `scrna_integration`. Always pass
  `--rerun-triggers=mtime` and explicit targets.
- **Item 4 is not bit-exact for malignancy.** The streamed reference mean accumulates in float64 and
  divides, where the old code called `.mean()` on a float32 block. This is strictly more accurate but
  may shift `cnv_threshold` in the last decimals and flip a handful of borderline cells. It matches
  the pattern the no-reference branch in the same script already used. `malignancy_labeled.h5ad` is
  not in `baseline.pinned_artifacts`, so nothing frozen depends on it.
- Items 8-10 (HVG-before-scVI, `.raw` after `subset=True`) and 11-12 (peak-RSS logging, polars port
  of the five table workloads) remain unimplemented.
- `scrna.n_top_genes: 3000` is still **unused by scVI** — `scrna_integration.py` does no HVG
  selection, so the run below trains a full 24,048-gene decoder. This is item 9, deferred as a
  science change.

**FAIR Notes:** No provenance JSON was rewritten. The pinned v1.3.0 reference is untouched — it is
read as the concordance comparator only. The new cohort namespaces every artifact under
`gbm_cellxgene_56c4912d_full/`, so `baseline.pinned: true` and its 10 undefined rules stay in force.

**Invocation for the full-cohort run (NOT executed):**

This machine is configured to idle-sleep after **1 minute** (`pmset -g custom`: `sleep 1`,
`disksleep 10`). macOS idle sleep keys off user input, not CPU load, so an unattended multi-hour
job WILL be suspended mid-run unless a power assertion is held. `ttyskeepawake 1` covers it only
while a tty stays connected. Hence `caffeinate -ims` (idle + disk + system sleep; display is left
free to sleep). Assertions are held only for the lifetime of the wrapped command.

`/usr/bin/time -l` captures `maximum resident set size` across the largest reaped child — the only
way to get a peak-RSS number until item 11 (per-rule RSS logging) is implemented. Expect < 20 GB
if items 2-4 worked; the pre-fix code peaked around 32 GB on a 36 GB machine.

```
caffeinate -ims /usr/bin/time -l \
snakemake --use-conda --cores all --rerun-triggers=mtime \
  results/tables/gbm_cellxgene_56c4912d_full/cohort_concordance_summary.json \
  results/tables/gbm_cellxgene_56c4912d_full/nerve_cluster_sample_purity_v2.csv \
  results/tables/gbm_cellxgene_56c4912d_full/nerve_celltype_label_summary.csv \
  results/figures/gbm_cellxgene_56c4912d_full/nerve_scanvi_training_curves.png \
  results/figures/gbm_cellxgene_56c4912d_full/05_census_nerve_immune_explorer.html \
  2>&1 | tee logs/full_cohort_run_$(date +%Y%m%d_%H%M%S).log
```

Run it inside `tmux`/`screen` (or under `nohup`) — closing the terminal SIGHUPs `caffeinate`, which
kills the wrapped snakemake and drops the assertions. Verify the assertion is live during the run
with `pmset -g assertions | grep -i caffeinate`.

Note on `--cores all` (14): `resources: mem_mb` is ignored by the scheduler unless
`--resources mem_mb=N` is also passed. It is deliberately omitted here — every rule declares the
28000 default, so supplying it would serialize all 170 ingest/QC jobs (each of which is a small
per-donor backed read) and cost far more wall-clock than it protects. The heavy rules
(`ds_scrna_integration`, `ds_scrna_malignancy`, `ds_nerve_scanvi_retrain`) are singletons and run
alone regardless.

---

## [2026-08-05] Phase 0-1 — Compartment integrity fix: plan review + baseline audit (Test Oracle)

**Phase:** Compartment integrity fix (`markdowns/plan_compartment_integrity_fix.md`), Phases 0-1 of 6.

### Action — Phase 0: plan reviewed against the code

Verified every structural claim in the plan against `workflow/rules/datasets.smk`, the scripts, and
the on-disk artifacts. **All seven defects (D1-D7) are real and correctly located.** Four
corrections (C1-C4) were written into the plan; they change *how* the fix is built, not *what* it
fixes:

- **C1** — `dataset_gene_symbol_map.py:35` writes `chromosome = "unknown"` for all 61,497 genes, so
  D4's "the symbol map already carries a chromosome column" is not executable. Added a new
  `ds_gene_positions` rule (Ensembl 113 GTF, SHA-pinned like `download_msigdb_gmt`) feeding
  `infercnvpy`. **This is new work Phase 3 did not previously account for.**
- **C2** — `nerve_cells.markers` is passed to `ds_scrna_qc` (`datasets.smk:165`) across 170 donors x
  2 arms. Putting the new TME panels there would invalidate every QC artifact and therefore
  `ds_scrna_integration` — the 11h13m scVI train the plan exists to preserve. New panels go in a
  separate `annotation_markers:` block.
- **C3** — `immune_cell_subset.py:55` matches a single `source_label` string exactly. Adding a
  macrophage panel drops those cells out of the immune compartment entirely, and the "purity >=95%"
  gate would still pass because purity is not a size check. `source_label` -> `source_labels` list,
  plus a size gate.
- **C4** — `nerve_scanvi.labels_key` and `nerve_cells.leiden_resolution` are global keys the
  reference cohort also reads; both need per-arm scoping.

### Action — Phase 1: `ds_compartment_audit` built and run on existing artifacts

New `workflow/scripts/compartment_audit.py` + `ds_compartment_audit` rule + `compartment_audit:`
config block. Reads **obs only** (never `.X`) from artifacts that already exist; both arms complete
in ~13 s. Registered as a terminal target in `rule all`.

Useful discovery: **`cell_type` survives in `nerve_cells.h5ad` and `immune_cells_labeled.h5ad`.**
D6's overwrite only hits `nerve_cells_v2.h5ad` (the scANVI side-branch), so the audit needs no
join back to the 1M-cell parent object.

### Outcome — baseline reproduces the S1PR1 report exactly

| metric | full arm | capped arm | gate |
|---|---|---|---|
| nerve neural fraction | **0.1108** | 0.1155 | >=0.85 FAIL |
| nerve malignant / myeloid / vascular | 0.5888 / 0.2754 / 0.0187 | 0.6155 / 0.2469 / 0.0149 | — |
| tumor compartment malignant fraction | **0.4790** | 0.4559 | >=0.85 FAIL |
| malignancy recall | **0.1799** | 0.2482 | >=0.80 FAIL |
| max nerve-cluster endothelial fraction | **0.9289** (c24) | 0.4526 (c16) | <=0.20 FAIL |
| immune purity | 0.9946 | 0.9978 | >=0.95 PASS |

Every figure in `markdowns/plan_compartment_integrity_fix.md` and
`markdowns/s1pr1_localization_report.md` is reproduced to the stated precision (11.1% / 58.9% /
27.5% / 47.9% / 18.0% / 99.5% / 92.9%). The tumor-compartment contamination breakdown also matches:
25.5% macrophage, 10.6% microglia, 5.4% oligodendrocyte, 4.3% monocyte. 5 of 40 nerve clusters are
neural-dominant; 23 are malignant-dominant and 9 myeloid-dominant.

The finding is now **reproducible on disk** rather than existing only as prose — three CSVs plus a
machine-readable gate table per arm, with provenance.

Per-arm size gates measured and written into config (C4): immune baseline n = 329,608 / 169,617;
nerve bands 45-75k (full, Census neuroglial total 62,632) and 30-52k (capped, total 43,857).

### Tool failure worth recording

Running the new rule **without `--use-conda`** fails with
`/bin/bash: /Users/jarrettevans/Documents/Biomedical: No such file or directory`. Snakemake invokes
`sys.executable` by absolute path and does not quote it, so the space in "Biomedical Data Science"
splits the command. The pipeline works only because it is always run with `--use-conda` (which
activates an env and calls `python` from PATH). **`--use-conda` is not optional on this machine.**

Also re-confirmed the plan's warning empirically: a dry-run without `--rerun-triggers=mtime`
schedules **698 jobs** (340 ingest + 340 QC + both scVI trains), triggered by pre-existing code and
software-environment drift on `ds_ingest_dataset` / `ds_scrna_qc` / `ds_gene_symbol_map`. With
`--rerun-triggers=mtime --allowed-rules ds_compartment_audit` it is exactly 2 jobs. Note `--quiet`
has the same `nargs='+'` target-swallowing behaviour as `--allowed-rules` — pass targets first.

### Verification

- `flake8 workflow/scripts/compartment_audit.py` exit 0.
- Provenance JSON written per arm; artifacts non-empty.
- `results/pinned_reference_verification.json` still reads `pass: true, 37/37` — v1.3.0 untouched.

### Open / next

`compartment_audit.enforce` is **false** (baseline mode); flips to true in Phase 4. Phases 2-4
(annotate fixes, CNV rebuild on infercnvpy, D5/D6/D7 + conda-env enforcement + numba pin) are
approved to run through. **Phase 5 — the 6-10h full-arm and 4-7h capped-arm re-run — requires
researcher approval before it starts.**

---

## [2026-08-05] Phases 2-4 — Compartment integrity fix: annotation, CNV, masks, environment

**Phase:** Compartment integrity fix, Phases 2-4 of 6. Phase 5 (the re-run) awaits approval.

### Phase 2 — D1/D2/D3, `scrna_annotate.py`

- **D1.** `sc.tl.score_genes` assumes log-normalized input and does not check; the Census arms
  carry raw UMIs, so panels were compared on absolute count scales. `.X` is now normalized +
  log1p'd **in place** for scoring and restored to counts before writing — a normalized copy is
  16.5 GB on top of 16.5 GB and does not fit in 36 GB. Round-trip verified bit-exact (20k x 3k
  matrix incl. a 53,027-count outlier), guarded at runtime by a sum-drift check.
- **D2.** Argmax now runs over per-panel **z-scores**, so panels compete on relative enrichment.
  Annotate resolution 1.0 -> 2.0 on its own config key. Clusters whose top-vs-second margin is
  below a floor **and cross a compartment boundary** are labelled `ambiguous`; same-compartment
  ties (macrophage vs neutrophil) keep the top label. That refinement cut ambiguous from 17.7% to
  7.5% and recovered ~6,150 correctly-placed immune cells.
- **D3.** Astrocyte panel drops VIM for SLC1A2/SLC1A3/ALDH1L1/GJA1; seven lineages that had **no
  panel at all** (macrophage, mural, mast, B/plasma, NK, neutrophil, DC) now have one.
- **C2.** All of this lives in a new `annotation_markers:` block, NOT `nerve_cells.markers`, which
  `scrna_qc` also reads per-sample and which would have cascaded into the 11 h scVI train.
- **C3.** `immune_cells.source_label` (str) -> `source_labels` (list).

### Phase 3 — D4, `scrna_malignancy.py` + new `download_gene_positions` rule

Three stacked defects: gene order came from the Ensembl **accession counter**, not coordinates;
library-size normalization was missing entirely; and the "normal" reference was drawn from the same
annotation this rule should be independent of (degrading to T-cells-only in the capped arm).

No coordinates existed anywhere in the project (`dataset_gene_symbol_map.py` writes the literal
string `"unknown"`), so a new SHA-pinned Ensembl 113 GTF rule supplies them — 100% of both arms'
genes map. Windows now run **within** each contig.

That reached only 0.62/0.40. The remaining problem was the score itself: measured by genome-wide
spread, **normal neural cells scored HIGHER than malignant ones** (0.063 vs 0.057), because an
oligodendrocyte differs from an immune reference across whole chromosomes for reasons unrelated to
dosage. Replaced with a **chr7-gain minus chr10-loss contrast**, measured inside a single cell so
the cell-type baseline cancels. +7/-10 is a WHO 2021 IDH-wildtype GBM criterion and also fell out
of this cohort unprompted as the most-up and most-down contig (+0.043 / -0.038).

    AUC        0.859 -> 0.958
    precision  0.479 -> 0.934   (gate >= 0.85)  PASS
    recall     0.180 -> 0.894   (gate >= 0.80)  PASS

### Phase 4 — D5/D6/D7, compartment definition, environment

- **D5.** Exact canonical-label matching replaces substring matching, plus a hard fail when a
  configured `cell_types` entry matches no observed label. The old bidirectional substring test is
  what dropped 70,881 `opc` cells in silence.
- **D6.** `nerve_celltype_labels.py` writes `cell_type_marker_label`; **`cell_type` is never
  touched**. `nerve_scanvi.labels_key` follows, per-arm overridable.
- **D7.** Both remaining placeholder-on-empty paths hard-fail.
- Nerve compartment split into `glia` / `neuron` (`nerve_subcompartment`); neurons minted as one
  `nerve_neuron` LIANA group, glia per cluster. The `nerve_` prefix is kept so the concordance pair
  key stays comparable with the pinned v1.3.0 reference.

### Researcher decision — nerve compartment narrowed (2026-08-05)

Post-fix the nerve compartment reached 61.9% neural (from 11.1%), short of the 85% gate. The
residual was concentrated in two labels:

    oligodendrocyte    90.2% truly neural       excluded: opc        18.9%
    excitatory_neuron  84.2%                              astrocyte  10.0%
    neuron             61.5%

These are the AC-like/OPC-like malignant states sharing GFAP/PTPRZ1/SLC1A3 with normal glia; Census
reports **347 astrocytes in 1,006,344 cells (0.03%)** against the ~1.5% the panel calls. They are
CNV-quiet on +7/-10, so the caller cannot remove them, and tightening does not help — at 0 SD the
compartment is still only 78% neural and has lost 58% of its true neural cells.

**Decision: drop `astrocyte` and `opc` from `nerve_cells.cell_types`; gate relaxed 0.85 -> 0.80.**
Measured result 82.9% neural, retaining 69.2% of true neural cells.

> **KNOWN COST — must be reported, not treated as a finding.** OPCs are a real neural population
> (21,460 cells in the full arm) and relevant to neuron-glia-tumor crosstalk. They remain annotated
> and present in every artifact; they are excluded from THIS COMPARTMENT only. Their absence from
> the interaction tables is a masking decision.

### Defect 1 (conda env enforcement) — CLOSED, root cause was neither candidate approach

`markdowns/task_conda_env_enforcement.md` proposed reinstalling Snakemake outside the venv or
rebuilding `claude_science`. **Neither was needed.** Root cause: the venv is *active in the calling
shell*, so its exported `VIRTUAL_ENV` and PATH entry re-shadow `conda activate` inside every job
subshell. Unsetting `VIRTUAL_ENV` and stripping the venv from PATH before launching Snakemake makes
all four acceptance criteria pass. `CONDA_PREFIX` must be **left set** — unsetting it makes conda's
own deactivate-script lookup raise inside `posixpath.join`.

New `scripts/run_snakemake.sh` does this; new `conda_env_smoke_test` rule asserts it (14/14 pass:
`sys.executable` in the conda env, igraph 0.11.8, torch 2.12.0, liana 1.7.1, infercnvpy 0.4.3, all
resolving inside the env). **Phase 5 must be launched via the wrapper or enforcement silently
reverts.**

`numba==0.65.0` + `llvmlite==0.47.0` pinned in `scrna.yaml` (previously transitive and unpinned;
numba's threading layer caused both SIGSEGVs). This changed the env hash; the env was rebuilt and
re-verified now, deliberately, rather than mid-run.

### Tool failures worth recording

- **Snakemake deletes a failed job's declared outputs** — and the compartment audit is the one rule
  where that is exactly backwards, since a failing gate is when its cross-tabs matter most. Hit
  live: enabling `enforce` deleted the full arm's baseline tables. The audit now also writes an
  **undeclared sidecar** (`results/compartment_audit_snapshots/<arm>/`) that nothing in the DAG can
  remove. Baseline regenerated via a `--configfile` overlay with `enforce: false`.
- **`from __future__ import annotations` breaks under Snakemake `script:`** — Snakemake prepends its
  preamble, so the future import is no longer first and raises SyntaxError. Unnecessary on 3.12.
- **Ensembl's HTTPS mirror stalls mid-transfer** on the 64 MB GTF; the downloader now resumes by
  byte range and treats HTTP 416 as "already complete".

### Verification

- `flake8` exit 0 on all 13 edited scripts.
- `conda_env_smoke_test`: 14/14.
- Enforcing audit **correctly fails** on the pre-fix artifacts (7 gates), and the sidecar survives.
- `results/pinned_reference_verification.json` still `pass: true, 37/37`.

### Deviation from plan — infercnvpy not used

The plan preferred replacing the hand-rolled smoother with `infercnvpy`. It appeared absent from
the environment, which is precisely what the venv-shadowing defect looked like; it is in fact
installed (0.4.3). By the time that was clear the in-place rebuild already met the gate at
0.934/0.894, and swapping in a library unproven at 1M cells — with an `X_cnv` matrix in the ~10 GB
range — would risk memory on a 36 GB machine for no measured gain. Recorded, not quietly dropped.

### Next

**Phase 5 — the 6-10 h full-arm and 4-7 h capped-arm re-run — requires researcher approval.**
Launch via `scripts/run_snakemake.sh`. Re-entry at `ds_scrna_annotate`; the 11 h 13 m scVI train is
preserved because `X_scVI` is label-free.

---

## [2026-08-06] Phase 5 — full-arm re-run: all 10 gates pass; two corrections to the record

**Phase:** Compartment integrity fix, Phase 5. Full arm essentially complete; capped arm pending.

### Result — the fix works

`ds_compartment_audit` ran ENFORCING against the rebuilt full arm and passed **10/10 gates**:

| gate | before | after | required |
|---|---|---|---|
| nerve compartment neural | 11.1% | **95.39%** | >=80% |
| tumor compartment malignant | 47.9% | **92.96%** | >=85% |
| malignancy precision | 0.479 | **0.9296** | >=0.85 |
| malignancy recall | 0.180 | **0.8940** | >=0.80 |
| immune purity | 99.5% | **99.63%** | >=95% |
| max nerve-cluster endothelial | 92.9% | **0.0000** | <=20% |
| nerve compartment n | 377,343 | **37,945** | 25k-50k |
| neuron group n | 0 | **4,275** | >=1,500 |

`nerve_c24` — the 92.9%-endothelial cluster that made S1PR1 look nerve-side and started this whole
investigation — is gone as a failure mode.

scANVI retrained from scratch on the corrected 37,945-cell subset: classifier accuracy 0.9828.
Read with care: the compartment now has TWO anchoring labels (neuron, oligodendrocyte) where it had
five, so the task is easier and this is NOT comparable to the previous 0.8446.

### CORRECTION 1 — LIANA was never failing; I killed a healthy run

An earlier commit message (8ca071e) claims `rank_aggregate` was thrashing swap and "would never
have finished", and introduced a per-group cell cap on that basis. **That diagnosis was wrong.**

Both LIANA rules had ALREADY COMPLETED, uncapped, before I intervened:

    nerve_tumor_immune_interactions.csv   44,139 rows   13:52   provenance written
    nerve_tumor_interactions.csv          11,280 rows   13:55   provenance written

The provenance JSON is written last, so its presence is proof of completion. The three-way ran on
all 952,087 labelled cells x 24,135 genes in ~10 minutes.

What I actually measured (swap 12.5/13.3 GB, ~1s CPU per 20s wall) was `nerve_cell_heterogeneity`
— GSEA, which legitimately ran 1h50m at 11.5 GB and completed normally at 15:53. My process filter
matched a LIANA process that was already exiting, and I attributed another rule's memory pressure
to it. `liana.max_cells_per_group` is therefore set to **0 (disabled)**; the mechanism is retained
and documented for a cohort that genuinely does not fit, but nothing here needed it.

The other fixes from that episode stand on their own evidence and are unaffected: parent-sourcing
cured a real silent scale defect and cut peak RSS from 30+ GB to ~9 GB.

### CORRECTION 2 — D8, a silent scale defect (pre-existing, affected the original run too)

The interaction rules concatenated compartments on DIFFERENT scales and then normalized the result
as one object:

    malignancy_labeled.h5ad    raw counts
    nerve_cells.h5ad           already log1p  (normalize_total + log1p in *_cell_subset)
    immune_cells_labeled.h5ad  already log1p

So malignant cells were transformed ONCE and nerve/immune cells TWICE. Every cross-compartment
ligand-receptor comparison — in this run and in the original — was between differently-transformed
data. It stayed hidden because the combined X.max() is dominated by the raw malignant cells, so
every scale check saw "counts".

Fixed by taking EXPRESSION from the shared parent (all three compartments are subsets of it) and
using the compartment files only for LABELS, with assertions that the three are pairwise disjoint.
Verified: combined matrix now X.max() 58,860 -> 9.110 on one consistent scale, versus 8,291 -> 8.780
before.

### Defects found while getting here (all fixed)

- Four consumers still assumed the `nerve_c{N}` label shape and broke on the pooled `nerve_neuron`
  group. One of them (`_nerve_cluster`) would have SILENTLY blanked every neuron row rather than
  raising. Now routed through shared `fair_utils.nerve_group_key` / `nerve_group_sort_key`.
- `nerve_batch_qc` now computes purity over the same groups the interaction rules mint, so the
  neuron group has a row to join onto (23 -> 21 rows; three clusters were entirely neuronal).
- scANVI anchoring labels were hardcoded to the v1.x compartment and demanded >=1% OPC from a
  compartment that deliberately excludes OPC. Now derived from `nerve_cells.cell_types`.
- Undeclared crash-safe checkpoints survived a legitimate input change and were loaded against the
  wrong roster. Now fingerprinted (n_obs, n_vars, obs/var hashes, keys, label set, seed) and
  discarded on mismatch. Verified firing in production.
- `max_cells_per_group` was indented at 8 spaces inside a rule nested under `if not
  BASELINE_PINNED:` (12-space params), making it a rule keyword and breaking parsing for the WHOLE
  workflow — a 1-second failure that masqueraded as a run failure.

### Operational

- macOS has no `setsid(1)`, and a harness-tracked background task gets reaped mid-run. Long runs
  now launch via `scratchpad/daemonize.py` (double-fork + `os.setsid`), which orphans them to init
  so nothing upstream can kill them.
- `--quiet` has the same `nargs='+'` target-swallowing trap as `--allowed-rules`.

### Next

Full arm: 5 light jobs remain (cluster annotations, batch_qc_v2, annotate_cluster_qc, concordance,
notebook). Capped arm: 18 jobs from `ds_scrna_annotate`, not yet started.

### [2026-08-06 16:25] Full arm COMPLETE — rc=0, all verification met

`gbm_cellxgene_56c4912d_full` finished all 14 jobs. Goal-backward verification:

1. **ZERO nerve-side S1PR1 rows** in `nerve_tumor_immune_top_pairs_with_qc.csv` (zero S1PR1 rows at
   all). The finding that started this investigation — SPP1->S1PR1 called immune->nerve and
   tumor->nerve — was an artifact of `nerve_c24` being 92.9% endothelial. That cluster's
   contamination is now 0.0000 and the rows are gone.
2. **10/10 compartment-integrity gates PASS** (audit runs enforcing, so this is self-checking):
   nerve 95.39% neural, tumor 92.96% malignant, malignancy precision 0.9296 / recall 0.8940,
   immune 99.63% pure, max nerve-cluster endothelial 0.0000, nerve n=37,945, neuron group n=4,275.
3. `results/pinned_reference_verification.json` still **pass: true, 37/37** — the v1.3.0 reference
   was never touched.
4. scANVI-v2 retrained on the corrected subset; `nerve_cells_v2.h5ad` present; `cell_type`
   preserved throughout (defect D6 closed).

**Concordance against the pinned v1.3.0 reference moved, as predicted:**

    reference significant pairs   3,368
    dataset significant pairs     2,283
    shared                        1,724
    Jaccard overlap               0.439
    Spearman rho (shared)         0.6182

This is NOT like-for-like and must not be read as a replication failure. The reference's own nerve
compartment was built by the same uncorrected logic this work removed, it carries no author
annotation, and it was scored against differently-normalized data (defect D8). Its status is
*unknown*, not *cleared*. Demote it from headline until the reference is itself re-derived.

**Interaction tables:** 44,139 rows, all nerve-involving; 2,609 rows carry the pooled
`nerve_neuron` group; 15,504 rows carry `batch_qc_pass=False` (dominated by the neuron group's
donor-dominance failure, recorded not exempted).

Capped arm `gbm_cellxgene_56c4912d` started automatically at 16:25:13 (18 jobs from
`ds_scrna_annotate`).

Defects fixed during this final stretch, both the same species — a consumer that parsed the nerve
group label by hand instead of via the shared helper:
- `nerve_batch_qc` emitted only the interaction keying, dropping the entirely-neuronal clusters
  (11, 14, 19) that the marker/enrichment tables key on. Purity now emits per-cluster rows AND a
  pooled neuron row (24 rows).
- `annotate_cluster_qc` derives the nerve key in TWO places; only one had been routed through
  `nerve_group_key`. The three-way path still used a bare `removeprefix("nerve_c")`.

### [2026-08-06 18:10] PHASE 5 COMPLETE — both arms, rc=0

`gbm_cellxgene_56c4912d` finished 18/18 at 18:10:07. Both Census arms are now rebuilt on the
corrected pipeline and both pass every gate.

| | FULL (1.0M) | CAPPED (615k) |
|---|---|---|
| compartment gates | **10/10 PASS** | **10/10 PASS** |
| nerve neural fraction | 95.39% (was 11.1%) | 94.76% (was 11.6%) |
| tumor malignant fraction | 92.96% (was 47.9%) | 93.47% (was 45.6%) |
| malignancy precision / recall | 0.930 / 0.894 | 0.935 / 0.875 |
| immune purity | 99.63% | 99.63% |
| max nerve-cluster endothelial | 0.0000 (was 0.929) | 0.0000 (was 0.453) |
| nerve compartment n | 37,945 | 28,936 |
| Census neuroglial truth | 37,561 | 27,860 |
| neuron group n | 4,275 | 4,240 |
| **nerve-side S1PR1 rows** | **0** | **0** |

The two arms were corrected independently and each landed within ~1-4% of its OWN Census
neuroglial count. Two cohorts of different sequencing depth converging on their own ground truth
is the strongest available evidence that the fix is real rather than tuned to one dataset. The
depth comparison — the reason both arms exist — is preserved.

Concordance against the pinned v1.3.0 reference:

    FULL    jaccard 0.4390   spearman 0.6182   shared 1,724
    CAPPED  jaccard 0.4739   spearman 0.6620   shared 1,914

Both moved down from the pre-fix figures. **This is not a replication failure and must not be
reported as one.** The reference's own nerve compartment was built by the uncorrected logic this
work removed, it carries no author annotation, and it was scored against differently-normalized
data (defect D8). Its status is *unknown*, not *cleared*. Notably the two corrected arms agree with
each other far better than either agrees with the reference — consistent with the reference being
the outlier.

Guards still green after the full re-run:
  pinned v1.3.0 reference   pass 37/37 unchanged
  conda env enforcement     pass 14/14

### Open items for the researcher

1. **Re-derive the 40-axis shortlist** (`results/tables/nerve_crosstalk_lead_targets.csv`). It
   predates the fix. Note S1PR1, CXCR4 and LRP1 were never on it — they appeared only in raw LIANA
   tables, and whatever shortlist named them was produced outside this repo.
2. **The neuron group fails batch purity** (dominant_sample_fraction 0.5032; one donor supplies
   half of all ~4,275 neurons across 7 contributing samples). Recorded, NOT exempted. Any
   neuron-side axis is substantially one patient's biology.
3. **The immune compartment doubled** (1.80x full, 2.00x capped) because T/NK/B cells are in it for
   the first time. Every pre-fix "immune" result was myeloid-only, so immune-side findings change
   in character, not just magnitude.
4. **OPC, astrocyte, generic neuron and ependymal are excluded from the nerve compartment** by
   researcher decision. They remain annotated and present in every artifact — their absence from
   the interaction tables is a masking decision, not a biological finding.
5. **The v1.3.0 reference should be re-derived** on the corrected pipeline before concordance is
   used as evidence either way.

---

### [2026-08-07] | Phase: Post-compartment-fix notebook audit | Status: COMPLETE

**Action:** Researcher asked whether `notebooks/05_census_nerve_immune_explorer.py` references the
newly generated data after the 2026-08-06 re-run. Audited it against the artifacts on disk.

**Outcome:** Mechanically yes, semantically no. Every path the notebook resolves points at a table
rebuilt on 2026-08-06, but the logic reading those tables was written against the pre-fix
compartments. Four defects, one of them silent:

1. **`nerve_neuron` join failure (correctness).** `_prepare` derived its join key with
   `str.replace("nerve_c", "")`. Neurons are carried as one pooled LIANA group, `nerve_neuron`,
   which contains no `nerve_c` substring, so the id passed through unchanged and matched no
   annotation row. **2,679 of 48,581 nerve-side rows** lost their cell type, dropped out of the
   Panel D interface matrix entirely, and rendered as a bare id in Panel C — no error raised. This
   is the exact failure `fair_utils.nerve_group_key` exists to prevent and that
   `annotate_cluster_qc.py` already guards against; the notebook was the one consumer that never
   adopted the helper. Worse than its row count: the pooled neuron group is precisely the
   donor-dominated population that needs to be *visible and caveated*, not invisible.
2. **Panel A's compartment mapping was a hardcoded copy of the old definitions.** It still listed
   astrocyte/opc/generic-neuron/ependymal as nerve and `microglia` alone as immune. Measured
   against the live `annotation_summary.csv` for the capped arm, it reported nerve **194,080** and
   immune **125,671** where config gives **36,712** and **343,280** — nerve overstated 5x, immune
   understated 3x — and named seven labels as excluded that are now inside the immune compartment.
3. **Three false narrative claims.** "resolves only 5 cell types … no opc/ependymal/tumor_gbm"
   (there are now 18 labels including all three); "`mean_confidence` 5.96-15.29" (live range
   0.95-28.11); and a pointer to `blocker_census_annotation_scoring.md` as an *unfixed* defect —
   that is defect D1 and `scrna_annotate.py:102-135` now normalizes to log1p before `score_genes`
   and reverts after. Footer's "11 of 35 nerve clusters fail batch QC" is now 9 of 26.
4. **Panel F framed the v1.3.0 concordance as an agreement result** ("both rose once the scales
   matched"), which is exactly the reading the 2026-08-06 entry above forbids.

**Fixes applied** (researcher-directed scope: correctness + stale prose; remove Panel E; §2.5
option (a) for concordance):

- Import `nerve_group_key` / `nerve_group_sort_key` from `fair_utils` instead of parsing ids by
  hand; synthesize the pooled-neuron annotation row the Leiden-only annotation rule cannot emit;
  add a `[FAIR-ALERT]` guard that raises if any nerve group lands without a label, mirroring the
  pipeline-side check.
- Panel A now reads `nerve_cells.cell_types` and `immune_cells.source_labels` from config at run
  time. Its callout states the §2.4 masking decision with the measured per-label purity, and that
  the tumor compartment is CNV-derived so it cuts across every label in the chart.
- Panel B sorts nerve groups numerically-then-named and surfaces the neuron group's donor
  dominance inline.
- **Panel E removed** — it scored the current tables against a 2026-07-15 shortlist derived when
  "nerve" was 59% malignant. Replaced with a stub recording why, so it is not silently re-added.
- Panel F demoted to a diagnostic with the three reasons the reference is not a valid comparator.

**Artifacts:** `notebooks/05_census_nerve_immune_explorer.py`,
`notebooks/__marimo__/session/05_census_nerve_immune_explorer.py.json`,
`.claude/plans/plan_notebook05_post_compartment_fix.md`,
`markdowns/post_compartment_fix_next_steps.md` (researcher's review doc, now tracked — the notebook
cites it in six places).

**Verification:** `flake8` clean. Headless `marimo export html` run on **both** arms, exit 0, no
error cells. Join completeness confirmed — the only unmatched nerve key on either arm is `neuron`,
the group the notebook now synthesizes. Panel A renders 36,712 nerve / 343,280 immune on the capped
arm. Panel D regains 90 (capped) / 97 (full) neuron rows at the default thresholds that it
previously discarded. Purity callout renders 9 of 26 groups failing (capped) and 9 of 24 (full).

**Tool Versions:** marimo 0.23.1 (matches the `workflow/envs/notebooks.yaml` pin), pandas 2.x.

**Open Issues:**
- `markdowns/blocker_census_annotation_scoring.md` still reads `Status: OPEN` although its defect
  (D1) is fixed. One-line docs correction, not done here.
- Deferred by researcher decision, all still open: a compartment-audit panel reading
  `compartment_audit_gates.csv` / `malignancy_confusion.csv` / `nerve_compartment_cluster_audit.csv`
  (no notebook reads them today); the §3.1 direct S1PR1/CXCR4/LRP1 x `cell_type` cross-tab; §2.1
  shortlist re-derivation; a cross-arm concordance artifact (§2.5 option (c)).
- `notebooks/04_tme_nerve_immune_explorer.py` was NOT audited. Its cohort was not rebuilt, so it
  does not carry defects 1-2, but its concordance framing is worth the same pass.

**FAIR Notes:** No artifact changed — this run modified no table and required no pipeline
re-execution. The notebook's export sidecar now hashes 11 inputs rather than 12, the dropped one
being the withdrawn pre-fix shortlist.

---

### [2026-08-07] | Phase: Post-compartment-fix notebook audit (04) | Status: COMPLETE

**Action:** Same pass over `notebooks/04_tme_nerve_immune_explorer.py` as the entry above did for
05. This notebook reads the **pinned v1.3.0 reference** tables at `results/tables/` (root), which
the compartment fix did **not** rebuild — so the audit question is different: not "does it read the
new data" but "does it say what its own data is worth".

**Outcome:** Two of notebook 05's four defects do not apply; the two that remain are the ones that
matter most here.

*Does not apply, verified rather than assumed:*

1. **No `nerve_neuron` join failure.** The reference cohort's groups are uniformly `nerve_c{N}`;
   the notebook's prefix slice is exactly equivalent to `fair_utils.nerve_group_key` on this table
   and leaves **0** rows unmatched. Switched to the helper anyway, plus the same `[FAIR-ALERT]`
   orphan guard, so the notebook cannot silently lose a group if it is ever pointed at a rebuilt
   cohort.
2. **Panel A's hardcoded compartment map is correct here — and must stay hardcoded.** This is the
   deliberate *opposite* of the fix applied to notebook 05. Config has moved past these pinned
   artifacts: today's `nerve_cells.cell_types` selects **18,405** labelled cells against the
   **106,603** actually modelled in the pinned LR run, so reading config at run time would
   understate the nerve compartment ~6x. The literals are now commented as the definition in force
   at freeze time, with an explicit DO-NOT-"FIX" note.

*Applies, and was the real gap:*

3. **Nothing on the page said this cohort is unaudited.** Its nerve compartment was built by the
   logic the fix removed, and unlike both Census arms it carries **no author annotation**, so there
   is no oracle to check it against. Where that logic could be measured it came out 59% malignant /
   11% neural; here the contamination is **unknown, not measured, and not cleared**. Added as a
   top-of-notebook banner and as footer limit 1, which now dominates the others.
4. **Stale prose.** The header and footer still excluded the Census cohort on the grounds of the
   raw-counts LIANA defect — RESOLVED 2026-07-26, and that cohort has had its own notebook since.
   Panel A still described the OPC drop as an open naming miss (defect D5, fixed; OPCs are now
   excluded by explicit decision) and described `t_cell` as sitting outside the immune subset
   without noting that this makes every immune result on the page a **myeloid** result.

**Shortlist provenance settled, and a correction to yesterday's work.** Notebook 05's Panel E stub
asserted the curated shortlist was "derived from tables in which the nerve compartment was 59%
malignant and 27% myeloid". That figure is the Census full arm's, and the shortlist is not from
there: **34 of its 40 `best_mag` values reproduce against the reference interaction table to 1e-6,
and 0 against either Census arm.** It is reference-derived. Corrected the wording in 05 to the
defensible claim — the shortlist comes from a cohort that has never been audited, so its basis is
unknown rather than measured-bad. Also found the shortlist's `n_rows` no longer matches the live
reference table (median gap 186 rows across 40 axes), because it is dated 2026-07-15 and the
`_with_qc` table was regenerated 2026-07-21; Panel E of 04 now says so.

**Panel E kept in 04, removed in 05 — deliberately asymmetric.** In 04 it traces an axis back to
the rows that support it within one consistent cohort, which is a like-for-like trace. In 05 it
scored *corrected* tables against the same shortlist, which reads as a replication test and is not
one.

**Artifacts:** `notebooks/04_tme_nerve_immune_explorer.py`,
`notebooks/05_census_nerve_immune_explorer.py` (Panel E wording correction),
both `notebooks/__marimo__/session/*.json`.

**Verification:** `flake8` clean on both. Headless `marimo export html` on 04, exit 0, no error
cells; rendered Panel A shows 106,603 modelled nerve against 136,587 labelled and a 2,786-cell
coverage gap, all three matching hand computation against `annotation_summary.csv` and the pinned
provenance. Notebook 05 re-exported after its edit, exit 0.

**Open Issues:**
- `duckdb` is imported by notebook 04 and never used. Left alone — removing it churns the cell
  return signature for no functional gain.
- `markdowns/blocker_census_annotation_scoring.md` still reads `Status: OPEN` though defect D1 is
  fixed (carried over from the previous entry).
- Notebook 03 (`03_nerve_tumor_immune_explorer.py`) has not had this pass. It is the pure LR view
  over the same pinned reference tables, so it very likely carries defect 3 above.

**FAIR Notes:** No artifact changed; no pipeline re-execution. The pinned v1.3.0 tables were read
only.

---

### [2026-08-07] | Phase: Post-fix lead-axes panel (notebook 05) | Status: COMPLETE

**Action:** Researcher generated `results/tables/nerve_immune_lead_axes_postfix.csv` (183 rows x 15
cols) — the §2.1 re-derivation — and asked for a notebook 05 panel pointing at it. It lands in the
slot where Panel E was removed on 2026-08-07.

**Provenance established before writing anything.** This shortlist is Census-derived across **both**
arms: all 183 axes are present in both Census interaction tables and **0** of their magnitudes match
the pinned v1.3.0 reference. That is the exact inverse of the withdrawn 2026-07-15 shortlist
(34/40 matching the reference, 0 matching Census), and it is what makes this one usable where its
predecessor was not.

**Structural facts that drove the design:**

- **The primary key is `(axis, nerve_side)`, not `axis`** — 144 distinct axes over 183 rows, 39 of
  them ranked separately on both nerve sides (`APP|CD74` is rank 1 on each), with ranks restarting
  per side. A panel keyed on `axis` would have silently collapsed those pairs. Guarded by a
  `[FAIR-ALERT]` raise on duplicate keys.
- **`n_rows` is a FULL-ARM count with no capped counterpart.** Measured under the reconstructed
  filter it agrees with the full arm on 118/183 axes but with the capped arm on only **47/183**.
  Comparing it while the capped arm is active would have rendered ~3/4 of the table as
  "disagreement" that was really just the wrong arm. It is now shown for reference and compared
  only on the full arm. **This was caught by running both arms, not by reading the file.**
- **The file's two magnitude columns were not built the same way.** Under the reconstructed rule
  (QC-passing, `magnitude_rank <= 0.05`, restricted to each axis's own interfaces),
  `capped_best_mag` reproduces **exactly, 128/128**, while `best_mag` reaches only **128/167**.
  Worth resolving at the source.
- **55 axes have a null `capped_best_mag` while still being present in the capped table** — they
  cleared the bar in the full arm and not the capped one. Surfaced as a `cross_arm` column rather
  than as missing data; this is the §2.5(c) cross-arm signal, and it is better evidence than the
  v1.3.0 concordance in Panel F.

**Panel design** (researcher-directed): curated values rendered beside live recomputed ones with
disagreement flagged; arm-aware via an explicit `dataset` -> column map (never a
`endswith("_full")` heuristic, which would mislabel a future arm); the Panel E stub replaced,
keeping a compressed note on why the old shortlist was withdrawn. Local widgets (tier, nerve side,
disagreements-only) kept separate from the Panel C/D filter block so the two cannot interfere.

**Two columns deliberately not rendered.** `withdrawn` is empty in all 183 rows, and `min_pval` is
`0.0` in all 183 (the permutation floor at `n_perms=1000`). A column that is entirely null or
entirely constant carries no information and invites a reader to infer meaning from it.

**Artifacts:** `notebooks/05_census_nerve_immune_explorer.py`,
`results/tables/nerve_immune_lead_axes_postfix.csv` (NOT tracked — see below),
`notebooks/__marimo__/session/05_census_nerve_immune_explorer.py.json`.

**Verification:** `flake8` clean. Headless `marimo export html` on both arms, exit 0, no error
cells. Rendered counters match the offline audit exactly — FULL: 183/183 shown, magnitude 128/167,
`n_rows` 118/183, 0 capped-only; CAPPED: magnitude 128/128, `n_rows` comparison correctly suppressed,
55 full-only. Both arms report 55 pooled-neuron axes and CXCR4 3 / LRP1 9 / S1PR1 0. Export sidecar
now hashes 12 inputs. The disagreements-only filter branch was exercised separately (the checkbox
defaults off, so the export never reached it) and is index-safe under a tier subset.

**Open Issues:**
- **`withdrawn` is empty in all 183 rows** although `tier_v2` marks 14 axes
  `1b_approved_withdrawn_only` and `agents_flagged` carries `[WITHDRAWN]` tags inline. Looks like a
  bug in the generating script. Not fixed here — that means regenerating the CSV.
- **The generating script is not in this repository.** Until it is, the panel's live recomputation
  is the only reproducible statement of what these axes are, and the curated/live gap cannot be
  closed. This is the same condition §2.1 flagged about the previous shortlist.
- Notebook 03 still has not had the post-compartment-fix pass (carried from the previous entry).

**FAIR Notes:** [FAIR-ALERT] No artifact changed and no pipeline re-execution, but the new input is
a FAIR gap on three counts at once: **no producing Snakemake rule**, **no generating script in the
repository**, and **`results/` is gitignored** (0 results files are tracked — this is the project's
standing convention, not an oversight, because results are normally derived artifacts). It was
therefore NOT force-added; doing so would override a deliberate .gitignore for a file that policy
says should be regenerable. The consequence is that this input can be neither rebuilt nor restored
from git, so a clean checkout cannot render Panel E. The notebook's missing-artifact callout
special-cases it rather than sending the reader after a Snakemake rule that does not exist, and
Panel E carries the same alert. **A rule under `workflow/rules/` that emits this table would close
the reproducibility gap, the versioning gap and the curated/live gap together** — recommended as the
next step, and it is also what CLAUDE.md's "never run one-off scripts for analysis steps" requires.

---

### [2026-08-19] | Phase: Lead-axes generator brought into the pipeline | Status: COMPLETE

**Action:** Close the reproducibility gap flagged in the 2026-08-07 entry above — give
`results/tables/nerve_immune_lead_axes_postfix.csv` a producing Snakemake rule. The generator
arrived as four Claude Science REPL transcripts in `claude_science/code_export/` (gitignored, since
`claude_science/` is a venv).

**Outcome:** New `rule nerve_immune_lead_axes` reproduces the 2026-08-07 artifact **byte-for-byte**
— sha256 `f30d33ed71c8038f8cc989a3d29d48ade873ad2837e47570d1d067bbe00f6286` before and after. That
exact match is the proof the reimplementation is faithful; it was the acceptance gate for the work.

Root cause of the `best_mag` discrepancy Panel E reports, now understood and documented: the CSV is
a `pd.concat` of **two halves built by different rules**, not one table.

| | oligodendrocyte (128 rows) | neuron (55 rows) |
|---|---|---|
| source arm | full | capped |
| batch QC | required | **not applied** |
| clusters | all | `nerve_neuron` only |
| sort key | `['best_mag','min_pval']` | `['best_mag']` only |
| dense rank | `rank_full` 1..128 | `rank_capped` 1..55 |
| `capped_best_mag` | populated | **column never created → NaN** |

So `best_mag` means "full arm, QC-passing" for 128 rows and "capped arm, no QC, neuron-only" for 55,
and `n_rows` is arm-inconsistent the same way. That is why a single-arm live recomputation reproduces
`capped_best_mag` 128/128 but `best_mag` only 128/167. **Reproduced deliberately, not fixed** — the
fix changes the science and is scoped to a separate commit so it can be reviewed and reverted alone.

**The four drug columns cannot be recomputed and are now vendored.** `tier_v2`, `agents_flagged`,
`glioma_trials`, `withdrawn` came from ChEMBL and ClinicalTrials.gov via a Claude Science MCP
connector. None of the generator's `handoff/*.json` response caches survived, and one input is
permanently lost: a prior-session artifact addressed only as UUID
`34d720d0-5c16-46bb-92e6-d8d6b6a02ec3`, which supplied fallback agents for genes with no ChEMBL
mechanism record. Vendoring is lossless because all four columns are a pure function of the `axis`
string — verified: zero disagreement across the 39 axes appearing on both nerve sides — so a 144-row
axis-keyed snapshot reproduces all 183 rows.

**Artifacts:**
- `workflow/rules/leads.smk` (new), `workflow/scripts/nerve_immune_lead_axes.py` (new)
- `reference/drug_annotation/` (new tracked dir — `data/` and `results/` are both gitignored, so
  neither could host a committed input): pinned snapshot `..._2026-08-07.csv` (144 axes, sha256
  `a80b1de3c98bd0ed5191df2101315c05f4fd8e2e72d6b3b6021fd4b615781acc`), `curated_agent_map_2026-08-07.json`,
  `MANIFEST.json`
- `config/config.yaml`: new `lead_axes:` block; `Snakefile`: include + `rule all` target
- `provenance/nerve_immune_lead_axes_provenance.json` (new; FAIR validation now 68 records)

**Tool Versions:** pandas 2.3.3 (`workflow/envs/scrna.yaml`), Snakemake 9.x

**Verification:** dry run schedules exactly 1 job (no upstream LIANA rebuild); byte-identity gate
passed; 183 rows / 15 columns; 128 + 55 split; `capped_best_mag` non-null on exactly the oligo rows;
oligo `rank_full` == 1..128; neuron `rank_capped` == 1..55; tier counts match §8 of
`markdowns/GBM_TME_Crosstalk_Analysis.md` (46/14/50/38/35); `fair_validate_metadata` passes;
`flake8` clean. Notebook 05 was **not** re-rendered — its input is byte-identical, so Panel E cannot
have changed.

**Open Issues:**
1. **`withdrawn` is empty in all 183 rows — a real bug, carried forward on purpose.** The generator
   sourced it from `DRUG.withdrawn_agents`, populated by the bulk chembl_id lookup, which per the
   export's own README does not carry withdrawal status; the corrected set came from a later by-name
   sweep and was never written back. Recoverable — the data survives inline as `[WITHDRAWN]` tags in
   `agents_flagged` and as the 14 axes tiered `1b_approved_withdrawn_only`.
2. Arm-inconsistent `best_mag` / `n_rows` (see table above) — fix alongside notebook 05 Panel E and
   §8 of the crosstalk markdown, whose 128/128 and 128/167 counts go stale the moment it changes.
3. `rule refresh_drug_annotation` (opt-in ChEMBL + ClinicalTrials.gov REST re-query, writing a NEW
   dated snapshot) is **specified but not yet implemented**. Until it exists the drug annotation is
   frozen, not refreshable.

**FAIR Notes:** The provenance and versioning gaps are closed — the table is now rebuildable offline
from committed inputs, with a provenance sidecar. The *recomputability* gap is NOT closed and cannot
be: `reference/drug_annotation/MANIFEST.json` states this plainly rather than implying otherwise,
including that the ChEMBL release behind the snapshot was never recorded and that glioma-trial
coverage is bounded by a hand-curated 20-agent list (absence of a trial ≠ no trial exists).

---

### [2026-08-19] | Phase: Lead-axes open items (1/3) — `withdrawn` recovered | Status: COMPLETE

**Action:** Fix open item 1 from the entry above — `withdrawn` was empty in all 183 rows of
`nerve_immune_lead_axes_postfix.csv`, a generator bug rather than a property of the data.

**Outcome:** New `rule derive_withdrawn_agents` recovers the column offline from the inline
`[WITHDRAWN]` tags that `agents_flagged` already carried, writing a new dated snapshot
`reference/drug_annotation/nerve_immune_axis_drug_annotation_2026-08-19.csv` (sha256
`eed0995c…`). The 2026-08-07 file is kept untouched as the audit record — rewriting a dated
artifact would falsify what it contained on that date.

**13 of 144 axes** populated, covering **3 agents: BENZIODARONE, PRENYLAMINE, PROBUCOL** —
matching the generator's own report text verbatim. In the 183-row table that is 20 rows (some
axes appear on both nerve sides), 14 of them tier `1b_approved_withdrawn_only`. Zero 1b axes
were left empty, so tier and agent tags now corroborate each other instead of contradicting.

Three of the 13 are tiered `1_approved_available`, not 1b — `CALM1|INSR`, `CALM1|PDE1A`,
`CALM1|PDE1C`. That is correct, not a leak: axis tier is the *best* tier among constituent
genes, so a withdrawn CALM1 agent does not demote an axis whose partner gene has a clean
approved one. These are exactly the "CALM1 axes are a trap" cases the generator's report calls
out — the agent a reader would reach for is withdrawn even though the axis reads tier 1.

**Verification:** the rebuilt table differs from the previous build in **exactly one column**
(`withdrawn`) — every other column compared equal, and the rule carries its own guards that
raise `[FAIR-ALERT]` if any passthrough column moves, if the row count changes, if the parsed
count disagrees with the tag count, or if a 1b axis ends up empty. Table sha256 moves
`f30d33ed…` → `62e66e0d…`. `flake8` clean.

**Artifacts:** `workflow/scripts/derive_withdrawn_agents.py` (new), `rule derive_withdrawn_agents`
in `workflow/rules/leads.smk`, the 2026-08-19 snapshot,
`provenance/derive_withdrawn_agents_provenance.json`, `config/config.yaml` (repointed),
`reference/drug_annotation/MANIFEST.json` (`corrections` section added; the `withdrawn` defect
moved to RESOLVED).

**Tool Versions:** pandas 2.3.3 (`workflow/envs/scrna.yaml`)

**Open Issues:** items 2 (arm-inconsistent `best_mag`/`n_rows`) and 3 (`refresh_drug_annotation`)
still open — next two commits. Notebook 05 Panel E still says `withdrawn` "is empty in all 183
rows … looks like a bug in the generating script"; that text goes stale with this commit and is
corrected in the item-2 commit, which is where the notebook is touched.

**FAIR Notes:** The correction is itself a provenance-stamped Snakemake rule, not a manual edit,
so the fix is as reproducible as the artifact it fixes. No network access — this recovers
misplaced data, it does not introduce new data, so nothing here depends on ChEMBL's current state.

---

### [2026-08-19] | Phase: Lead-axes open items (2/3) — arm-explicit magnitudes | Status: COMPLETE

**Action:** Fix open item 2 — `best_mag` and `n_rows` were arm-inconsistent: full-arm
QC-passing for the 128 oligodendrocyte rows, capped-arm non-QC pooled-neuron for the 55
neuron rows, with nothing in the row saying which.

**Outcome:** Both replaced by four arm-labelled columns — `full_best_mag`,
`capped_best_mag`, `full_n_rows`, `capped_n_rows` — **all 183/183 populated**. No NaNs,
because each half is already a cross-arm intersection, so every axis has a counterpart in
the other arm. The 55 `capped_best_mag` NaNs disappear as a side effect. Table is now
183×16; sha256 `62e66e0d…` → `9869b698…`.

Verified as a pure relabelling, not a recomputation: the oligodendrocyte half's
`full_best_mag`/`full_n_rows` equal the old `best_mag`/`n_rows` exactly, and the neuron
half's `capped_best_mag`/`capped_n_rows` equal its old `best_mag`/`n_rows` exactly. Ranks
still dense 1..128 and 1..55; tier counts unchanged at 46/14/50/38/35.

**What this does NOT fix, stated in the script docstring, the rule docstring, Panel E and
§8:** the *selection rule* still differs between halves — oligodendrocyte is QC-passing
all-cluster, neuron is no-QC pooled-neuron-only. That is what defines the two halves and
is labelled by `compartment_side`. The arm is now explicit; the filter is not. Comparing
a `full_best_mag` across the two sides is still comparing two different filters.

**Panel E reframed from caveat to regression check.** Two corrections were needed to make
the oracle actually correct, both found by testing rather than assumed:
1. It applied the oligodendrocyte rule to every row. The neuron half is built without the
   QC filter and restricted to `nerve_neuron`, so all 55 neuron rows would have read as
   disagreements. The oracle now branches on `compartment_side`, as the rule does.
2. It restricted matches to each axis's own `interfaces`. The rule *derives* that column
   from the matched rows rather than filtering on it, and each half's value is computed on
   one arm while the recomputation runs on both — so the restriction silently undercut
   `n_rows` by 5/183 on the capped arm and 1/183 on the full one. Removed.

With both fixed the recomputation agrees **183/183 on magnitude and 183/183 on row counts
in both arms** — confirmed in the rendered HTML, not just in a harness. Panel E now turns
the callout red and prints a REGRESSION banner if either falls below 183.

**Also fixed:** `rule ds_census_nerve_immune_notebook` declared the *withdrawn predecessor*
`nerve_crosstalk_lead_targets.csv` as an input but never `nerve_immune_lead_axes_postfix.csv`,
the file Panel E actually reads — so rebuilding the shortlist did not re-render the notebook
and a stale render could hide a regression. Now declared.

**Artifacts:** `workflow/scripts/nerve_immune_lead_axes.py`, `workflow/rules/notebooks.smk`,
`notebooks/05_census_nerve_immune_explorer.py` (Panel E compute + view cells, the
`_unbuildable`/FAIR-ALERT machinery removed from the loader),
`markdowns/GBM_TME_Crosstalk_Analysis.md` §8 (now tracked in git for the first time).

**Verification:** 183 rows / 16 cols; four arm columns 183/183 non-null; `best_mag` and
`n_rows` absent; relabelling equalities all True; notebook rendered for both arms with zero
tracebacks and 183/183 agreement; `flake8` clean.

**Open Issues:** §8 of the markdown carried a pre-existing arithmetic error — "125 on the
oligodendrocyte side and 55 on the pooled-neuron side" sums to 180, not 183. The real split
is 128 + 55. Corrected. Item 3 (`refresh_drug_annotation`) remains open.

**FAIR Notes:** The notebook's oracle deliberately duplicates the rule's selection logic —
an oracle sharing an implementation with the thing it checks cannot catch anything. The cost
is that the two must be edited together; that is stated in the cell so the next editor knows.

---

### [2026-08-19] | Phase: Lead-axes open items (3/3) — refresh rule | Status: COMPLETE

**Action:** Fix open item 3 — implement `rule refresh_drug_annotation` so the drug annotation
can be re-derived from public APIs instead of only being pinned.

**Outcome:** Implemented and, contrary to plan, **fully exercised**. ChEMBL was returning HTTP
500s and timeouts throughout planning, so the approved plan was to ship the ChEMBL half
unexercised; EBI recovered mid-session (now serving **ChEMBL_37**, released 2026-05-01), so the
whole rule was run for real. That mattered — it caught two bugs that code review did not, both
of which would have shipped as confidently wrong output.

**Bug 1 — `/target/search` is the wrong endpoint.** It rejects `target_type`/`organism` filters
with HTTP 400, and 400s outright on short symbols (`q=C3` fails even bare). Replaced with an
exact gene-symbol match on the non-search `/target` endpoint via
`target_components__target_component_synonyms__component_synonym__iexact`. Also strictly more
precise: exactly one human target per gene instead of a fuzzy ranked list, so a near-miss cannot
be mistaken for a hit. `NLGN1` correctly resolves to nothing at any target type, matching the
original generator's recorded finding.

**Bug 2 — the dangerous one.** ChEMBL_37 serialises `max_phase` as the **string** `'4.0'`, so
`phase == 4` was never true and every approved drug fell through to `3_chembl_target_no_agent`.
The run exited 0 and wrote a clean-looking 144-axis snapshot asserting **zero approved or
clinical agents across all 156 genes** — including CXCR4, whose plerixafor record had been
verified by hand minutes earlier. Not a crash: a plausible table making a false scientific claim.
Fixed by coercion, plus `1 <= phase <= 3` for the clinical band (which also catches ChEMBL's
fractional phases). A **plausibility guard** now aborts if no axis reaches an approved or
clinical tier — drift is expected, total collapse of the top tiers is a parse failure.

**Fail-loud semantics verified under real failure**, not simulated: the first run aborted on gene
`C3` after 4 backed-off retries and wrote no output file. A partial snapshot would have demoted
every unresolved gene to `4_no_chembl_target`.

**Drift vs the pinned snapshot** (same 144 axes; `tier_v2` unchanged on **120/144, 83%**):

| tier | pinned | refreshed | Δ |
|---|---:|---:|---:|
| `1_approved_available` | 38 | 29 | −9 |
| `1b_approved_withdrawn_only` | 10 | 12 | +2 |
| `2_clinical` | 40 | 46 | +6 |
| `3_chembl_target_no_agent` | 33 | 46 | +13 |
| `4_no_chembl_target` | 23 | 11 | −12 |

The 12 axes leaving `4_no_chembl_target` are an improvement — exact symbol resolution finds
targets the original full-text search missed. The 9 leaving `1_approved_available` are the
expected cost of the permanently lost prior-pass artifact, which supplied agents for genes with
no ChEMBL mechanism record; without it those genes tier on ChEMBL evidence alone.
`glioma_trials` matches exactly on 130/144, the 14 differences being trials registered since
2026-08-07.

**Cross-validation worth recording:** the `withdrawn` column re-queried straight from ChEMBL's
`withdrawn_flag` matches the column derived in commit 1/3 by parsing inline `[WITHDRAWN]` tags on
**13/13 axes, character for character**. Two fully independent derivations agreeing is the
strongest evidence available that the tag-parsing recovery was correct.

**NOT ADOPTED.** `config.lead_axes.drug_annotation` still points at the derived 2026-08-19
snapshot. The refreshed file is committed as evidence the rule works and as the basis for a
future adoption decision, which is the researcher's call, not a build step.

**Artifacts:** `workflow/scripts/refresh_drug_annotation.py` (new), `rule refresh_drug_annotation`
in `workflow/rules/leads.smk`, `reference/drug_annotation/nerve_immune_axis_drug_annotation_refreshed_2026-08-19.csv`,
`refresh_manifest_2026-08-19.json`, `provenance/refresh_drug_annotation_2026-08-19.json`,
MANIFEST.json (refresh section, drift table, gotchas).

**Tool Versions:** pandas 2.3.3, requests 2.32.3 (`workflow/envs/scrna.yaml`); ChEMBL_37;
ClinicalTrials.gov API v2

**Open Issues:** None of the three items from the 2026-08-19 lead-axes entry remain open. Whether
to adopt the refreshed snapshot is an open *decision*, not an open defect.

**FAIR Notes:** The refreshed snapshot records its ChEMBL release — the pin the 2026-08-07 original
never had, and the single biggest reason that one is unreproducible. Refreshed snapshots use a
distinct `_refreshed_` filename infix so the rule's date wildcard can never resolve to a committed
pinned/derived file and shadow it.

---

### [2026-08-19] | Phase: Adopt the refreshed drug annotation | Status: COMPLETE

**Action:** Researcher accepted the ChEMBL-only answer, including its cost. Pointed
`config.lead_axes.drug_annotation` at the refreshed 2026-08-19 snapshot
(sha256 `0280e1f8…`, ChEMBL_37) and propagated it.

**Outcome:** The drug annotation is no longer inherited from an artifact nobody can
inspect — every value now traces to a recorded database release. Table sha256
`9869b698…` → `1ca78122…`.

**Adoption swaps only the annotation.** Verified on a `(axis, compartment_side)` key
rather than positionally, because rows sort by `tier_v2` and reordering otherwise masks
as a diff: `rank_full`, `rank_capped`, `full_best_mag`, `capped_best_mag`, `full_n_rows`,
`capped_n_rows`, `nerve_side`, `interfaces`, `immune` and `min_pval` are **all
bit-identical**. Changed: `tier_v2` 33/183, `agents_flagged` 82/183, `glioma_trials`
15/183, and `withdrawn` **0/183**.

That last number is the interesting one. The `withdrawn` column re-queried from ChEMBL's
`withdrawn_flag` is character-identical to the one derived days earlier by parsing inline
`[WITHDRAWN]` tags, so adoption does not disturb it at all.

**Tier movement in the 183-row table:**

| tier | before | after | Δ |
|---|---:|---:|---:|
| `1_approved_available` | 46 | 34 | −12 |
| `1b_approved_withdrawn_only` | 14 | 18 | +4 |
| `2_clinical` | 50 | 57 | +7 |
| `3_chembl_target_no_agent` | 38 | 57 | +19 |
| `4_no_chembl_target` | 35 | 17 | −18 |

The 12 tier-1 losses concentrate in the ITGB1 axes (`ITGB1|VCAN`, `ITGB1|LGALS1`,
`ITGB1|LGALS3BP`, the `ITGA*_ITGB1|SPP1` family, `CD14|ITGB1`) and the CALM/PDE1 axes —
exactly the genes whose tier-1 status came from the lost prior-pass artifact rather than
from ChEMBL. The 18 leaving `4_no_chembl_target` are the gain from exact gene-symbol
resolution.

**`CALM1|PDE1A` and `CALM1|PDE1C` moving to `1b_approved_withdrawn_only` is a correction,
not a loss.** These are the "CALM1 axes are a trap" cases the original generator's own
report warned about: the only approved agents on the CALM1 end are benziodarone and
prenylamine, both market-withdrawn. Under the old annotation a clean partner agent held
the axis at tier 1 and hid that; the ChEMBL-only view surfaces it.

**Verification:** notebook 05 re-rendered for both arms — Panel E still **183/183** on
magnitude and row counts, no regression banner, zero tracebacks (expected: the oracle
checks geometry, which adoption does not touch). `config` ↔ `MANIFEST._active_snapshot`
↔ table tier counts all cross-checked consistent. All three snapshots present on disk.

**Artifacts:** `config/config.yaml`, `reference/drug_annotation/MANIFEST.json` (rewritten
— it had accreted across three commits into claiming "nothing in this repository can
regenerate it", which the refresh rule had made false), `markdowns/GBM_TME_Crosstalk_Analysis.md`
§8 (tier table + provenance paragraph), both rendered notebook HTMLs.

**Open Issues:** None. All three lead-axes items are closed and the adoption decision is
made.

**FAIR Notes:** Both superseded snapshots stay committed and are marked "audit record, do
not delete or edit" in the manifest — each records what was true at its date. The manifest
now separates what is reproducible (the active snapshot, re-runnable against ChEMBL_37)
from what never will be (the 2026-08-07 original), instead of applying the original's
limits to all three.

---

### [2026-08-19] | Phase: Remove Panel F from notebook 05 | Status: COMPLETE

**Action:** Researcher asked for Panel F (whole-table overlap with the pinned v1.3.0
reference) to be removed — the reference is no longer a valid comparator, so the panel
has no value.

**Outcome:** Removed. The panel's own callout already argued at length that the
reference is not a comparator (its nerve compartment was built by the logic this
pipeline removed; no author annotation so never audited; scored against
differently-normalized data, defect D8; structurally unreproducible). A panel whose
callout tells the reader not to use its numbers is an invitation to use them anyway.
Cross-arm agreement — every Panel E axis clears the bar in **both** Census arms — is the
replication evidence this cohort actually has, and it does not need v1.3.0 to stand up.

Policy unchanged (§2.5 option (a)): v1.3.0 stays pinned and is simply no longer compared
against. `ds_cohort_concordance` still runs and still writes its outputs; the notebook
just no longer reads them. A removal note in the notebook records what was there and why
it went, so this does not read later as an accidental deletion.

**Removed with it:** the `concordance` and `shared_pairs` loads, their `paths` entries,
and their declarations on `ds_census_nerve_immune_notebook`.

**Also removed — a spurious dependency found while checking.** The notebook rule declared
`lead_targets = results/tables/nerve_crosstalk_lead_targets.csv` (the withdrawn 2026-07-15
reference shortlist), but notebook 05 never opens that file — its only mention is prose
inside a Panel E callout. Only notebook 04 reads it. The declaration was tying the Census
notebook to the pinned reference cohort for no reason.

**[BUG FOUND] Cyclic dependency, mine, latent since the adoption commit.** The full-DAG
dry run raised `CyclicGraphException on rule nerve_immune_lead_axes`. Adopting the
refreshed snapshot pointed `config.lead_axes.drug_annotation` at
`refresh_drug_annotation`'s own output, while that rule read the lead-axes table:
`nerve_immune_lead_axes -> refreshed snapshot -> refresh_drug_annotation -> lead-axes
table -> nerve_immune_lead_axes`.

It went unnoticed because every check after adoption was scoped with `--allowed-rules`,
which never builds the whole DAG — a real gap in that commit's verification, not a
harmless oversight: `rule all` would have failed for anyone running the pipeline plainly.
Fixed by sourcing the refresh rule's axis list from the committed, static 2026-08-07
snapshot instead of from the table it feeds. That breaks the cycle without changing what
gets queried, and a stale axis list is caught loudly rather than silently, because
`nerve_immune_lead_axes` already raises `[FAIR-ALERT]` when an axis has no annotation row.

**Also fixed:** `_cp` (compartment_pair array) was left dead in the Panel E cell when the
spurious `interfaces` restriction was removed earlier today — caught by `flake8 --select=F`
on the notebook, which had not been run against it before.

**Verification:** notebook parses; `flake8 --select=F` clean; rendered for **both** arms
with zero tracebacks; panels present are A–E with F absent and no "Jaccard overlap" text
remaining; Panel E still **183/183** on magnitude and row counts in both arms. Full
`rule all` dry run resolves with **0 cyclic errors**; the opt-in refresh target still
resolves to exactly 1 job. HTML shrank ~964K→764K (full arm) and ~1.0M→804K (capped).

**Artifacts:** `notebooks/05_census_nerve_immune_explorer.py`, `workflow/rules/notebooks.smk`,
`workflow/rules/leads.smk`, `workflow/scripts/refresh_drug_annotation.py`,
`markdowns/GBM_TME_Crosstalk_Analysis.md` (§3 DAG + §7 panel description), both rendered HTMLs.

**Open Issues:** None.

**FAIR Notes:** Removing a panel removes a claim from the record, so the notebook keeps an
in-place comment explaining what Panel F showed, why it was withdrawn, and that the
underlying rule still runs — future readers should not have to reconstruct that from git.

---

### [2026-08-19] | Phase: Demote the reference-cohort notebooks | Status: COMPLETE

**Action:** Researcher asked whether the non-05 notebooks still make sense, given they
read the pinned v1.3.0 reference. Investigated, then demoted all five.

**Outcome:** `01_explore_gbm_data`, `02_nerve_enrichment_explorer`,
`03_nerve_tumor_immune_explorer`, `04_tme_nerve_immune_explorer` and
`nerve_tumor_exploration` are no longer built by `rule all`. Their rules remain in
`workflow/rules/notebooks.smk` and can be invoked explicitly; the `.py` files stay.

**The finding that made this urgent:** `rule all` was actively rebuilding all five. They
are gated on `SAMPLES` (17 configured) and a dry run showed every one *pending*, so a
plain pipeline run would regenerate five HTMLs into `results/figures/` alongside notebook
05's, looking equally current — while every nerve-side number in them derives from the
compartment `nerve_cell_subset.py` was fixed on 2026-08-06, which measured 59% malignant
and 11% neural. This was not a dormant wart; it was a live source of void deliverables.

Every existing render predated the fix (May 23 – Jul 26). They cannot be corrected:
`data/processed/nerve_cells.h5ad` was deleted on 2026-07-21 and v1.3.0 is pinned and
structurally unreproducible.

**Demoted rather than deleted**, because the two things at issue are separable. The
rendered HTMLs and the automatic rebuild were the hazard; the notebook *code* is not.
If a v1.4.0 baseline is ever rebuilt through the corrected pipeline — which the removed
Panel F callout explicitly contemplated — 02/03/04 become usable again immediately. The
project already treats superseded artifacts this way: v1.3.0 stays pinned, the 2026-08-07
annotation snapshot is kept as an audit record, feature branches are retained.

**Banners added in two places per notebook**, because they serve different readers: the
module docstring for anyone opening the `.py`, and a `kind="danger"` marimo cell placed
second (right after imports) so it is the first thing in any rendered HTML. Both state
what the defect was, that nerve-side claims are void, that it cannot be corrected, why the
file is kept, and where the current analysis lives.

**Stale HTMLs removed** — but only after verifying every one is re-renderable: all five
rules had 100% of their declared inputs present on disk (2, 3, 3, 3 and 12 inputs
respectively, 0 missing). Deleting an unreproducible record would have destroyed the only
copy; deleting a reproducible one just removes a stale artifact.

**Verification:** all six notebooks parse; `02` re-rendered end-to-end as a smoke test —
0 tracebacks, banner present in the HTML, "59% malignant" text confirmed rendered; `rule
all` dry run now schedules `ds_census_nerve_immune_notebook` (×2 arms) and **no** reference
notebook rule; `results/figures/` root holds no notebook HTMLs, Census renders living
under `results/figures/<dataset>/`.

**Artifacts:** `Snakefile` (rule all block replaced with an explanatory comment carrying
the explicit re-render command), the five notebook `.py` files.

**Open Issues:** Notebook 05's docstring still names `notebooks/04_tme_nerve_immune_explorer.py`
as the reference-cohort sibling. That remains accurate — 04 still exists and still covers
that cohort — so it was left alone.

**FAIR Notes:** Demotion is recorded in three places a reader might look: the `Snakefile`
comment explaining why the targets are absent and how to render them anyway, the notebook
docstrings, and the rendered banner. Removing the targets silently would have looked like
an oversight and invited someone to add them back.

---

### [2026-08-19] | Phase: Census notebook set replaces the archived reference set | Status: COMPLETE

**Action:** Archive the five reference-cohort notebooks and replace them, where the data
supports it, with Census-cohort notebooks.

**Outcome:** `notebooks/` now holds four Census notebooks, all built by `rule all` for
both arms. The five reference notebooks moved to `notebooks/archive/` with their rules
deleted (see the previous entry).

| archived | replacement | reason |
|---|---|---|
| `01_explore_gbm_data` | `01_census_cohort_qc` | 170 donors vs 17 |
| `02_nerve_enrichment_explorer` | `02_census_nerve_enrichment` | compartment now 95.4% neural, was 11% |
| `nerve_tumor_exploration` | — | declined; overlaps Panels C/D of notebook 05 |
| `03_nerve_tumor_immune_explorer` | — | notebook 05 already *is* its Census equivalent |
| `04_tme_nerve_immune_explorer` | — | same |
| — | `06_census_compartment_audit` | **new, no ancestor possible** |

**Notebook 06 is the one that could not have existed for the reference cohort** and is
the reason the asymmetry matters: it cross-tabulates every compartment against the
CELLxGENE author annotation, and the TCGA reference carries no author annotation at all.
Five panels — the 10 enforcing gates with margin-to-threshold, compartment composition
against the oracle, the malignancy confusion matrix with its false-positive breakdown,
per-cluster nerve purity, and donor purity v1 vs scANVI-v2.

**[FINDING] It immediately surfaced something no existing artifact shows.** The nerve
compartment passes at 95.4% neural in aggregate, but **5 of 23 clusters in the full arm
sit below the 0.80 bar** — `c12` (1,314 cells, 16% neural / 83% malignant), `c19`, `c21`,
`c17`, `c18` — together 1,677 of 37,945 cells (4.4%), three of them majority-malignant.
The capped arm shows 4 of 25, 1,042 of 28,936 (3.6%). This is not a regression and does
not contradict the gate: it is the residue the mask could not separate, small enough that
the aggregate stays well clear. But it is exactly the failure mode a per-cluster view
exists to catch, and `batch_qc_pass` does **not** cover it — that flag is donor
dominance, not cell identity. The panel names the clusters and says to cross-reference
the interaction table's `nerve_cluster` column before reporting an axis resting on them.

**Second finding, from notebook 01:** the capped arm needs **63 of 170 donors** to reach
half its cells against **35 of 170** in the full arm — direct evidence that
`subsample_per_donor` did what it was designed to do, measured rather than assumed.

**Deliberately not built:** a nerve↔tumour LR notebook. `ds_nerve_tumor_interaction`
output exists, but Panels C and D of notebook 05 already browse the three-way table that
subsumes it; a fourth notebook would have duplicated it.

**Not portable, and stated rather than faked:** `nerve_clinical_association` and
`nerve_leiden_resolution_sweep` have no `ds_` twin, and the Census `gdc_clinical.tsv` is a
generated stub. Clinical association is genuinely lost with the reference cohort.

**Fixed while building:** `nerve_cluster_sample_purity.csv` writes `cluster` as str while
the v2 table writes int64, so the join raised. Normalised to str and kept the join outer —
v1 has 24 clusters and v2 20, and a cluster in only one is a real difference between
clusterings, not missing data. The panel now says the two are different partitions so the
scatter is not misread as paired measurements.

**Verification:** all four notebooks parse, `flake8 --select=F` clean, and render for both
arms with **0 cell errors**. Full `rule all` DAG resolves with no cyclic or missing-rule
errors. Notebook 01 reports 170 donors / 1,020,902 cells (full) and 624,688 (capped);
notebook 02 reports 160/460 and 180/500 enrichment rows from donor-dominated clusters;
notebook 06 reports 10/10 gates on both arms.

**Artifacts:** `notebooks/01_census_cohort_qc.py`, `notebooks/02_census_nerve_enrichment.py`,
`notebooks/06_census_compartment_audit.py`, three new rules in `workflow/rules/notebooks.smk`,
`Snakefile` (`rule all` now lists four Census notebooks per arm).

**Open Issues:** The 5 sub-0.80-neural nerve clusters are now visible but not acted on. Whether
to exclude them from LR aggregation, or flag axes that depend on them, is a scientific
decision for the researcher.

**FAIR Notes:** Notebook 01 globs the 170 per-donor QC files rather than declaring 340 inputs;
`annotation_summary.csv` is the declared edge because `ds_scrna_annotate` sits downstream of
`ds_scrna_qc` for every sample and so cannot exist before them. Stated in the rule docstring so
the shortcut is visible rather than implicit.

---

### [2026-08-19] | Phase: Notebook prose audit after the archive/replace work | Status: COMPLETE

**Action:** Researcher asked whether notebook 05's language had gone stale. Audited it and
the two notebooks written earlier today.

**Outcome:** Three classes of problem, one of them a factual error of mine.

**[BUG, mine] Notebook 01 described its own data backwards.** It called the population
"cells post-QC" and its limits section claimed "cells here are already post-QC ... this is
the surviving population, not the raw one." That is inverted: `scrna_qc.py` writes
`*_qc_metrics.csv` at line 37, *before* the filters at lines 93-95, so the file is the
**pre-filter** population. Introduced this morning when the notebook was written and
caught only by cross-checking its 1,020,902 against notebook 05's donor count.

Fixing it made the notebook better rather than merely correct: because the metrics are
pre-filter, the thresholds can be replayed from config to show what filtering removed.
Doing so reproduces the documented post-QC size **exactly** — 1,020,902 ingested →
1,006,344 passing (98.6%), matching `markdowns/GBM_TME_Crosstalk_Analysis.md` line 78.
Panel A now carries that funnel and names the hardest-hit donors (`LB3771T` 26%,
`MGH143` 25%, `LB4130T` 18%), which is a QC signal that was previously invisible.

**[FACT] Notebook 05 said 169 donors; the cohort has 170.** Five occurrences, all wrong,
contradicting the crosstalk markdown (which says 170 in four places), the 170
`*_qc_metrics.csv` files and the 170 distinct `sample_id`s in `*_gene_presence.csv`.
Corrected. The `169617` in `config.baseline.immune_baseline_n` is an unrelated cell count
and is probably where the digit came from.

**[STALE] References to the now-archived notebook 04.** Three survived the archive commit:
the footer still pointed at `notebooks/04_tme_nerve_immune_explorer.py` (wrong path), and
two prose mentions read as though it were live. Corrected, and the framing changed — 04 was
described as a "Companion", which it no longer is.

**Also updated in notebook 05:**
- Header and docstring now name the sibling notebooks (01, 02, 06) rather than a single
  archived predecessor, so a reader landing here knows where the cohort QC, enrichment and
  audit views are.
- The "both arms pass 10/10 compartment gates" claim now points at notebook 06, where it is
  rendered — **and states the finding notebook 06 surfaced**, that 5 of 23 nerve clusters
  sit below the 0.80 neural bar despite the aggregate passing at 95.4%. A page that quotes
  the aggregate should not omit that.
- The drug-annotation limit was written when the annotation was an un-refreshable snapshot.
  It now says the columns are re-derived from **ChEMBL_37** via `refresh_drug_annotation`
  and adopted 2026-08-19, while keeping the two limits that still hold (snapshot not live
  query; trial coverage bounded by a curated 20-agent list).

**Verification:** all three notebooks parse, `flake8 --select=F` clean, render on both arms
with 0 cell errors. Notebook 05 Panel E still 183/183 on both arms. Rendered HTML confirmed
to contain "170 donors", the sibling-notebook list and `ChEMBL_37`, and to contain neither
"169 donors" nor the stale `notebooks/04_tme` path.

**Artifacts:** `notebooks/01_census_cohort_qc.py`, `notebooks/05_census_nerve_immune_explorer.py`,
`markdowns/GBM_TME_Crosstalk_Analysis.md` (footer now distinguishes ingested from post-QC),
six re-rendered HTMLs.

**Open Issues:** None from this pass.

**FAIR Notes:** The QC thresholds used for the funnel are read from `config.scrna` rather
than hardcoded — if they change and the notebook is not re-run, the funnel would otherwise
silently describe a filter that no longer exists.

---

## [2026-08-27] Disk reclamation Tier 1 — git garbage pack + SCP393 cohort removed

**Phase:** Housekeeping ahead of public release. No analytical artifact touched; no rule re-run.

**Trigger:** The project directory measured **176 GB**, crowding out other work on the machine.
`data/` accounts for 154 GB of it. This entry covers only Tier 1 — items carrying zero
scientific risk. Tiers 2 and 3 are proposed but **NOT executed**, pending researcher approval.

### Removed 1 — `.git/objects/pack/tmp_pack_xP24MQ` (2.5 GB)

A leftover temporary pack from a `git repack`/`gc` interrupted on 2026-04-18. Git itself
classified it as unreachable garbage, not as a pack in use:

    before:  count: 1071  size: 70396 KB  in-pack: 0  packs: 0  size-pack: 0
             garbage: 1   size-garbage: 2599654 KB

Every reachable object in the repository was loose (~70 MB); the 2.5 GB file was referenced
by nothing.

**Verification before deletion (goal-backward):** full `git fsck --no-progress` reported only
`dangling` entries — no missing or broken objects. The file was then *renamed out of*
`.git/objects/` rather than deleted outright, and the repository was re-verified in that state:
`git fsck` clean, `git count-objects -v` reporting `garbage: 0`, `HEAD` resolving to `c08aa4c`,
`git log` and `git status` normal, and **all five branches resolving** (`main`, plus the four
feature branches, which are retained as audit trail and were not touched). Only after that did
the quarantined copy get deleted.

    after `git gc --prune=now`:  count: 0  in-pack: 1055  packs: 1  size-pack: 12.78 MiB

`.git/` is now **13 MB**, down from 2.5 GB. No history, ref, or branch was altered.

### Removed 2 — `broad_data/SCP393/` (1.3 GB, 13 files)

The Neftel et al. IDHwt-GBM replication cohort, triaged 2026-07-16 and **superseded that same
day** by the CELLxGENE Census pull (`cellxgene_data/gbm_10x_raw.h5ad`, 1.29M cells), which was
assessed HIGH suitability where SCP393's nerve arm was too thin. See
`markdowns/assessment_cellxgene_gbm_cohort.md`.

**Verification before deletion:** untracked by git (`git ls-files broad_data/` → 0), and
unreferenced across `workflow/`, `config/`, `scripts/`, `Snakefile`, and `notebooks/` by a
recursive grep for `broad_data` and `SCP393` over all `.py`/`.smk`/`.yaml`/`.yml`/`.sh` files.
No file had been modified since the 2026-07-17 triage.

**Manifest of removed files** (re-downloadable from the Broad Single Cell Portal, accession
SCP393 — this is published third-party data, not a project-generated artifact):

    569M  other/IDHwtGBM.processed.SS2.CNA.txt
    416M  expression/5e16df92771a5b0eb30ca010/IDHwtGBM.processed.10X.counts.mtx
    296M  expression/IDHwtGBM.processed.SS2.logTPM.txt.gz
     26M  expression/5e16dae3771a5b0eb30c9ff8/IDHwtGBM.processed.10X.counts.2.mtx
    1.4M  metadata/IDHwt.GBM.Metadata.SS2.txt
    plus 8 files under 1 MB (genes/cells TSVs, cluster, documentation, supplemental info)

### Outcome

**3.8 GB reclaimed. 176 GB → 172 GB.** No analytical artifact, provenance record, or pinned
reference was affected; `results/`, `provenance/`, `data/`, and `.snakemake/conda/` untouched.

**Open Issues:**

- **Tier 2 (~69 GB, awaiting researcher go-ahead).** The archived TCGA reference cohort:
  67 GB of root-level `data/processed/*.h5ad` (17 sample pairs at 42 GB + its scVI chain at
  24 GB) plus `data/raw/gdc_extract/` at 2.3 GB. `ds_cohort_concordance` reads only
  `results/tables/nerve_tumor_immune_interactions_with_qc.csv`, so no `.h5ad` is required by
  any surviving deliverable.
  **BLOCKER — must be fixed first:** `Snakefile:54` requests
  `data/processed/integrated_latent.h5ad` whenever `samples:` is non-empty and is **not**
  wrapped in `pinned_target()`, unlike the ~15 lines below it. Deleting the h5ads without
  first setting `samples: []` or adding that guard would schedule an ~11 h scVI retrain of
  the very cohort `baseline.pinned` exists to protect. Same exposure on `qc_summary.csv`,
  `annotation_summary.csv`, `cnv_heatmap.png`, and the `immune_*` root targets.
- **Tier 3 (~37 GB, not scheduled).** The 680 per-donor `{donor}.h5ad` / `{donor}_qc.h5ad`
  files across both live Census arms — deterministic derivatives of the 7.1 GB Census pull via
  `ds_ingest_dataset` → `ds_scrna_qc`, no scVI involved. Verify with a dry run before and
  after; sequence after Tier 2 so two changes are never in flight at once.

**FAIR Notes:** Both removals are of **non-project-generated** data — one a git internal
temp file, one a published third-party cohort retrievable from its original accession. No
project-generated intermediate artifact was deleted, so no provenance record was orphaned.
The SCP393 manifest is recorded above so the removal is auditable and the cohort is
re-obtainable at its stated accession.

---

## [2026-08-27] Disk reclamation Tier 2, Step 0 — v1.3.0 reference archived; a pinned artifact found inside the proposed deletion set

**Phase:** Housekeeping. Step 0 of the Tier 2 plan (see previous entry). **No deletion performed.**

### Why this step existed

The surviving v1.3.0 reference output is explicitly unreproducible (`baseline.pinned: true`, the
producing rules are not defined, `data/processed/nerve_cells.h5ad` destroyed 2026-07-21). It lives
under `results/` and `provenance/*`, both of which are **gitignored** — so it existed in exactly one
place on one disk, with 65 GB of deletions about to happen around it.

### FINDING — a pinned artifact was inside the Wave B deletion set

`provenance/pinned_reference_v1.3.0.json` defines the reference as **37 artifacts** under sha256
drift detection. Reading it rather than trusting the directory layout revealed that

    data/processed/nerve_cells_counts.h5ad   1,544,544,351 bytes

is artifact #1 of those 37 — and it sits in the root `data/processed/` chain that the Tier 2 plan
proposed deleting in Wave B. Deleting it would have flipped
`results/pinned_reference_verification.json` from `pass: true, n_missing: 0` to a failing state and
destroyed a pinned, unreproducible artifact.

**Correction to the Tier 2 plan: `data/processed/nerve_cells_counts.h5ad` is now EXCLUDED from
deletion.** It is the only root-level `.h5ad` appearing in the pinned manifest; the 17 sample pairs,
`integrated_latent`, `annotated`, `malignancy_labeled`, `immune_cells{,_labeled}`, `nerve_cells_v2`
and `nerve_cells_counts_labeled` are **not** in it and remain in scope.

### Integrity verified before and after

All 37 pinned artifacts re-hashed against the 2026-07-29 freeze: **37/37 byte-identical, zero
drift**, both before archiving and again afterwards. The reference has not moved since it was frozen.

### Archive created

`v1.3.0_reference_archive_2026-08-27.tar.gz` — 83 MB raw, **26 MB compressed**, 163 files:

- all root-level `results/tables/*.csv|json` (64 files, incl. the 19 MB
  `nerve_tumor_immune_interactions_with_qc.csv` that `ds_cohort_concordance` reads)
- all root-level `results/figures/*.png`
- all 70 `provenance/*.json`, including `pinned_reference_v1.3.0.json` itself
- `results/{fair_validation_report,conda_env_smoke_test,pinned_reference_verification}.json`
  and `snakemake_report.html`
- `data/raw/gdc_extract/MANIFEST.txt` — ids + md5s for all 17 GDC looms, which is what makes the
  2.3 GB of raw loom data re-downloadable rather than lost
- `SHA256SUMS.txt` covering all 163

**Deliberately excluded:** `data/processed/nerve_cells_counts.h5ad` (1.54 GB). It is no longer being
deleted, so it needs no deletion-protection, and including it would have bloated a 26 MB archive to
1.6 GB for no gain. It remains single-copy on disk — see Open Issues.

**Verification (goal-backward, not assumed):** the tarball was extracted to a clean directory and
every file checked against its recorded hash — **163/163 OK, 0 failed**. An unverified archive is
not a backup.

    sha256(v1.3.0_reference_archive_2026-08-27.tar.gz)
      = f5f788d6f79a88bcf52db7b0c0be762ec83bd78c377ec4b611ee9f733557cfc7

**Open Issues:**

- **The archive is still on the same disk it protects against.** It guards the deletion operation,
  not drive failure. 26 MB — it should go to the Zenodo deposit already needed for the public
  notebook release, or any off-machine backup.
- **`data/processed/nerve_cells_counts.h5ad` (1.54 GB) remains single-copy and unreproducible.**
  Not at risk from Tier 2 any more, but it belongs in a real off-machine backup alongside the
  archive.
- Tier 2 Steps 1-7 remain **not started**, pending researcher go-ahead. Revised recovery:
  ~63 GB (Wave A 42 GB + Wave B 20.9 GB), down from the 65 GB estimated before this finding.

**FAIR Notes:** The archive scope was derived from `pinned_reference_v1.3.0.json` — the project's
own declaration of what the reference *is* — rather than from directory structure. That is what
surfaced the misplaced 1.54 GB artifact; a scope guessed from `du` output would have deleted it.

### Addendum — Steps 1-2 (guard + dry run). Two further findings; HOLDING before any deletion.

**Step 1 applied:** `Snakefile:54` now wraps `integrated_latent.h5ad` in `pinned_target()`, matching
the ~15 pinned targets below it. Confirmed **inert while the file exists** — `pinned_target` returns
the path unchanged, and the dry run still lists it in `rule all: input:`. It only takes effect once
the file is gone.

**FINDING A — the reference cohort is ALREADY dirty, independent of this cleanup.** A dry run at
the project's canonical invocation (`--use-conda --rerun-triggers mtime`) schedules **11 jobs**, four
of which are reference-cohort rules:

| rule | why | consequence |
|---|---|---|
| `scrna_malignancy` | `data/external/ensembl/ensembl113_gene_positions.tsv` is newer than `malignancy_labeled.h5ad` | would **re-run reference CNV** |
| `scrna_qc_report` | 4 of 17 `_qc_metrics.csv` newer than `qc_summary.csv` (2026-05-24) | rebuilds `qc_summary.csv` |
| `immune_cell_subset`, `immune_cluster_annotations` | downstream of the above | rebuild immune artifacts |

The gene-positions file was rebuilt by the 2026-08-05/06 fix (defect D4, coordinate-ordered genes).
So **a bare `scripts/run_snakemake.sh --use-conda --cores all` today would already attempt to re-run
CNV on the pinned reference** — a pre-existing landmine, not something this cleanup introduced.
(A bare dry run *without* `--rerun-triggers mtime` schedules **778** jobs on code/env/param triggers,
which is why every documented invocation in this project pins that flag.)

**FINDING B — the one-line guard from Step 1 is NOT sufficient for Wave B.** Because
`scrna_malignancy` is already scheduled, deleting its input `annotated.h5ad` makes it unsatisfiable
and it schedules `scrna_annotate` -> `scrna_integration` (~11 h scVI) -> `scrna_qc` x17 -> the GDC
looms. `pinned_target` only suppresses a *missing* file that `rule all` requests directly; it cannot
suppress a dirty rule whose other outputs (`cnv_heatmap.png`, `annotation_summary.csv`) are still
demanded unguarded and still exist.

**Measured fix — `samples: []`.** Re-running the dry run with `--config 'samples=[]'` collapses every
`if SAMPLES else []` guard and drops the reference from `rule all` entirely:

    default (samples populated):  11 jobs  — incl. scrna_malignancy, scrna_qc_report,
                                             immune_cell_subset, immune_cluster_annotations
    with samples=[]:               7 jobs  — all four reference rules GONE; the remaining 7 are
                                             Census-arm/cross-arm work unrelated to this cleanup

**Revised Step 1: set `samples: []` in `config/config.yaml`** (the `pinned_target` guard is retained
as correct-by-symmetry, and protects the case where `samples:` is ever restored). This both unblocks
Wave B and defuses Finding A. It is a config change only, fully reversible, and matches the fact
that the reference is already retired per README and `baseline.pinned`.

**Status: HOLDING.** Steps 0-2 complete, nothing deleted. Steps 3+ (Wave A, 42 GB) await researcher
go-ahead on the revised Step 1.

---

## [2026-08-27] Disk reclamation Tier 2, Steps 1-3 — reference retired from `rule all`; Wave A executed

**Phase:** Housekeeping. Revised Step 1 + Step 3 (Wave A). Wave B **not** executed.

### Step 1 — `samples: []` in `config/config.yaml`

The 17 GDC accessions are **commented out in place, not deleted**, with a restore note; every
reference target in `rule all` is written `... if SAMPLES else []`, so an empty list drops the whole
reference block from the default target and restoring the list restores prior behaviour exactly.

Rationale recorded inline in `config.yaml`: the cohort was already archived in substance
(`baseline.pinned: true`, producing rules undefined, `nerve_cells.h5ad` destroyed 2026-07-21,
notebooks in `notebooks/archive/`, README calls it not a valid comparator) — only `rule all` still
treated it as live — **and** it was silently dirty against post-fix gene positions.

`Snakefile:54`'s `pinned_target()` guard from the previous entry is retained as correct-by-symmetry;
it becomes the active protection if `samples:` is ever restored.

**Verified:** config parses, `samples == []`, both Census datasets present, `baseline.pinned` still
`True`. Dry run **11 jobs -> 7 jobs**, with `scrna_malignancy`, `scrna_qc_report`,
`immune_cell_subset` and `immune_cluster_annotations` gone. The 7 remaining are Census-arm and
cross-arm jobs that were already pending and are unrelated to this cleanup.

### Step 3 — Wave A: 17 reference sample pairs deleted (45.5 GB)

**Deletion set built programmatically from the captured accession list, not from a shell glob**, and
gated on four assertions that all had to pass before any `rm`:

1. exactly 34 files (17 accessions x `{,_qc}.h5ad`)
2. **none present in `pinned_reference_v1.3.0.json`** — the check that caught `nerve_cells_counts.h5ad`
   in Step 0
3. all 34 exist
4. every path is root-level `data/processed/` — structurally cannot touch a namespaced Census arm

    data/processed: 151G -> 109G     34/34 removed, none remaining

### Post-deletion verification (goal-backward)

| check | result |
|---|---|
| pinned reference re-hashed | **37/37 unchanged**, 0 drifted, 0 missing |
| dry run job count | **7 — identical to pre-deletion** |
| `scrna_qc` / `scrna_integration` scheduled? | **no** — no cascade, the ~11 h scVI retrain did not arm |
| Census arms | untouched, both still resolve |

**Outcome: 45.5 GB reclaimed. Project 172 GB -> 130 GB** (176 GB at session start).

**Open Issues:**

- **Wave B not executed** (~20.9 GB), pending researcher go-ahead. Remaining root `data/processed`:
  `malignancy_labeled` 5.7G, `annotated` 5.7G, `integrated_latent` 5.6G, `immune_cells_labeled` 1.6G,
  `immune_cells` 1.6G — plus `data/raw/gdc_extract/` 2.3 GB (keep `MANIFEST.txt`).
  **Excluded from Wave B by decision:** `nerve_cells_counts.h5ad` (1.4G, sha256-pinned) and, held
  back as higher-risk, `nerve_cells_v2.h5ad` (1.5G) + `nerve_cells_counts_labeled.h5ad` (1.4G) — the
  scANVI-v2 branch is deliberately NOT pinned, so its rules exist and could schedule against the
  already-destroyed `nerve_cells.h5ad`.
- The 26 MB archive still sits on the disk it protects; it belongs off-machine.
- `CHANGELOG.md`, `Snakefile`, `config/config.yaml` modified and uncommitted.

**FAIR Notes:** No pinned artifact, provenance record, or `results/` deliverable was touched. The
reference's citable output survives intact on disk *and* in the verified archive; what was removed
was per-sample intermediate input, regenerable in principle from the GDC looms via the retained
`MANIFEST.txt` md5s.

---

## [2026-08-27] Disk reclamation Tier 2, Step 5 — Wave B executed. Tier 2 complete.

**Phase:** Housekeeping, final Tier 2 step. 176 GB -> **107 GB** across the session.

### Removed (22 files, 24.1 GB)

Five root-level reference chain artifacts — `malignancy_labeled.h5ad`, `annotated.h5ad`,
`integrated_latent.h5ad`, `immune_cells.h5ad`, `immune_cells_labeled.h5ad` — plus the 17 GDC
`.seurat.1000x1000.loom` files under `data/raw/gdc_extract/`.

**Held back by decision (4.3 GB), as planned:**

| file | reason |
|---|---|
| `nerve_cells_counts.h5ad` | sha256-pinned artifact #1 of the 37 |
| `nerve_cells_v2.h5ad` | scANVI-v2 branch, deliberately NOT pinned — its rules exist |
| `nerve_cells_counts_labeled.h5ad` | same branch |

### The looms are recoverable, and that was verified rather than assumed

`data/raw/gdc_extract/MANIFEST.txt` was **kept** (and is inside the archive). Before deleting, all
17 looms were hashed against it: **17/17 md5 and byte-size identical to their GDC records**, so the
2.3 GB is re-downloadable with `gdc-client download -m MANIFEST.txt`. Deleting the looms without
first proving the manifest matched would have made the loss silent and permanent.

    data/  111G -> 89G          22/22 removed, MANIFEST.txt retained

### Post-deletion verification

| check | result |
|---|---|
| pinned reference re-hashed | **37/37 unchanged**, 0 drifted, 0 missing |
| dry run | **7 jobs — unchanged from before Wave A** |
| `scrna_*` cascade armed? | **no** |
| Census arms | 46 GB + 38 GB, both intact |
| reference deliverables | 64 root `results/tables` files, 70 `provenance/*.json`, and the 19 MB table `ds_cohort_concordance` reads — all present |
| held-back files | all 3 present |

### Session outcome

    176 GB -> 107 GB      69 GB reclaimed

    Tier 1  git garbage pack + SCP393            3.8 GB
    Wave A  17 reference sample pairs           45.5 GB
    Wave B  reference chain + GDC looms         24.1 GB

`data/` is now 89 GB, of which 84 GB is the two live Census arms — legitimate working data.

**Open Issues:**

- **Tier 3 (~37 GB) remains available and unstarted** — the 680 per-donor `{donor}{,_qc}.h5ad`
  files across both Census arms, deterministic derivatives of the 7.1 GB Census pull via
  `ds_ingest_dataset` -> `ds_scrna_qc`, no scVI involved. Would land the project near 70 GB.
- **The 26 MB archive is still on the disk it protects.** It guards against the deletions performed
  here, not drive failure. `nerve_cells_counts.h5ad` (1.4 GB, pinned, unreproducible) is likewise
  single-copy. Both belong in the Zenodo deposit needed for the public notebook release.
- The reference cohort is now **retired in fact as well as in documentation**. Anything that wants
  to rebuild it needs the GDC looms re-downloaded and `samples:` restored — and per
  `markdowns/plan_pin_v1_3_0_reference.md` it still would not reproduce v1.3.0.

**FAIR Notes:** Every removed byte was either (a) a per-sample or chain *intermediate*, or (b) raw
input carrying a verified md5 in a retained manifest. No pinned artifact, no `provenance/` record,
no `results/` deliverable, and no Census-arm file was touched. The reference's citable output
survives on disk and in a checksum-verified archive.

---

## [2026-08-27] Disk reclamation Tier 3 — per-donor Census intermediates removed. 176 GB -> 71 GB.

**Phase:** Housekeeping, final tier. Both **live** Census arms — a different risk class from Tier 2.

### Why this needed a different method

Tier 2 operated on a retired cohort whose producing rules are undefined (`baseline.pinned`). Tier 3
touches the two **active** arms, whose `ds_ingest_dataset` and `ds_scrna_qc` rules exist and are
reachable from `rule all`. A missing declared output normally schedules its producer, so the danger
was cascading into `ds_scrna_integration` — an ~11 h scVI train, twice.

**Files were therefore MOVED to a staging directory on the same filesystem, not deleted**, and the
DAG was interrogated in that state before anything became irreversible. An instant rename is a free
experiment; a deletion is not.

    stage capped arm (340 files)   -> dry run: 7 jobs, unchanged
    stage full arm too (340 files) -> dry run: 7 jobs, unchanged
    => no ds_scrna_qc, no ds_ingest_dataset, no ds_scrna_integration scheduled
    => only then: rm -rf staging

Snakemake does not regenerate an input whose downstream output is already up to date, so the arms
remain satisfied from their chain artifacts alone. This was **verified, not assumed** — the reasoning
above is exactly the kind that is right in principle and wrong in practice often enough to test.

### Removed (680 files, 37 GB)

Per-donor `{donor}.h5ad` and `{donor}_qc.h5ad` for all 170 donors in each arm:

    gbm_cellxgene_56c4912d        340 files, 26.1 GB
    gbm_cellxgene_56c4912d_full   340 files, 13.6 GB
    data/  89G -> 52G

Both arms retain their full chain — `integrated_latent`, `annotated`, `malignancy_labeled`,
`immune_cells{,_labeled}`, `nerve_cells{,_v2,_counts,_counts_labeled}` — plus `gene_symbol_map.tsv`.

**Regenerable:** these are deterministic derivatives of `cellxgene_data/gbm_10x_raw.h5ad` (7.1 GB,
retained) via `ds_ingest_dataset` -> `ds_scrna_qc`. No scVI, no MPS training.

### Verification

| check | result |
|---|---|
| dry run | **7 jobs — unchanged across all three tiers** |
| pinned reference | **37/37 unchanged**, 0 drifted, 0 missing |
| arm chains | 22 GB + 26 GB, complete |
| notebook inputs | 170 `_qc_metrics.csv` + 170 `_gene_presence.csv` per arm, all present |
| `results/tables/{arm}/` | 48 MB + 57 MB, intact |

### Session outcome — all tiers

    176 GB -> 71 GB      105 GB reclaimed (60%)

    Tier 1  git garbage pack + SCP393            3.8 GB
    Tier 2  Wave A, 17 reference sample pairs   45.5 GB
    Tier 2  Wave B, reference chain + looms     24.1 GB
    Tier 3  680 per-donor Census intermediates  37.0 GB

Remaining 71 GB: `data/` 52 GB (both live arm chains), `.snakemake/conda/` 8.2 GB (path-keyed,
never rebuild), `cellxgene_data/` 7.1 GB (the Census source), `claude_science/` 2.3 GB,
`results/` 1.1 GB.

**Open Issues:**

- **Nothing further is safe to delete without losing something expensive or unreproducible.** The
  remaining `data/` is scVI output; `.snakemake/conda/` is explicitly protected by CLAUDE.md.
- **Off-machine backup is now the outstanding risk, not disk space.** Single-copy and unreproducible:
  `v1.3.0_reference_archive_2026-08-27.tar.gz` (26 MB) and `data/processed/nerve_cells_counts.h5ad`
  (1.4 GB, sha256-pinned). Both belong in the Zenodo deposit required for the public notebook release.

**FAIR Notes:** No pinned artifact, provenance record, `results/` deliverable, or arm chain artifact
was touched in any tier. Everything removed was a regenerable intermediate or raw input with a
retained, verified checksum manifest.
