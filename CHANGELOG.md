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
Phase:          Annotation — Gene ID Namespace Fix
Last Updated:   2026-04-20
Active Agent:   lead-researcher
Current Task:   MyGene.info Ensembl→symbol cache wired into loom_to_h5ad and scrna_annotate
Blocked On:     —
Next Action:    Re-run pipeline end-to-end (requires network for first MyGene.info fetch); verify nerve_cell_subset now returns neurons/astrocytes/oligodendrocytes
```

---

## Session Log

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
