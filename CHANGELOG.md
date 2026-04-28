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
