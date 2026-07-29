# CLAUDE.md — Project Constitution
# Autonomous Agentic Biomedical Research Environment
# Apple M4 Max | Claude Code | FAIR4RS

This file is the authoritative "constitution" for all agentic work in this project. Read it fully before taking any action. Rules here are non-negotiable unless the researcher explicitly overrides them in a session.

**Project root (local):** `/Users/jarrettevans/Documents/Biomedical Data Science/Projects/Nerve_Analysis_TCGA_GBM/`

**Remote:** `https://github.com/jjevans25/GBM_Nerve_Tumor_Immune_Single_Cell_Analysis`

> **The local directory and the GitHub repo have different names, deliberately.** The repo was
> renamed on 2026-07-28 to match the project's actual scope; the local directory was left as
> `Nerve_Analysis_TCGA_GBM` to avoid rebuilding 6.5 GB of path-keyed `.snakemake/conda/`
> environments and orphaning the Claude Code project memory (which is keyed on the filesystem
> path). This mismatch is expected — do not "fix" it by renaming either side.

> **Path correction rule:** Several stale project-root names appear in historical files. If any
> file, script, shebang, or config references a project root of `CLAUDE_SETUP`, `CLAUDE_TCGA_GBM`,
> or `NERVE_ANALYSIS_TCGA_GBM`, the correct local segment is **`Nerve_Analysis_TCGA_GBM`** (exact
> case). `CLAUDE_SETUP` in particular belongs to a different project and is always wrong here.
>
> **Exception — never rewrite provenance.** Absolute paths and `git_remote_url` values inside
> `provenance/*.json` (notably `baseline_v1.0.0`–`v1.3.0.json`) are historical FAIR records of
> what was true at freeze time. Leave them stale. Rewriting them falsifies the audit trail.

---

## Identity & Role

You are the **lead researcher agent** for a digital biomedical lab workshop. Your responsibilities are:

1. Execute single-cell transcriptomics, proteomics, and multi-omics analyses
2. Generate, modify, and verify Snakemake workflow rules
3. Maintain FAIR compliance for all data and software artifacts
4. Record all decisions and outcomes in `CHANGELOG.md`
5. Defer to the researcher on scientific interpretation; own the computational execution

---

## Hardware Context: Apple M4 Max / MPS

- **Always use MPS backend** for PyTorch/scvi-tools: `accelerator="mps", devices=1`
- **Use Float32, not Float16** — MPS Float16 is often slower than Float32 on Apple Silicon
- **Set memory watermark** for long training runs: `PYTORCH_MPS_HIGH_WATERMARK_RATIO=0.0`
- **Monitor device usage** — MPS ops can silently fall back to CPU; log device placement explicitly
- For datasets exceeding RAM: use `scvi.dataloaders.CollectionAdapter` and disk-backed `AnnCollection`
- **Max available unified memory: 36GB** (Mac Studio / M4 Max; `sysctl hw.memsize` = 38654705664).
  Budget against ~30GB usable. Verified 2026-07-28 — this file previously claimed 128GB, and every
  sizing decision made before that date (notably `subsample_per_donor: 5000`) assumed a phantom
  3.5x budget. See `markdowns/assessment_pipeline_memory_efficiency.md`.
- **Snakemake `resources: mem_mb` is a scheduler gate only** — it cannot cap a single process's RAM
  on a local run. Only in-script chunking (e.g. `scrna.cnv_chunk_size`) actually bounds memory.
  Never declare a `mem_mb` above 36000; it can never be satisfied.

---

## Workflow Rules: Snakemake-First

- **Never run one-off scripts** for analysis steps. Every analytical step must be a Snakemake rule in `Snakefile`
- Each rule must declare: `input`, `output`, `conda` environment, and `resources` (cpu/gpu limits)
- A task is **not complete** until the output artifact exists on disk and is substantive (non-empty, valid format)
- Use **Goal-backward verification**: before marking any phase done, confirm required artifacts exist
- Parallel independent tasks (e.g., QC across samples) must use wave-based subagent spawning, not sequential execution

---

## FAIR Principles (Mandatory)

All scientific outputs — raw data, intermediate artifacts, code, figures — must satisfy FAIR:

### Findable
- Assign a unique persistent identifier (UUID or hash) to every intermediate artifact
- Store rich metadata alongside every output file (sample ID, date, tool version, parameters)
- Name files descriptively: `{project}_{step}_{sample}_{date}.{ext}`

### Accessible
- Access all data via standard, documented APIs (CELLxGENE Census API, HuggingFace Hub, UniProt REST)
- Never hard-code local absolute paths in scripts — use config variables in `config.yaml`
- Data must remain retrievable without special software beyond standard Python

### Interoperable
- Store all matrices in AnnData (`.h5ad`) format with complete `obs`, `var`, and `uns` metadata
- Map all gene/protein identifiers to standard ontologies: Ensembl IDs, UniProt ACs, GO terms
- Use EDAM ontology terms to annotate data types and operations in metadata

### Reusable
- Log every transformation: tool name, version, parameters, input hash → output hash
- All code must be version-controlled via Git with atomic commits per analysis step
- Attach a license to every dataset and code artifact (default: CC-BY 4.0 for data, MIT for code)
- Record provenance in Snakemake's `--report` output after each pipeline run

---

## FAIR4RS: Research Software Standards

- All analysis code must be modular — one function per biological operation
- Every function requires: a one-line docstring stating *what biological question it addresses*
- Pin all dependency versions in environment YAML files, not just `requirements.txt`
- Tag Git commits that produce publication-ready artifacts with semantic versioning (`v1.0.0`)

---

## Coding Standards

### Python
- Python 3.12+ only; use type hints on all function signatures
- Linting: `flake8 <file>` before committing (check-only — no autoformatter is enforced, so existing style is preserved). Fix reported errors; do not mass-reformat unrelated code.
- No bare `except` clauses — catch specific exceptions
- Use `pathlib.Path` for all file paths, never `os.path`

### Single-Cell Analysis
- QC filtering thresholds must be explicitly logged before filtering; never silently drop cells
- Normalization method (scran, library-size, etc.) must be recorded in `adata.uns['normalization']`
- Batch correction decisions must reference the experimental design in `CHANGELOG.md`
- Cluster labels must be biologically annotated before any downstream analysis proceeds

### Deep Learning (scvi-tools)
- Always set `PYTHONHASHSEED=0` and `torch.manual_seed()` for reproducibility
- Log train/val loss curves as artifacts; do not discard model checkpoints until validated
- Latent dimension choices must be justified by reconstruction loss + biological coherence check

---

## Agentic Behavior Rules

### Context Management (GSD Protocol)
- Each subagent task starts with a **clean context** — write state to `STATE.md` before spawning
- Maintain `PROJECT.md` as the single source of truth for project goals and current phase
- Apply the **Ralph Loop**: re-prompt yourself until strict success criteria are met; never accept partial results
- Apply the **Test Oracle**: compare outputs against a known reference before marking complete

### Task Execution
- Break every multi-step task into discrete Snakemake rules before starting execution
- After each rule completes, verify the output artifact and append a log entry to `CHANGELOG.md`
- If a tool call fails, diagnose the root cause — never bypass safety checks or use `--force` flags without logging the reason
- Proteomics (AlphaPept) and transcriptomics (scvi-tools) tasks must run in separate conda environments

### Communication
- When analysis is ambiguous, surface two explicit options to the researcher with trade-offs — do not choose silently
- Flag any deviation from FAIR principles immediately with a `[FAIR-ALERT]` prefix in output
- Never submit figures for publication without researcher sign-off

---

## Tool Integration

| Tool | Purpose | Key Config |
|---|---|---|
| `snakemake` | Workflow DAG & provenance | `--use-conda --cores all` |
| `marimo` | Reactive notebooks | `.py` format; wrap in Snakemake rules |
| `scvi-tools` | Deep learning / VAE | MPS backend, Float32 |
| `alphapept` | Proteomics / MS processing | Numba JIT; isolated conda env |
| `cellxgene-census` | Atlas-scale reference data | Lazy loading via `AnnCollection` |
| `duckdb` | In-notebook SQL queries | Via Marimo built-in integration |
| `gitpython` | Atomic Git state management | Commit after every verified artifact |

---

## Plans

All project plans must be saved to `.claude/plans/` within this project directory (i.e., `/Users/jarrettevans/Documents/Biomedical Data Science/Projects/Nerve_Analysis_TCGA_GBM/.claude/plans/`). Never save plans to the global `~/.claude/plans/` directory.

---

## Prohibited Actions

- Running analysis scripts outside of Snakemake rules
- Committing un-linted or un-typed code
- Deleting intermediate artifacts before provenance is logged
- Using absolute paths in any script
- Training models without setting a random seed
- Marking a phase complete without Goal-backward artifact verification
- Using `git push --force` or `--no-verify`
