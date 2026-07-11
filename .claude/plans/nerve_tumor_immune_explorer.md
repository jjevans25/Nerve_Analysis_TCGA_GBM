# Plan: Nerve–Tumor–Immune Interaction Analysis + Marimo Explorer

> **CLAUDE.md compliance note:** Per project constitution, plans live in
> `.claude/plans/` inside the project. This file is the plan-mode working copy;
> the **first execution step** is to copy it to
> `NERVE_ANALYSIS_TCGA_GBM/.claude/plans/nerve_tumor_immune_explorer.md`.

## Context

The pipeline already infers **nerve↔tumor** ligand–receptor communication
(LIANA+ consensus, `nerve_tumor_interaction.py` → 1,207 significant L-R pairs
between the malignant compartment and 27 nerve Leiden clusters). The **immune
compartment is annotated but never entered any interaction analysis**: the
whole-dataset annotation lumps all myeloid/lymphoid cells into a single
`cell_type_predicted == "microglia"` group (46,777 cells) because its Leiden
resolution is coarse. Verified in-object: the AnnData carry the full
**20,420-gene** matrix with `gene_symbol` in `var`, and immune subtype markers
are present (microglia 6/7, macrophage/TAM 6/6, T-cell 4/6, NK 2/3, DC 2/4), so
subtypes are resolvable by re-clustering the immune subset.

**Goal:** resolve immune subtypes, compute three-way nerve–tumor–immune L-R
communication, and ship an interactive marimo notebook to explore it.

**Design decisions (confirmed with researcher):**
- **Immune compartment:** full subclustering → microglia / TAM / T-cell / NK / DC.
- **Nerve compartment:** kept at 27 Leiden clusters (`nerve_c0..c27`), matching
  the existing nerve↔tumor analysis.
- **Method:** LIANA+ `rank_aggregate`, consensus resource — the house method
  (consistency with the existing interaction script).

## Approach (Snakemake-first, FAIR-compliant)

Three stages. Every new script reuses `workflow/scripts/fair_utils.py`
(`log_transformation`, `stamp_artifact`, `verify_artifact`, `write_provenance`)
and every rule emits a provenance JSON, per CLAUDE.md.

### Stage 1 — Immune subclustering (new upstream analysis)

Mirrors the nerve-subset workflow. **No model retraining** — reuse the existing
`X_scVI` latent already stored in `malignancy_labeled.h5ad`.

**New config block** in `config/config.yaml` (pattern-copy `nerve_cells:`):
```yaml
immune_cells:
  leiden_resolution: 1.0            # tune upward if subtypes under-split
  source_label: "microglia"        # cell_type_predicted value = immune blob
  subtype_markers:
    microglia:  [P2RY12, CX3CR1, TMEM119, AIF1, C1QA, C1QB, C1QC]
    tam:        [CD68, CD163, MRC1, MARCO, LYVE1, F13A1]
    t_cell:     [CD2, CD3D, CD3E, CD8A, CD4, IL7R]
    nk_cell:    [GNLY, KLRD1, NKG7]
    dendritic:  [FLT3, LAMP3, CLEC9A, FCER1A]
  batch_qc:                         # reuse thresholds/pattern from nerve_cells
    dominant_fraction_max: 0.5
    min_contributing_fraction: 0.01
    min_contributing_samples: 3
    exclude_clusters: []            # populate after inspecting composition
```

**New rule `immune_cell_subset`** (add to a new `workflow/rules/immune.smk`,
`include:` it in `Snakefile`) + **`workflow/scripts/immune_cell_subset.py`**
(adapt `nerve_cell_subset.py`):
- Input `data/processed/malignancy_labeled.h5ad`.
- Select `cell_type_predicted == source_label` & `~is_malignant`.
- `sc.pp.neighbors(use_rep="X_scVI")` → `sc.tl.leiden` (→ `immune_leiden`) →
  `sc.tl.umap`.
- Outputs: `data/processed/immune_cells.h5ad`,
  `results/figures/immune_cells_umap.png`,
  `results/tables/immune_cluster_composition.csv`, provenance JSON.

**New rule `immune_cluster_annotations`** + **`workflow/scripts/immune_cluster_annotations.py`**
(adapt the marker-scoring half of `scrna_annotate.py`):
- Score the five `subtype_markers` panels per `immune_leiden` cluster
  (`sc.tl.score_genes`), assign dominant `immune_subtype` per cluster.
- Outputs: `data/processed/immune_cells_labeled.h5ad`,
  `results/tables/immune_cluster_annotations.csv`, provenance JSON.
- Small/low-confidence clusters flagged; genuinely artifactual ones added to
  `immune_cells.batch_qc.exclude_clusters` (GBM is immune-cold → expect small
  T/NK/DC clusters; handle with the existing exclusion convention).

### Stage 2 — Three-way interaction (new)

**New rule `nerve_tumor_immune_interaction`** (in `immune.smk`) +
**`workflow/scripts/nerve_tumor_immune_interaction.py`** (adapt
`nerve_tumor_interaction.py`):
- Inputs: `malignancy_labeled.h5ad` (malignant), `nerve_cells.h5ad`
  (`nerve_c{N}` from `nerve_leiden`), `immune_cells_labeled.h5ad`
  (`immune_{subtype}` from `immune_subtype`).
- `_to_symbol_index` + inner-join concat on shared HGNC symbols (reuse the
  existing helper logic); sanity-check `X.max()` for log1p (`use_raw=False`).
- Build `cell_label` = compartment-tagged group; restrict `groupby_pairs` to
  **cross-compartment** directional pairs only (nerve↔tumor, nerve↔immune,
  immune↔tumor — both directions), exclude within-compartment pairs.
- `li.mt.rank_aggregate(..., resource_name="consensus", n_perms=1000,
  expr_prop=0.10, seed=random_seed, use_raw=False)`.
- Add `compartment_pair` and `direction` columns; sort stably.
- Outputs: `results/tables/nerve_tumor_immune_interactions.csv`,
  `nerve_tumor_immune_top_pairs.csv`,
  `results/figures/nerve_tumor_immune_sig_heatmap.png`,
  `nerve_tumor_immune_dotplot.png`, provenance JSON.

**QC-flag variant:** extend `workflow/scripts/annotate_cluster_qc.py` (and its
rule) to left-join batch-QC verdicts and emit
`nerve_tumor_immune_interactions_with_qc.csv` +
`nerve_tumor_immune_top_pairs_with_qc.csv` (same `_with_qc` convention the
notebooks consume; applies nerve-cluster QC flags and drops `exclude_clusters`).

### Stage 3 — Marimo notebook (new)

**New `notebooks/03_nerve_tumor_immune_explorer.py`** (adapt
`notebooks/nerve_tumor_exploration.py`, the newest/most sophisticated notebook).
Follow all conventions: module docstring (`Marimo reactive notebook:` /
`Addresses:` / `Source rule:`), `__generated_with = "0.23.1"`,
`app = marimo.App(width="wide", app_title=...)`, named `@app.cell` functions,
`_load_config` resolving paths from `config["dirs"]` via
`Path(__file__).parent.parent` (no hardcoded paths), DuckDB for CSV ingestion,
`mo.stop` input guards naming the Snakemake rule to run first, matplotlib +
seaborn plots, `_with_qc` batch-QC flagging (warn callout + `*` ticks), trailing
provenance callout, and an export cell writing a filtered CSV + `.provenance.txt`.

Consumes `nerve_tumor_immune_interactions_with_qc.csv` +
`nerve_tumor_immune_top_pairs_with_qc.csv`. Views:
1. **Compartment-pair selector** (`mo.ui.radio`/`multiselect`): nerve↔tumor,
   nerve↔immune, immune↔tumor — filters all downstream views.
2. Filtered L-R table + significance heatmap (compartment/cluster × direction).
3. Top-K L-R dotplot per selected pairing (matplotlib scatter, magnitude vs
   specificity rank), with L-R and cluster search boxes + threshold sliders.
4. **Relay-circuit view** (the novel three-way angle): surface signaling triads
   — e.g. tumor→immune ligand paired with immune→nerve ligand — to reveal
   candidate tumor↦immune↦nerve relays.

**New rule `nerve_tumor_immune_notebook`** in `workflow/rules/notebooks.smk`
(copy `nerve_tumor_exploration_notebook`), using the shared
`workflow/scripts/run_notebook_export.py`. Output
`results/figures/03_nerve_tumor_immune_explorer.html` + provenance JSON.

### Wiring

- `Snakefile`: `include: "workflow/rules/immune.smk"`; add the new `_with_qc`
  tables + the notebook HTML to `rule all` (gated on `SAMPLES`, matching the
  existing notebook-target pattern).

## Files created / modified

**Create**
- `workflow/rules/immune.smk` (immune subset, annotations, 3-way interaction rules)
- `workflow/scripts/immune_cell_subset.py`
- `workflow/scripts/immune_cluster_annotations.py`
- `workflow/scripts/nerve_tumor_immune_interaction.py`
- `notebooks/03_nerve_tumor_immune_explorer.py`

**Modify**
- `config/config.yaml` — add `immune_cells:` block
- `workflow/scripts/annotate_cluster_qc.py` + its rule in `nerve_cells.smk` —
  emit the new `_with_qc` interaction variants
- `workflow/rules/notebooks.smk` — add `nerve_tumor_immune_notebook`
- `Snakefile` — `include` immune.smk; extend `rule all`

**Reuse (no change)**
- `workflow/scripts/fair_utils.py`, `workflow/scripts/run_notebook_export.py`,
  `workflow/envs/scrna.yaml` (LIANA/scvi), `workflow/envs/notebooks.yaml` (marimo)

## Coding standards (CLAUDE.md)

Python 3.12 + type hints; `pathlib.Path` (no `os.path` in new logic);
`ruff format` + `ruff check` before commit; no bare `except`; set
`PYTHONHASHSEED` + seeds for determinism; config-driven paths only; atomic git
commit per verified artifact.

## Verification (goal-backward)

Run from project root using the project venv `claude_science/` (per
`execution_instructions.md`).

1. **DAG dry-run:** `claude_science/bin/python3 -m snakemake --use-conda --cores 8 -n`
   — confirm the new immune + interaction + notebook jobs appear, no cycles.
2. **Build immune subset & inspect:**
   `snakemake --use-conda --cores 8 data/processed/immune_cells_labeled.h5ad`
   then check `results/tables/immune_cluster_annotations.csv` — confirm
   microglia/TAM and (if present) T/NK/DC subtypes separate; tune
   `immune_cells.leiden_resolution` if under-split; add tiny/artifact clusters
   to `exclude_clusters`.
3. **Run interaction:**
   `snakemake --use-conda --cores 8 results/tables/nerve_tumor_immune_interactions_with_qc.csv`
   — verify non-empty, `compartment_pair` covers all three pairings, and
   grepping `immune_`/`microglia`/`t_cell` returns hits (the gap this closes).
4. **Export notebook:**
   `snakemake --use-conda --cores 1 results/figures/03_nerve_tumor_immune_explorer.html`
   — confirm HTML renders; then drive it interactively with
   `marimo edit notebooks/03_nerve_tumor_immune_explorer.py` to confirm the
   compartment-pair selector, dotplot, and relay-circuit view react correctly.
5. **FAIR check:** confirm a provenance JSON exists for every new rule output
   (`ls provenance/*immune*.json`) and each carries input hashes + tool versions.
6. Append a `CHANGELOG.md` entry (new science) once artifacts verify.

## Notes / risks

- **Immune-cold tumor:** T/NK/DC clusters may be small; expected — use the
  `batch_qc.exclude_clusters` convention rather than forcing splits.
- **Leiden resolution** is the main tunable; `immune_cluster_composition.csv` is
  the inspection point before locking it.
- LIANA permutations are the only heavy compute (no VAE retrain); runtime is
  comparable to the existing nerve↔tumor interaction run.
