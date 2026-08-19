# Bring `nerve_immune_lead_axes_postfix.csv` into the Snakemake pipeline

## Context

`results/tables/nerve_immune_lead_axes_postfix.csv` (183 axes) is the only input to
notebook 05 Panel E that no rule produces. It was generated out of band in a Claude
Science session; `results/` is gitignored, so the file can neither be rebuilt nor
restored from git. Notebook 05 raises a `[FAIR-ALERT]` and shows curated-vs-live
numbers side by side instead of trusting it, and `markdowns/GBM_TME_Crosstalk_Analysis.md`
§8 closes with "Committing the generator would close that gap permanently."

`claude_science/code_export/` now supplies that generator as four Claude Science REPL
transcripts. This plan turns the reproducible part into a real rule and vendors the
part that cannot be recomputed.

**Outcome:** the table gets a producing rule, a provenance sidecar, and a committed
input for every column — so it is rebuildable offline and its drug annotation is
version-controlled rather than orphaned.

### What the generator actually does (verified, not inferred)

The CSV is a `pd.concat` of **two independently-built halves**, which is the root of
the `best_mag` discrepancy the notebook reports:

| | oligodendrocyte half (128 rows) | neuron half (55 rows) |
|---|---|---|
| source arm | `..._full` | `..._56c4912d` (capped) |
| row filter | `cellphone_pvals<=0.05 & magnitude_rank<=0.05 & batch_qc_pass` | same **minus `batch_qc_pass`**, and `nerve_cluster=='nerve_neuron'` only |
| sort key | `['best_mag','min_pval']` | `['best_mag']` only |
| `rank_full` | dense `1..128` (`index+1`) | sparse: position in unfiltered full-arm table |
| `rank_capped` | sparse: position in unfiltered capped table | dense `1..55` (`index+1`) |
| `capped_best_mag` | capped-arm value | **column never created → 55 NaNs** |

So `best_mag` means "full arm, QC-passing" for 128 rows and "capped arm, no QC,
neuron-only" for 55. That is why a single-arm live recomputation reproduces
`capped_best_mag` 128/128 but `best_mag` only 128/167. `n_rows` is arm-inconsistent
the same way. Both explorations independently rebuilt both halves from the on-disk
interaction tables and reproduced the CSV exactly (128/128 and 55/55 on every
numeric column).

### What cannot be recomputed

`tier_v2`, `agents_flagged`, `glioma_trials`, `withdrawn` came from
`host.mcp("chembl"/"clinical-trials", ...)` calls that run only inside Claude Science.
**No `handoff/*.json` cache survived** — I searched the repo; the directory does not
exist. One input is permanently lost: a prior-session artifact referenced only by UUID
`34d720d0-5c16-46bb-92e6-d8d6b6a02ec3`, which supplied fallback agents for genes with
no ChEMBL mechanism record (`evidence_source == 'prior_reference_pass'`).

Verified property that makes vendoring safe: those four columns are a **pure function
of the `axis` string** — across the 39 axes that appear on both nerve sides, zero
disagreements. So a 144-row axis-keyed snapshot reproduces all 183 rows losslessly.

### Decisions taken

1. **Vendor the drug annotation as a committed snapshot, plus a separate opt-in
   refresh rule** that re-queries public ChEMBL + ClinicalTrials.gov REST into a *new*
   dated file. Default builds stay offline, deterministic, and exact; refresh never
   mutates the pinned snapshot silently.
2. **Reproduce byte-identically first, fix the defects in a second commit.** Exact
   reproduction is the only available proof that the reimplementation is correct.
3. **Scope: `nerve_immune_lead_axes_postfix.csv` only.** No sibling tables, figure,
   or report.

---

## Commit 1 — rule that reproduces the table byte-for-byte

### 1a. Preserve the irreplaceable artifact first

`results/` is gitignored, so the existing CSV has **no git safety net** and the new
rule will overwrite it. Before anything else:

```bash
mkdir -p reference/drug_annotation
cp results/tables/nerve_immune_lead_axes_postfix.csv \
   "$SCRATCH/nerve_immune_lead_axes_postfix.ORIGINAL.csv"
shasum -a 256 results/tables/nerve_immune_lead_axes_postfix.csv > "$SCRATCH/original.sha256"
```

This copy is the acceptance oracle for the whole commit. Do not skip it.

### 1b. Vendor the drug annotation (new tracked directory)

`data/` and `results/` are both gitignored, so neither can host a committed input.
Create a tracked top-level `reference/` directory — semantically honest (external
reference data, not a provenance record) and needs no `.gitignore` negation, unlike
the `!provenance/cl15_split_*.csv` pattern.

- `reference/drug_annotation/nerve_immune_axis_drug_annotation_2026-08-07.csv`
  — 144 rows, columns `axis, tier_v2, agents_flagged, glioma_trials, withdrawn`,
  extracted from the preserved original via
  `df[['axis',...]].drop_duplicates('axis').sort_values('axis')`.
- `reference/drug_annotation/curated_agent_map_2026-08-07.json`
  — the two hardcoded literals the refresh rule will need:
  `AGENTS` (20 names, `02_chembl_druggability_pipeline.py:226-228`) and `AG2GENE`
  (20 entries, same file `:250-253`).
- `reference/drug_annotation/MANIFEST.json` — must record, plainly:
  source = ChEMBL + ClinicalTrials.gov via Claude Science MCP on 2026-08-07;
  **ChEMBL release version unknown/unpinned**; the lost prior-pass artifact UUID and
  which genes depend on it; sha256 of the snapshot; and the known `withdrawn`-column
  bug (below).

### 1c. Config block

Add to `config/config.yaml` (rules read config and pass values via `params:`; scripts
never open config themselves):

```yaml
lead_axes:
  pval_max: 0.05
  magnitude_rank_max: 0.05
  full_arm:   "gbm_cellxgene_56c4912d_full"
  capped_arm: "gbm_cellxgene_56c4912d"
  neuron_group: "nerve_neuron"
  drug_annotation: "reference/drug_annotation/nerve_immune_axis_drug_annotation_2026-08-07.csv"
  compartment_side_oligo:  "oligodendrocyte (QC-passing)"
  compartment_side_neuron: "neuron (pooled; donor-caveated)"
```

### 1d. `workflow/scripts/nerve_immune_lead_axes.py`

New script. Follow `workflow/scripts/cohort_concordance.py` exactly for boilerplate —
module docstring stating the biological question, `snakemake.input/output/params`
(never argparse), `sys.path.insert(0, "workflow/scripts")` then
`from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance`,
`log = snakemake.log[0]`, `pathlib.Path`, type hints on every function.

**Gotcha:** do *not* add `from __future__ import annotations` — Snakemake's `script:`
directive prepends a preamble, so it stops being the first statement and raises
`SyntaxError`. This is documented at `workflow/scripts/download_gene_positions.py:22-25`.

Reuse `fair_utils.nerve_group_key` (`fair_utils.py:166`) for the cluster key — do not
re-slice the `nerve_c{N}` prefix by hand.

Functions, mirroring `01_build_corrected_axis_tables.py`:

- `_prep(df, ann)` — coerce the three `*_batch_qc_pass` columns through the
  `{True:True,'True':True,False:False,'False':False}` map then
  `.fillna(True).astype(bool)`; derive `cluster_key` via `nerve_group_key`; merge
  `nerve_cell_type` from the annotation `label` split on `|` taking `parts[1].strip()`;
  force `neuron(pooled)` where `cluster_key == 'neuron'`; build the undirected axis as
  `'|'.join(sorted([ligand_complex, receptor_complex]))`.
- `_axes_table(df, qc_only=True)` — the oligo aggregation; sort
  `['best_mag','min_pval']`.
- `_pooled_neuron_axes(df)` — `nerve_cluster == neuron_group`, **no QC filter**;
  sort `['best_mag']` only.
- `_assemble(...)` — intersect each half across arms, take rows from the correct arm,
  assign dense vs sparse ranks per the table above, map `capped_best_mag` onto the
  oligo half only, join the vendored annotation on `axis`, set
  `nerve_side` (oligo: `nerve_types.replace('', '—')`; neuron: literal
  `'neuron(pooled)'`) and `compartment_side`, then emit the 15 columns in file order
  and `sort_values(['tier_v2','rank_full'])` / `(['tier_v2','rank_capped'])` per half
  before `pd.concat(..., ignore_index=True)`.

Three faithful-reproduction details that are easy to lose:
- The column list in `03_lead_targets_report_postfix.py:2` contains `qc_pass`, which
  exists in neither half and is silently dropped by the `[c for c in cols if c in
  df.columns]` filter → the file has **15 columns, not 16**. Keep it dropped.
- The differing sort keys between halves are load-bearing: pandas' stable mergesort
  means ties resolve by prior order, so using `['best_mag','min_pval']` for the neuron
  half would shift ranks.
- Row order is the two sorted halves concatenated — **not** global rank order.

Close with `verify_artifact`, `log_transformation`, `stamp_artifact` +
`write_provenance`, matching the call site at `cohort_concordance.py:143-158`. Use
`ontology_operation="operation:3501"`.

### 1e. `workflow/rules/leads.smk`

New file; `include:` it from `Snakefile` after `notebooks.smk`. The rule is
**un-wildcarded** — the output sits at `results/tables/` root and consumes both arms in
one table, so it cannot use the `{dataset}` wildcard. Model it on
`ds_cohort_concordance` (`workflow/rules/datasets.smk:710-733`), which likewise mixes a
namespaced input with a root-level one:

```python
rule nerve_immune_lead_axes:
    """Cross-arm reproducible nerve–immune signalling axes, tiered by druggability."""
    input:
        full_lr    = os.path.join(config["dirs"]["tables"], config["lead_axes"]["full_arm"],
                                  "nerve_tumor_immune_interactions_with_qc.csv"),
        capped_lr  = os.path.join(config["dirs"]["tables"], config["lead_axes"]["capped_arm"],
                                  "nerve_tumor_immune_interactions_with_qc.csv"),
        full_ann   = ...  # nerve_cluster_annotations.csv, full arm
        capped_ann = ...  # nerve_cluster_annotations.csv, capped arm
        drug       = config["lead_axes"]["drug_annotation"],
    output:
        table      = os.path.join(config["dirs"]["tables"], "nerve_immune_lead_axes_postfix.csv"),
        provenance = os.path.join(config["dirs"]["provenance"], "nerve_immune_lead_axes_provenance.json"),
    log:    os.path.join(config["dirs"]["logs"], "nerve_immune_lead_axes.log")
    conda:  "../envs/scrna.yaml"
    resources: mem_mb = 8000, threads = 1
    params: ...  # the lead_axes config block
    script: "../scripts/nerve_immune_lead_axes.py"
```

Use `../envs/scrna.yaml` — every existing pure-pandas rule does, and
`config.env_enforcement.expected_versions` is written against it. Add the table to
`rule all`. Do **not** gate on `BASELINE_PINNED`; this is a census-arm product, not a
frozen v1.3.0 artifact.

### 1f. Acceptance gate

The commit is not done until this passes:

```bash
diff <(sort "$SCRATCH/nerve_immune_lead_axes_postfix.ORIGINAL.csv") \
     <(sort results/tables/nerve_immune_lead_axes_postfix.csv) && echo "BYTE-IDENTICAL"
```

Float formatting is the likely failure mode. If only float repr differs, compare with
`pandas.testing.assert_frame_equal(..., rtol=0, atol=0)` and record in CHANGELOG that
the match is exact-valued rather than exact-bytes — do not paper over a *value*
difference, which would mean the reconstruction is wrong.

### 1g. Refresh rule (opt-in, never in `rule all`)

`rule refresh_drug_annotation` in the same `.smk`, writing
`reference/drug_annotation/nerve_immune_axis_drug_annotation_{date}.csv` — a **new**
file, so the pinned snapshot is never overwritten. Model the network handling on
`workflow/scripts/download_gene_positions.py` (retry + backoff, manifest sidecar).

- ChEMBL REST (`https://www.ebi.ac.uk/chembl/api/data/`): `target` search by
  `pref_name`/gene synonym filtered to `Homo sapiens` + `SINGLE PROTEIN`, then
  `PROTEIN COMPLEX`/`PROTEIN FAMILY` for genes with no single-protein hit; `mechanism`
  by `target_chembl_id`; `molecule` by name for `max_phase`, `withdrawn_flag`,
  `black_box_warning`. No API key.
- ClinicalTrials.gov v2 (`https://clinicaltrials.gov/api/v2/studies`):
  `query.cond=glioma OR glioblastoma`, `query.intv=<agent>`, over the vendored
  `AGENTS` list; map back to genes via vendored `AG2GENE`.
- Tiering logic from `02_chembl_druggability_pipeline.py:489-503` — gene tier
  `1_approved_available` / `1b_approved_withdrawn_only` / `2_clinical` / `3_` / `4_`,
  then axis roll-up as the minimum over constituent genes with
  `T2 = {'1_':1, '1b':1.5, '2_':2, '3_':3, '4_':4}` keyed on `t[:2]`.
  Genes are split `axis.split('|')` then `part.split('_')` to expand LIANA subunits.
- **Must record the ChEMBL release in its manifest** — the missing release pin is the
  main reason the original is unreproducible.

The rule must state in its docstring that output **will** differ from the 2026-08-07
snapshot (both databases have moved) and that adopting it is a deliberate config bump,
reviewed by diff.

---

## Commit 2 — fix the defects, with the notebook

Separate and independently revertible, because this changes the science.

1. **`withdrawn` is empty in all 183 rows — a real bug.** `axis_drug` sources it from
   `DRUG.withdrawn_agents`, populated from the bulk chembl_id lookup, which per the
   export's own README "does not carry withdrawal status". The corrected set came from
   the later by-name sweep and was never written back. Yet `tier_v2` marks 14 axes
   `1b_approved_withdrawn_only` and `agents_flagged` carries inline `[WITHDRAWN]` tags —
   so the data exists and is recoverable by parsing `agents_flagged`. Already tracked as
   open in `CHANGELOG.md` (2026-08-07 entry).
2. **Arm-inconsistent `best_mag` / `n_rows`.** Add `capped_best_mag` for the neuron
   rows, add an explicit `full_best_mag`, and label `n_rows` by arm. Keep the existing
   column names populated as they are today so the notebook does not break silently.
3. **Update notebook 05 Panel E** in the same commit: drop the `_unbuildable` entry and
   the `[FAIR-ALERT]` (`notebooks/05_census_nerve_immune_explorer.py:117-149`), and
   update the curated-vs-live comparison at `:816-821` and `:965-967`, whose 128/128 and
   128/167 counts become stale the moment the columns are fixed.
4. Update `markdowns/GBM_TME_Crosstalk_Analysis.md` §8 — the "not produced by this
   pipeline" paragraph and the closing sentence this work exists to retire.

---

## Verification

Run through `scripts/run_snakemake.sh` — never bare `snakemake`. An active venv in the
calling shell silently reverts the conda pins, and without `--use-conda` the unquoted
interpreter path splits on the space in "Biomedical Data Science", giving an empty log
and a tracebackless `CalledProcessError`.

Scope the run with `--allowed-rules` so a stale mtime on the 11–13 MB interaction
tables cannot trigger an expensive upstream LIANA rebuild:

```bash
scripts/run_snakemake.sh results/tables/nerve_immune_lead_axes_postfix.csv \
  --use-conda --cores 1 --allowed-rules nerve_immune_lead_axes
```

1. **Dry run first** (`-n -r`) — confirm only `nerve_immune_lead_axes` is scheduled.
2. **Byte-identity gate** (§1f) — the primary acceptance test.
3. **Row/column shape**: 183 rows, 15 columns; 128 `oligodendrocyte (QC-passing)` +
   55 `neuron (pooled; donor-caveated)`; `capped_best_mag` non-null on exactly the 128
   oligo rows; oligo `rank_full` exactly `1..128`; neuron `rank_capped` exactly `1..55`.
4. **Tier counts** match §8 of the markdown: `1_approved_available` 46,
   `1b_approved_withdrawn_only` 14, `2_clinical` 50, `3_chembl_target_no_agent` 38,
   `4_no_chembl_target` 35.
5. **Provenance**: `provenance/nerve_immune_lead_axes_provenance.json` exists and
   `rule fair_validate_metadata` (first entry of `rule all`) passes.
6. **Notebook**: re-render notebook 05 and confirm Panel E loads and its live-vs-curated
   panel still agrees (commit 1) / reflects the fix (commit 2).
7. `flake8 workflow/scripts/nerve_immune_lead_axes.py` clean before committing.
8. Append a `CHANGELOG.md` entry per the project constitution; commit on the current
   branch `fix/notebook04-post-compartment-fix` or a new one — do not delete branches.

## Risks

- **The existing CSV is overwritten by the first successful run and is not in git.**
  §1a is the only protection.
- **The drug annotation is frozen, not recomputed.** Commit 1 closes the provenance and
  versioning gaps; it does not make the ChEMBL/trials annotation recomputable, and the
  lost prior-pass artifact means some genes' tiers can never be re-derived exactly. The
  MANIFEST must say so plainly rather than implying the table is fully reproducible.
- The refresh rule's output will not match the snapshot, by construction.
