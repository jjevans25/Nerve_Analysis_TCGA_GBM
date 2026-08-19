# Execution Instructions — GBM Nerve–Tumor–Immune Pipeline

> ## ⚠️ Read first — two invocation rules that are not optional
>
> **1. Always launch through `scripts/run_snakemake.sh`.**
>
> ```bash
> scripts/run_snakemake.sh <targets…> --use-conda --cores all
> ```
>
> Do **not** run `snakemake` directly, and do **not** `source claude_science/bin/activate` first.
> `claude_science` is a venv whose `bin/python` symlinks into `/opt/anaconda3`. When it is active
> in the calling shell, `VIRTUAL_ENV` is exported and `claude_science/bin` sits first on `PATH`;
> Snakemake's job shells inherit that, `conda activate` loses, and **every pin in
> `workflow/envs/*.yaml` is silently un-enforced**. The run succeeds and the results are computed
> against the wrong package versions. `run_snakemake.sh` scrubs the environment before handing
> off.
>
> **2. `--use-conda` is not optional.** Without it, Snakemake's unquoted interpreter path splits
> on the space in `Biomedical Data Science` and the run dies with an **empty log and no
> traceback**.
>
> **3. Targets go FIRST.** `--allowed-rules`, `--forcerun` and `--quiet` all take `nargs='+'` and
> will swallow any target placed after them.

> ## ⚠️ The v1.3.0 reference cohort is pinned and archived
>
> `config/config.yaml` sets `baseline.pinned: true`. The 10 rules that would rewrite reference
> nerve/immune artifacts are **not defined at parse time**, so no invocation can schedule them.
>
> This is now a historical concern rather than an operational one: the reference cohort was
> **archived on 2026-08-19**. Its nerve compartment was built by the logic the 2026-08-06
> compartment fix removed (59 % malignant / 11 % neural), and it cannot be corrected — its
> `nerve_cell_subset` run is insulated by a v1.1.0 frozen barcode list that bypasses the fixed
> mask entirely, so re-running it would faithfully reproduce the defect. Its notebooks are in
> `notebooks/archive/` and are not built.
>
> - `--forcerun nerve_cell_subset` fails with "no such rule". That is the pin working.
> - `verify_pinned_reference` is the **only** drift detector for those files — while the pin is
>   on, their producing rules do not exist, so Snakemake cannot notice a modified or deleted
>   pinned table on its own.
> - See README "Read this first" for what reviving it would actually take.

## Prerequisites

All commands run from the project root:

```
/Users/jarrettevans/Documents/Biomedical Data Science/Projects/Nerve_Analysis_TCGA_GBM/
```

`scripts/run_snakemake.sh` finds its own Python. You do **not** need to activate anything — see
rule 1 above.

For long runs:

```bash
tmux new -s gbmfull
caffeinate -ims scripts/run_snakemake.sh …
```

macOS idle-sleep keys off user input, not CPU load, so an overnight run is suspended without
`caffeinate`.

---

## Step 1 — Validate the DAG (dry run)

Always first.

```bash
scripts/run_snakemake.sh --use-conda --cores 8 -n
```

What to check:

- the DAG resolves with **no errors** (in particular no `CyclicGraphException`)
- only `ds_*` rules, `nerve_immune_lead_axes` and the four notebook rules appear
- **no reference-cohort rule you did not ask for** — if a nerve cascade appears, add
  `--allowed-rules` before running for real

Job count varies with which arms are enabled and what already exists on disk, so there is no
fixed expected number.

---

## Step 2 — Smoke test one donor

Verifies ingest + QC and builds the `scrna` conda env (~5–10 min the first time).

```bash
scripts/run_snakemake.sh \
  results/tables/gbm_cellxgene_56c4912d_full/3182_qc_metrics.csv \
  --use-conda --cores 4
```

Then check marker detectability before committing to the full run:

```bash
cat results/tables/gbm_cellxgene_56c4912d_full/3182_gene_presence.csv
```

`01_census_cohort_qc` Panel D shows this across all 170 donors once the pipeline has run — a
marker absent from most donors cannot carry a cell-type call.

---

## Step 3 — Full pipeline

```bash
caffeinate -ims scripts/run_snakemake.sh --use-conda --cores all --rerun-triggers mtime
```

Approximate timings on an M4 Max, **per arm**:

| Stage | Approx. time |
|---|---|
| 170× `ds_ingest_dataset` | 20–40 min |
| 170× `ds_scrna_qc` | 30–60 min |
| `ds_scrna_integration` (scVI VAE, MPS) | 2–5 h |
| `ds_scrna_annotate` | 15–30 min |
| `ds_scrna_malignancy` (CNV, chunked) | 30–60 min |
| `ds_nerve_cell_subset` / `ds_immune_cell_subset` | 10–25 min |
| `ds_compartment_audit` | 2–5 min |
| `ds_nerve_cell_heterogeneity` + GSEA | 15–40 min |
| `ds_nerve_tumor_immune_interaction` (LIANA) | 10–30 min |
| `ds_nerve_scanvi_retrain` | 45–120 min (MPS) |
| `nerve_immune_lead_axes` | seconds |
| 4× notebook export | 1–2 min |

Budget overnight-plus for a cold full-arm run. The 170 ingest and QC jobs parallelise; the
singletons serialise, with scVI integration the long pole.

Reference-cohort stages are skipped — those outputs are frozen tables.

---

## Step 3b — One arm only

```bash
scripts/run_snakemake.sh \
  results/tables/gbm_cellxgene_56c4912d_full/compartment_audit_gates.csv \
  results/figures/gbm_cellxgene_56c4912d_full/05_census_nerve_immune_explorer.html \
  --use-conda --cores all --rerun-triggers mtime
```

The two arms are **separate namespaces, not a config edit**. `gbm_cellxgene_56c4912d` is capped
at `subsample_per_donor: 5000` (624,688 cells); `gbm_cellxgene_56c4912d_full` is uncapped
(1,020,902 → 1,006,344 post-QC). Raising the cap in place would overwrite the capped arm and
destroy it as a comparison. Cross-arm agreement is this project's replication evidence, so both
arms must exist.

---

## Step 4 — Goal-backward verification

### 4a. The compartment audit is the gate — check it first

```bash
column -s, -t results/tables/gbm_cellxgene_56c4912d_full/compartment_audit_gates.csv
```

Expect **10/10 PASS**. `ds_compartment_audit` fails the build on any FAIL, so a green run
already implies this — but read the margins. `nerve_neural_fraction` (0.954) and
`tumor_malignant_fraction` (0.930) are the two that were catastrophically wrong before the fix.

```bash
# Per-cluster purity — the aggregate can pass while hiding a contaminated cluster
python3 -c "
import pandas as pd
c = pd.read_csv('results/tables/gbm_cellxgene_56c4912d_full/nerve_compartment_cluster_audit.csv')
print(c[c.frac_neural < 0.8][['nerve_leiden','n_cells','frac_neural','frac_malignant']])"
```

5 of 23 clusters are currently below the 0.80 bar (4.4 % of cells, three majority-malignant).
That is expected residue, not a regression — but an axis resting on them is suspect, and
`batch_qc_pass` does **not** flag it (that column is donor dominance, not cell identity).

### 4b. Per-arm artifacts

```bash
ARM=gbm_cellxgene_56c4912d_full
ls -lh results/tables/$ARM/{compartment_audit_gates,malignancy_confusion,annotation_summary}.csv \
       results/tables/$ARM/nerve_tumor_immune_interactions_with_qc.csv \
       results/tables/$ARM/nerve_enrichment_with_qc.csv \
       results/tables/$ARM/nerve_cluster_sample_purity{,_v2}.csv
ls -lh results/figures/$ARM/*.html
```

### 4c. Cross-arm shortlist

```bash
ls -lh results/tables/nerve_immune_lead_axes_postfix.csv
python3 -c "
import pandas as pd
d = pd.read_csv('results/tables/nerve_immune_lead_axes_postfix.csv')
print(len(d), 'rows /', d.axis.nunique(), 'axes'); print(d.tier_v2.value_counts())"
```

Expect **183 rows / 144 axes**, tiers 34 / 18 / 57 / 57 / 17.

### 4d. LIANA sanity

```bash
head -3 results/tables/gbm_cellxgene_56c4912d_full/nerve_tumor_immune_interactions_with_qc.csv
```

`specificity_rank` must not be empty. Both interaction scripts hard-fail on an off-scale `.X`,
but an empty rank column is the fingerprint of the pre-2026-07-26 raw-counts defect.

### 4e. Pinned reference drift

```bash
scripts/run_snakemake.sh verify_pinned_reference --use-conda --cores 1
# -> results/pinned_reference_verification.json
```

> `data/processed/nerve_cells.h5ad` and `results/figures/nerve_cells_umap.png` are **absent by
> design** for the reference cohort — lost 2026-07-21. Do not treat their absence as a failure.
> Both Census arms have their own under `{cohort}/`.

---

## Step 5 — FAIR validation and report

```bash
scripts/run_snakemake.sh results/fair_validation_report.json \
  --use-conda --cores 1 --allowed-rules fair_validate_metadata
python3 -c "import json;print(json.load(open('results/fair_validation_report.json'))['total_provenance_records'])"
# 70 as of 2026-08-19

scripts/run_snakemake.sh --use-conda --cores all --report results/snakemake_report.html
```

---

## Step 6 — Interactive exploration (marimo)

Four notebooks, each built per arm. Set `GBM_DATASET` to choose the arm; it defaults to the
first entry in `config.datasets`.

```bash
GBM_DATASET=gbm_cellxgene_56c4912d_full \
  claude_science/bin/python3 -m marimo edit notebooks/06_census_compartment_audit.py
```

```bash
claude_science/bin/python3 -m marimo run notebooks/05_census_nerve_immune_explorer.py
```

| Notebook | Scope |
|---|---|
| `01_census_cohort_qc.py` | Cohort scale, per-donor QC, filtering funnel, composition, marker detectability |
| `02_census_nerve_enrichment.py` | Nerve-cluster GSEA, shared-profile overlap, marker genes |
| `05_census_nerve_immune_explorer.py` | Three-way LR browser, interface matrix, cross-arm lead axes |
| `06_census_compartment_audit.py` | **The Test Oracle** — read this one first |

All four are exported to `results/figures/{arm}/*.html` by the pipeline; editing is only needed
for interactive work.

`notebooks/archive/` holds five superseded reference-cohort notebooks. They have **no rules** and
cannot be built — their nerve-side numbers are void. They carry a `SUPERSEDED` banner explaining
why.

---

## Targeted re-runs

> ⚠️ `--forcerun nerve_cell_subset` does not work and should not be made to work. While
> `baseline.pinned: true` that rule is not defined. Tune resolution on a Census arm instead,
> where the `ds_*` equivalents are live.

Always dry-run (`-n`) first and read the job list before dropping it.

**Re-run an arm's nerve chain** (e.g. after changing `leiden_resolution`):

```bash
scripts/run_snakemake.sh \
  results/tables/gbm_cellxgene_56c4912d_full/compartment_audit_gates.csv \
  --use-conda --cores all --forcerun ds_nerve_cell_subset
```

**Re-run only the LIANA rules** (~7 min, reference untouched):

```bash
scripts/run_snakemake.sh --use-conda --cores all --rerun-triggers mtime \
  --forcerun ds_nerve_tumor_interaction ds_nerve_tumor_immune_interaction
```

**Rebuild the cross-arm shortlist** without touching anything upstream:

```bash
scripts/run_snakemake.sh results/tables/nerve_immune_lead_axes_postfix.csv \
  --use-conda --cores 1 --allowed-rules nerve_immune_lead_axes
```

**Re-render one notebook**:

```bash
scripts/run_snakemake.sh \
  results/figures/gbm_cellxgene_56c4912d_full/06_census_compartment_audit.html \
  --use-conda --cores 1 --allowed-rules ds_census_compartment_audit_notebook --force
```

**Refresh the drug annotation** (opt-in, needs network, ~10 min). Writes a **new dated file** and
never touches the pinned one; adopting it means bumping `config.lead_axes.drug_annotation`
deliberately after a diff:

```bash
scripts/run_snakemake.sh \
  reference/drug_annotation/nerve_immune_axis_drug_annotation_refreshed_$(date +%F).csv \
  --use-conda --cores 1 --allowed-rules refresh_drug_annotation
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Run succeeds but package versions are wrong | A venv was active in the calling shell — conda pins un-enforced | Never `source claude_science/bin/activate`; always use `scripts/run_snakemake.sh` |
| Empty log, `CalledProcessError`, no traceback | `--use-conda` omitted; interpreter path split on the space in "Biomedical Data Science" | Add `--use-conda` |
| A target is ignored / "no targets" | Target placed after `--allowed-rules`/`--forcerun` (both `nargs='+'`) | Put targets first |
| `SIGSEGV` in `sc.pp.neighbors` | numba's OpenMP pool colliding with torch's libomp after MPS init | `NUMBA_THREADING_LAYER=workqueue` must be set **before** importing scanpy — not a memory problem |
| `CyclicGraphException on nerve_immune_lead_axes` | A refreshed drug snapshot was adopted while `refresh_drug_annotation` still read the lead-axes table | Fixed 2026-08-19 — the refresh rule sources its axis list from the static 2026-08-07 snapshot |
| `There is no rule named nerve_cell_subset` | `baseline.pinned: true` — rule not defined | Expected. Use the `ds_*` equivalent |
| Snakemake plans a nerve cascade you didn't ask for | Missing `--allowed-rules` on a rule that reads pinned tables | Add `--allowed-rules <rule> --rerun-triggers mtime`; verify with `-n` |
| A compartment gate FAILS | A compartment mask regressed | **Stop.** Do not read any interaction table from that arm. Check `compartment_audit.csv` for which class contaminated it |
| LIANA rows all have empty `specificity_rank` | `.X` on the wrong scale at the LIANA call | Check `normalize_counts` for that arm; the scripts now hard-fail rather than proceed |
| MPS silent CPU fallback | Op unsupported on MPS | Check logs for `[MPS-ALERT]`; `hardware.device: cpu` as fallback |
| GSEA timeout / empty enrichment | gseapy Enrichr API unreachable | Needs internet; pipeline continues with an empty CSV |
| Census ingest exhausts memory | Full 1.29M-cell file loaded | Confirm `filter.dataset_id` and `subsample_per_donor` in the `datasets:` entry. `resources: mem_mb` is a scheduler gate, **not** a memory cap |
| Pinned artifact hash mismatch | A frozen v1.3.0 table was modified or deleted | `verify_pinned_reference` reports which; restore from git — do **not** re-freeze to absorb drift |
| `refresh_drug_annotation` writes no output | It aborted rather than emit a partial snapshot | Intended. An unresolved gene would silently become `4_no_chembl_target`, which reads as a finding rather than a network error |

---

## Key output files

**Per arm** — `results/tables/{arm}/`, `results/figures/{arm}/`:

| File | Description |
|---|---|
| `compartment_audit_gates.csv` | The 10 enforcing gates — **check first** |
| `compartment_audit.csv` | Compartment composition vs the Census oracle |
| `malignancy_confusion.csv` | CNV confusion matrix + false positives by cell type |
| `nerve_compartment_cluster_audit.csv` | Per-cluster neural purity |
| `annotation_summary.csv` | Cell-type counts and mean confidence |
| `nerve_cluster_markers_with_qc.csv` | Wilcoxon DE markers per nerve cluster |
| `nerve_enrichment_with_qc.csv` | GO BP/MF enrichment per cluster |
| `nerve_tumor_immune_interactions_with_qc.csv` | Three-way LR pairs — **interpret from `_with_qc`** |
| `nerve_tumor_interactions_with_qc.csv` | Two-way tumour→nerve LR pairs |
| `nerve_cluster_sample_purity.csv` / `_v2.csv` | Donor purity, Leiden and scANVI-v2 (full arm: 14/20 pass v2; capped: 16/20) |
| `cnv_heatmap.png`, `nerve_cells_umap.png`, `nerve_dotplot.png` | Figures |
| `01/02/05/06_*.html` | Rendered notebooks |

**Cross-arm** — `results/tables/`:

| File | Description |
|---|---|
| `nerve_immune_lead_axes_postfix.csv` | 183 rows / 144 axes significant in **both** arms, tiered by druggability |

**Committed inputs and provenance:**

| File | Description |
|---|---|
| `reference/drug_annotation/MANIFEST.json` | Active snapshot, drift table, and the annotation's limits |
| `provenance/pinned_reference_v1.3.0.json` | SHA-256 manifest of every pinned artifact |
| `provenance/{arm}/*.json` | FAIR provenance per rule output |
| `results/fair_validation_report.json` | Provenance completeness check |
| `results/snakemake_report.html` | DAG, rule stats, provenance |

---

## Picking up where you left off

```bash
claude --continue   # or: claude -c
```

`CHANGELOG.md` is the persistent lab notebook — the most recent entries record what was done,
what is open, and why. Read the last few before resuming.
