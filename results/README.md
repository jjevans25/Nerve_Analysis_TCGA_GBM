# `results/` — what is published here, and what is deliberately not

`results/` is a generated directory and is gitignored as a whole. A curated slice is
force-added so that the analysis is checkable without running the pipeline: 92 files,
~49 MB. Everything here was produced by the Snakemake workflow in this repository and
carries a matching record under `provenance/`.

Two cohort arms are published, and comparing them is the point — an axis is only
treated as a lead if it clears significance in **both**:

| arm | cells | note |
|---|---|---|
| `gbm_cellxgene_56c4912d_full` | 1,006,344 post-QC | the full cohort; the arm quoted in the write-ups |
| `gbm_cellxgene_56c4912d` | 614,951 | 5,000 cells per donor, retained as an independent depth |

## Start here

| file | what it is |
|---|---|
| `tables/nerve_immune_lead_axes_postfix.csv` | The headline result. 183 rows / 144 unique axes significant in both arms, tiered by druggability. Read `selection_arm` and `qc_filter_applied` before comparing two rows — the two halves are selected under different rules, on purpose. |
| `tables/<arm>/compartment_audit_gates.csv` | The ten enforcing gates that license every other number here. If these do not pass, nothing downstream is trustworthy. |
| `tables/<arm>/nerve_compartment_cluster_audit.csv` | Per-cluster purity. Read this alongside the gates: the 95.4 % aggregate hides five clusters below the 0.80 neural bar, and `batch_qc_pass` does not cover them. |
| `tables/<arm>/nerve_tumor_immune_interactions_with_qc.csv` | Every tested ligand–receptor interaction (44,139 on the full arm), QC flags joined on. Rows that fail QC are **retained and flagged**, never dropped. |
| `figures/<arm>/05_census_nerve_immune_explorer.html` | The rendered dashboard. Note this is a static export — the filters are frozen. |
| `fair_validation_report.json` | Provenance audit: 797 records checked against the declared environment pins. |

## Reading these honestly

- **`cellphone_pvals = 0` does not mean p = 0.** Permutation p-values at `n_perms=1000`;
  the floor is 0.001. 59.9 % of rows sit exactly there, tied, so the column cannot rank
  anything. Rank on `magnitude_rank` / `lrscore`.
- **QC failures are flagged, not removed.** `batch_qc_pass = False` marks donor-dominated
  groups. They are kept because the dominance test cannot distinguish real biology
  preserved in one donor's tissue from a patient-driven artifact.
- **Purity and donor diversity are different failure modes.** `batch_qc_pass` tests only
  the second. Cluster 12 is 83.5 % malignant against the external annotation and passes it.
- **Masked labels are a decision, not a finding.** Astrocytes, OPCs and ependymal cells
  were excluded from the nerve compartment; their absence from the interaction tables says
  nothing about whether they participate in crosstalk.

## Deliberately not published

**The retired v1.3.0 reference cohort.** Its ~40 root-level tables are excluded. That
cohort's nerve compartment measured **11 % neural / 59 % malignant** — every nerve claim
it carries is void, and it cannot be corrected, because it is pinned and its
`nerve_cells.h5ad` was deleted. Publishing those tables next to the corrected ones would
hand a reader void numbers with nothing in the filename to warn them. They remain in git
history and in a checksum-verified archive.

**Model checkpoints** (`models/`, 740 MB) and the **pre-fix snapshot**
(`prephase5_snapshot_*`, 153 MB) — too large for git, and earmarked for a Zenodo deposit.

**Per-sample QC metrics** (680 files, 51 MB) and **pre-QC-join duplicates** — the latter
are strictly redundant: each `*_with_qc.csv` has the same rows and a superset of the
columns of its plain sibling.

## Regenerating

Always launch through `scripts/run_snakemake.sh` — a bare `snakemake` inherits an active
venv from the calling shell and silently un-enforces every pin in `workflow/envs/*.yaml`.
Targets go first; several flags take `nargs='+'` and will swallow them.

```bash
scripts/run_snakemake.sh results/tables/nerve_immune_lead_axes_postfix.csv \
  --use-conda --cores 1 --rerun-triggers mtime \
  --forcerun nerve_immune_lead_axes --allowed-rules nerve_immune_lead_axes
```

Pin `--allowed-rules`. `--forcerun` alone re-evaluates inputs and can reach
`refresh_drug_annotation`, which performs network I/O and overwrites the pinned drug
annotation snapshot.

Note that the per-sample intermediates under `data/processed/<arm>/` were reclaimed, so a
bare `rule all` now schedules a full cohort rebuild (365 jobs) rather than a no-op.
