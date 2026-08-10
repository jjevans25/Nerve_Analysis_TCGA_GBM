# Post-Compartment-Fix — Open Items and Next Steps

**Written:** 2026-08-06 · **Status:** for researcher review, nothing here has been executed
**Trigger:** completion of the compartment integrity fix (`markdowns/plan_compartment_integrity_fix.md`)
**Arms rebuilt:** `gbm_cellxgene_56c4912d_full` (1,006,344 cells) + `gbm_cellxgene_56c4912d` (614,951 cells)

---

## 1. Where things stand

Both Census arms were rebuilt from `ds_scrna_annotate` down and both pass all ten
compartment-integrity gates, with the audit running in enforcing mode (so the result is
self-checking, not asserted).

| | FULL (1.0M) | CAPPED (615k) |
|---|---|---|
| compartment gates | 10/10 PASS | 10/10 PASS |
| nerve neural fraction | 95.39% (was 11.1%) | 94.76% (was 11.6%) |
| tumor malignant fraction | 92.96% (was 47.9%) | 93.47% (was 45.6%) |
| malignancy precision / recall | 0.930 / 0.894 | 0.935 / 0.875 |
| immune purity | 99.63% | 99.63% |
| worst nerve-cluster endothelial | 0.0000 (was 0.929) | 0.0000 (was 0.453) |
| nerve compartment n | 37,945 | 28,936 |
| Census neuroglial truth | 37,561 | 27,860 |
| **nerve-side S1PR1 rows** | **0** | **0** |

The strongest evidence is not any single gate. Each arm was corrected independently and each landed
within 1–4% of **its own** Census neuroglial count. Two cohorts of different sequencing depth
converging separately on their own ground truth is much harder to obtain by tuning than one arm
clearing a threshold.

### What is now trustworthy

- The **tumor** compartment (93% malignant, was 48%) and the **immune** compartment (99.6% pure).
- The **malignancy call** itself — precision 0.93, recall 0.88–0.89, AUC 0.958.
- Cross-compartment LR tables, which for the first time compare data on a **single consistent
  normalization** (defect D8 below).
- `cell_type`, the CELLxGENE author annotation, is preserved in every artifact and can be used as
  an external check at any point.

### What is not

- **Anything computed before 2026-08-05.** This includes the 40-axis shortlist, the TME notebooks'
  interpretations, and any figure or table quoting LR magnitudes.
- **Concordance against the pinned v1.3.0 reference.** See §2.5.
- **Anything about OPCs, astrocytes, or ependymal cells in the nerve compartment** — they are
  deliberately excluded (§2.4).

---

## 2. Open items requiring a decision

### 2.1 Re-derive the 40-axis shortlist

`results/tables/nerve_crosstalk_lead_targets.csv` predates the fix and was derived from tables in
which the "nerve" compartment was 59% malignant and 27% myeloid. It should be regenerated from the
new `nerve_tumor_immune_interactions_with_qc.csv`.

Worth knowing before re-deriving: **S1PR1, CXCR4 and LRP1 were never on that shortlist.** They
appear only in the raw LIANA tables, so whichever document named them as leads was produced outside
this repository. Whatever selection logic that was, it is not reproducible here and should be
restated explicitly as part of the re-derivation.

### 2.2 The neuron group is donor-dominated

`dominant_sample_fraction = 0.5032` against a 0.5 threshold, with only **7 contributing samples**.
One donor supplies roughly half of all ~4,275 neurons. This is recorded, not exempted — 15,504 rows
in the full arm's `_with_qc` table carry `batch_qc_pass = False`, mostly from this.

**Any neuron-side ligand-receptor axis is substantially one patient's biology.** The plan
anticipated a purity flag but attributed it to small *n*; the real cause is donor dominance, which
is a stronger caveat. A neuron-side result should either be reported with this stated plainly, or
restricted to axes that survive a leave-that-donor-out check.

### 2.3 The immune compartment roughly doubled — and changed in kind

1.80× on the full arm, 2.00× on the capped arm. The cause is not a threshold change: `t_cell` was
its own annotation label and was **never** in `immune_cells.source_label`, so the pre-fix immune
compartment was 96.4% myeloid / 3.1% lymphoid.

**Every pre-fix "immune" finding was a myeloid finding.** Immune-side results now change in
character, not merely in magnitude. Any prior statement of the form "immune cells signal to X"
should be re-read as "myeloid cells signal to X" and then re-tested.

### 2.4 Four cell types are excluded from the nerve compartment by decision, not by biology

`astrocyte`, `opc`, generic `neuron`, and `ependymal` are excluded from `nerve_cells.cell_types`.
Measured purity against Census, on the labels' non-malignant cells:

| label | % truly neural | in compartment? |
|---|---|---|
| inhibitory_neuron | 99.5% | yes |
| oligodendrocyte | 95.8% | yes |
| excitatory_neuron | 83.9% | yes |
| generic `neuron` | 40.6% | **no** |
| ependymal | 31.6% | **no** |
| opc | 18.9% | **no** |
| astrocyte | 10.0% | **no** |

All four remain annotated and present in every artifact. **Their absence from the interaction
tables is a masking decision and must be reported as such** — it is not evidence that they do not
participate in crosstalk. See §3.2 and §3.3 for how to bring OPC and astrocyte back properly.

### 2.5 The v1.3.0 reference should be re-derived before concordance means anything

|  | Jaccard | Spearman ρ | shared pairs |
|---|---|---|---|
| FULL vs reference | 0.4390 | 0.6182 | 1,724 |
| CAPPED vs reference | 0.4739 | 0.6620 | 1,914 |

Both moved down. **This is not a replication failure.** The reference's own nerve compartment was
built by the exact logic this work removed, it carries no author annotation (so it has never been
audited), and it was scored against differently-normalized data.

The diagnostic detail: **the two corrected arms agree with each other markedly better than either
agrees with the reference.** That pattern points at the reference as the outlier.

The reference is also structurally unreproducible — `data/processed/nerve_cells.h5ad` was deleted
by a failed job on 2026-07-21 and the retrained latent re-clusters, so the frozen cl15 split no
longer maps. Three options:

- **(a)** Keep v1.3.0 pinned and stop using it for concordance; report cross-arm agreement instead.
- **(b)** Rebuild the TCGA reference arm through the corrected pipeline as a new **v1.4.0**
  baseline, and re-freeze. Highest cost, restores a real reference.
- **(c)** Use the capped arm as the reference for the full arm and vice versa — already available,
  and arguably the more honest comparison since both are audited.

Recommendation: **(a) now, (b) when a publication-facing baseline is needed.** Do not report the
current concordance numbers as a replication result under any option.

---

## 3. Recommended next analyses

Ordered by value-per-effort. Effort is rough wall-clock on the current machine.

### 3.1 Answer the original S1PR1 question directly — *small effort, high value*

The investigation began with "is S1PR1 on microglia or on malignant cells?" and the answer was
"neither — it's endothelial," which turned out to be a compartment artifact. That question has
never actually been answered on clean data.

Do it **directly**, not through LIANA: cross-tabulate S1PR1 expression (and CXCR4, LRP1) against
`cell_type` (Census) and `cell_type_predicted`, per arm, on `malignancy_labeled.h5ad`. This needs
no re-run — the artifacts exist — and it closes the loop on the question that motivated everything.

Note the corrected tables now contain **zero** S1PR1 rows at all, which is itself informative: the
axis did not merely move compartments, it failed to reach significance anywhere once the
compartments were clean.

### 3.2 Bring OPCs back as their own audited compartment — *medium effort, high value*

OPCs are 21,460 cells in the full arm and are central to the project's actual biological question:
OPC–neuron synapses and OPC-like malignant states are the best-characterised neuron–glioma
interface in the literature. Excluding them is the largest scientific cost of the current
compartment definition.

The reason for exclusion was that the `opc` label is only 18.9% truly neural — the panel
(PDGFRA/CSPG4/SOX10/OLIG1) cannot separate normal OPCs from OPC-like malignant cells. But the CNV
caller now works (precision 0.93), so a defensible construction is available:

> high-confidence OPC = `cell_type_predicted == opc` **AND** CNV-negative at a stricter threshold
> (`is_malignant_suspected` at a lower `cnv_exclusion_sd`) **AND** high `cell_type_confidence`

Then run it as a **separate compartment** with its own audit gate, rather than folding it back into
"nerve." The `cnv_exclusion_sd` knob already exists in config for exactly this. Expect a much
smaller but far cleaner OPC population; report its size honestly.

### 3.3 Resolve the astrocyte question — *medium effort, medium value*

Census reports **347 astrocytes in 1,006,344 cells (0.03%)**; the marker panel calls ~1.5%. Both
are probably wrong. GBM study authors label AC-like cells "malignant cell" by convention, so the
oracle almost certainly under-calls; the panel over-calls because GFAP/S100B/AQP4 are shared with
AC-like malignant states.

A dedicated pass — CNV-negative **and** AQP4/SLC1A2/ALDH1L1-high **and** lacking a tumor module
score — would give a defensible normal-astrocyte set, and the disagreement with Census would then
be a *result* rather than a nuisance. This one is genuinely open science, not cleanup.

### 3.4 Validate the +7/−10 CNV contrast independently — *small effort, high assurance*

The malignancy call now rests on a chr7-gain minus chr10-loss contrast. It is well-founded (a WHO
2021 diagnostic criterion for IDH-wildtype GBM) and it emerged from the data unprompted as the
single most-up and most-down contig. But it is disease-specific by construction, and it is worth
knowing how much it is load-bearing.

Cheap checks: (i) per-donor malignant fraction — donors should vary but none should be ~0% or
~100% unless clinically expected; (ii) confirm the flagged cells show the expected *other* GBM
events (chr9p/CDKN2A loss, chr13/14 loss) rather than only +7/−10, which would indicate the score
is finding real aneuploidy rather than a two-chromosome quirk; (iii) IDH-mutant or paediatric
donors, if any, should be visible as +7/−10-negative outliers.

### 3.5 Re-audit every downstream interpretation for the D8 scale defect — *small effort, high value*

Defect D8 (below) means that **all** cross-compartment LR magnitudes computed before this fix were
between differently-normalized data. The shortlist is the obvious casualty, but the same applies to
anything in `notebooks/`, `markdowns/project_writeup.*`, and any figure quoting `lr_means` or
`specificity_rank`. Worth a systematic sweep rather than fixing artefacts one at a time as they are
noticed.

### 3.6 Decide whether the neuronal hypothesis is answerable in this cohort at all — *no compute, high value*

The cohort contains **3,401 neurons (0.34%, ~20 per donor)**, and half of the recovered neurons come
from one donor. This is a hard ceiling that no amount of pipeline correction can lift.

Options: accept that neuron-side claims will be exploratory and say so; restrict neuronal claims to
the oligodendrocyte/glial interface, which is well-powered (33,670 cells, 95.8% pure); or source a
cohort with neuronal enrichment for the neuron-specific question. This is a scoping decision worth
making explicitly before more analysis is built on ~3.4k cells.

### 3.7 Sensitivity analysis on the compartment definition — *medium effort, medium value*

The nerve compartment definition changed three times during this work, each time on measured
purity. That was the right process, but it means the top LR axes have never been tested for
stability against it. Re-running the interaction step under two or three defensible definitions
(e.g. with/without OPC; strict vs standard CNV exclusion) and reporting which axes survive all of
them would convert a judgement call into a quantified robustness statement.

---

## 4. Traps to avoid

These are all things that have already bitten during this work.

- **Always launch via `scripts/run_snakemake.sh`.** A bare `snakemake` silently reverts to the
  venv's packages and the `workflow/envs/*.yaml` pins go unenforced. Running *without*
  `--use-conda` fails differently and confusingly (unquoted interpreter path splits on the space in
  "Biomedical Data Science", producing an empty log and no traceback).
- **Targets go FIRST on the command line.** `--allowed-rules`, `--forcerun` and `--quiet` all take
  `nargs='+'` and will swallow any target placed after them.
- **`--rerun-triggers=mtime` does not fire on params or code changes.** A config-only change needs
  an explicit `--forcerun <rule>`.
- **Snakemake deletes a failed job's declared outputs.** The compartment audit writes an undeclared
  sidecar at `results/compartment_audit_snapshots/<arm>/` precisely so a *failing* gate still
  leaves its evidence. Look there first when something stops near the end.
- **Do not read the absence of OPC/astrocyte/ependymal from the interaction tables as biology.**
- **Do not report the v1.3.0 concordance as a replication result.**
- **The scANVI classifier accuracy of 0.9828 is not comparable to the previous 0.8446** — the
  compartment now has two anchoring labels where it had five, so the task is easier.

---

## 5. Defects found and fixed, for reference

Recorded here because several were silent and could recur in adjacent code.

| # | Defect | Nature |
|---|---|---|
| D1 | `score_genes` run on raw UMI counts | wrong scale, silent |
| D2 | per-cluster argmax over uncalibrated scores at res 1.0 | wrong labels |
| D3 | no marker panel for 7 lineages; VIM in the astrocyte panel | cells forced into wrong labels |
| D4 | gene order from Ensembl accession number, not coordinates; reference set drawn from the same broken annotation; no library-size normalization | CNV score ≈ noise |
| D5 | substring label matching dropped 70,881 `opc` cells | silent |
| D6 | `cell_type` (the external oracle) overwritten | destroyed ground truth |
| D7 | placeholder-on-empty exited 0 | laundered failure as success |
| **D8** | **compartments concatenated on different scales, then normalized as one — tumor transformed once, nerve/immune twice** | **silent; affected every cross-compartment LR comparison, including the original run** |

Two further classes worth noting because they are structural rather than one-off:

- **A guard written against an old configuration will fire on a correct new state.** The scANVI
  label floor demanded ≥1% OPC from a compartment that deliberately excludes OPC. Fixed by deriving
  the labels from the compartment definition so the two cannot drift apart.
- **Crash-safe checkpoints at undeclared paths survive legitimate input changes too.** A scVI
  baseline trained on the pre-fix 377,343-cell roster was loaded against the corrected 37,945-cell
  one. Now fingerprinted and discarded on mismatch.

---

## 6. Suggested sequencing

1. **§3.1** — answer the S1PR1/CXCR4/LRP1 question directly. No re-run, closes the original loop.
2. **§2.1** — re-derive the shortlist, stating the selection logic explicitly this time.
3. **§3.4** — validate the CNV contrast. Cheap, and everything downstream rests on it.
4. **§2.5** — decide the reference policy (recommend (a) now).
5. **§3.6** — decide the scope of the neuronal claim before building further on ~3.4k cells.
6. **§3.2 / §3.3** — OPC and astrocyte as their own audited compartments. This is where the
   remaining science is.
7. **§3.5 / §3.7** — sweep older interpretations; quantify robustness.

Items 1–4 are days; 6 is the substantive next piece of work.
