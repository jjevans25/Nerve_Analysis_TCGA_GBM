# Where S1PR1 actually sits in the full GBM cohort — and a compartment-integrity problem it exposed

**Question asked:** is S1PR1 expression in this cohort on microglia (making an S1P modulator act on the
immune relay) or on the malignant compartment (making it act on the tumor directly)?

**Answer: neither.** S1PR1 in this cohort is a vascular gene. The binary in the question does not hold,
and the reason the pipeline's LIANA output placed S1PR1 on a "nerve" cluster is that the cluster is not neural.

---

## 1. S1PR1 is endothelial

Detection rate (percent of cells with non-zero S1PR1), computed on the matched per-cell normalized basis:

| cell type (census author annotation) | % cells S1PR1+ |
|---|---|
| endothelial cell | **47.0** |
| mature T cell | 5.8 |
| monocyte | 4.3 |
| macrophage | 3.0 |
| malignant cell | 1.39 |
| microglial cell | 1.31 |

Endothelium carries S1PR1 at **36x** the microglial rate and **34x** the malignant rate.
Neither candidate in the original binary is a meaningful S1PR1 source.

The marker cross-check confirms this is not a normalization artifact — S1PR1 co-varies with the endothelial
panel and not with the microglial panel:

| gene | endothelial | microglial | malignant |
|---|---|---|---|
| *S1PR1* | 47.0 | 1.3 | 1.4 |
| *PECAM1* | 78.1 | 22.0 | 1.0 |
| *CLDN5* | 92.0 | 2.0 | 2.6 |
| *VWF* | 77.9 | 0.9 | 0.8 |
| *P2RY12* | 2.8 | 32.0 | 1.2 |
| *CX3CR1* | 6.4 | 49.1 | 0.9 |
| *CSF1R* | 11.8 | 66.3 | 2.2 |
| *EGFR* | 21.8 | 5.7 | 55.9 |

## 2. Why the pipeline reported S1PR1 on a nerve cluster

LIANA assigned the S1PR1 receptor role to two clusters in the nerve compartment. Both are non-neural:

- **nerve_c24** — 3,264 cells, **93% endothelial cell** by census annotation, 35.2% S1PR1+.
  This is a blood vessel cluster sitting inside the nerve compartment.
- **nerve_c34** — 418 cells, 57% malignant cell, 14.8% S1PR1+.

The S1PR1 signal the interaction tables attribute to nerve is the vasculature.

## 3. The compartment-integrity finding

Auditing every nerve cluster against the census's own author annotation:

| true dominant identity | cells | share of "nerve" compartment |
|---|---|---|
| malignant | 231,639 | 61% |
| myeloid | 101,477 | 27% |
| **neural** | **33,354** | **8.8%** |
| mixed | 4,181 | 1% |
| lymphoid | 3,428 | 1% |
| endothelial | 3,264 | 1% |

Of 40 nerve clusters, only **5** are dominantly neural (nerve_c2, nerve_c25, nerve_c26, nerve_c28, nerve_c37);
23 are malignant, 8 myeloid, 1 endothelial, 1 lymphoid, 2 mixed.

**Root cause.** `workflow/scripts/nerve_cell_subset.py` builds the nerve mask from the pipeline's own
`cell_type_predicted` label while excluding cells flagged `is_malignant`. Both inputs fail here:
`cell_type_predicted` labels ~85% of the compartment "astrocyte", and the pipeline's CNV-based malignancy
caller did not flag the cells the census annotates as malignant. The exclusion therefore removed almost nothing.

**Scope.** Both census arms are affected almost identically:

| cohort | nerve cells | % neural | % malignant | % myeloid |
|---|---|---|---|---|
| full census | 377,343 | 11.1 | 58.9 | 27.2 |
| capped census | 270,520 | 11.6 | 61.6 | 24.4 |
| reference (GDC, 17 samples) | 106,603 | not testable | — | — |

(The percentages in this table use a broader neural definition than the per-cluster verdict in §3, which
requires a cluster to be *dominantly* neural; both are reported rather than reconciled to one number.)
The reference cohort carries no independent author annotation — only its own marker-score and model-predicted
labels — so the same cross-check cannot be constructed there. Whether it has the same problem is unknown, not ruled out.

## 4. What this means for the lead-target work

1. **S1PR1 does not survive as an immune-relay lead.** Its expression is vascular. An S1P modulator
   (fingolimod, siponimod, ozanimod) acting on this axis in GBM would be acting on endothelium — plausibly
   interesting for blood–brain-barrier permeability or vascular normalization, but that is a different
   hypothesis with different endpoints, not the neuro-immune relay the shortlist claimed.
2. **Every nerve-side axis in the full-cohort shortlist needs re-derivation.** With 61% of the compartment
   malignant and 27% myeloid, "nerve–tumor" interactions in these tables are substantially tumor–tumor and
   myeloid–tumor. This does not invalidate the axes as interactions; it invalidates the *compartment labels*
   on them, which is what the therapeutic rationale rested on.
3. **The 5 genuinely neural clusters are the honest starting point** for a nerve-side re-analysis:
   nerve_c2, nerve_c25, nerve_c26, nerve_c28, nerve_c37, totalling 33,354 cells.
4. **Not affected:** immune–tumor axes, which do not depend on the nerve mask.

## Caveats

- The census author annotation (`cell_type`) is treated as ground truth here. It is an external annotation,
  not a re-derivation from these counts.
- All full-cohort numbers inherit the environment-provenance caveat recorded in
  `markdowns/status_full_cohort_run_20260802.md`: the declared conda environment was not the one that
  executed, and the arm completed only after several relaunches.
- Detection rate depends on sequencing depth; endothelial cells are not systematically deeper here, but
  the percentages are not depth-corrected.

## Files

- `s1pr1_compartment_localization.png` — 4-panel figure
- `full_cohort_nerve_compartment_audit.csv` — per-cluster composition, S1PR1 stats, verdict (40 rows)
- `s1pr1_expression_by_celltype.csv` — S1PR1 by census cell type (17 rows)
- `s1pr1_panel_expression_matrix.csv` — gene x group detection/mean matrix (61 rows)
