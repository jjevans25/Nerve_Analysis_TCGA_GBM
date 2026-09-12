# The Drug Target That Wasn't There

*What a million glioblastoma cells say about nerve–immune crosstalk — and about
the question I started with.*

**Draft — science blog post, written 2026-09-11.** Companion to the technical
walkthrough in `GBM_TME_Crosstalk_Analysis.md`. Every number here resolves to an
artifact under `results/tables/gbm_cellxgene_56c4912d_full/`.

---

## A question with two answers, both wrong

I started this project with a narrow, practical question. Glioblastoma is wired
into the brain — neurons form functional synapses onto glioma cells, and that
electrical coupling drives proliferation. If you wanted to interrupt that, one
appealing handle was **S1PR1**, a sphingosine-1-phosphate receptor with drugs
already on the market. My pipeline's output put S1PR1 on a nerve cluster,
signalling to microglia.

So: is S1PR1 here on microglia, making an S1P modulator act on the immune relay?
Or on the malignant compartment, making it act on the tumor directly?

The answer was **neither**. S1PR1 in this cohort is a vascular gene:

| cell type (independent annotation) | % of cells S1PR1⁺ |
|---|---|
| **endothelial cell** | **47.0** |
| mature T cell | 5.8 |
| monocyte | 4.3 |
| macrophage | 3.0 |
| malignant cell | 1.4 |
| microglial cell | 1.3 |

Endothelium carries S1PR1 at 36× the microglial rate. It co-varies with *PECAM1*,
*CLDN5* and *VWF* — the endothelial panel — and not with *P2RY12* or *CX3CR1*.
This was not a normalization artifact. It was a **label** problem. The "nerve"
cluster S1PR1 sat on was `nerve_c24`: 3,264 cells, **93% endothelial** by the
independent annotation, 35.2% of them S1PR1⁺. The drug target was real. The
compartment I had put it in was not.

That finding is the reason this post exists, and it's worth being blunt about what
it means. The interesting result here is not an axis I found. It's an axis I lost.

---

## What was actually broken

Chasing the S1PR1 anomaly surfaced something systemic. The compartment
definitions — which cells count as "nerve," "tumor," "immune" — were wrong, and
wrong in ways that every downstream number silently inherited.

| | before | after |
|---|---|---|
| nerve compartment, fraction truly neural | **11.1 %** | **95.4 %** |
| tumor compartment, fraction truly malignant | 47.9 % | 93.0 % |
| worst nerve cluster's endothelial fraction | 0.929 | 0.000 |

An 11%-neural "nerve" compartment does not produce slightly noisy conclusions
about nerve biology. It produces conclusions about something else entirely,
expressed in the vocabulary of nerve biology.

The immune side had a quieter, worse version of the same problem. The immune
compartment was built by selecting cells with one label — `"microglia"` — so T
cells, NK cells, B/plasma cells and neutrophils never entered any interaction
analysis at all. **Every pre-fix "immune" finding was a myeloid finding.** Fixing
the definition to a list of six lineages grew the immune compartment by **1.80×**,
to 593,264 cells.

None of this announced itself. Nothing crashed. Every table was well-formed, and
every conclusion drawn from them was fluent and plausible.

---

## The fix, and how I know it worked

The cohort comes from the CZ CELLxGENE Census — 170 donors, 1,020,902 cells
ingested, 1,006,344 surviving QC. The Census ships an author cell-type annotation
harmonized to the Cell Ontology, and the single most useful decision in this
project was to **pull that annotation and never let it touch the pipeline.**

It is held back as an oracle. The pipeline labels cells from marker panels and a
CNV-based malignancy call, entirely independently, and then the two are compared
in a rule that **fails the build** if they disagree too much. Ten gates, all
enforcing:

| gate | observed | threshold |
|---|---|---|
| nerve neural fraction | 0.9539 | ≥ 0.80 |
| tumor malignant fraction | 0.9296 | ≥ 0.85 |
| malignancy precision | 0.9296 | ≥ 0.85 |
| malignancy recall | 0.8940 | ≥ 0.80 |
| immune purity | 0.9963 | ≥ 0.95 |
| worst nerve cluster endothelial | 0.0000 | ≤ 0.20 |

The malignancy call lands at F1 0.911 — 298,300 true positives against 22,578
false positives and 35,381 false negatives. Its largest error class is
**OPC-like cells** (10,264), which is exactly where a CNV-based caller should
struggle: OPC-like malignant states are the hardest to separate from real OPCs.

The evidence I find most convincing isn't any single gate. The analysis runs as
two arms — the full cohort and a 5,000-cells-per-donor subsample. Each was
corrected independently, and each landed within 1–4% of **its own** oracle count
(37,945 vs 37,561 neuroglial cells; 28,936 vs 27,860). Two cohorts of different
sequencing depth converging separately on their own ground truth is a much harder
thing to fake than one arm clearing a threshold.

![The corrected nerve compartment: 37,945 cells in 23 Leiden clusters (left), and
what they are (right).](../results/figures/gbm_cellxgene_56c4912d_full/nerve_cells_umap.png)

*The corrected nerve compartment — and a preview of its central limitation. Left:
37,945 cells, 23 clusters. Right: what they actually are. The green mass is
oligodendrocytes; the two small islands are excitatory and inhibitory neurons,
sitting cleanly apart from the glia rather than smeared into them, which is what a
compartment looks like when the labels are right. It also shows the asymmetry that
governs everything below — the glial side is an ocean and the neuron side is two
islands.*

---

## The pipeline re-derives biology it was never told about

With clean compartments, ligand–receptor inference ran across
tumor (320,878 cells) × nerve (37,945, in 21 groups) × immune (593,264, in 5
subtypes), producing 44,139 tested interactions, 2,673 of them significant. An
axis counts as a **lead** only if it clears significance in *both* arms
independently: **144 axes**, of which 34 have an approved drug available today.

Here is the part that made me trust the rest. Ranked purely by statistics, with no
literature input anywhere in the pipeline, the top nerve→tumor axes on the neuron
side are:

| rank | axis |
|---|---|
| 4 | **PTN\|PTPRZ1** |
| 5 | **NCAM1\|PTPRZ1** |
| 6 | **CNTN1\|PTPRZ1** |
| 7 | **NLGN1\|NRXN1** |

These are the four best-characterised neuron–glioma interactions in the
literature. PTN–PTPRZ1 is the canonical glioma stem cell axis. NLGN1–NRXN1 is the
synaptic adhesion pair at the center of the neuron–glioma synapse story. The
pipeline had no way to know any of that. It recovered them from expression
statistics alone.

That's a positive control I didn't design and couldn't have gamed, and it is the
main reason I believe the unfamiliar results sitting next to them.

It gets more pointed if you look at the nerve↔tumor interface on its own. Of the
11,280 interactions tested there, the **top 25 by effect size collapse onto just
nine distinct axes — and 19 of those 25 rows share a single receptor: PTPRZ1.**
Four different ligands are hitting it: `NCAM1`, `PTN`, `CNTN1` and `MDK`. All four
survive in both cohort arms.

![Top 25 nerve↔tumor ligand–receptor pairs. Each panel is one signalling source;
dot size is specificity, colour is effect
magnitude.](../results/figures/gbm_cellxgene_56c4912d_full/nerve_tumor_dotplot.png)

*The convergence, drawn. Each panel is a source — the leftmost is the malignant
compartment, the rest are individual nerve clusters and the pooled neuron group.
Read along the `NCAM1 -> PTPRZ1` row: **eleven separate nerve clusters** carry it,
and every one of them points at the same target, `malignant`. `PTN -> PTPRZ1`
accounts for six more rows. Eighteen of the 25 run nerve→tumor; the seven in the
leftmost panel are the tumor answering back.*

PTPRZ1 is the receptor tyrosine phosphatase that marks glioma stem-like cells, and
a funnel this narrow is the kind of thing that is either a real feature of the
tissue or an artifact of how the question was asked. I lean toward real, for one
reason: the eleven clusters carrying `NCAM1 -> PTPRZ1` were clustered
independently, and nothing in the pipeline encourages them to agree on a target.
But it is a hypothesis the data generated, not one it tested, and it is the
first thing I would want a wet-lab collaborator to push back on.

### What's next to them

The glial side — better powered, and where the novel material is — is topped by
druggable axes:

| rank | axis | immune partners | status |
|---|---|---|---|
| 1 | **APP\|CD74** | dendritic, microglia, NK, T, TAM | aducanumab / donanemab / lecanemab; milatuzumab |
| 10 | **FGFR2\|TIMP1** | dendritic, microglia, TAM | **9 named glioma trials** (pemigatinib, regorafenib) |
| 14 | EGFR\|S100A4 | dendritic, NK, T, TAM | large approved EGFR agent list |
| 19 | **CXCR4\|HMGB1** | dendritic, microglia, NK, T, TAM | plerixafor, 4 glioma trials |
| 23 | CD74_CXCR4\|MIF | dendritic, microglia, NK, T, TAM | milatuzumab, plerixafor |

`APP|CD74` ranks **first on both the neuron side and the glial side**,
independently. The immune–nerve interface both dominates the tested space (30,622
of 44,139 tested pairs) and survives it at the highest rate (2,118 of 2,673
significant).

![Significant ligand–receptor pairs per compartment pair and
direction.](../results/figures/gbm_cellxgene_56c4912d_full/nerve_tumor_immune_sig_heatmap.png)

*Where the crosstalk actually is. Each cell counts significant LR pairs
(`magnitude_rank < 0.05`) for one compartment pair in one direction. The
nerve↔immune axis carries 2,118 of the 2,673 significant pairs — 1,170 nerve→immune
and 948 immune→nerve — against 459 for nerve↔tumor and 96 for immune↔tumor. The
question this project set out to ask turned out to be the one the data had the most
to say about, which is not something you get to count on.*

---

## Four things I'd want a reader to hold onto

**The neuron side is substantially one patient.** The pooled neuron group is 4,275
cells from **7 donors**, one contributing 50.3%. Neurons are brutally rare in
dissociated tumor tissue, and no amount of pipeline correction lifts that ceiling.
Neuron-side claims — including the beautiful canonical axes above — deserve a
leave-that-donor-out check before anyone leans on them. The glial side, at 33,670
cells, does not have this problem.

**Purity and donor diversity are different failure modes, and only one has a
flag.** The 95.4% figure is an average. Five of 23 nerve clusters sit below the
0.80 neural bar — 1,677 cells, 4.4% of the compartment — and four of them reach
the interaction tables, carrying 12.8% of all rows. The worst, cluster 12, is
**83.5% malignant** against the oracle and **passes** the automated batch-QC check,
because 7 donors at a 0.44 dominant fraction is a perfectly healthy donor spread.
Concretely: 27 of the 144 lead axes draw some support from a purity-failing
cluster, and 3 rest on one entirely (`APOE|SCARB1`, `BCAN|EGFR`, `IGSF11|VSIR`).
`BCAN|EGFR` reaches tier 1 on a **single** interaction row.

**A p-value of 0 is not a p-value of 0.** These are permutation p-values at 1,000
shuffles, so the floor is 0.001. 59.9% of all rows — and 98.3% under the
dashboard's default filters — sit at exactly that floor, tied. The column cannot
rank anything. Ranking comes from effect magnitude.

**Absent is not the same as absent.** Astrocytes, OPCs and ependymal cells were
deliberately masked out of the nerve compartment. They are missing from these
tables because of a decision, not a finding — and OPCs in particular sit at the
best-characterised neuron–glioma interface there is. That exclusion is the
biggest open cost in this work.

---

## So, S1PR1

Zero rows. Across both arms, after the fix, S1PR1 does not appear in the
interaction tables at all. The axis that motivated the entire project failed to
reach significance anywhere once the compartments it depended on were correct.

I think that's the most useful thing here. A pipeline that only ever confirms is a
pipeline that can't tell you anything. This one was built so that its own founding
premise could fail, and then that premise failed — and the same machinery that
killed it handed back 144 axes that survive a test it was never able to pass.

---

## Look at it yourself

Everything above is browsable. The dashboard lets you filter by compartment
interface, nerve cell type, effect magnitude and druggability tier, search ligands
and receptors by name, and see the QC flag on every row rather than a cleaned-up
version of it.

→ **[explore the nerve–immune crosstalk dashboard]** *(link to be added)*

Three things worth doing once you're in there, in order:

1. **Filter to tier 1, then turn *off* "hide QC-failing rows."** The difference
   between those two views is the honest width of the result. Everything this post
   claims survives that toggle; not everything in the table does.
2. **Find `BCAN|EGFR`.** It sits in tier 1 — approved-drug-available, on EGFR's long
   agent list — on the strength of a single interaction row, from a cluster the
   external annotation calls 83.5% malignant. It is the cleanest example I have of
   why a druggability tier is a starting point for a conversation and not a result.
3. **Compare the two arms.** An axis present in one and absent in the other is not
   a finding. Cross-arm agreement is the only replication evidence here.

---

## The data

The analysis tables are in the repository — 92 files, both cohort arms: every
tested interaction with its QC flags, the compartment audit gates, per-cluster
purity, the lead-axes shortlist, and the rendered notebooks. `results/README.md`
documents what each file is and the four ways these tables are easy to misread.

Two things are deliberately **not** published. The retired v1.3.0 reference
cohort's tables are withheld: that cohort's nerve compartment measured 11% neural,
its nerve claims are void, and it can't be corrected, so putting it next to the
corrected tables with nothing in the filename to warn anyone seemed worse than
leaving it out. And the model checkpoints are too large for git.

One caveat on reproducibility, stated because the alternative is letting you
discover it: of 797 provenance records, **one** artifact — the full cohort's scVI
latent — was built under a torch version the environment spec doesn't pin (2.11.0
against a pinned 2.12.0), because for a period the declared conda environments
weren't the ones that actually executed. That's fixed, and a build gate now checks
every record against every pin. The one artifact was kept rather than rebuilt, as a
decision on the record: re-deriving it means retraining on a million cells and
renumbering every cluster in this post, and nobody knows whether it would change
anything. The second arm's latent is clean, and since the replication claim here is
cross-arm agreement, the result doesn't rest on that file. The engineering companion
post tells that story properly.

---

*Cohort: CELLxGENE Census GBM 10x, dataset `56c4912d`, 170 donors, 1,020,902 cells
ingested, 1,006,344 post-QC. Pipeline: Snakemake + scvi-tools + scanpy + LIANA +
marimo, on an Apple M4 Max. Compartment labels audited against the Census author
annotation, held back from the pipeline and used only as an oracle. The engineering
companion post covers how it was built and what broke along the way. Remaining open
items: `markdowns/post_compartment_fix_next_steps.md`.*
