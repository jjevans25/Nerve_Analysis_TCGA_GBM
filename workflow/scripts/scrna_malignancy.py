"""CNV-based tumor/normal separation; labels each cell with is_malignant.

Biological question: which cells in this tumor carry the large-scale copy-number
alterations that define malignancy (in IDH-wildtype GBM, canonically chromosome
7 gain and chromosome 10 loss), as opposed to the normal glia, immune and
vascular cells of the surrounding microenvironment?

The signal is a *run* of co-elevated or co-depressed expression across
physically adjacent loci, so it exists only when genes sit in genomic order and
windows never straddle a chromosome boundary.
"""

import os
import sys

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

matplotlib.use("Agg")

sys.path.insert(0, "workflow/scripts")
from counts_utils import is_log1p_scale, library_sizes
from fair_utils import H5AD_COMPRESSION, log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
os.environ["PYTHONHASHSEED"] = str(snakemake.params.random_seed)

log_transformation(log, "scrna_malignancy", f"Loading {snakemake.input.h5ad}")
adata = ad.read_h5ad(snakemake.input.h5ad)

# ---------------------------------------------------------------------------
# Step 1: Real chromosomal coordinates from the pinned Ensembl GTF
#
# This previously ordered genes by the numeric suffix of the Ensembl accession
# (ENSG00000141510 -> 141510). That is an accession counter, not a coordinate:
# it interleaves chromosomes arbitrarily, so every sliding window averaged
# unrelated loci and the CNV score was mostly noise (18% recall / 48% precision
# against the Census annotation). The Census API does not expose coordinates at
# all, which is why the proxy existed; `download_gene_positions` now supplies
# them properly.
# ---------------------------------------------------------------------------
positions = pd.read_csv(snakemake.input.gene_positions, sep="\t", dtype={"chromosome": str})
positions = positions.drop_duplicates(subset="ensembl_id", keep="first").set_index("ensembl_id")

gene_ids = (
    adata.var["feature_id"].astype(str)
    if "feature_id" in adata.var.columns
    else pd.Series(adata.var_names.astype(str), index=adata.var_names)
)
gene_ids = gene_ids.str.split(".").str[0]        # cohort IDs may be versioned
adata.var["ensembl_unversioned"] = gene_ids.values
joined = positions.reindex(gene_ids.values)
adata.var["chromosome"] = joined["chromosome"].values
adata.var["start"] = joined["start"].values
adata.var["genome_order"] = joined["genome_order"].values

placed = adata.var["genome_order"].notna().to_numpy()
n_placed = int(placed.sum())
log_transformation(log, "scrna_malignancy",
    f"Mapped {n_placed}/{adata.n_vars} genes to real chromosomal coordinates "
    f"({100 * n_placed / adata.n_vars:.1f}%) across "
    f"{adata.var.loc[placed, 'chromosome'].nunique()} contigs")

min_frac = float(snakemake.params.min_genes_placed_fraction)
if n_placed / adata.n_vars < min_frac:
    raise RuntimeError(
        f"[FAIR-ALERT] only {n_placed}/{adata.n_vars} genes "
        f"({n_placed / adata.n_vars:.1%}) carry coordinates, below the required "
        f"{min_frac:.0%}. CNV inference on an unplaced gene set is the exact "
        "failure this rule was rebuilt to remove — refusing to proceed."
    )

# Unplaced genes cannot participate in a positional window; drop them from the
# CNV pass rather than parking them at order -1 where they formed a phantom
# 'chromosome' at the head of the axis.
mt_mask = adata.var["chromosome"].fillna("").astype(str).str.upper().eq("MT").to_numpy()
cnv_gene_mask = placed & ~mt_mask          # MT has no meaningful copy number here
cnv_idx = np.flatnonzero(cnv_gene_mask)
order_within = np.argsort(adata.var["genome_order"].to_numpy()[cnv_idx], kind="stable")
gene_order_idx = cnv_idx[order_within]     # positions into adata.var, in genome order

ordered_contigs = adata.var["chromosome"].to_numpy()[gene_order_idx].astype(str)
# Half-open [start, stop) index ranges per contig, so no window straddles a
# chromosome boundary — a window spanning the chr1/chr2 junction would average
# loci that are megabases apart in reality.
contig_bounds: list[tuple[int, int]] = []
_edges = np.flatnonzero(ordered_contigs[1:] != ordered_contigs[:-1]) + 1
for lo_, hi_ in zip(np.r_[0, _edges], np.r_[_edges, len(ordered_contigs)]):
    contig_bounds.append((int(lo_), int(hi_)))

log_transformation(log, "scrna_malignancy",
    f"CNV gene axis: {len(gene_order_idx)} genes over {len(contig_bounds)} contigs "
    f"(dropped {adata.n_vars - len(gene_order_idx)}: unplaced or MT)")

# ---------------------------------------------------------------------------
# Step 1b: the diagnostic aneuploidy contrast
#
# A genome-wide "spread" score asks how far a cell's profile deviates from the
# reference *anywhere*. That conflates copy number with cell-type identity: an
# oligodendrocyte compared against an immune reference deviates across whole
# chromosomes purely because it runs a different expression programme. Measured
# on this cohort, normal neural cells scored HIGHER than malignant ones by that
# metric (0.063 vs 0.057), which is what capped it at AUC 0.86.
#
# Contrasting a gain against a loss removes that confound, because both terms
# are measured within the same cell and the cell-type baseline largely cancels.
# For IDH-wildtype glioblastoma the contrast is not arbitrary: +7/-10 is a WHO
# 2021 diagnostic criterion, and it emerges unprompted from this data as the
# single most-up and single most-down contig in Census-malignant cells
# (chr7 +0.043, chr10 -0.038). AUC 0.86 -> 0.96.
#
# Configurable, because it is disease-specific: a cohort without a canonical
# signature sets both lists empty and falls back to the genome-wide spread.
gain_contigs = {str(c) for c in snakemake.params.cnv_gain_contigs}
loss_contigs = {str(c) for c in snakemake.params.cnv_loss_contigs}
gain_cols = np.flatnonzero(np.isin(ordered_contigs, list(gain_contigs))) if gain_contigs else np.array([], int)
loss_cols = np.flatnonzero(np.isin(ordered_contigs, list(loss_contigs))) if loss_contigs else np.array([], int)
use_contrast = gain_cols.size > 0 and loss_cols.size > 0

if use_contrast:
    log_transformation(log, "scrna_malignancy",
        f"Aneuploidy contrast: gain {sorted(gain_contigs)} ({gain_cols.size} genes) minus "
        f"loss {sorted(loss_contigs)} ({loss_cols.size} genes)")
else:
    missing = (gain_contigs | loss_contigs) - set(ordered_contigs)
    log_transformation(log, "scrna_malignancy",
        f"No usable gain/loss contrast (configured gain={sorted(gain_contigs)}, "
        f"loss={sorted(loss_contigs)}; absent from the gene axis: {sorted(missing)}) — "
        "falling back to the genome-wide spread score, which cannot separate "
        "cell-type identity from copy number nearly as well.",
        status="WARNING")

# ---------------------------------------------------------------------------
# Step 2: Reference ("normal") cells
#
# Previously `cell_type_predicted in {t_cell, endothelial}` — drawn from the very
# annotation this rule is meant to be independent of, and in the capped arm
# `endothelial` was never assigned at all, silently degrading the reference to T
# cells alone. Now: confidently-called immune cells, which are the cleanest
# non-malignant population in a GBM dissociation (the compartment audit measures
# them at 99.5% pure), with a hard non-empty assertion.
# ---------------------------------------------------------------------------
reference_labels = {str(s).lower().strip() for s in snakemake.params.reference_labels}
_pred = adata.obs["cell_type_predicted"].astype(str).str.lower().str.strip()
is_ref = _pred.isin(reference_labels).to_numpy()

conf_q = float(snakemake.params.reference_confidence_quantile)
if conf_q > 0 and "cell_type_confidence" in adata.obs.columns and is_ref.any():
    conf = adata.obs["cell_type_confidence"].to_numpy(dtype=np.float64)
    cutoff = float(np.quantile(conf[is_ref], conf_q))
    is_ref = is_ref & (conf >= cutoff)
    log_transformation(log, "scrna_malignancy",
        f"Reference restricted to the top {100 * (1 - conf_q):.0f}% most confident "
        f"immune calls (confidence >= {cutoff:.4f})")

adata.obs["is_reference"] = is_ref
n_ref = int(is_ref.sum())
min_ref = int(snakemake.params.min_reference_cells)
log_transformation(log, "scrna_malignancy",
    f"Reference cells ({sorted(reference_labels)}): {n_ref} "
    f"({100 * n_ref / adata.n_obs:.1f}% of cohort)")

if n_ref < min_ref:
    raise RuntimeError(
        f"[FAIR-ALERT] only {n_ref} reference cells, below the required {min_ref}. "
        "A CNV baseline built from too few normal cells produces a threshold that "
        "is noise, which is how this rule reached 18% recall. Check that "
        "scrna_malignancy.reference_labels match the panel names emitted by "
        f"scrna_annotate — observed labels: {sorted(set(_pred.unique()))}."
    )

# ---------------------------------------------------------------------------
# Step 3: CNV scoring via sliding-window expression smoothing.
# Streamed over cells in row-blocks to bound peak memory — every operation
# here is row-wise and independent per cell; the only cross-cell quantity is
# the per-gene reference mean, computed once. Chunking is numerically identical
# to a whole-matrix pass but never materializes the full dense matrix (which is
# ~59 GB for a 615k×24k float32 cohort, and OOM-kills at whole-matrix scale).
# Cells with high CNV variance relative to reference → malignant.
# ---------------------------------------------------------------------------
n_obs = adata.n_obs
n_genes_ord = len(gene_order_idx)
ref_mask = adata.obs["is_reference"].values
chunk_size = int(snakemake.params.cnv_chunk_size)
clip = float(snakemake.params.cnv_clip)

if not sp.isspmatrix_csr(adata.X):
    adata.X = sp.csr_matrix(adata.X)

# Library-size normalization was missing entirely: log1p was applied straight to
# raw UMIs, so a cell with 20k counts looked systematically "amplified" relative
# to one with 2k regardless of copy number. Depth is a per-cell technical
# property; it has to be divided out before any expression is read as dosage.
_is_log1p_input = is_log1p_scale(adata.X)
lib = None if _is_log1p_input else library_sizes(adata.X)
if _is_log1p_input:
    log_transformation(log, "scrna_malignancy",
        "Input .X is already log1p-scaled; skipping library-size normalization",
        status="WARNING")

# Sliding-window size, applied WITHIN a contig.
window = min(int(snakemake.params.cnv_window), n_genes_ord)


def _load_ordered_block(lo: int, hi: int) -> np.ndarray:
    """Densify cells [lo:hi) in genome order, library-normalized and log1p'd."""
    blk = adata.X[lo:hi][:, gene_order_idx].toarray().astype(np.float32, copy=False)
    if lib is not None:
        scale = np.where(lib[lo:hi] > 0, 1e4 / np.maximum(lib[lo:hi], 1.0), 0.0)
        blk *= scale.astype(np.float32)[:, None]
        np.log1p(blk, out=blk)
    return blk


def _smooth_within_contigs(blk: np.ndarray) -> np.ndarray:
    """Running mean along the genome, restarting at every chromosome boundary.

    A window that straddles a contig junction averages loci megabases apart and
    manufactures a copy-number 'edge' where none exists.
    """
    out = np.empty_like(blk)
    for lo_, hi_ in contig_bounds:
        seg = blk[:, lo_:hi_]
        n_seg = hi_ - lo_
        w = min(window, n_seg)
        pad_l, pad_r = w // 2, w - w // 2 - 1
        padded = np.pad(seg, ((0, 0), (pad_l, pad_r)), mode="edge")
        cs = np.empty((padded.shape[0], padded.shape[1] + 1), dtype=np.float32)
        cs[:, 0] = 0.0
        np.cumsum(padded, axis=1, out=cs[:, 1:])
        out[:, lo_:hi_] = (cs[:, w:w + n_seg] - cs[:, :n_seg]) / w
    return out


# Per-gene reference baseline (mean of log-normalized reference cells).
# Streamed in row-blocks like the main loop below: the reference set is not
# small (55,477 cells in the last Census run = 5.3 GB dense in one allocation),
# and materializing it whole was the last unchunked densification in this script.
ref_pos = np.flatnonzero(ref_mask)
ref_sum = np.zeros((1, n_genes_ord), dtype=np.float64)
for lo in range(0, ref_pos.size, chunk_size):
    sel = ref_pos[lo:lo + chunk_size]
    ref_blk = adata.X[sel][:, gene_order_idx].toarray().astype(np.float32, copy=False)
    if lib is not None:
        _s = np.where(lib[sel] > 0, 1e4 / np.maximum(lib[sel], 1.0), 0.0)
        ref_blk *= _s.astype(np.float32)[:, None]
        np.log1p(ref_blk, out=ref_blk)
    ref_sum += ref_blk.sum(axis=0, keepdims=True)
    del ref_blk
ref_mean = (ref_sum / ref_pos.size).astype(np.float32)
del ref_sum, ref_pos

# Pre-select heatmap cells before the loop so only 500 smoothed rows are kept.
n_plot = min(500, n_obs)
rng = np.random.default_rng(snakemake.params.random_seed)
plot_idx = np.sort(rng.choice(n_obs, size=n_plot, replace=False))
plot_buffer = np.empty((n_plot, n_genes_ord), dtype=np.float32)

# Stream CNV scores in row-blocks.
cnv_scores = np.empty(n_obs, dtype=np.float32)
cnv_spread = np.empty(n_obs, dtype=np.float32)
for lo in range(0, n_obs, chunk_size):
    hi = min(lo + chunk_size, n_obs)
    blk = _load_ordered_block(lo, hi)
    blk -= ref_mean
    # Clip before smoothing: a handful of extreme genes (a burst-expressed
    # marker, a dropout artifact) would otherwise dominate the window mean and
    # be read as a copy-number event. Standard infercnv practice.
    np.clip(blk, -clip, clip, out=blk)

    smoothed = _smooth_within_contigs(blk)
    del blk

    # Re-center each cell on its own median across windows. Without this, a
    # cell's global expression offset (dissociation stress, cell size, cycle
    # phase) reads as genome-wide amplification.
    smoothed -= np.median(smoothed, axis=1, keepdims=True)

    # Spread of the smoothed profile: retained as a diagnostic, and used as the
    # score when no gain/loss contrast is configured.
    cnv_spread[lo:hi] = smoothed.std(axis=1)
    if use_contrast:
        cnv_scores[lo:hi] = (
            smoothed[:, gain_cols].mean(axis=1) - smoothed[:, loss_cols].mean(axis=1)
        )
    else:
        cnv_scores[lo:hi] = cnv_spread[lo:hi]

    # Retain smoothed rows for any pre-selected heatmap cells in this block.
    sel = plot_idx[(plot_idx >= lo) & (plot_idx < hi)]
    if sel.size:
        plot_buffer[np.searchsorted(plot_idx, sel)] = smoothed[sel - lo]
    del smoothed

adata.obs["cnv_score"] = cnv_scores
adata.obs["cnv_spread"] = cnv_spread

_score_definition = (
    f"gain({sorted(gain_contigs)}) minus loss({sorted(loss_contigs)}) mean smoothed deviation"
    if use_contrast else "genome-wide spread of smoothed deviation"
)

# Threshold from the reference distribution: a cell is called malignant when its
# genome-wide profile is more displaced than a confident normal cell's, by more
# than `n_sd` reference standard deviations. Reference cells are non-malignant by
# construction, so this is a false-positive rate control on known normals.
n_sd = float(snakemake.params.cnv_threshold_sd)
ref_cnv = cnv_scores[ref_mask]
threshold = float(ref_cnv.mean() + n_sd * ref_cnv.std())

adata.obs["is_malignant"] = (cnv_scores > threshold).astype(bool)
n_malignant = int(adata.obs["is_malignant"].sum())

# A second, deliberately more sensitive call. `is_malignant` is balanced: it is
# what BUILDS the tumor compartment, so its false positives become tumor cells.
# A compartment that must be *clean* — the nerve compartment — has the opposite
# need, and should exclude anything merely suspected. Consumers pick the flag
# that matches their error cost instead of sharing one threshold that can only
# be right for one of them.
excl_sd = float(snakemake.params.cnv_exclusion_sd)
excl_threshold = float(ref_cnv.mean() + excl_sd * ref_cnv.std())
adata.obs["is_malignant_suspected"] = (cnv_scores > excl_threshold).astype(bool)
n_suspected = int(adata.obs["is_malignant_suspected"].sum())
log_transformation(log, "scrna_malignancy",
    f"Sensitive exclusion flag at {excl_sd} SD (threshold={excl_threshold:.4f}): "
    f"{n_suspected}/{adata.n_obs} suspected ({100 * n_suspected / adata.n_obs:.1f}%)")

log_transformation(log, "scrna_malignancy",
    f"CNV threshold={threshold:.4f} (reference mean {ref_cnv.mean():.4f} + "
    f"{n_sd} x SD {ref_cnv.std():.4f}); "
    f"{n_malignant}/{adata.n_obs} cells flagged as malignant "
    f"({100 * n_malignant / adata.n_obs:.1f}%)")

# Concordance against the external oracle, when the cohort carries one. This is
# the number the Phase-3 gate is written against; recording it here means the
# CNV rule reports its own accuracy rather than deferring entirely to the audit.
if "cell_type" in adata.obs.columns:
    truth = (adata.obs["cell_type"].astype(str) == "malignant cell").to_numpy()
    pred = adata.obs["is_malignant"].to_numpy(dtype=bool)
    tp, fp, fn = int((pred & truth).sum()), int((pred & ~truth).sum()), int((~pred & truth).sum())
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    log_transformation(log, "scrna_malignancy",
        f"Against Census cell_type: precision={precision:.4f}, recall={recall:.4f} "
        f"(TP={tp}, FP={fp}, FN={fn})")

# ---------------------------------------------------------------------------
# Step 4: CNV heatmap (sampled cells × sorted genes), sorted by malignant flag
# ---------------------------------------------------------------------------
malignant_plot = adata.obs["is_malignant"].values[plot_idx]
order = np.argsort(malignant_plot.astype(int))
plot_buffer = plot_buffer[order]
n_normal_plot = int((~malignant_plot).sum())

fig, ax = plt.subplots(figsize=(14, 6))
im = ax.imshow(
    plot_buffer,
    aspect="auto",
    cmap="RdBu_r",
    vmin=-0.5,
    vmax=0.5,
    interpolation="nearest",
)
ax.axhline(n_normal_plot - 0.5, color="black", linewidth=1.5, linestyle="--")
ax.set_xlabel("Genes (genomic order)")
ax.set_ylabel(f"Cells (n={n_plot}; dashed = malignant boundary)")
ax.set_title("CNV Inference Heatmap — Smoothed Expression Signal (genes in true genomic order)")
plt.colorbar(im, ax=ax, label="Centered log-normalized expression")
fig.tight_layout()
fig.savefig(snakemake.output.cnv_plot, dpi=150, bbox_inches="tight")
plt.close(fig)

# --- Write output ------------------------------------------------------------
adata.write_h5ad(snakemake.output.h5ad, compression=H5AD_COMPRESSION)
verify_artifact(snakemake.output.h5ad, min_size_bytes=1024)

prov = stamp_artifact(
    output_path=snakemake.output.h5ad,
    rule_name="scrna_malignancy",
    input_paths=[snakemake.input.h5ad, snakemake.input.gene_positions],
    tool_versions={"scanpy": sc.__version__, "anndata": ad.__version__, "numpy": np.__version__},
    parameters={
        "cnv_score_definition":  _score_definition,
        "cnv_gain_contigs":      sorted(gain_contigs),
        "cnv_loss_contigs":      sorted(loss_contigs),
        "cnv_window_size":       window,
        "cnv_clip":              clip,
        "cnv_threshold":         float(threshold),
        "cnv_threshold_sd":      n_sd,
        "cnv_exclusion_sd":      excl_sd,
        "cnv_exclusion_threshold": excl_threshold,
        "n_malignant_suspected": n_suspected,
        "gene_order":            "Ensembl chromosome + start coordinate (per-contig windows)",
        "n_genes_placed":        int(n_placed),
        "n_genes_in_cnv_axis":   int(n_genes_ord),
        "n_contigs":             len(contig_bounds),
        "reference_labels":      sorted(reference_labels),
        "n_reference_cells":     int(n_ref),
        "n_malignant":           int(n_malignant),
        "n_cells":               adata.n_obs,
        "pct_malignant":         round(100 * n_malignant / adata.n_obs, 2),
    },
    description="CNV-scored AnnData with is_malignant label; per-chromosome sliding-window smoothing "
                "on library-normalized expression, centred on a confident-immune reference",
    ontology_operation="operation:3225",  # EDAM: Copy number variation detection
)
write_provenance(prov, snakemake.output.provenance)

log_transformation(log, "scrna_malignancy", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.h5ad, snakemake.output.cnv_plot])
