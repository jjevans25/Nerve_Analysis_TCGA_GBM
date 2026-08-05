"""Shared raw-count recovery (no snakemake side effects; safe to import)."""

import numpy as np
import scipy.sparse as sp


def recover_counts_from_log1p(X) -> sp.csr_matrix:
    """log1p(integer counts) -> integer counts via expm1 + nearest-int rounding.

    Baseline .X is Seurat SCT log1p (no library-size normalization upstream);
    scvi-tools needs raw integer counts. Returns a CSR int32 matrix.
    """
    if sp.issparse(X):
        data = np.rint(np.expm1(X.data)).astype(np.int32)
        out = sp.csr_matrix((data, X.indices, X.indptr), shape=X.shape, dtype=np.int32)
        out.eliminate_zeros()
        return out
    return sp.csr_matrix(np.rint(np.expm1(np.asarray(X))).astype(np.int32))


# A log1p-normalized expression matrix cannot plausibly exceed this. With a
# 1e4 target sum, a gene taking the entire library gives log1p(1e4) ≈ 9.2; the
# headroom to 50 tolerates unusual normalizations while still being orders of
# magnitude below raw UMI maxima (the Census cohort's raw .X peaks at 53,027).
LOG1P_MAX_PLAUSIBLE = 50.0


def is_log1p_scale(X, max_plausible: float = LOG1P_MAX_PLAUSIBLE) -> bool:
    """Is this matrix on a log1p scale, rather than raw counts?

    Guards the silent-wrong-answer failure mode that produced invalid Census
    ligand-receptor tables: LIANA and score_genes both assume log1p input and
    neither checks, so raw counts yield empty specificity ranks and infinite
    log-fold-changes instead of an error. See
    markdowns/blocker_census_liana_raw_counts.md.
    """
    return float(X.max()) <= max_plausible


def library_sizes(X: sp.csr_matrix, block: int = 100_000) -> np.ndarray:
    """Total UMIs per cell — the denominator for library-size normalization."""
    out = np.zeros(X.shape[0], dtype=np.float64)
    for lo in range(0, X.shape[0], block):
        hi = min(lo + block, X.shape[0])
        out[lo:hi] = np.asarray(X[lo:hi].sum(axis=1)).ravel()
    return out


def normalize_log1p_inplace(
    X: sp.csr_matrix,
    lib: np.ndarray,
    target_sum: float = 1e4,
    block: int = 100_000,
) -> None:
    """Library-size normalize + log1p a CSR counts matrix IN PLACE, allocating nothing.

    Marker-gene scoring assumes log1p input. A full normalized copy of the Census
    cohort is 16.5 GB on top of the 16.5 GB original — 33 GB on a 36 GB machine,
    which does not fit. Scaling `X.data` row-block-wise through `indptr` touches
    only the existing buffer, so peak memory is unchanged. Pair with
    `denormalize_log1p_inplace` to hand raw counts back to downstream consumers.
    """
    if not sp.isspmatrix_csr(X):
        raise TypeError(f"expected CSR, got {type(X)!r}")
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(lib > 0, target_sum / lib, 0.0).astype(np.float32)
    for lo in range(0, X.shape[0], block):
        hi = min(lo + block, X.shape[0])
        start, stop = X.indptr[lo], X.indptr[hi]
        counts = np.diff(X.indptr[lo:hi + 1])
        X.data[start:stop] *= np.repeat(scale[lo:hi], counts)
    np.log1p(X.data, out=X.data)


def denormalize_log1p_inplace(
    X: sp.csr_matrix,
    lib: np.ndarray,
    target_sum: float = 1e4,
    block: int = 100_000,
) -> None:
    """Invert `normalize_log1p_inplace`, restoring integer counts IN PLACE.

    Downstream consumers (scVI, the CNV caller, LIANA's own normalization step)
    all require raw counts in `.X`. Counts are integers, so `rint` absorbs the
    float32 round-trip error exactly for the magnitudes seen here.
    """
    np.expm1(X.data, out=X.data)
    with np.errstate(divide="ignore", invalid="ignore"):
        unscale = np.where(lib > 0, lib / target_sum, 0.0).astype(np.float32)
    for lo in range(0, X.shape[0], block):
        hi = min(lo + block, X.shape[0])
        start, stop = X.indptr[lo], X.indptr[hi]
        counts = np.diff(X.indptr[lo:hi + 1])
        X.data[start:stop] *= np.repeat(unscale[lo:hi], counts)
    np.rint(X.data, out=X.data)
