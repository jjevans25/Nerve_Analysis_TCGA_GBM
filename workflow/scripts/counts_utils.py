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
