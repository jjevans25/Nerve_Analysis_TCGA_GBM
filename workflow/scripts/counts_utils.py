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
