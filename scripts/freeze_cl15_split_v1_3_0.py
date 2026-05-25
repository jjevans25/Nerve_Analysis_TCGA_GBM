"""Freeze the cl15 sub-Leiden split assignment for the v1.3.0 surgical relabel.

The v1.2.0 rerun confirmed that ependymal-panel methodology alone cannot
resolve cl15's MIXED status (cluster-level argmax is robust to per-cell score
shifts). The v1.3.0 fix is surgical: split cl15's four sub-Leiden
subpopulations into separately-labeled nerve_leiden cluster IDs.

This script re-derives the cl15 sub-Leiden split with the *identical*
parameters used by scripts/diagnose_cl15_ependymal.py (resolution=0.5 on the
X_scVI latent, igraph flavor, random_state=0) and persists the per-barcode
assignment as a committed artifact. The pipeline (workflow/scripts/
nerve_cell_subset.py) then applies the relabel by barcode, so no live
re-clustering happens during the run — mirroring the freeze-artifact pattern
of scripts/freeze_nerve_subset_v1_1_0.py.

Target-ID assignment is by marker/score *signature*, not by sub-cluster
number (sub-Leiden numbering is not guaranteed stable across runs):

    sub with max mean score_ependymal               -> 24  ependymal
    sub with max mean score_microglia (immune tail)  -> 27  artifact (dropped)
    of the remaining two:
        higher mean score_excitatory_neuron          -> 25  neuron
        the other                                     -> 26  transitional

Falsification gate (the doc's NEW claim, DO_THIS_NEXT_post_v1.2.0_rerun.md
lines 90-94): the sub-population sizes must still match the v1.1.0 evidence
pack {773, 692, 378, 11}. If they have materially drifted, the frozen X_scVI
latent has leaked somewhere and the split must be re-diagnosed before use.

One-shot — not a Snakemake rule. Run once, commit the output CSV, record the
SHA256 in CHANGELOG.md under the v1.3.0 entry.

Usage:
    .snakemake/conda/<scrna-env>/bin/python scripts/freeze_cl15_split_v1_3_0.py
"""

import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data/processed/nerve_cells.h5ad"
OUT = ROOT / "provenance/cl15_split_v1_3_0.csv"

SOURCE_CLUSTER = "15"
# Expected sub-population sizes from results/tables/cl15_ependymal_diagnosis.json
# (v1.1.0 evidence pack). Used as the freeze-leak falsification gate.
EXPECTED_SORTED_SIZES = [11, 378, 692, 773]
SIZE_TOLERANCE = 0.05  # +/-5% per matched population

# Target nerve_leiden IDs + cell-type labels for each signature.
TARGET = {
    "ependymal": ("24", "ependymal"),
    "neuron": ("25", "neuron"),
    "transitional": ("26", "transitional"),
    "artifact": ("27", "artifact"),
}


def _check_sizes(sizes: list[int]) -> None:
    """Falsification gate: sub-population sizes must match the v1.1.0 pack."""
    got = sorted(sizes)
    if len(got) != len(EXPECTED_SORTED_SIZES):
        raise RuntimeError(
            f"[FAIR-ALERT] cl15 sub-Leiden produced {len(got)} sub-clusters "
            f"(sizes {got}); expected {len(EXPECTED_SORTED_SIZES)} "
            f"({EXPECTED_SORTED_SIZES}). The frozen X_scVI latent has drifted "
            "since the v1.1.0 evidence pack — diagnose the freeze leak before "
            "regenerating the split (see DO_THIS_NEXT_post_v1.2.0_rerun.md "
            "lines 90-94)."
        )
    for g, e in zip(got, EXPECTED_SORTED_SIZES):
        if abs(g - e) > max(1, round(e * SIZE_TOLERANCE)):
            raise RuntimeError(
                f"[FAIR-ALERT] cl15 sub-population size {g} deviates from the "
                f"expected {e} by more than {SIZE_TOLERANCE:.0%}. Sub-Leiden "
                f"sizes {got} vs expected {EXPECTED_SORTED_SIZES}. The freeze "
                "may be leaking — investigate before committing the split."
            )


def main() -> None:
    if not SRC.exists():
        raise FileNotFoundError(
            f"nerve_cells.h5ad missing at {SRC}. Run against the current "
            "frozen (v1.2.0) cut artifact."
        )

    # Backed read: obs/obsm load into memory; .X stays on disk (file is large).
    adata = ad.read_h5ad(SRC, backed="r")
    leiden = adata.obs["nerve_leiden"].astype(str)
    cl15_mask = (leiden == SOURCE_CLUSTER).to_numpy()
    n_cl15 = int(cl15_mask.sum())
    if n_cl15 == 0:
        raise RuntimeError(
            f"No cells in nerve_leiden == '{SOURCE_CLUSTER}'. Has the split "
            "already been applied, or did the cluster IDs change?"
        )
    print(f"Source            : {SRC}")
    print(f"cl{SOURCE_CLUSTER} cells          : {n_cl15}")

    cl15_pos = np.where(cl15_mask)[0]
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    obs_cl15 = adata.obs.iloc[cl15_pos][["sample_id", *score_cols]].copy()
    obs_cl15.index = adata.obs_names[cl15_pos]
    x_scvi = np.asarray(adata.obsm["X_scVI"])[cl15_pos]
    adata.file.close()

    # Re-derive the sub-Leiden split — identical params to
    # scripts/diagnose_cl15_ependymal.py (X_scVI latent only; no .X needed
    # since target assignment uses obs score signatures, not DE markers).
    sub = ad.AnnData(obs=obs_cl15, obsm={"X_scVI": x_scvi})
    sc.pp.neighbors(sub, use_rep="X_scVI", random_state=0)
    sc.tl.leiden(sub, resolution=0.5, key_added="cl15_subleiden", random_state=0,
                 flavor="igraph", n_iterations=2, directed=False)
    sub_id = sub.obs["cl15_subleiden"].astype(str)
    sizes = sub_id.value_counts().sort_index()
    print(f"sub-cluster sizes : {sizes.to_dict()}")
    _check_sizes(list(sizes.to_numpy()))
    print("falsification gate: PASS (sizes match v1.1.0 evidence pack)")

    # Mean score per sub-cluster (signature for target assignment).
    means = sub.obs.groupby("cl15_subleiden", observed=True)[score_cols].mean()

    ependymal_sub = means["score_ependymal"].idxmax()
    artifact_sub = means["score_microglia"].idxmax()
    remaining = [s for s in means.index if s not in {ependymal_sub, artifact_sub}]
    if len(remaining) != 2:
        raise RuntimeError(
            f"Expected exactly 2 sub-clusters after assigning ependymal "
            f"({ependymal_sub}) and artifact ({artifact_sub}); got {remaining}."
        )
    neuron_sub = means.loc[remaining, "score_excitatory_neuron"].idxmax()
    transitional_sub = [s for s in remaining if s != neuron_sub][0]

    sub_to_sig = {
        ependymal_sub: "ependymal",
        neuron_sub: "neuron",
        transitional_sub: "transitional",
        artifact_sub: "artifact",
    }
    print("\nsignature assignment (sub-cluster -> target):")
    for s in sizes.index:
        sig = sub_to_sig[s]
        tid, tlabel = TARGET[sig]
        print(f"  sub {s} (n={sizes[s]:>4})  ep={means.loc[s,'score_ependymal']:+.3f} "
              f"exc_neu={means.loc[s,'score_excitatory_neuron']:+.3f} "
              f"micro={means.loc[s,'score_microglia']:+.3f}  -> cl{tid} {tlabel}")

    # Per-barcode assignment table.
    rows = []
    for barcode, s in sub_id.items():
        sig = sub_to_sig[s]
        tid, tlabel = TARGET[sig]
        rows.append((barcode, s, tid, tlabel))
    rows.sort(key=lambda r: (r[2], r[0]))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        fh.write("barcode,sub_id,target_cluster,target_cell_type\n")
        for barcode, s, tid, tlabel in rows:
            fh.write(f"{barcode},{s},{tid},{tlabel}\n")
    sha = hashlib.sha256(OUT.read_bytes()).hexdigest()

    print(f"\nOutput            : {OUT}")
    print(f"N barcodes        : {len(rows)}")
    print("target cluster counts:")
    for sig, (tid, tlabel) in TARGET.items():
        n = sum(1 for r in rows if r[2] == tid)
        print(f"  cl{tid} {tlabel:<13}: {n}")
    print(f"SHA256            : {sha}")
    print("\nRecord this SHA256 in CHANGELOG.md under the v1.3.0 entry.")


if __name__ == "__main__":
    main()
