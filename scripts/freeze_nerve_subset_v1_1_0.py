"""Freeze the v1.1.0 nerve-subset cell list for panel-tightening insulation.

Captures the obs_names from data/processed/nerve_cells.h5ad as committed at
the v1.1.0 tag. The v1.2.0 rerun reads this list via
config.nerve_cells.frozen_subset_file to ensure the nerve-subset cell
composition (and therefore nerve_leiden cluster IDs 0–23) survives the
ependymal-panel change in config.

One-shot — not a Snakemake rule. Run once, commit the output, retire when
v1.3.0 does a full rerun and accepts the cluster renumber.

Usage:
    .snakemake/conda/<scrna-env>/bin/python scripts/freeze_nerve_subset_v1_1_0.py
"""

import hashlib
from pathlib import Path

import anndata as ad

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data/processed/nerve_cells.h5ad"
OUT = ROOT / "provenance/nerve_subset_v1_1_0.txt"


def main() -> None:
    if not SRC.exists():
        raise FileNotFoundError(
            f"v1.1.0 nerve_cells.h5ad missing at {SRC}. The freeze script must "
            "run against the v1.1.0 cut artifact."
        )

    adata = ad.read_h5ad(SRC, backed="r")
    barcodes = list(adata.obs_names)
    n = len(barcodes)
    adata.file.close()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(barcodes) + "\n")
    sha = hashlib.sha256(OUT.read_bytes()).hexdigest()

    print(f"Source        : {SRC}")
    print(f"Output        : {OUT}")
    print(f"N barcodes    : {n}")
    print(f"SHA256        : {sha}")
    print()
    print("Record this SHA256 in CHANGELOG.md under the v1.2.0 entry.")


if __name__ == "__main__":
    main()
