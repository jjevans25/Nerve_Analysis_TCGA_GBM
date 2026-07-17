"""Derive the per-donor sample sheet + MANIFEST for a CELLxGENE h5ad replication cohort.

SETUP utility (run once, outside the Snakemake DAG — like building MANIFEST.txt),
not an analytical step. The sample list must exist at Snakefile *parse* time so
`expand()` over donors can build the DAG, hence it is a committed sample sheet
rather than a rule output.

Reads the source h5ad `.obs` in backed mode (metadata only — never loads the
matrix), restricts to the configured Census study, and writes:
  data/raw/<dataset>/samples.txt   one donor_id per line (the {sample} wildcard values)
  data/raw/<dataset>/MANIFEST.txt  tab-separated donor_id / n_cells / source_file / study

Usage:
  python scripts/derive_dataset_sample_sheet.py <dataset_key> [config/config.yaml]
"""

import sys
from pathlib import Path

import anndata as ad
import yaml


def main() -> None:
    dataset_key = sys.argv[1] if len(sys.argv) > 1 else "gbm_cellxgene_56c4912d"
    config_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("config/config.yaml")

    cfg = yaml.safe_load(config_path.read_text())
    entry = cfg["datasets"][dataset_key]

    raw_file = Path(entry["raw_file"])
    sample_key = entry.get("sample_key", "donor_id")
    filt = entry.get("filter", {}) or {}

    print(f"[derive] dataset={dataset_key} source={raw_file} sample_key={sample_key} filter={filt}")
    adata = ad.read_h5ad(raw_file, backed="r")
    obs = adata.obs
    print(f"[derive] full obs: {obs.shape[0]:,} rows")

    mask = obs.index == obs.index  # all True
    for col, val in filt.items():
        mask = mask & (obs[col].astype(str) == str(val))
    sub = obs.loc[mask]
    print(f"[derive] after filter: {sub.shape[0]:,} rows")

    counts = sub[sample_key].astype(str).value_counts()
    donors = sorted(counts.index.tolist())
    print(f"[derive] {len(donors)} unique {sample_key} values")

    out_dir = Path(cfg["dirs"]["data_raw"]) / dataset_key
    out_dir.mkdir(parents=True, exist_ok=True)

    samples_txt = out_dir / "samples.txt"
    samples_txt.write_text("\n".join(donors) + "\n")

    study = next(iter(filt.values()), "NA") if filt else "NA"
    manifest = out_dir / "MANIFEST.txt"
    lines = ["donor_id\tn_cells\tsource_file\tstudy"]
    lines += [f"{d}\t{int(counts[d])}\t{raw_file.name}\t{study}" for d in donors]
    manifest.write_text("\n".join(lines) + "\n")

    print(f"[derive] wrote {samples_txt} ({len(donors)} donors)")
    print(f"[derive] wrote {manifest}")


if __name__ == "__main__":
    main()
