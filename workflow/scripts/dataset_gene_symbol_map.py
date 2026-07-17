"""Build a dataset-scoped Ensembl→symbol map from the cohort's own .var (Stage A).

The reference cohort's symbol map is a MyGene.info cache keyed on *versioned*
loom Ensembl IDs. An external cohort (CELLxGENE Census, unversioned Ensembl) is
not covered by that cache, so scrna_annotate's fallback re-join would yield an
all-NaN gene_symbol and silently collapse annotation. The Census .var already
carries the authoritative mapping (feature_id → feature_name), so we materialize
it here as the same 3-column TSV (ensembl_id, gene_symbol, chromosome) the
reused scripts expect. Chromosome is set to 'unknown' — MT flagging downstream
uses the symbol prefix, not this column.
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
raw_file = snakemake.params.raw_file

log_transformation(log, "dataset_gene_symbol_map", f"Reading var from {raw_file} (backed)")
var = ad.read_h5ad(raw_file, backed="r").var

if {"feature_id", "feature_name"} <= set(var.columns):
    ensembl = var["feature_id"].astype(str).values
    symbol = var["feature_name"].astype(str).values
else:
    ensembl = var.index.astype(str).values
    symbol = var.index.astype(str).values

out = pd.DataFrame({"ensembl_id": ensembl, "gene_symbol": symbol, "chromosome": "unknown"})
out = out.drop_duplicates(subset="ensembl_id", keep="first")

Path(snakemake.output.cache).parent.mkdir(parents=True, exist_ok=True)
out.to_csv(snakemake.output.cache, sep="\t", index=False)
verify_artifact(snakemake.output.cache, min_size_bytes=16)
log_transformation(log, "dataset_gene_symbol_map",
    f"Wrote {len(out)} ensembl→symbol rows")

prov = stamp_artifact(
    output_path=snakemake.output.cache,
    rule_name="dataset_gene_symbol_map",
    input_paths=[raw_file],
    tool_versions={"anndata": ad.__version__, "pandas": pd.__version__},
    parameters={"dataset": snakemake.params.dataset, "n_genes": int(len(out))},
    description="Dataset-scoped Ensembl→symbol map extracted from cohort .var",
    ontology_operation="operation:2409",
)
write_provenance(prov, snakemake.output.provenance)
log_transformation(log, "dataset_gene_symbol_map", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.cache])
