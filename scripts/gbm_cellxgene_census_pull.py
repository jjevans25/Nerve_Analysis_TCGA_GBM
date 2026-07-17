import cellxgene_census
import scanpy as sc
import pathlib

# Open a connection to the default Census build
with cellxgene_census.open_soma() as census:
    print("Querying CELLxGENE Census for Glioblastoma data...")

    # Extract the AnnData slice
    # is_primary_data=True filters out duplicate cells from meta-analyses
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="homo_sapiens",
        obs_value_filter="disease == 'glioblastoma' and is_primary_data == True",
        column_names={
            "obs": [
                "dataset_id", 
                "donor_id",        # Maps to your patient/sample ID for batch correction
                "assay",           # To filter for 10x Genomics
                "cell_type",       # Provided ontology labels 
                "tissue_general"
            ]
        }
    )

# Filter for 10x Genomics assays locally
# The Census tracks various 10x chemistries (e.g., "10x 3' v3", "10x 5' v2")
# String matching captures all of them while dropping Smart-seq2, Drop-seq, etc.
is_10x = adata.obs["assay"].str.contains("10x", na=False)
adata_10x = adata[is_10x].copy()

# Clean up the object to match your ingest pipeline
# Map their donor_id to the sample_id your pipeline expects
adata_10x.obs["sample_id"] = adata_10x.obs["donor_id"]

print(f"Successfully retrieved {adata_10x.n_obs} 10x GBM cells across {adata_10x.obs['sample_id'].nunique()} unique samples.")

# Define the target directory and ensure it exists
output_dir = pathlib.Path("cellxgene_data")
output_dir.mkdir(parents=True, exist_ok=True)

# Define the exact file path
file_path = output_dir / "gbm_10x_raw.h5ad"

# Save the AnnData object to disk
# Using compression saves significant SSD space with minimal read/write overhead
adata_10x.write_h5ad(file_path, compression="gzip")

print(f"Successfully saved {adata_10x.n_obs} cells to {file_path}")