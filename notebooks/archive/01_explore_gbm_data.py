"""
Marimo reactive notebook: Interactive exploration of TCGA GBM loom files.
Addresses: What is the cell/gene composition and QC landscape of each GBM sample?

SUPERSEDED 2026-08-19 — reads the pinned v1.3.0 reference cohort.

This notebook is retained as a record of the pre-fix analysis and is NOT built by
`rule all`. Its nerve compartment was assembled by the logic removed from
`workflow/scripts/nerve_cell_subset.py` on 2026-08-06, which measured 59% malignant
and 11% neural, so every nerve-side claim rendered here is void. It cannot be
corrected: `data/processed/nerve_cells.h5ad` was deleted by a failed job on
2026-07-21 and the v1.3.0 reference is pinned and unreproducible.

The current analysis is `notebooks/05_census_nerve_immune_explorer.py`, over the
CELLxGENE Census arms. This file is kept because it becomes usable again if a
v1.4.0 baseline is ever rebuilt through the corrected pipeline; render it
explicitly via its rule in `workflow/rules/notebooks.smk` if you need the
historical view.
"""

import marimo

__generated_with = "0.23.1"
app = marimo.App(width="wide", app_title="TCGA GBM — Sample Explorer")


@app.cell
def _imports():
    import sys
    import uuid
    from datetime import datetime
    from pathlib import Path

    import anndata as ad
    import duckdb
    import loompy
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import yaml

    return Path, ad, datetime, duckdb, loompy, mo, np, pd, plt, sc, uuid, yaml


@app.cell
def _superseded_banner(mo):
    mo.callout(
        mo.md(
            "**SUPERSEDED — do not read the nerve-side numbers on this page.**\n\n"
            "This notebook reads the pinned **v1.3.0 reference cohort**, whose nerve "
            "compartment was built by the logic removed from `nerve_cell_subset.py` on "
            "2026-08-06. That compartment measured **59% malignant and 11% neural**: "
            "the cells it labelled 'nerve' were largely tumour and myeloid, so every "
            "nerve-side enrichment, interaction and lead axis rendered below is an "
            "artefact of the defect rather than a finding.\n\n"
            "It **cannot be corrected**. `data/processed/nerve_cells.h5ad` was deleted "
            "by a failed job on 2026-07-21 and the v1.3.0 reference is pinned and "
            "structurally unreproducible, so there is no path to re-rendering this "
            "against clean compartments.\n\n"
            "Retained as a record of what was believed before the fix, and because the "
            "code becomes reusable if a v1.4.0 baseline is ever rebuilt through the "
            "corrected pipeline. It is **not** built by `rule all`.\n\n"
            "Current analysis: `notebooks/05_census_nerve_immune_explorer.py` "
            "(CELLxGENE Census arms, compartments audited at 95.4% neural against an "
            "external oracle)."
        ),
        kind="danger",
    )
    return



@app.cell
def _load_config(Path, yaml):
    """Load project config — all paths sourced here, no hardcoded strings."""
    _config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    with open(_config_path) as _f:
        config = yaml.safe_load(_f)
    loom_dir = Path(config["loom_dir"])
    manifest_path = Path(config["loom_manifest"])
    qc_params = config["scrna"]
    sample_ids: list[str] = config["samples"]
    return config, loom_dir, manifest_path, qc_params, sample_ids


@app.cell
def _header(mo):
    mo.md("""
    # TCGA GBM — Interactive Sample Explorer
    **Project:** TCGA Glioblastoma Multiforme | Single-Cell RNA-seq
    **Format:** Seurat 1000×1000 sparse loom files
    **Purpose:** Exploratory QC prior to Snakemake pipeline execution
    """)
    return


@app.cell
def _manifest_table(duckdb, manifest_path, mo):
    """DuckDB query of MANIFEST.txt — fast metadata scan across all 17 samples."""
    _conn = duckdb.connect()
    manifest_df = _conn.execute(
        f"""
        SELECT
            id,
            regexp_extract(filename, '[^/]+$') AS loom_file,
            round(size / 1e6, 1)               AS size_mb,
            state
        FROM read_csv('{manifest_path}', delim='\t', header=true)
        ORDER BY size_mb DESC
        """
    ).df()
    _conn.close()

    mo.md("## Sample Manifest (17 GDC files, ordered by size)")
    return (manifest_df,)


@app.cell
def _show_manifest(manifest_df, mo):
    mo.ui.table(manifest_df, selection=None)
    return


@app.cell
def _sample_selector(mo, sample_ids: list[str]):
    sample_dropdown = mo.ui.dropdown(
        options=sample_ids,
        value=sample_ids[0],
        label="Select GDC sample UUID",
    )
    return (sample_dropdown,)


@app.cell
def _show_selector(mo, qc_params, sample_dropdown):
    mo.hstack(
        [
            sample_dropdown,
            mo.callout(
                mo.md(
                    f"**QC thresholds (from config.yaml):** "
                    f"min_genes={qc_params['min_genes']} | "
                    f"max_genes={qc_params['max_genes']} | "
                    f"max_pct_mito={qc_params['max_pct_mito']}%"
                ),
                kind="info",
            ),
        ],
        gap=1,
    )
    return


@app.cell
def _load_loom(loom_dir, loompy, mo, sample_dropdown):
    """Load the selected loom file in read-only mode."""
    _sample_id = sample_dropdown.value
    _loom_files = list((loom_dir / _sample_id).glob("*.loom"))
    if not _loom_files:
        mo.stop(True, mo.callout(mo.md(f"No loom file found for {_sample_id}"), kind="danger"))

    loom_path = _loom_files[0]
    loom_conn = loompy.connect(str(loom_path), mode="r")
    n_genes, n_cells = loom_conn.shape
    return loom_conn, loom_path, n_cells, n_genes


@app.cell
def _loom_summary(loom_path, mo, n_cells, n_genes, sample_dropdown):
    mo.md(f"""
    ## Sample: `{sample_dropdown.value}`
    **File:** `{loom_path.name}`
    **Dimensions:** {n_genes:,} genes × {n_cells:,} cells
    """)
    return


@app.cell
def _build_adata(ad, loom_conn, mo, np):
    """Convert loom to AnnData and compute QC metrics."""
    _matrix = loom_conn[:, :]  # genes × cells; transpose to cells × genes
    adata = ad.AnnData(X=_matrix.T.astype(np.float32))
    adata.var_names = [str(g) for g in loom_conn.ra.get("Gene", range(loom_conn.shape[0]))]
    adata.obs_names = [str(c) for c in loom_conn.ca.get("CellID", range(loom_conn.shape[1]))]

    # Mitochondrial gene flag (MT- prefix after Ensembl → gene symbol mapping)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    mo.md("*Computing QC metrics (n_genes_by_counts, total_counts, pct_counts_mt)…*")
    return (adata,)


@app.cell
def _qc_metrics(adata, sc):
    """Log QC threshold parameters before any filtering step."""
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return


@app.cell
def _qc_plots(adata, plt, qc_params):
    """Reactive QC histograms — updates when sample selection changes."""
    _fig, _axes = plt.subplots(1, 3, figsize=(14, 4))

    _axes[0].hist(adata.obs["n_genes_by_counts"], bins=60, color="#4C72B0", edgecolor="white", linewidth=0.3)
    _axes[0].axvline(qc_params["min_genes"], color="red", linestyle="--", label=f"min={qc_params['min_genes']}")
    _axes[0].axvline(qc_params["max_genes"], color="orange", linestyle="--", label=f"max={qc_params['max_genes']}")
    _axes[0].set_xlabel("Genes per cell")
    _axes[0].set_ylabel("Cell count")
    _axes[0].set_title("n_genes_by_counts")
    _axes[0].legend(fontsize=8)

    _axes[1].hist(adata.obs["total_counts"], bins=60, color="#55A868", edgecolor="white", linewidth=0.3)
    _axes[1].set_xlabel("UMI counts per cell")
    _axes[1].set_ylabel("Cell count")
    _axes[1].set_title("total_counts")

    _mt_pct = adata.obs["pct_counts_mt"]
    _axes[2].hist(_mt_pct, bins=60, color="#C44E52", edgecolor="white", linewidth=0.3)
    _axes[2].axvline(qc_params["max_pct_mito"], color="red", linestyle="--", label=f"max={qc_params['max_pct_mito']}%")
    _axes[2].set_xlabel("% mitochondrial counts")
    _axes[2].set_ylabel("Cell count")
    _axes[2].set_title("pct_counts_mt")
    _axes[2].legend(fontsize=8)

    _fig.suptitle("QC Distributions — thresholds from config.yaml", fontsize=11)
    _fig.tight_layout()
    return


@app.cell
def _qc_summary_table(adata, mo, pd, qc_params):
    """Tabular QC summary with pass/fail cell counts relative to config thresholds."""
    _obs = adata.obs
    _n_total = len(_obs)
    _pass_min_genes  = (_obs["n_genes_by_counts"] >= qc_params["min_genes"]).sum()
    _pass_max_genes  = (_obs["n_genes_by_counts"] <= qc_params["max_genes"]).sum()
    _pass_mito       = (_obs["pct_counts_mt"]     <= qc_params["max_pct_mito"]).sum()
    _pass_all        = (
        (_obs["n_genes_by_counts"] >= qc_params["min_genes"]) &
        (_obs["n_genes_by_counts"] <= qc_params["max_genes"]) &
        (_obs["pct_counts_mt"]     <= qc_params["max_pct_mito"])
    ).sum()

    _summary = pd.DataFrame(
        {
            "Metric":       ["Total cells", "≥ min_genes", "≤ max_genes", "≤ max_pct_mito", "Passing all filters"],
            "Count":        [_n_total, int(_pass_min_genes), int(_pass_max_genes), int(_pass_mito), int(_pass_all)],
            "Pct of total": [
                "100%",
                f"{100 * _pass_min_genes  / _n_total:.1f}%",
                f"{100 * _pass_max_genes  / _n_total:.1f}%",
                f"{100 * _pass_mito       / _n_total:.1f}%",
                f"{100 * _pass_all        / _n_total:.1f}%",
            ],
        }
    )

    mo.vstack([
        mo.md("### QC Filter Summary"),
        mo.ui.table(_summary, selection=None),
    ])
    return


@app.cell
def _provenance(datetime, mo, sample_dropdown, uuid):
    """FAIR provenance record for this notebook execution."""
    _run_id    = str(uuid.uuid4())
    _timestamp = datetime.utcnow().isoformat() + "Z"
    _prov = {
        "run_id":    _run_id,
        "timestamp": _timestamp,
        "notebook":  "notebooks/01_explore_gbm_data.py",
        "sample":    sample_dropdown.value,
        "tool":      "marimo==0.23.1",
    }
    mo.callout(
        mo.md(
            f"**FAIR Provenance**  \n"
            f"`run_id:` {_run_id}  \n"
            f"`timestamp:` {_timestamp}  \n"
            f"`sample:` {sample_dropdown.value}"
        ),
        kind="success",
    )
    return


@app.cell
def _nerve_cell_header(mo):
    mo.md("""
    ---
    ## Nerve Cell Heterogeneity Explorer
    Visualizes results from the nerve-cell subset and heterogeneity Snakemake rules.
    Requires the pipeline to have been run through the `nerve_cell_heterogeneity` rule.

    *Per-cluster readers (`nerve_enrichment.csv`, `nerve_cluster_markers.csv`,
    `nerve_tumor_interactions.csv`, `nerve_tumor_top_pairs.csv`) now have
    `_with_qc.csv` companions produced by the `annotate_cluster_qc` rule —
    each row carries `batch_qc_pass` and purity context so QC-failing clusters
    can be flagged in any downstream view.*
    """)
    return


@app.cell
def _nerve_artifacts(Path, config, mo):
    """Check which nerve-cell result artifacts exist."""
    results_dir = Path(config["dirs"]["results"])
    _artifacts = {
        "nerve_cells.h5ad":               Path(config["dirs"]["data_processed"]) / "nerve_cells.h5ad",
        "nerve_cluster_markers.csv":       results_dir / "tables" / "nerve_cluster_markers.csv",
        "nerve_enrichment.csv":            results_dir / "tables" / "nerve_enrichment.csv",
        "nerve_cells_umap.png":            results_dir / "figures" / "nerve_cells_umap.png",
        "nerve_dotplot.png":               results_dir / "figures" / "nerve_dotplot.png",
        "nerve_abundance_heatmap.png":     results_dir / "figures" / "nerve_abundance_heatmap.png",
        "gdc_clinical.tsv":                Path(config["dirs"]["data_external"]) / "gdc_clinical.tsv",
    }
    _status = {k: "✓" if v.exists() else "✗ (not yet produced)" for k, v in _artifacts.items()}
    nerve_artifacts = _artifacts
    mo.callout(
        mo.md("**Nerve-cell pipeline artifact status:**\n\n" +
              "\n".join(f"- `{k}`: {s}" for k, s in _status.items())),
        kind="info",
    )
    return (nerve_artifacts,)


@app.cell
def _nerve_umap(mo, nerve_artifacts):
    """Display nerve-cell UMAP if available."""
    _path = nerve_artifacts["nerve_cells_umap.png"]
    if _path.exists():
        mo.image(str(_path), alt="Nerve-cell UMAP", width="100%")
    else:
        mo.callout(mo.md("UMAP not yet generated — run `nerve_cell_subset` rule first."), kind="warn")
    return


@app.cell
def _nerve_markers_table(duckdb, mo, nerve_artifacts):
    """Top DE markers per nerve-cell cluster via DuckDB."""
    _path = nerve_artifacts["nerve_cluster_markers.csv"]
    _has_data = _path.exists() and sum(
        1 for line in _path.read_text().splitlines() if line.strip()
    ) > 1
    if not _has_data:
        mo.stop(True, mo.callout(
            mo.md("Marker CSV not yet generated or contains no data rows."), kind="warn",
        ))

    _conn = duckdb.connect()
    nerve_markers_df = _conn.execute(
        f"""
        SELECT cluster,
               names      AS gene,
               round(logfoldchanges, 3) AS log2fc,
               round(scores, 3)         AS score,
               round(pvals_adj, 4)      AS padj
        FROM read_csv('{_path}', header=true)
        WHERE pvals_adj < 0.05
        ORDER BY cluster, score DESC
        LIMIT 200
        """
    ).df()
    _conn.close()

    if nerve_markers_df.empty:
        mo.stop(True, mo.callout(
            mo.md("No nerve-cell markers available yet — upstream heterogeneity step found 0 cells."),
            kind="warn",
        ))

    mo.vstack([
        mo.md("### Top Differential Markers per Nerve-Cell Cluster"),
        mo.ui.table(nerve_markers_df, selection=None),
    ])
    return


@app.cell
def _nerve_enrichment_table(duckdb, mo, nerve_artifacts):
    """Top GSEA enrichment terms per nerve cluster."""
    _path = nerve_artifacts["nerve_enrichment.csv"]
    _has_data = _path.exists() and sum(
        1 for line in _path.read_text().splitlines() if line.strip()
    ) > 1
    if not _has_data:
        mo.stop(True, mo.callout(
            mo.md("Enrichment CSV not yet generated or contains no data rows."), kind="warn",
        ))

    _conn = duckdb.connect()
    enr_df = _conn.execute(
        f"""
        SELECT cluster,
               gene_set_library,
               Term,
               round(CAST("Adjusted P-value" AS DOUBLE), 4) AS padj,
               "Overlap"
        FROM read_csv('{_path}', header=true)
        WHERE CAST("Adjusted P-value" AS DOUBLE) < 0.05
        ORDER BY cluster, padj
        LIMIT 100
        """
    ).df()
    _conn.close()

    if enr_df.empty:
        mo.stop(True, mo.callout(
            mo.md("No GSEA enrichment results available yet."), kind="warn",
        ))

    mo.vstack([
        mo.md("### GSEA Enrichment — Top Terms per Nerve-Cell Cluster"),
        mo.ui.table(enr_df, selection=None),
    ])
    return


@app.cell
def _nerve_dotplot(mo, nerve_artifacts):
    """Canonical marker dot plot."""
    _path = nerve_artifacts["nerve_dotplot.png"]
    if _path.exists():
        mo.image(str(_path), alt="Nerve-cell marker dot plot", width="100%")
    else:
        mo.callout(mo.md("Dot plot not yet generated."), kind="warn")
    return


@app.cell
def _nerve_abundance(mo, nerve_artifacts):
    """Per-sample cluster abundance heatmap."""
    _path = nerve_artifacts["nerve_abundance_heatmap.png"]
    if _path.exists():
        mo.image(str(_path), alt="Nerve-cell cluster abundance heatmap", width="100%")
    else:
        mo.callout(mo.md("Abundance heatmap not yet generated."), kind="warn")
    return


@app.cell
def _clinical_table(duckdb, mo, nerve_artifacts):
    """GDC clinical metadata summary."""
    _path = nerve_artifacts["gdc_clinical.tsv"]
    if not _path.exists():
        mo.stop(True, mo.callout(mo.md("Clinical metadata not yet fetched."), kind="warn"))

    _conn = duckdb.connect()
    clinical_df = _conn.execute(
        f"""
        SELECT file_uuid, case_id, primary_diagnosis,
               tumor_grade, tissue_type, gender, age_at_index
        FROM read_csv('{_path}', delim='\t', header=true)
        ORDER BY primary_diagnosis
        """
    ).df()
    _conn.close()

    mo.vstack([
        mo.md("### GDC Clinical Metadata"),
        mo.ui.table(clinical_df, selection=None),
    ])
    return


if __name__ == "__main__":
    app.run()
