"""
Marimo reactive notebook: nerve-cluster GSEA enrichment for the CELLxGENE Census arms.
Addresses: which biological themes characterise each nerve-cell cluster, and which
clusters share an enrichment profile?

Census replacement for the archived notebooks/archive/02_nerve_enrichment_explorer.py.
Same question, materially different answer: the archived version ran over the pinned
v1.3.0 reference, whose nerve compartment was 59% malignant and 11% neural, so its
"nerve cluster" enrichment largely described tumour and myeloid biology. The
compartment behind this page is 95.4% neural against a held-out oracle (notebook 06).

The QC caveat is not incidental and is surfaced everywhere below: a third of the
enrichment rows here come from clusters that FAIL the donor-dominance test. They are
retained rather than dropped — the test cannot separate "real biology preserved in
one donor's tissue" from "one patient's artifact" — but a term enriched only in a
single-donor cluster is a statement about that donor.

Source rules: workflow/rules/datasets.smk -> ds_nerve_cell_heterogeneity (GSEA over
MSigDB GO sets), ds_annotate_cluster_qc (joins donor purity + batch_qc_pass).
"""

import marimo

__generated_with = "0.23.1"
app = marimo.App(width="wide", app_title="GBM Census — Nerve Enrichment")


@app.cell
def _imports():
    import os
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import yaml

    return Path, mo, np, os, pd, plt, yaml


@app.cell
def _load_config(Path, mo, os, yaml):
    """Resolve cohort-namespaced paths. Cohort is overridable via GBM_DATASET."""
    project_root = Path(__file__).parent.parent
    with open(project_root / "config" / "config.yaml") as _f:
        config = yaml.safe_load(_f)

    _datasets = list(config.get("datasets", {}))
    dataset = os.environ.get("GBM_DATASET") or (
        _datasets[0] if _datasets else "gbm_cellxgene_56c4912d"
    )
    ds_tables = project_root / config["dirs"]["tables"] / dataset

    paths = {
        "enrichment": ds_tables / "nerve_enrichment_with_qc.csv",
        "markers": ds_tables / "nerve_cluster_markers_with_qc.csv",
        "annotations": ds_tables / "nerve_cluster_annotations.csv",
    }
    _missing = {k: p for k, p in paths.items() if not p.exists()}
    if _missing:
        mo.stop(True, mo.callout(mo.md(
            f"Required artifacts not found for cohort `{dataset}`:\n\n"
            + "\n".join(f"- `{p}`" for p in _missing.values())
            + "\n\n```\nscripts/run_snakemake.sh \\\n"
            f"  results/tables/{dataset}/nerve_enrichment_with_qc.csv \\\n"
            "  --use-conda --cores all --rerun-triggers mtime\n```"
        ), kind="danger"))
    return dataset, paths


@app.cell
def _load(pd, paths):
    enrich_df = pd.read_csv(paths["enrichment"])
    markers_df = pd.read_csv(paths["markers"])
    ann_df = pd.read_csv(paths["annotations"])
    # `label` is written as "c{cluster} | {dominant cell type} | {marker argmax}".
    _parts = ann_df["label"].astype(str).str.split("|", expand=True)
    ann_df = ann_df.copy()
    ann_df["cell_type"] = _parts[1].str.strip() if _parts.shape[1] > 1 else ""
    ann_df["cluster_key"] = ann_df["cluster"].astype(str).str.strip()
    for _d in (enrich_df, markers_df):
        _d["cluster_key"] = _d["cluster"].astype(str).str.strip()
        _d["batch_qc_pass"] = _d["batch_qc_pass"].astype(str).str.lower().isin(
            ["true", "1"])
    _lut = dict(zip(ann_df["cluster_key"], ann_df["cell_type"]))
    for _d in (enrich_df, markers_df):
        _d["cell_type"] = _d["cluster_key"].map(_lut).fillna("(unlabelled)")
    return ann_df, enrich_df, markers_df


@app.cell
def _header(dataset, enrich_df, mo):
    _n = len(enrich_df)
    _pass = int(enrich_df["batch_qc_pass"].sum())
    _nclust = enrich_df["cluster_key"].nunique()
    mo.md(f"""
    # GBM Census — Nerve Cluster Enrichment

    **Cohort: `{dataset}`** — **{_n:,} enriched terms** across **{_nclust} nerve
    clusters**, of which **{_pass:,} ({_pass / _n:.0%})** come from clusters that pass
    the donor-dominance test.

    Census replacement for the archived `02_nerve_enrichment_explorer.py`. The
    question is the same; the answer is not. The archived version ran over the pinned
    v1.3.0 reference, whose "nerve" compartment was **59% malignant and 11% neural** —
    so its cluster enrichment was largely describing tumour and myeloid biology under
    a nerve label. The compartment behind this page measures **95.4% neural** against
    a held-out annotation (notebook 06).

    > **{_n - _pass} of {_n} rows come from donor-dominated clusters.** They are shown,
    > not dropped, and flagged throughout. A GO term enriched only in a cluster that is
    > 83% one donor is a statement about that donor.
    """)
    return


@app.cell
def _filters(enrich_df, mo):
    _libs = sorted(enrich_df["gene_set_library"].dropna().unique())
    _types = sorted(enrich_df["cell_type"].dropna().unique())
    lib_select = mo.ui.multiselect(options=_libs, value=_libs, label="Gene-set library")
    type_select = mo.ui.multiselect(options=_types, value=_types, label="Nerve cell type")
    qc_only = mo.ui.checkbox(value=True, label="QC-passing clusters only")
    max_padj = mo.ui.slider(0.001, 0.25, value=0.05, step=0.001,
                            label="Max adjusted p-value")
    return lib_select, max_padj, qc_only, type_select


@app.cell
def _browser(enrich_df, lib_select, max_padj, mo, qc_only, type_select):
    _d = enrich_df[
        enrich_df["gene_set_library"].isin(lib_select.value or [])
        & enrich_df["cell_type"].isin(type_select.value or [])
        & (enrich_df["Adjusted P-value"] <= max_padj.value)
    ]
    if qc_only.value:
        _d = _d[_d["batch_qc_pass"]]
    _cols = ["cluster", "cell_type", "Term", "Adjusted P-value", "Overlap",
             "dominant_sample_fraction", "n_contributing_samples", "batch_qc_pass"]
    mo.vstack([
        mo.md("---\n## Panel A — Term browser"),
        mo.hstack([lib_select, type_select], gap=2),
        mo.hstack([max_padj, qc_only], gap=2),
        mo.md(f"### Showing **{len(_d):,}** of {len(enrich_df):,} terms"),
        mo.ui.table(_d[[c for c in _cols if c in _d.columns]].sort_values(
            "Adjusted P-value"), selection=None, page_size=20),
        mo.callout(mo.md(
            "**`Overlap` is the honest column here.** A term like `3/17` means three "
            "of the cluster's marker genes fall in a 17-gene set — significant by the "
            "hypergeometric test, and still three genes. Adjusted p-value orders terms "
            "within a cluster; it is not an effect size and does not make a 3-gene "
            "overlap a strong claim.\n\n"
            "**Enrichment is computed over each cluster's marker genes**, so it "
            "inherits whatever the differential test produced. A cluster with few "
            "distinctive markers yields few terms — absence of enrichment is not "
            "absence of biology."
        ), kind="info"),
    ])
    return


@app.cell
def _profile_header(mo):
    mo.md("""
    ---
    ## Panel B — Which clusters share a profile
    Jaccard overlap of significant term sets. Clusters sharing terms are often the
    same cell type split by depth or donor rather than by biology.
    """)
    return


@app.cell
def _profile_panel(enrich_df, max_padj, mo, np, plt, qc_only):
    _d = enrich_df[enrich_df["Adjusted P-value"] <= max_padj.value]
    if qc_only.value:
        _d = _d[_d["batch_qc_pass"]]
    _sets = {k: set(g["Term"]) for k, g in _d.groupby("cluster_key") if len(g)}
    _keys = sorted(_sets, key=lambda x: (len(x), x))

    if len(_keys) < 2:
        mo.stop(True, mo.callout(mo.md(
            "Fewer than two clusters survive the current filters — relax the "
            "p-value threshold or untick *QC-passing only* to compare profiles."
        ), kind="warn"))

    _m = np.zeros((len(_keys), len(_keys)))
    for _i, _a in enumerate(_keys):
        for _j, _b in enumerate(_keys):
            _u = _sets[_a] | _sets[_b]
            _m[_i, _j] = (len(_sets[_a] & _sets[_b]) / len(_u)) if _u else 0.0

    _fig, _ax = plt.subplots(figsize=(1.6 + 0.42 * len(_keys),
                                      1.4 + 0.38 * len(_keys)))
    _im = _ax.imshow(_m, cmap="viridis", vmin=0, vmax=1)
    _labels = [f"c{k}" for k in _keys]
    _ax.set_xticks(range(len(_keys)), _labels, rotation=90, fontsize=7)
    _ax.set_yticks(range(len(_keys)), _labels, fontsize=7)
    _ax.set_title("Jaccard overlap of significant term sets")
    _fig.colorbar(_im, ax=_ax, shrink=0.75)
    _fig.tight_layout()

    # Most-similar off-diagonal pairs.
    _pairs = [(_keys[i], _keys[j], _m[i, j])
              for i in range(len(_keys)) for j in range(i + 1, len(_keys))]
    _pairs.sort(key=lambda t: -t[2])
    _top = _pairs[:5]
    mo.vstack([
        mo.center(_fig),
        mo.md("**Most similar cluster pairs**\n\n" + "\n".join(
            f"- `c{a}` ↔ `c{b}` — Jaccard {v:.2f} "
            f"({len(_sets[a] & _sets[b])} shared terms)" for a, b, v in _top)
            if _top else "*No pairs to compare.*"),
        mo.callout(mo.md(
            "**High overlap is a question, not a conclusion.** Two clusters sharing "
            "most of their terms may be one population the clustering split — by "
            "sequencing depth, by donor, or by resolution — rather than two biological "
            "states. Check their donor composition in notebook 06 Panel E before "
            "treating them as distinct.\n\n"
            "**Jaccard on term sets ignores significance and direction**: a term "
            "scraping in at p=0.049 counts the same as one at p=1e-12, and GO sets "
            "overlap heavily by construction, which inflates similarity between any "
            "two clusters enriched for broadly neural themes."
        ), kind="info"),
    ])
    return


@app.cell
def _markers_header(mo):
    mo.md("""
    ---
    ## Panel C — The marker genes behind the terms
    Enrichment is downstream of these. If a cluster's markers look wrong, its terms
    are wrong regardless of their p-values.
    """)
    return


@app.cell
def _markers_panel(markers_df, mo, qc_only):
    _d = markers_df.copy()
    if qc_only.value:
        _d = _d[_d["batch_qc_pass"]]
    _top = (_d.sort_values(["cluster_key", "scores"], ascending=[True, False])
            .groupby("cluster_key").head(8))
    _cols = ["cluster", "cell_type", "gene_symbol", "scores", "logfoldchanges",
             "pvals_adj", "dominant_sample_fraction", "batch_qc_pass"]
    mo.vstack([
        mo.ui.table(_top[[c for c in _cols if c in _top.columns]].round(4),
                    selection=None, page_size=20),
        mo.callout(mo.md(
            "**Top 8 markers per cluster by score.** `pvals_adj` of exactly `0.0` is a "
            "floating-point underflow, not a p-value of zero — with tens of thousands "
            "of cells the Wilcoxon statistic saturates, so ranking has to come from "
            "`scores` and `logfoldchanges`.\n\n"
            "These are one-vs-rest markers **within the nerve compartment**, not "
            "against the whole cohort. A gene marking a cluster here is distinctive "
            "among nerve cells; it may be unremarkable cohort-wide."
        ), kind="info"),
    ])
    return


@app.cell
def _footer(dataset, enrich_df, mo):
    _fail = int((~enrich_df["batch_qc_pass"]).sum())
    mo.md(f"""
    ---
    **Cohort.** `{dataset}`. Set `GBM_DATASET` to switch arms.

    **Limits.**
    1. **{_fail} of {len(enrich_df)} enrichment rows come from donor-dominated
       clusters.** Shown and flagged, never dropped — but a term from a cluster that is
       one patient is a fact about that patient.
    2. Enrichment runs on marker genes, so it inherits the differential test's
       decisions. Few markers means few terms, which is not the same as no biology.
    3. GO gene sets overlap heavily; the Jaccard view in Panel B is inflated by that
       and should be read as a prompt to inspect, not a clustering.
    4. `astrocyte`, `opc`, generic `neuron` and `ependymal` are masked out of the nerve
       compartment by decision, so no cluster here can be enriched for them. Their
       absence is a choice made upstream, not a result.

    **FAIR.** Inputs from `ds_nerve_cell_heterogeneity` (GSEA over MSigDB GO sets) and
    `ds_annotate_cluster_qc` (donor purity and `batch_qc_pass` joined on).
    """)
    return


if __name__ == "__main__":
    app.run()
