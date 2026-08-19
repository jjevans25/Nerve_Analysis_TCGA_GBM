"""
Marimo reactive notebook: cohort QC and composition for the CELLxGENE Census GBM arms.
Addresses: how much data is there, how uneven is it across donors, and does any
donor's QC profile differ enough to distort a cohort-level result?

Census replacement for the archived notebooks/archive/01_explore_gbm_data.py, which
asked the same questions of the 17-sample TCGA reference. The scale difference is
the point: 170 donors against 17, with an author annotation to check composition
against. The reference version is archived because its downstream compartments were
defective and unfixable — see CHANGELOG 2026-08-19.

Per-donor rather than per-cell by design. The arm carries ~1.02M cells across 170
`*_qc_metrics.csv` files; a million-point scatter would render slowly and say less
than 170 summarised donors. Distribution shape within a donor is deliberately
reduced to median plus IQR — this notebook is about *between-donor* variation,
which is what a cohort-level analysis is exposed to.

Source rules: workflow/rules/datasets.smk -> ds_scrna_qc, ds_scrna_annotate.
"""

import marimo

__generated_with = "0.23.1"
app = marimo.App(width="wide", app_title="GBM Census — Cohort QC")


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

    qc_files = sorted(ds_tables.glob("*_qc_metrics.csv"))
    presence_files = sorted(ds_tables.glob("*_gene_presence.csv"))
    annotation_path = ds_tables / "annotation_summary.csv"

    if not qc_files or not annotation_path.exists():
        mo.stop(True, mo.callout(mo.md(
            f"QC artifacts not found for cohort `{dataset}` under `{ds_tables}`.\n\n"
            "```\nscripts/run_snakemake.sh \\\n"
            f"  results/tables/{dataset}/annotation_summary.csv \\\n"
            "  --use-conda --cores all --rerun-triggers mtime\n```"
        ), kind="danger"))
    return annotation_path, dataset, presence_files, qc_files


@app.cell
def _load_qc(np, pd, qc_files):
    """Summarise each donor's per-cell QC into one row. ~1.02M cells -> 170 rows."""
    _rows = []
    for _p in qc_files:
        _sample = _p.name.replace("_qc_metrics.csv", "")
        _d = pd.read_csv(_p)
        if not len(_d):
            continue
        _rows.append({
            "sample_id": _sample,
            "n_cells": len(_d),
            "median_genes": float(_d["n_genes_by_counts"].median()),
            "q25_genes": float(_d["n_genes_by_counts"].quantile(0.25)),
            "q75_genes": float(_d["n_genes_by_counts"].quantile(0.75)),
            "median_counts": float(_d["total_counts"].median()),
            "median_pct_mt": float(_d["pct_counts_mt"].median()),
            "max_pct_mt": float(_d["pct_counts_mt"].max()),
        })
    qc_df = pd.DataFrame(_rows).sort_values("n_cells", ascending=False).reset_index(drop=True)
    # Robust outlier flag: median absolute deviation, not standard deviation — a few
    # extreme donors would inflate an SD and hide themselves inside it.
    for _col in ("median_genes", "median_counts", "median_pct_mt"):
        _med = qc_df[_col].median()
        _mad = float(np.median(np.abs(qc_df[_col] - _med))) or float("nan")
        qc_df[f"{_col}_mad_z"] = (qc_df[_col] - _med) / (1.4826 * _mad)
    return (qc_df,)


@app.cell
def _header(dataset, mo, qc_df):
    mo.md(f"""
    # GBM Census — Cohort QC and Composition

    **Cohort: `{dataset}`** — **{len(qc_df)} donors**, **{int(qc_df['n_cells'].sum()):,}
    cells** post-QC.

    Census replacement for the archived `01_explore_gbm_data.py`, which asked the same
    questions of the 17-sample TCGA reference. Ten times the donors, and — unlike the
    reference — an independent author annotation to check composition against
    (notebook 06).

    Everything here is **per-donor**. The arm has ~1M cells; a million-point scatter
    would render slowly and tell you less than {len(qc_df)} summarised donors. The
    question a cohort analysis is actually exposed to is *between-donor* variation,
    not the shape of any one donor's distribution.
    """)
    return


@app.cell
def _scale_header(mo):
    mo.md("""
    ---
    ## Panel A — Cohort scale and how unevenly it is distributed
    A cohort of N donors is not N equal contributions. The concentration here bounds
    every downstream per-donor QC test.
    """)
    return


@app.cell
def _scale_panel(mo, np, plt, qc_df):
    _n = qc_df["n_cells"].to_numpy()
    _sorted = np.sort(_n)[::-1]
    _cum = np.cumsum(_sorted) / _sorted.sum()

    _fig, _ax = plt.subplots(1, 2, figsize=(11.5, 3.8))
    _ax[0].bar(range(len(_sorted)), _sorted, width=1.0, color="#2980b9")
    _ax[0].set_xlabel("donor (rank by cell count)")
    _ax[0].set_ylabel("cells")
    _ax[0].set_title("Cells per donor")
    _ax[1].plot(range(1, len(_cum) + 1), _cum, color="#c0392b", lw=1.8)
    _ax[1].axhline(0.5, ls=":", lw=1, color="black")
    _ax[1].set_xlabel("donors, largest first")
    _ax[1].set_ylabel("cumulative share of cells")
    _ax[1].set_ylim(0, 1)
    _ax[1].set_title("Concentration")
    _fig.tight_layout()

    _half = int(np.searchsorted(_cum, 0.5) + 1)
    _top10 = float(_cum[min(9, len(_cum) - 1)])
    mo.vstack([
        mo.center(_fig),
        mo.md(
            f"**{_half} of {len(_n)} donors supply half the cells.** The largest 10 "
            f"supply **{_top10:.0%}**. Median {int(np.median(_n)):,} cells per donor; "
            f"range {int(_n.min()):,}–{int(_n.max()):,}."
        ),
        mo.callout(mo.md(
            "**This is the ceiling on every donor-level QC claim downstream.** The "
            "batch-QC test used throughout this project asks whether a cluster draws "
            "from enough distinct donors; when the cohort itself is this concentrated, "
            "a cluster can look donor-diverse and still be dominated by a handful of "
            "large contributors.\n\n"
            "It is also why the capped arm exists. `subsample_per_donor` caps the "
            "largest donors so they cannot dominate the latent space, at the cost of "
            "discarding real cells — which is exactly why both arms are carried and "
            "every lead axis has to clear the bar in both."
        ), kind="info"),
    ])
    return


@app.cell
def _qc_header(mo):
    mo.md("""
    ---
    ## Panel B — Per-donor QC profile
    Flagged donors are **not** excluded anywhere in this pipeline. The flag is a
    prompt to check whether a result leans on them.
    """)
    return


@app.cell
def _qc_filters(mo):
    mad_threshold = mo.ui.slider(
        2.0, 6.0, value=3.5, step=0.5,
        label="Outlier threshold (robust MAD z-score)")
    return (mad_threshold,)


@app.cell
def _qc_panel(mad_threshold, mo, plt, qc_df):
    _t = mad_threshold.value
    _d = qc_df.copy()
    _d["is_outlier"] = (
        (_d["median_genes_mad_z"].abs() > _t)
        | (_d["median_counts_mad_z"].abs() > _t)
        | (_d["median_pct_mt_mad_z"].abs() > _t)
    )

    _fig, _ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for _i, (_col, _label) in enumerate([
        ("median_genes", "median genes / cell"),
        ("median_counts", "median UMIs / cell"),
        ("median_pct_mt", "median % mitochondrial"),
    ]):
        _ax[_i].scatter(_d["n_cells"], _d[_col], s=18, alpha=0.75,
                        c=["#c0392b" if o else "#2980b9" for o in _d["is_outlier"]],
                        edgecolor="white", linewidth=0.4)
        _ax[_i].set_xscale("log")
        _ax[_i].set_xlabel("cells in donor (log)")
        _ax[_i].set_ylabel(_label)
    _ax[0].set_title(f"Red = |MAD z| > {_t} on any metric")
    _fig.tight_layout()

    _out = _d[_d["is_outlier"]].sort_values("n_cells", ascending=False)
    _cols = ["sample_id", "n_cells", "median_genes", "median_counts",
             "median_pct_mt", "max_pct_mt"]
    mo.vstack([
        mo.center(_fig),
        mad_threshold,
        mo.md(f"**{len(_out)} of {len(_d)} donors flagged** at |MAD z| > {_t}."),
        mo.ui.table(_out[_cols].round(2) if len(_out) else _d[_cols].head(0),
                    selection=None, page_size=10),
        mo.callout(mo.md(
            "**Robust z-scores, not standard deviations.** A handful of extreme donors "
            "inflate an SD and then hide inside it; median absolute deviation does not "
            "move when the tail does.\n\n"
            "**Nothing here filters anything.** These donors are present in every "
            "compartment and every interaction table. A high mitochondrial fraction is "
            "as often a tissue-handling property of one sample as a quality failure, "
            "and dropping donors on that basis would silently narrow the cohort. Use "
            "this to ask whether a specific result leans on a flagged donor — Panel B "
            "of notebook 05 asks that question per cluster."
        ), kind="info"),
    ])
    return


@app.cell
def _composition_header(mo):
    mo.md("""
    ---
    ## Panel C — What the cohort is annotated as
    Pipeline `cell_type_predicted`, with its mean confidence. This is the pipeline's
    own call — notebook 06 checks it against the held-out Census annotation.
    """)
    return


@app.cell
def _composition_panel(annotation_path, mo, pd, plt):
    _a = pd.read_csv(annotation_path).sort_values("n_cells", ascending=False)
    _a["share"] = _a["n_cells"] / _a["n_cells"].sum()

    _fig, _ax = plt.subplots(figsize=(9, max(2.6, 0.34 * len(_a))))
    _c = _ax.barh(_a["cell_type_predicted"][::-1], _a["n_cells"][::-1],
                  color="#34495e", height=0.62)
    _ax.set_xlabel("cells")
    _ax.set_title("Predicted cell-type composition")
    for _bar, _conf in zip(_c, _a["mean_confidence"][::-1]):
        _ax.text(_bar.get_width() * 1.01, _bar.get_y() + _bar.get_height() / 2,
                 f"conf {_conf:.2f}", va="center", fontsize=6.5, color="#555")
    _ax.margins(x=0.16)
    _fig.tight_layout()

    _low = _a[_a["mean_confidence"] < 0.8]
    mo.vstack([
        mo.center(_fig),
        mo.ui.table(_a.round(4), selection=None, page_size=12),
        mo.callout(mo.md(
            "**These are the pipeline's own labels, so they cannot validate "
            "themselves.** `mean_confidence` is the marker-score margin behind each "
            "call, not a probability of being right. Notebook 06 does the actual check "
            "— cross-tabulating each compartment against the CELLxGENE author "
            "annotation, which no pipeline rule is allowed to read.\n\n"
            + (f"**{len(_low)} label(s) sit below 0.80 mean confidence** "
               f"({', '.join(_low['cell_type_predicted'].astype(str))}). Treat "
               f"compartment membership drawn from those as provisional."
               if len(_low) else
               "Every label carries ≥0.80 mean confidence.")
            + "\n\nComposition is also **not** compartment membership: compartments "
            "additionally drop CNV-malignant cells and mask several neural labels by "
            "decision, so these counts are upper bounds."
        ), kind="info"),
    ])
    return


@app.cell
def _presence_header(mo):
    mo.md("""
    ---
    ## Panel D — Marker detectability across donors
    Whether each canonical marker is detected at all, per donor. A marker missing in
    many donors cannot support a cell-type call made with it.
    """)
    return


@app.cell
def _presence_panel(mo, pd, plt, presence_files):
    if not presence_files:
        mo.stop(True, mo.md("*No `*_gene_presence.csv` files in this cohort.*"))
    _p = pd.concat([pd.read_csv(f) for f in presence_files], ignore_index=True)
    _p["present"] = _p["present"].astype(str).str.lower().isin(["true", "1"])
    _rate = (_p.groupby(["cell_type", "gene"])["present"]
             .mean().reset_index(name="detection_rate"))
    _pivot = _rate.pivot(index="gene", columns="cell_type", values="detection_rate")

    _fig, _ax = plt.subplots(figsize=(1.4 + 1.05 * _pivot.shape[1],
                                      max(2.6, 0.3 * _pivot.shape[0])))
    _im = _ax.imshow(_pivot.fillna(0).values, cmap="RdYlGn", vmin=0, vmax=1,
                     aspect="auto")
    _ax.set_xticks(range(_pivot.shape[1]), _pivot.columns, rotation=45, ha="right")
    _ax.set_yticks(range(_pivot.shape[0]), _pivot.index, fontsize=7)
    _ax.set_title(f"Fraction of {_p['sample_id'].nunique()} donors detecting each marker")
    _fig.colorbar(_im, ax=_ax, shrink=0.7, label="detection rate")
    _fig.tight_layout()

    _weak = _rate[_rate["detection_rate"] < 0.5].sort_values("detection_rate")
    mo.vstack([
        mo.center(_fig),
        mo.ui.table(_rate.round(3).sort_values("detection_rate"),
                    selection=None, page_size=12),
        mo.callout(mo.md(
            (f"**{len(_weak)} marker/cell-type pair(s) are detected in fewer than half "
             f"the donors** — lowest: "
             + ", ".join(f"`{r.gene}` ({r.cell_type}, {r.detection_rate:.0%})"
                         for r in _weak.head(4).itertuples())
             + ". A marker absent from most donors cannot carry a cell-type call: "
             "wherever it appears in a marker panel, the call rests on its "
             "companions."
             if len(_weak) else
             "**Every marker is detected in at least half the donors.**")
            + "\n\nDetection here means *non-zero in that donor*, which is a floor, "
            "not an expression level. A gene can be universally detected and still be "
            "useless for discriminating one cell type from another."
        ), kind="info" if _weak.empty else "warn"),
    ])
    return


@app.cell
def _footer(dataset, mo, qc_df):
    mo.md(f"""
    ---
    **Cohort.** `{dataset}` — {len(qc_df)} donors, {int(qc_df['n_cells'].sum()):,} cells.
    Set `GBM_DATASET` to switch arms.

    **Limits.**
    1. Per-donor summaries hide within-donor structure by construction. A donor whose
       cells are bimodal in depth shows one median here.
    2. QC flags are advisory and filter nothing. No donor is excluded anywhere in this
       pipeline on the basis of this page.
    3. `cell_type_predicted` composition is the pipeline's own output. The independent
       check is notebook 06.
    4. Cells here are already post-QC — `ds_scrna_qc` applied `min_genes`/`max_genes`
       and mitochondrial thresholds upstream, so this is the surviving population, not
       the raw one. The filtering those thresholds did is not visible from this page.

    **FAIR.** Inputs from `ds_scrna_qc` (per-donor `*_qc_metrics.csv`,
    `*_gene_presence.csv`) and `ds_scrna_annotate` (`annotation_summary.csv`).
    """)
    return


if __name__ == "__main__":
    app.run()
