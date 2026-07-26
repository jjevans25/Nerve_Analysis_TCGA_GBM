"""
Marimo reactive notebook: TME x nerve x immune ligand-receptor exploration for the
CELLxGENE Census GBM replication cohort (gbm_cellxgene_56c4912d).
Addresses: Which nerve cell types interface with which immune subtypes in an
independent 169-donor cohort, and do the curated lead axes from the 17-sample
reference cohort replicate there?
Source rules: workflow/rules/datasets.smk -> ds_nerve_tumor_immune_interaction,
ds_annotate_cluster_qc, ds_nerve_cluster_annotations, ds_cohort_concordance.

Sibling of notebooks/04_tme_nerve_immune_explorer.py, which covers the pinned
v1.3.0 reference cohort. Deliberately a separate file rather than a cohort switch:
the two cohorts differ in what exists (no clinical metadata here, no curated
target list, 169 donors instead of 17), so the panels are not the same.

REQUIRES the 2026-07-26 normalization fix. Before it, this cohort's LIANA tables
were computed on raw UMI counts and every row had an empty specificity_rank.
See markdowns/blocker_census_liana_raw_counts.md.
"""

import marimo

__generated_with = "0.23.1"
app = marimo.App(
    width="wide",
    app_title="GBM Census - Nerve x Immune Replication Explorer",
)


@app.cell
def _imports():
    import hashlib
    import json
    import os
    from datetime import datetime
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import yaml

    return Path, datetime, hashlib, json, mo, np, os, pd, plt, sns, yaml


@app.cell
def _hint_helper(mo):
    """Reusable click-to-expand inline definition -- uses native <details>."""

    def hint(label: str, tip: str):
        return mo.md(
            f"<details style='display:inline-block; margin:0 1em 0 0;'>"
            f"<summary style='cursor:pointer; user-select:none;'>{label} ⓘ</summary>"
            f"<div style='font-size:0.9em; color:#555; margin-top:0.25em;'>{tip}</div>"
            f"</details>"
        )

    return (hint,)


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

    tables_dir = project_root / config["dirs"]["tables"]
    ds_tables = tables_dir / dataset
    ds_prov = project_root / config["dirs"]["provenance"] / dataset

    paths = {
        "interactions": ds_tables / "nerve_tumor_immune_interactions_with_qc.csv",
        "top_pairs": ds_tables / "nerve_tumor_immune_top_pairs_with_qc.csv",
        "nerve_annotations": ds_tables / "nerve_cluster_annotations.csv",
        "immune_annotations": ds_tables / "immune_cluster_annotations.csv",
        "nerve_purity": ds_tables / "nerve_cluster_sample_purity.csv",
        "immune_cluster_purity": ds_tables / "immune_cluster_sample_purity.csv",
        "immune_purity": ds_tables / "immune_subtype_sample_purity.csv",
        "annotation_summary": ds_tables / "annotation_summary.csv",
        "concordance": ds_tables / "cohort_concordance_summary.json",
        "shared_pairs": ds_tables / "cohort_concordance_shared_pairs.csv",
        "lr_provenance": ds_prov / "nerve_tumor_immune_interaction_provenance.json",
        # Reference-cohort input: the curated shortlist has no per-cohort twin.
        "reference_lead_targets": tables_dir / "nerve_crosstalk_lead_targets.csv",
    }

    _missing = [p for p in paths.values() if not p.exists()]
    if _missing:
        mo.stop(
            True,
            mo.callout(
                mo.md(
                    f"Required artifacts not found for cohort `{dataset}`:\n\n"
                    + "\n".join(f"- `{p}`" for p in _missing)
                    + "\n\nBuild with:\n\n```\nsnakemake --use-conda --cores all "
                    "--rerun-triggers mtime -- \\\n"
                    f"  results/tables/{dataset}/nerve_cluster_annotations.csv \\\n"
                    f"  results/tables/{dataset}/cohort_concordance_summary.json\n```"
                ),
                kind="danger",
            ),
        )
    return dataset, ds_tables, paths


@app.cell
def _header(dataset, mo):
    mo.md(f"""
    # GBM Census — Nerve × Immune Replication Explorer

    **Cohort: `{dataset}`** — the CELLxGENE Census GBM 10x replication cohort
    (169 donors), independent of the 17-sample TCGA reference.

    Companion to `notebooks/04_tme_nerve_immune_explorer.py`, which covers the
    reference cohort. Kept separate because the two cohorts do not offer the same
    evidence: this one has **no usable clinical metadata** (its `gdc_clinical.tsv`
    is a generated stub) and **no curated target list** of its own.

    > **This notebook depends on the 2026-07-26 normalization fix.** Before it, the
    > LIANA tables here were computed on raw UMI counts — `X.max() = 53027` where
    > log1p data peaks near 9 — and every one of the 71,189 rows had an empty
    > `specificity_rank`, with 13,492 infinite log-fold-changes. Those tables were
    > void. The fix normalizes to counts-per-10k + log1p before the LIANA call and
    > adds a hard guard that refuses to run on an unnormalized matrix.
    """)
    return


@app.cell
def _load_tables(json, paths, pd):
    """Read every input and join cluster identity onto the LR rows."""
    _read = dict(low_memory=False)
    interactions_df = pd.read_csv(paths["interactions"], **_read)
    top_pairs_df = pd.read_csv(paths["top_pairs"], **_read)
    nerve_ann_df = pd.read_csv(paths["nerve_annotations"])
    immune_ann_df = pd.read_csv(paths["immune_annotations"])
    nerve_purity_df = pd.read_csv(paths["nerve_purity"])
    immune_cluster_purity_df = pd.read_csv(paths["immune_cluster_purity"])
    immune_purity_df = pd.read_csv(paths["immune_purity"])
    annotation_df = pd.read_csv(paths["annotation_summary"])
    shared_pairs_df = pd.read_csv(paths["shared_pairs"])
    lead_targets_df = pd.read_csv(paths["reference_lead_targets"])
    with open(paths["concordance"]) as _f:
        concordance = json.load(_f)
    with open(paths["lr_provenance"]) as _f:
        lr_prov = json.load(_f)

    # `label` is written by nerve_cluster_annotations.py as
    # f"c{cluster} | {dominant cell type} | {argmax canonical score}".
    _parts = nerve_ann_df["label"].astype(str).str.split("|", expand=True)
    nerve_ann_df = nerve_ann_df.copy()
    nerve_ann_df["nerve_cell_type"] = _parts[1].str.strip()
    nerve_ann_df["nerve_score_type"] = _parts[2].str.strip()
    nerve_ann_df["nerve_type_agrees"] = (
        nerve_ann_df["nerve_cell_type"] == nerve_ann_df["nerve_score_type"]
    )
    nerve_ann_df["cluster_key"] = nerve_ann_df["cluster"].astype(str)

    _bool_map = {True: True, "True": True, False: False, "False": False}
    _ann_cols = [
        "cluster_key", "nerve_cell_type", "nerve_score_type", "nerve_type_agrees",
        "label", "interpretation", "top_markers",
    ]

    def _prepare(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["nerve_cluster"] = out["nerve_cluster"].fillna("").astype(str)
        out["immune_subtype"] = out["immune_subtype"].fillna("").astype(str)
        for c in ("batch_qc_pass", "nerve_batch_qc_pass", "immune_batch_qc_pass"):
            if c in out.columns:
                out[c] = out[c].map(_bool_map)
        out["batch_qc_pass"] = out["batch_qc_pass"].fillna(True).astype(bool)
        out["cluster_key"] = (
            out["nerve_cluster"].str.replace("nerve_c", "", regex=False).str.strip()
        )
        out = out.merge(nerve_ann_df[_ann_cols], on="cluster_key", how="left")
        out["nerve_cell_type"] = out["nerve_cell_type"].fillna("")
        out["nerve_label"] = out["label"].fillna("")
        out["nerve_interpretation"] = out["interpretation"].fillna("")
        out["pair_unit"] = [
            (nl or nc) if cp == "nerve-tumor"
            else (isub if cp == "immune-tumor" else f"{nl or nc} × {isub}")
            for cp, nc, nl, isub in zip(
                out["compartment_pair"].astype(str), out["nerve_cluster"],
                out["nerve_label"], out["immune_subtype"],
            )
        ]
        return out.drop(columns=["label", "interpretation"])

    interactions_df = _prepare(interactions_df)
    top_pairs_df = _prepare(top_pairs_df)
    return (
        annotation_df,
        concordance,
        immune_ann_df,
        immune_cluster_purity_df,
        immune_purity_df,
        interactions_df,
        lead_targets_df,
        lr_prov,
        nerve_ann_df,
        nerve_purity_df,
        shared_pairs_df,
        top_pairs_df,
    )


@app.cell
def _integrity_banner(interactions_df, lr_prov, mo):
    """Confirm at read time that the normalization fix is present in this table."""
    _p = lr_prov["parameters"]
    _spec_missing = int(interactions_df["specificity_rank"].isna().sum())
    _inf = int(
        interactions_df["lr_logfc"].isin([float("inf"), float("-inf")]).sum()
    )
    _n = len(interactions_df)

    _meta = mo.md(
        f"- LIANA run: **{str(lr_prov.get('created_at', ''))[:10]}** · "
        f"tumor **{int(_p['n_tumor_cells']):,}** · nerve **{int(_p['n_nerve_cells']):,}** "
        f"({int(_p['n_nerve_clusters'])} clusters) · immune "
        f"**{int(_p['n_immune_cells']):,}** ({int(_p['n_immune_subtypes'])} subtypes)\n"
        f"- `n_perms={_p['n_perms']}`, `expr_prop={_p['expr_prop']}`, "
        f"resource `{_p['resource_name']}`"
    )

    if _spec_missing == 0 and _inf == 0:
        _view = mo.callout(
            mo.md(
                f"**Normalization verified.** All {_n:,} rows carry a "
                f"`specificity_rank` and none has an infinite `lr_logfc` — the "
                f"signature of the raw-counts defect is absent."
            ),
            kind="success",
        )
    else:
        _view = mo.callout(
            mo.md(
                f"**Do not trust these results.** {_spec_missing:,} of {_n:,} rows have "
                f"no `specificity_rank` and {_inf:,} have an infinite `lr_logfc`. That is "
                f"the raw-counts signature described in "
                f"`markdowns/blocker_census_liana_raw_counts.md`. Re-run with the fix:\n\n"
                f"```\nsnakemake --use-conda --cores all --rerun-triggers mtime \\\n"
                f"  --forcerun ds_nerve_tumor_immune_interaction -- <concordance target>\n```"
            ),
            kind="danger",
        )
    mo.vstack([mo.md("## Provenance and integrity"), _meta, _view])
    return


@app.cell
def _census_header(mo):
    mo.md("""
    ---
    ## Panel A — What this cohort is made of, and what was modelled
    """)
    return


@app.cell
def _compartment_census(annotation_df, interactions_df, lr_prov, mo, pd, plt):
    """Cohort annotation census vs the compartments the LR analysis actually saw."""
    _p = lr_prov["parameters"]
    _groups = set(interactions_df["source"].astype(str)) | set(
        interactions_df["target"].astype(str)
    )
    _immune_groups = {g.removeprefix("immune_") for g in _groups if g.startswith("immune_")}

    _nerve_feeders = {
        "neuron", "excitatory_neuron", "inhibitory_neuron",
        "astrocyte", "oligodendrocyte", "ependymal", "opc",
    }

    def _compartment(ct: str) -> str:
        if ct in _nerve_feeders:
            return "nerve"
        if ct == "microglia":
            return "immune"
        return "not modelled"

    _c = annotation_df.copy()
    _c["compartment"] = _c["cell_type_predicted"].map(_compartment)
    _c["modelled"] = _c["compartment"] != "not modelled"
    _c["pct_of_cohort"] = (100 * _c["n_cells"] / _c["n_cells"].sum()).round(2)
    census_df = _c[["cell_type_predicted", "n_cells", "pct_of_cohort",
                    "compartment", "modelled"]].sort_values(
        "n_cells", ascending=False).reset_index(drop=True)

    _total = int(census_df["n_cells"].sum())
    _unmod = census_df.loc[~census_df["modelled"]]
    _n_unmod = int(_unmod["n_cells"].sum())

    _fig, _ax = plt.subplots(figsize=(8, 3.4))
    _colors = {"nerve": "#4C72B0", "immune": "#DD8452", "not modelled": "#BBBBBB"}
    _ax.barh(census_df["cell_type_predicted"][::-1], census_df["n_cells"][::-1],
             color=[_colors[c] for c in census_df["compartment"][::-1]])
    _ax.set_xlabel("cells (marker-argmax annotation)")
    _ax.set_title("Census cohort annotation census — grey = never enters the LR analysis")
    _fig.tight_layout()

    _compartment_tbl = pd.DataFrame({
        "compartment": ["tumor (CNV-malignant)", "nerve", "immune", "— total modelled"],
        "cells_in_LR_run": [
            int(_p["n_tumor_cells"]), int(_p["n_nerve_cells"]),
            int(_p["n_immune_cells"]),
            int(_p["n_tumor_cells"]) + int(_p["n_nerve_cells"]) + int(_p["n_immune_cells"]),
        ],
        "groups_in_LR_table": [
            1, int(_p["n_nerve_clusters"]), len(_immune_groups), "",
        ],
    })

    mo.vstack([
        mo.center(_fig),
        mo.hstack([
            mo.vstack([mo.md("**Annotation census**"),
                       mo.ui.table(census_df, selection=None)]),
            mo.vstack([mo.md("**Compartments the LR analysis saw**"),
                       mo.ui.table(_compartment_tbl, selection=None)]),
        ], widths=[3, 2], gap=2),
        mo.callout(
            mo.md(
                f"**{_n_unmod:,} of {_total:,} annotated cells "
                f"({100 * _n_unmod / _total:.1f}%) never enter any interaction result** — "
                f"a much larger gap than the reference cohort's 1.5%.\n\n"
                + "\n".join(f"- **{r.cell_type_predicted}** — {r.n_cells:,} cells"
                            for r in _unmod.itertuples())
                + "\n\nTwo things differ sharply from the reference cohort and are worth "
                "knowing before reading any result below:\n\n"
                "- **The excluded `t_cell` population is large.** `immune_cell_subset` takes "
                "only cells labelled `microglia`, then sub-clusters T cells back out of that "
                "blob. So the T cells in the LR analysis are the ones found *inside* the "
                "myeloid blob, while a separately-annotated T-cell population of this size "
                "sits outside it entirely. Treat `immune_t_cell` results as a subset of the "
                "cohort's T cells, not all of them.\n"
                "- **This cohort's annotation resolves only 5 cell types** against the "
                "reference's 10 — no `endothelial`, `opc`, `ependymal` or `tumor_gbm` "
                "category exists here at all. That is itself a suspected artifact of "
                "marker scoring on raw counts (`mean_confidence` 5.96–15.29 here vs "
                "0.11–0.73 in the reference); it is a **separate, unfixed defect** from the "
                "LIANA one — see `markdowns/blocker_census_annotation_scoring.md`."
            ),
            kind="warn",
        ),
    ])
    return (census_df,)


@app.cell
def _purity_header(mo):
    mo.md("""
    ---
    ## Panel B — Patient diversity behind each group
    """)
    return


@app.cell
def _patient_purity(immune_cluster_purity_df, mo, nerve_purity_df, np, plt):
    """Per-cluster patient-dominance for both compartments."""
    def _panel(ax, df, title):
        d = df.sort_values("cluster").copy()
        x = np.arange(len(d))
        colors = ["#C44E52" if not bool(p) else "#4C72B0" for p in d["pass_overall"]]
        ax.bar(x, d["dominant_sample_fraction"], color=colors)
        ax.axhline(0.5, ls="--", lw=1, color="black")
        ax.set_xticks(x)
        ax.set_xticklabels([str(c) for c in d["cluster"]], rotation=90, fontsize=7)
        ax.set_ylabel("dominant sample fraction")
        ax.set_title(title)
        ax2 = ax.twinx()
        ax2.plot(x, d["n_contributing_samples"], color="#55A868", marker="o", ms=2.5, lw=1)
        ax2.set_ylabel("# contributing samples", color="#55A868")

    _fig, _axes = plt.subplots(2, 1, figsize=(11, 7))
    _panel(_axes[0], nerve_purity_df,
           "Nerve clusters — patient dominance (red = fails batch QC)")
    _panel(_axes[1], immune_cluster_purity_df,
           "Immune Leiden clusters — patient dominance")
    _fig.tight_layout()

    _n_fail = int((~nerve_purity_df["pass_overall"].astype(bool)).sum())
    mo.vstack([
        mo.center(_fig),
        mo.callout(
            mo.md(
                f"**{_n_fail} of {len(nerve_purity_df)} nerve clusters fail batch QC** "
                f"(dominant-sample fraction ≥ 0.5 or too few contributing donors). "
                f"Their rows are marked `*` and can be dropped in Panel C.\n\n"
                "*Why this is a bar chart and not the `cluster × patient` heatmap used in "
                "notebook 04: this cohort has 169 donors, so that matrix is 27 × 169 and "
                "unreadable. The purity summary carries the same verdict in a legible form.*"
            ),
            kind="info",
        ),
    ])
    return


@app.cell
def _browser_header(mo):
    mo.md("""
    ---
    ## Panel C — Interaction browser, labelled by biology
    """)
    return


@app.cell
def _filters(interactions_df, mo):
    """Reactive filter widgets -- drive Panels C and D."""
    _pairs = sorted(interactions_df["compartment_pair"].astype(str).unique())
    _dirs = sorted(interactions_df["direction"].astype(str).unique())
    _types = sorted(t for t in interactions_df["nerve_cell_type"].unique() if t)

    source_toggle = mo.ui.switch(value=True, label="Use top_pairs (off = full interactions)")
    pair_select = mo.ui.radio(options=_pairs, value=_pairs[0], label="Compartment interface")
    direction_select = mo.ui.multiselect(options=_dirs, value=_dirs, label="Direction(s)")
    celltype_select = mo.ui.multiselect(
        options=_types, value=_types, label="Nerve cell type(s)")
    magnitude_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="magnitude_rank ≤")
    pval_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="cellphone_pvals ≤")
    spec_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=1.0, step=0.01, label="specificity_rank ≤")
    ligand_search = mo.ui.text(label="Ligand contains", placeholder="e.g. NLGN1")
    receptor_search = mo.ui.text(label="Receptor contains", placeholder="e.g. NRXN")
    hide_qc_fail = mo.ui.checkbox(value=False, label="Hide QC-failing rows")
    return (
        celltype_select, direction_select, hide_qc_fail, ligand_search,
        magnitude_slider, pair_select, pval_slider, receptor_search,
        source_toggle, spec_slider,
    )


@app.cell
def _show_filters(
    celltype_select, direction_select, hide_qc_fail, hint, ligand_search,
    magnitude_slider, mo, pair_select, pval_slider, receptor_search,
    source_toggle, spec_slider,
):
    mo.vstack([
        mo.md("### Filters"),
        mo.hstack([source_toggle, hide_qc_fail], gap=2),
        pair_select,
        direction_select,
        celltype_select,
        mo.hstack([magnitude_slider, pval_slider, spec_slider], gap=2),
        mo.hstack([
            hint("magnitude_rank",
                 "LIANA aggregate magnitude rank across methods (RRA). Lower = stronger."),
            hint("cellphone_pvals",
                 "CellPhoneDB permutation p-value (~1000 permutations). Lower = more specific."),
            hint("specificity_rank",
                 "LIANA aggregate specificity rank. Lower = more group-specific. "
                 "Empty for every row before the 2026-07-26 normalization fix — "
                 "default 1.0 keeps all rows."),
        ], gap=2),
        mo.hstack([ligand_search, receptor_search], gap=2),
    ], gap=1)
    return


@app.cell
def _filtered_view(
    celltype_select, direction_select, hide_qc_fail, interactions_df, ligand_search,
    magnitude_slider, pair_select, pval_slider, receptor_search, source_toggle,
    spec_slider, top_pairs_df,
):
    """Apply every filter reactively to the chosen source dataframe."""
    _src = top_pairs_df if source_toggle.value else interactions_df
    _df = _src[_src["compartment_pair"].astype(str) == pair_select.value].copy()

    _dirs = direction_select.value or list(_df["direction"].astype(str).unique())
    _df = _df[_df["direction"].astype(str).isin(_dirs)]

    if celltype_select.value is not None:
        _keep = set(celltype_select.value)
        _df = _df[(_df["nerve_cell_type"] == "") | (_df["nerve_cell_type"].isin(_keep))]

    if hide_qc_fail.value:
        _df = _df[_df["batch_qc_pass"].astype(bool)]

    _df = _df[
        (_df["magnitude_rank"] <= magnitude_slider.value)
        & (_df["cellphone_pvals"] <= pval_slider.value)
        & (_df["specificity_rank"].fillna(1.0) <= spec_slider.value)
    ]

    if ligand_search.value:
        _df = _df[_df["ligand_complex"].str.contains(
            ligand_search.value, case=False, na=False)]
    if receptor_search.value:
        _df = _df[_df["receptor_complex"].str.contains(
            receptor_search.value, case=False, na=False)]

    filtered_df = _df.sort_values("magnitude_rank").reset_index(drop=True)
    return (filtered_df,)


@app.cell
def _show_filtered(filtered_df, mo, pair_select):
    _cols = [
        "nerve_label", "nerve_cell_type", "immune_subtype", "direction",
        "ligand_complex", "receptor_complex", "lrscore", "magnitude_rank",
        "specificity_rank", "cellphone_pvals", "batch_qc_pass", "nerve_type_agrees",
    ]
    mo.vstack([
        mo.md(f"### `{pair_select.value}` — filtered rows: **{len(filtered_df):,}**"),
        mo.ui.table(filtered_df[[c for c in _cols if c in filtered_df.columns]],
                    page_size=25, selection=None),
    ])
    return


@app.cell
def _matrix_header(mo):
    mo.md("""
    ---
    ## Panel D — Cell type × immune subtype interface matrix
    """)
    return


@app.cell
def _celltype_interface_matrix(filtered_df, mo, plt, sns):
    """Significant LR-pair counts by nerve cell type x immune subtype."""
    _df = filtered_df[
        (filtered_df["nerve_cell_type"] != "") & (filtered_df["immune_subtype"] != "")
    ]
    mo.stop(
        _df.empty,
        mo.callout(
            mo.md("*Needs rows with both a nerve and an immune side — select the "
                  "`immune-nerve` interface in Panel C.*"),
            kind="info",
        ),
    )

    _dirs = sorted(_df["direction"].astype(str).unique())
    _fig, _axes = plt.subplots(1, len(_dirs), figsize=(1 + 4.2 * len(_dirs), 3.4),
                               squeeze=False)
    for _i, _d in enumerate(_dirs):
        _counts = (
            _df[_df["direction"].astype(str) == _d]
            .groupby(["nerve_cell_type", "immune_subtype"], observed=True)
            .size().unstack("immune_subtype", fill_value=0)
        )
        _ax = _axes[0][_i]
        sns.heatmap(_counts, annot=True, fmt="d", cmap="mako_r", ax=_ax,
                    cbar_kws={"label": "# LR pairs"})
        _ax.set_title(_d)
        _ax.set_xlabel("immune subtype")
        _ax.set_ylabel("nerve cell type" if _i == 0 else "")
    _fig.tight_layout()

    _amb = sorted(_df.loc[_df["nerve_type_agrees"] == False, "nerve_label"].unique())  # noqa: E712
    mo.vstack([
        mo.center(_fig),
        mo.callout(
            mo.md(
                "**How to read:** pair counts surviving the Panel C thresholds, rolled up "
                "from Leiden clusters to their dominant cell type. Counts, not effect "
                "sizes — a cell type spread over more clusters accumulates more pairs."
                + ("\n\n**Ambiguous clusters in this rollup** (dominant cell type "
                   "disagrees with marker-score argmax): "
                   + ", ".join(f"`{a}`" for a in _amb) if _amb else "")
            ),
            kind="info",
        ),
    ])
    return


@app.cell
def _replication_header(mo):
    mo.md("""
    ---
    ## Panel E — Do the reference cohort's curated leads replicate here?
    The 40 hand-curated axes in `nerve_crosstalk_lead_targets.csv` were derived from
    the 17-sample reference cohort. This panel looks each one up in *this* cohort's
    LR table — the single most useful cross-cohort question available.
    """)
    return


@app.cell
def _replication_table(interactions_df, lead_targets_df, mo, pd):
    """Look up every curated reference axis in the Census LR table."""
    _lig = interactions_df["ligand_complex"].astype(str)
    _rec = interactions_df["receptor_complex"].astype(str)

    _rows = []
    for _r in lead_targets_df.itertuples():
        _a, _b = str(_r.mol_a), str(_r.mol_b)
        _hits = interactions_df[((_lig == _a) & (_rec == _b)) | ((_lig == _b) & (_rec == _a))]
        _sig = _hits[_hits["magnitude_rank"] <= 0.05]
        _rows.append({
            "axis": _r.axis,
            "target_class": _r.target_class,
            "ref_best_mag": _r.best_mag,
            "census_rows": len(_hits),
            "census_sig_rows": len(_sig),
            "census_best_mag": _hits["magnitude_rank"].min() if len(_hits) else None,
            "census_best_lrscore": _hits["lrscore"].max() if len(_hits) else None,
            "census_interfaces": ",".join(sorted(_hits["compartment_pair"].astype(str).unique())),
            "replicates": "yes" if len(_sig) else ("present" if len(_hits) else "absent"),
        })
    replication_df = pd.DataFrame(_rows)

    _n_yes = int((replication_df["replicates"] == "yes").sum())
    _n_present = int((replication_df["replicates"] == "present").sum())
    _n_absent = int((replication_df["replicates"] == "absent").sum())

    mo.vstack([
        mo.callout(
            mo.md(
                f"**{_n_yes} of {len(replication_df)} curated axes replicate** at "
                f"`magnitude_rank ≤ 0.05` in this cohort · {_n_present} present but not "
                f"significant · {_n_absent} absent entirely.\n\n"
                "`absent` usually means the ligand or receptor did not survive this "
                "cohort's gene filtering or `expr_prop` threshold, not that the biology "
                "is contradicted. Check `census_rows = 0` against the gene lists before "
                "concluding non-replication."
            ),
            kind="info",
        ),
        mo.ui.table(replication_df.sort_values(
            ["replicates", "census_best_mag"]).reset_index(drop=True),
            selection=None, page_size=20),
    ])
    return (replication_df,)


@app.cell
def _concordance_header(mo):
    mo.md("""
    ---
    ## Panel F — Whole-table concordance with the reference cohort
    """)
    return


@app.cell
def _concordance_panel(concordance, mo, np, plt, shared_pairs_df):
    """Cohort-level agreement, from ds_cohort_concordance."""
    _c = concordance
    _fig, _ax = plt.subplots(figsize=(5, 4.6))
    _ax.scatter(shared_pairs_df["reference_magnitude"],
                shared_pairs_df["dataset_magnitude"],
                s=6, alpha=0.35, edgecolors="none")
    _lim = [0, 1]
    _ax.plot(_lim, _lim, ls="--", lw=1, color="black")
    _ax.set_xlim(_lim)
    _ax.set_ylim(_lim)
    _ax.set_xlabel("reference magnitude_rank")
    _ax.set_ylabel("Census magnitude_rank")
    _ax.set_title(f"Shared significant LR pairs (n={len(shared_pairs_df):,})")
    _fig.tight_layout()

    mo.vstack([
        mo.hstack([
            mo.center(_fig),
            mo.md(f"""
**Jaccard overlap** {_c['jaccard_overlap']:.4f}
&nbsp;

**Spearman ρ (shared)** {_c['spearman_rho_shared']:.4f}
&nbsp;

- reference significant pairs: **{_c['n_reference_significant_pairs']:,}**
- Census significant pairs: **{_c['n_dataset_significant_pairs']:,}**
- shared: **{_c['n_shared_pairs']:,}** / union **{_c['n_union_pairs']:,}**
- `pval_max` = {_c['pval_max']}
            """),
        ], widths=[3, 2], gap=2),
        mo.callout(
            mo.md(
                "**These numbers superseded a void set.** Computed before the "
                "normalization fix they were Jaccard **0.4248** / ρ **0.5357**, comparing "
                "log1p-scored reference magnitudes against raw-count-scored Census ones. "
                "Both rose once the scales matched. Points near the diagonal are pairs "
                "the two cohorts rank alike; the axes are ranks, so **lower-left is "
                "stronger in both**."
            ),
            kind="info",
        ),
    ])
    return


@app.cell
def _export_button(mo):
    export_button = mo.ui.run_button(label="Export current filtered view → results/tables/")
    mo.center(export_button)
    return (export_button,)


@app.cell
def _export(
    celltype_select, dataset, datetime, direction_select, ds_tables, export_button,
    filtered_df, hashlib, hide_qc_fail, ligand_search, magnitude_slider, mo,
    pair_select, paths, pval_slider, receptor_search, source_toggle, spec_slider,
):
    """Persist the filtered view with FAIR provenance -- hashes every input."""
    mo.stop(not export_button.value, mo.md("*Click the button above to export.*"))

    _ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    _out = ds_tables / f"census_nerve_immune_exploration_{_ts}.csv"
    _prov = ds_tables / f"census_nerve_immune_exploration_{_ts}.provenance.txt"

    filtered_df.to_csv(_out, index=False)
    _hashes = {k: hashlib.sha256(p.read_bytes()).hexdigest()[:12] for k, p in paths.items()}
    _prov.write_text(
        f"timestamp: {_ts}\n"
        f"notebook: notebooks/05_census_nerve_immune_explorer.py\n"
        f"cohort: {dataset}\n"
        f"n_rows_exported: {len(filtered_df)}\n"
        "inputs_sha256_12:\n"
        + "".join(f"  {k}: {h}  ({paths[k].name})\n" for k, h in sorted(_hashes.items()))
        + "filters:\n"
        f"  source_toggle_top_pairs: {source_toggle.value}\n"
        f"  compartment_pair: {pair_select.value}\n"
        f"  directions: {direction_select.value}\n"
        f"  nerve_cell_types: {celltype_select.value}\n"
        f"  hide_qc_failing: {hide_qc_fail.value}\n"
        f"  magnitude_rank_le: {magnitude_slider.value}\n"
        f"  cellphone_pvals_le: {pval_slider.value}\n"
        f"  specificity_rank_le: {spec_slider.value}\n"
        f"  ligand_contains: {ligand_search.value!r}\n"
        f"  receptor_contains: {receptor_search.value!r}\n"
    )
    mo.callout(
        mo.md(f"**Exported** `{_out.name}` ({len(filtered_df):,} rows)\n\n"
              f"**Provenance** `{_prov.name}` — sha256[:12] of all {len(_hashes)} inputs "
              f"plus every filter value."),
        kind="success",
    )
    return


@app.cell
def _footer(dataset, mo):
    mo.md(f"""
    ---
    **Cohort.** `{dataset}` — CELLxGENE Census GBM 10x, 169 donors. The reference
    cohort lives in `notebooks/04_tme_nerve_immune_explorer.py`; run both and compare.

    **Limits carried by this cohort, all surfaced above.**
    1. No clinical metadata — `gdc_clinical.tsv` here is a generated stub, so the
       clinical-association panel from notebook 04 has no counterpart.
    2. The curated lead axes in Panel E come from the *reference* cohort; there is no
       per-cohort curation.
    3. Annotation resolves only 5 cell types (no endothelial / opc / ependymal), and
       a large separately-annotated T-cell population sits outside the immune
       compartment. Suspected raw-counts marker-scoring defect — unfixed, see
       `markdowns/blocker_census_annotation_scoring.md`.
    4. 11 of 35 nerve clusters fail batch QC in this cohort.

    **FAIR.** Inputs from `ds_nerve_tumor_immune_interaction`, `ds_annotate_cluster_qc`,
    `ds_nerve_cluster_annotations` and `ds_cohort_concordance`. Exports carry a
    `.provenance.txt` sidecar hashing every input.
    """)
    return


if __name__ == "__main__":
    app.run()
