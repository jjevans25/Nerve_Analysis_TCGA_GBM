"""
Marimo reactive notebook: TME context layer over the three-way nerve-tumor-immune
ligand-receptor analysis in TCGA GBM (reference cohort, v1.3.0).
Addresses: What is this tumor microenvironment actually made of, which parts of it
did the interaction analysis model, and which *cell types* -- not opaque Leiden ids
-- are talking to which immune subtypes?
Source rules: workflow/rules/immune.smk -> nerve_tumor_immune_interaction.py,
nerve_cells.smk -> {annotate_cluster_qc, nerve_cluster_annotations,
nerve_clinical_association}, annotation.smk -> scrna_annotate.py.

Companion to notebooks/03_nerve_tumor_immune_explorer.py, which stays the pure
ligand-receptor view. This notebook adds the biological and provenance context
that 03 structurally cannot show.
"""

import marimo

__generated_with = "0.23.1"
app = marimo.App(
    width="wide",
    app_title="GBM - TME x Nerve x Immune Context Explorer",
)


@app.cell
def _imports():
    import hashlib
    import json
    from datetime import datetime
    from pathlib import Path

    import duckdb
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import yaml

    return Path, datetime, duckdb, hashlib, json, mo, np, pd, plt, sns, yaml


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
def _load_config(Path, mo, yaml):
    """Resolve every input path from config.yaml -- no hardcoded absolute paths."""
    project_root = Path(__file__).parent.parent
    _config_path = project_root / "config" / "config.yaml"
    with open(_config_path) as _f:
        config = yaml.safe_load(_f)

    tables_dir = project_root / config["dirs"]["tables"]
    prov_dir = project_root / config["dirs"]["provenance"]

    paths = {
        "interactions": tables_dir / "nerve_tumor_immune_interactions_with_qc.csv",
        "top_pairs": tables_dir / "nerve_tumor_immune_top_pairs_with_qc.csv",
        "nerve_annotations": tables_dir / "nerve_cluster_annotations.csv",
        "immune_annotations": tables_dir / "immune_cluster_annotations.csv",
        "immune_composition": tables_dir / "immune_cluster_composition.csv",
        "nerve_purity": tables_dir / "nerve_cluster_sample_purity.csv",
        "immune_purity": tables_dir / "immune_subtype_sample_purity.csv",
        "annotation_summary": tables_dir / "annotation_summary.csv",
        "lead_targets": tables_dir / "nerve_crosstalk_lead_targets.csv",
        "clinical": tables_dir / "nerve_clinical_association.csv",
        "lr_provenance": prov_dir / "nerve_tumor_immune_interaction_provenance.json",
    }

    _missing = [p for p in paths.values() if not p.exists()]
    if _missing:
        mo.stop(
            True,
            mo.callout(
                mo.md(
                    "Required artifacts not found:\n\n"
                    + "\n".join(f"- `{p}`" for p in _missing)
                    + "\n\nThese are **pinned v1.3.0 artifacts** — the rules that "
                    "produced them are deliberately not defined while "
                    "`baseline.pinned` is set in `config/config.yaml`, because "
                    "`data/processed/nerve_cells.h5ad` was deleted and is "
                    "unreproducible.\n\n"
                    "If one is missing it has been moved or deleted, and it cannot "
                    "be regenerated: restore it from git or a backup. Run "
                    "`snakemake --use-conda verify_pinned_reference` to see exactly "
                    "what drifted."
                ),
                kind="danger",
            ),
        )
    return config, paths, tables_dir


@app.cell
def _header(mo):
    mo.md("""
    # GBM - TME × Nerve × Immune Context Explorer

    **Scope: the 17-sample TCGA-GBM reference cohort (v1.3.0) only.**

    This notebook is the *context* companion to
    `notebooks/03_nerve_tumor_immune_explorer.py`. 03 answers "which ligand-receptor
    pairs are strong?"; this one answers the questions that come before and after it:

    - **What is this TME made of, and how much of it did we actually model?**
    - **Which nerve *cell types* -- not Leiden ids -- interface with which immune subtypes?**
    - **How stale or patient-driven is the evidence behind any given pair?**

    > **The replication cohort is deliberately excluded.** The CELLxGENE Census cohort's
    > LIANA tables were computed on raw UMI counts rather than log1p-normalised data
    > (its run logged `X.max() = 53027.0`; the reference logged `7.762`). Every one of
    > its 71,189 rows has an empty `specificity_rank` and 13,492 have `lr_logfc = inf`.
    > See `markdowns/blocker_census_liana_raw_counts.md`.
    """)
    return


@app.cell
def _glossary(mo):
    """Column definitions -- collapsed by default."""
    mo.accordion(
        {
            "Glossary - columns added by this notebook": mo.md("""
    03's glossary covers the raw LIANA columns (`lrscore`, `magnitude_rank`,
    `specificity_rank`, `cellphone_pvals`, `batch_qc_pass`, ...). This notebook joins
    biological context on top of them:

    **`nerve_cell_type`** — the majority `cell_type_predicted` among the cluster's
    cells ("dominant"). The cluster's identity by per-cell vote; this is what the
    cell-type filter and the Panel D rollup use.

    **`nerve_score_type`** — argmax of the canonical marker-module scores for the
    cluster. A second, independent identity signal.

    **`nerve_type_agrees`** — whether the two above agree. `False` ⇒ ambiguous
    cluster identity; treat its cell-type rollup with care.

    **`nerve_label`** — `c{N} | dominant | score-argmax`, straight from
    `nerve_cluster_annotations.csv`. The human-readable cluster name.

    **`nerve_interpretation`** — one-sentence cluster summary including dominance
    percentage, top markers and best canonical score.

    **`modelled`** (Panel A) — whether an annotated cell type entered the LR
    analysis at all. `False` ⇒ present in the cohort but invisible to every
    interaction result.

    *Cluster identity comes from `nerve_cluster_annotations.py`, which builds
    `label = f"c{cluster} | {dominant} | {best_module}"`.*
    """)
        }
    )
    return


@app.cell
def _load_tables(json, paths, pd):
    """Read every input CSV and join biological identity onto the LR rows."""
    _read = dict(low_memory=False)
    interactions_df = pd.read_csv(paths["interactions"], **_read)
    top_pairs_df = pd.read_csv(paths["top_pairs"], **_read)
    nerve_ann_df = pd.read_csv(paths["nerve_annotations"])
    immune_ann_df = pd.read_csv(paths["immune_annotations"])
    immune_comp_df = pd.read_csv(paths["immune_composition"])
    nerve_purity_df = pd.read_csv(paths["nerve_purity"])
    immune_purity_df = pd.read_csv(paths["immune_purity"])
    annotation_df = pd.read_csv(paths["annotation_summary"])
    lead_targets_df = pd.read_csv(paths["lead_targets"])
    clinical_df = pd.read_csv(paths["clinical"])
    with open(paths["lr_provenance"]) as _f:
        lr_prov = json.load(_f)

    # --- Parse cluster identity out of the `label` column ----------------------
    # nerve_cluster_annotations.py writes: f"c{cluster} | {dominant} | {best_module}"
    _parts = nerve_ann_df["label"].astype(str).str.split("|", expand=True)
    nerve_ann_df = nerve_ann_df.copy()
    nerve_ann_df["nerve_cell_type"] = _parts[1].str.strip()
    nerve_ann_df["nerve_score_type"] = _parts[2].str.strip()
    nerve_ann_df["nerve_type_agrees"] = (
        nerve_ann_df["nerve_cell_type"] == nerve_ann_df["nerve_score_type"]
    )
    nerve_ann_df["cluster_key"] = nerve_ann_df["cluster"].astype(str)

    # --- Normalise the LR tables and join identity on -------------------------
    _bool_map = {True: True, "True": True, False: False, "False": False}
    _ann_cols = [
        "cluster_key",
        "nerve_cell_type",
        "nerve_score_type",
        "nerve_type_agrees",
        "label",
        "interpretation",
        "top_markers",
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
        # Left join: 25 clusters survive into the _with_qc table vs 27 annotated
        # (21 and 27 are dropped as artifacts) -- unmatched rows are expected.
        out = out.merge(nerve_ann_df[_ann_cols], on="cluster_key", how="left")
        out["nerve_cell_type"] = out["nerve_cell_type"].fillna("")
        out["nerve_label"] = out["label"].fillna("")
        out["nerve_interpretation"] = out["interpretation"].fillna("")
        # Biology-readable interacting unit, per interface.
        out["pair_unit"] = [
            (nl or nc) if cp == "nerve-tumor"
            else (isub if cp == "immune-tumor" else f"{nl or nc} × {isub}")
            for cp, nc, nl, isub in zip(
                out["compartment_pair"].astype(str),
                out["nerve_cluster"],
                out["nerve_label"],
                out["immune_subtype"],
            )
        ]
        return out.drop(columns=["label", "interpretation"])

    interactions_df = _prepare(interactions_df)
    top_pairs_df = _prepare(top_pairs_df)
    return (
        annotation_df,
        clinical_df,
        immune_ann_df,
        immune_comp_df,
        immune_purity_df,
        interactions_df,
        lead_targets_df,
        lr_prov,
        nerve_ann_df,
        nerve_purity_df,
        top_pairs_df,
    )


@app.cell
def _provenance_banner(immune_purity_df, lr_prov, mo):
    """Surface the staleness gap between the LR run and the current immune tables."""
    _p = lr_prov["parameters"]
    _lr_immune = int(_p["n_immune_cells"])
    _now_immune = int(immune_purity_df["n_cells"].sum())
    _created = str(lr_prov.get("created_at", ""))[:10]

    _rows = mo.md(
        f"- LIANA run date: **{_created}** &nbsp;·&nbsp; "
        f"tumor **{int(_p['n_tumor_cells']):,}** · "
        f"nerve **{int(_p['n_nerve_cells']):,}** ({int(_p['n_nerve_clusters'])} clusters) · "
        f"immune **{_lr_immune:,}** ({int(_p['n_immune_subtypes'])} subtypes)\n"
        f"- LIANA settings: `n_perms={_p['n_perms']}`, `expr_prop={_p['expr_prop']}`, "
        f"resource `{_p['resource_name']}` · "
        f"**{int(_p['n_lr_rows_total']):,}** LR rows, "
        f"**{int(_p['n_lr_rows_sig']):,}** significant"
    )

    if _lr_immune == _now_immune:
        _view = mo.vstack([_rows, mo.callout(
            mo.md("Immune compartment is consistent between the LR run and the "
                  "current annotation tables."), kind="success")])
    else:
        _view = mo.vstack([_rows, mo.callout(
            mo.md(
                f"**Provenance caveat -- the immune compartment was re-run after "
                f"the LR scores were computed.** The ligand-receptor scores below "
                f"derive from an immune subset of **{_lr_immune:,}** cells, but "
                f"`immune_cluster_annotations.csv` / `immune_subtype_sample_purity.csv` "
                f"(and Panel B) describe a later re-run of **{_now_immune:,}** cells "
                f"with a materially different subtype mix.\n\n"
                f"The `_with_qc` join is still valid because it keys on **subtype "
                f"name**, which is stable across both runs, and all subtypes pass "
                f"batch QC in both. But do **not** read a subtype's cell count in "
                f"Panel B as the sample size behind its LR score. Re-running "
                f"`nerve_tumor_immune_interaction` would resolve this -- it is "
                f"blocked by the deleted `nerve_cells.h5ad`."
            ),
            kind="warn")])
    mo.vstack([mo.md("## Provenance of the interaction evidence"), _view])
    return


@app.cell
def _compartment_census_header(mo):
    mo.md("""
    ---
    ## Panel A - What is this TME made of, and what did we model?
    """)
    return


@app.cell
def _compartment_census(annotation_df, interactions_df, lr_prov, mo, pd, plt):
    """Whole-cohort annotation census vs the three compartments the LR analysis saw."""
    _p = lr_prov["parameters"]

    # Which group labels actually appear anywhere in the interaction table?
    _groups = set(interactions_df["source"].astype(str)) | set(
        interactions_df["target"].astype(str)
    )
    _has_tumor = "malignant" in _groups
    _immune_groups = {g.removeprefix("immune_") for g in _groups if g.startswith("immune_")}

    # Map each annotated cell type to the compartment it feeds, if any.
    # nerve_cell_subset.py matches config nerve_cells.cell_types by substring after
    # lower-casing and replacing "_" with " "; immune_cell_subset.py takes exactly
    # cell_type_predicted == "microglia".
    _nerve_feeders = {
        "neuron", "excitatory_neuron", "inhibitory_neuron",
        "astrocyte", "oligodendrocyte", "ependymal",
    }
    _census = annotation_df.copy()

    def _compartment(ct: str) -> str:
        if ct in _nerve_feeders:
            return "nerve"
        if ct == "microglia":
            return "immune"
        return "not modelled"

    _census["compartment"] = _census["cell_type_predicted"].map(_compartment)
    _census["modelled"] = _census["compartment"] != "not modelled"
    _census["pct_of_cohort"] = (
        100 * _census["n_cells"] / _census["n_cells"].sum()
    ).round(2)
    census_df = _census[
        ["cell_type_predicted", "n_cells", "pct_of_cohort", "compartment", "modelled"]
    ].sort_values("n_cells", ascending=False).reset_index(drop=True)

    _total = int(census_df["n_cells"].sum())
    _unmodelled = census_df.loc[~census_df["modelled"]]
    _n_unmodelled = int(_unmodelled["n_cells"].sum())

    _modelled_cells = (
        int(_p["n_tumor_cells"]) + int(_p["n_nerve_cells"]) + int(_p["n_immune_cells"])
    )

    _fig, _ax = plt.subplots(figsize=(8, 4.2))
    _colors = {"nerve": "#4C72B0", "immune": "#DD8452", "not modelled": "#BBBBBB"}
    _ax.barh(
        census_df["cell_type_predicted"][::-1],
        census_df["n_cells"][::-1],
        color=[_colors[c] for c in census_df["compartment"][::-1]],
    )
    _ax.set_xlabel("cells (marker-argmax annotation)")
    _ax.set_title("Cohort annotation census — grey = never enters the LR analysis")
    _fig.tight_layout()

    _compartment_tbl = pd.DataFrame(
        {
            "compartment": ["tumor (CNV-malignant)", "nerve", "immune", "— total modelled"],
            "cells_in_LR_run": [
                int(_p["n_tumor_cells"]),
                int(_p["n_nerve_cells"]),
                int(_p["n_immune_cells"]),
                _modelled_cells,
            ],
            "groups_in_LR_table": [
                1 if _has_tumor else 0,
                int(_p["n_nerve_clusters"]),
                len(_immune_groups),
                "",
            ],
        }
    )

    mo.vstack(
        [
            mo.center(_fig),
            mo.hstack(
                [
                    mo.vstack([mo.md("**Annotation census**"),
                               mo.ui.table(census_df, selection=None, page_size=12)]),
                    mo.vstack([mo.md("**Compartments the LR analysis actually saw**"),
                               mo.ui.table(_compartment_tbl, selection=None)]),
                ],
                widths=[3, 2],
                gap=2,
            ),
            mo.callout(
                mo.md(
                    f"**Coverage gap: {_n_unmodelled:,} of {_total:,} annotated cells "
                    f"({100 * _n_unmodelled / _total:.1f}%) never enter any interaction "
                    f"result.** Verified directly against the interaction table: no "
                    f"`endothelial` or `opc` group appears as a source or target anywhere "
                    f"in its {len(interactions_df):,} rows.\n\n"
                    + "\n".join(
                        f"- **{r.cell_type_predicted}** — {r.n_cells:,} cells"
                        for r in _unmodelled.itertuples()
                    )
                    + "\n\nTwo of these are worth a decision:\n\n"
                    "- **`endothelial` (vasculature)** is a canonical TME compartment and "
                    "is annotated here, but no rule ever subsets it. Perivascular niche "
                    "signalling is therefore entirely absent from every result in this "
                    "project.\n"
                    "- **`opc`** looks like a naming miss, not a choice. "
                    "`config.yaml nerve_cells.cell_types` lists "
                    "`\"oligodendrocyte precursor cell\"`, but `scrna_annotate.py` emits "
                    "the label `opc`. `nerve_cell_subset.py` matches by substring, and "
                    "neither string contains the other — so OPCs are silently dropped "
                    "from the nerve compartment.\n"
                    "- **`t_cell`** here is the *marker-argmax* label on the full cohort; "
                    "it is separate from the T cells inside the immune subset, which are "
                    "sub-clustered out of the `microglia` blob."
                ),
                kind="warn",
            ),
        ]
    )
    return (census_df,)


@app.cell
def _patient_composition_header(mo):
    mo.md("""
    ---
    ## Panel B - Who contributed these cells?
    Batch context for every group behind the interaction results.
    """)
    return


@app.cell
def _patient_composition(immune_comp_df, mo, nerve_purity_df, np, plt, sns):
    """Immune per-patient composition matrix + nerve per-cluster purity summary."""
    # --- Immune: a real cluster x sample matrix exists -------------------------
    _mat = (
        immune_comp_df.pivot_table(
            index="immune_leiden", columns="sample_id", values="n_cells",
            aggfunc="sum", fill_value=0,
        )
    )
    _frac = _mat.div(_mat.sum(axis=1), axis=0)
    _fig1, _ax1 = plt.subplots(
        figsize=(max(6, 0.45 * _frac.shape[1] + 3), max(3, 0.3 * _frac.shape[0] + 2))
    )
    sns.heatmap(
        _frac, cmap="rocket_r", ax=_ax1,
        cbar_kws={"label": "fraction of cluster's cells"},
        xticklabels=[s[:8] for s in _frac.columns], yticklabels=_frac.index,
    )
    _ax1.set_title("Immune Leiden cluster × patient composition (row-normalised)")
    _ax1.set_xlabel("sample (first 8 chars of UUID)")
    _ax1.set_ylabel("immune_leiden")
    _fig1.tight_layout()

    # --- Nerve: only the purity summary survives -------------------------------
    _np_df = nerve_purity_df.sort_values("cluster").copy()
    _x = np.arange(len(_np_df))
    _fig2, _ax2 = plt.subplots(figsize=(max(6, 0.4 * len(_np_df) + 3), 3.6))
    _bar_colors = [
        "#C44E52" if not bool(p) else "#4C72B0"
        for p in _np_df["pass_overall"]
    ]
    _ax2.bar(_x, _np_df["dominant_sample_fraction"], color=_bar_colors)
    _ax2.axhline(0.5, ls="--", lw=1, color="black")
    _ax2.set_xticks(_x)
    _ax2.set_xticklabels([f"c{c}" for c in _np_df["cluster"]], rotation=90, fontsize=8)
    _ax2.set_ylabel("dominant sample fraction")
    _ax2.set_title(
        "Nerve cluster patient-dominance (red = fails batch QC; dashed = 0.5 threshold)"
    )
    _ax2b = _ax2.twinx()
    _ax2b.plot(_x, _np_df["n_contributing_samples"], color="#55A868", marker="o", ms=3, lw=1)
    _ax2b.set_ylabel("# contributing samples", color="#55A868")
    _fig2.tight_layout()

    mo.vstack(
        [
            mo.center(_fig1),
            mo.center(_fig2),
            mo.callout(
                mo.md(
                    "**Why the two sides look different.** The immune side has a full "
                    "`cluster × patient` matrix (`immune_cluster_composition.csv`). The "
                    "nerve side does not: `nerve_cluster_composition.csv` was never "
                    "produced for the reference cohort, and its producing rule "
                    "(`nerve_cell_subset`) cannot be re-run — "
                    "`data/processed/nerve_cells.h5ad` was deleted and is documented as "
                    "unreproducible, and re-running it would break the frozen "
                    "`cl15_split_v1_3_0` guard and overwrite the pinned v1.3.0 tables. "
                    "So the nerve panel is the per-cluster purity summary instead: "
                    "dominance fraction plus contributing-sample count."
                ),
                kind="info",
            ),
        ]
    )
    return


@app.cell
def _browser_header(mo):
    mo.md("""
    ---
    ## Panel C - Interaction browser, labelled by biology
    Same evidence as notebook 03, but every nerve cluster carries its cell-type
    identity, and you can filter by cell type directly.
    """)
    return


@app.cell
def _filters(interactions_df, mo):
    """Reactive filter widgets -- drive Panels C and D."""
    _pairs = sorted(interactions_df["compartment_pair"].astype(str).unique())
    _directions = sorted(interactions_df["direction"].astype(str).unique())
    _cell_types = sorted(t for t in interactions_df["nerve_cell_type"].unique() if t)

    source_toggle = mo.ui.switch(
        value=True, label="Use top_pairs (off = full interactions)"
    )
    pair_select = mo.ui.radio(
        options=_pairs, value=_pairs[0], label="Compartment interface"
    )
    direction_select = mo.ui.multiselect(
        options=_directions, value=_directions, label="Direction(s)"
    )
    celltype_select = mo.ui.multiselect(
        options=_cell_types, value=_cell_types, label="Nerve cell type(s)"
    )
    magnitude_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="magnitude_rank ≤"
    )
    pval_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="cellphone_pvals ≤"
    )
    ligand_search = mo.ui.text(label="Ligand contains", placeholder="e.g. NLGN1")
    receptor_search = mo.ui.text(label="Receptor contains", placeholder="e.g. NRXN")
    hide_qc_fail = mo.ui.checkbox(value=False, label="Hide QC-failing rows")
    ambiguous_only = mo.ui.checkbox(
        value=False, label="Only clusters with ambiguous identity"
    )
    return (
        ambiguous_only,
        celltype_select,
        direction_select,
        hide_qc_fail,
        ligand_search,
        magnitude_slider,
        pair_select,
        pval_slider,
        receptor_search,
        source_toggle,
    )


@app.cell
def _show_filters(
    ambiguous_only,
    celltype_select,
    direction_select,
    hide_qc_fail,
    hint,
    ligand_search,
    magnitude_slider,
    mo,
    pair_select,
    pval_slider,
    receptor_search,
    source_toggle,
):
    mo.vstack(
        [
            mo.md("### Filters"),
            mo.md(
                "*Pick a **compartment interface** first — Panels C and D are scoped "
                "to it. The cell-type filter only bites on interfaces that involve "
                "nerve.*"
            ),
            mo.hstack([source_toggle, hide_qc_fail, ambiguous_only], gap=2),
            pair_select,
            direction_select,
            celltype_select,
            mo.hstack([magnitude_slider, pval_slider], gap=2),
            mo.hstack(
                [
                    hint(
                        "magnitude_rank",
                        "LIANA aggregate magnitude rank across methods (RRA). Lower = stronger.",
                    ),
                    hint(
                        "cellphone_pvals",
                        "CellPhoneDB permutation p-value (~1000 permutations). Lower = more specific.",
                    ),
                    hint(
                        "ambiguous identity",
                        "Clusters where the majority cell_type_predicted disagrees with the "
                        "argmax canonical marker score.",
                    ),
                ],
                gap=2,
            ),
            mo.hstack([ligand_search, receptor_search], gap=2),
        ],
        gap=1,
    )
    return


@app.cell
def _filtered_view(
    ambiguous_only,
    celltype_select,
    direction_select,
    hide_qc_fail,
    interactions_df,
    ligand_search,
    magnitude_slider,
    pair_select,
    pval_slider,
    receptor_search,
    source_toggle,
    top_pairs_df,
):
    """Apply every filter reactively to the chosen source dataframe."""
    _src = top_pairs_df if source_toggle.value else interactions_df
    _df = _src[_src["compartment_pair"].astype(str) == pair_select.value].copy()

    _dirs = direction_select.value or list(_df["direction"].astype(str).unique())
    _df = _df[_df["direction"].astype(str).isin(_dirs)]

    # Cell-type filter applies only where a nerve side exists; rows with no nerve
    # compartment (immune-tumor) carry an empty cell type and must not be dropped.
    if celltype_select.value is not None:
        _keep = set(celltype_select.value)
        _df = _df[
            (_df["nerve_cell_type"] == "") | (_df["nerve_cell_type"].isin(_keep))
        ]

    if ambiguous_only.value:
        _df = _df[_df["nerve_type_agrees"] == False]  # noqa: E712

    if hide_qc_fail.value and "batch_qc_pass" in _df.columns:
        _df = _df[_df["batch_qc_pass"].astype(bool)]

    _df = _df[
        (_df["magnitude_rank"] <= magnitude_slider.value)
        & (_df["cellphone_pvals"] <= pval_slider.value)
    ]

    if ligand_search.value:
        _df = _df[
            _df["ligand_complex"].str.contains(ligand_search.value, case=False, na=False)
        ]
    if receptor_search.value:
        _df = _df[
            _df["receptor_complex"].str.contains(
                receptor_search.value, case=False, na=False
            )
        ]

    filtered_df = _df.sort_values("magnitude_rank").reset_index(drop=True)
    return (filtered_df,)


@app.cell
def _show_filtered(filtered_df, mo, pair_select):
    _cols = [
        "nerve_label",
        "nerve_cell_type",
        "immune_subtype",
        "direction",
        "ligand_complex",
        "receptor_complex",
        "lrscore",
        "magnitude_rank",
        "specificity_rank",
        "cellphone_pvals",
        "batch_qc_pass",
        "nerve_type_agrees",
    ]
    _view = filtered_df[[c for c in _cols if c in filtered_df.columns]]
    mo.vstack(
        [
            mo.md(f"### `{pair_select.value}` — filtered rows: **{len(filtered_df):,}**"),
            mo.ui.table(_view, page_size=25, selection=None),
        ]
    )
    return


@app.cell
def _celltype_matrix_header(mo):
    mo.md("""
    ---
    ## Panel D - Cell type × immune subtype interface matrix
    The rollup notebook 03 cannot produce: it only knows Leiden ids.
    """)
    return


@app.cell
def _celltype_interface_matrix(filtered_df, mo, plt, sns):
    """Significant LR-pair counts aggregated by nerve cell type x immune subtype."""
    _df = filtered_df[
        (filtered_df["nerve_cell_type"] != "") & (filtered_df["immune_subtype"] != "")
    ]
    mo.stop(
        _df.empty,
        mo.callout(
            mo.md(
                "*This matrix needs rows with **both** a nerve and an immune side. "
                "Select the `immune-nerve` interface in Panel C.*"
            ),
            kind="info",
        ),
    )

    _dirs = sorted(_df["direction"].astype(str).unique())
    _fig, _axes = plt.subplots(
        1, len(_dirs), figsize=(1 + 4.2 * len(_dirs), 3.6), squeeze=False
    )
    for _i, _d in enumerate(_dirs):
        _sub = _df[_df["direction"].astype(str) == _d]
        _counts = (
            _sub.groupby(["nerve_cell_type", "immune_subtype"], observed=True)
            .size()
            .unstack("immune_subtype", fill_value=0)
        )
        _ax = _axes[0][_i]
        sns.heatmap(
            _counts, annot=True, fmt="d", cmap="mako_r", ax=_ax,
            cbar_kws={"label": "# LR pairs"},
        )
        _ax.set_title(_d)
        _ax.set_xlabel("immune subtype")
        _ax.set_ylabel("nerve cell type" if _i == 0 else "")
    _fig.tight_layout()

    _amb = sorted(
        _df.loc[_df["nerve_type_agrees"] == False, "nerve_label"].unique()  # noqa: E712
    )
    mo.vstack(
        [
            mo.center(_fig),
            mo.callout(
                mo.md(
                    "**How to read:** each cell counts ligand-receptor pairs surviving "
                    "the Panel C thresholds, rolled up from Leiden clusters to their "
                    "dominant cell type. Counts are **pair counts, not effect sizes** — "
                    "a cell type spread over more clusters accumulates more pairs, so "
                    "compare within a column before comparing across rows."
                    + (
                        "\n\n**Ambiguous clusters folded into this rollup:** "
                        + ", ".join(f"`{a}`" for a in _amb)
                        + ". Their dominant cell type disagrees with their marker-score "
                        "argmax, so their contribution to a row is uncertain."
                        if _amb
                        else ""
                    )
                ),
                kind="info",
            ),
        ]
    )
    return


@app.cell
def _lead_targets_header(mo):
    mo.md("""
    ---
    ## Panel E - Curated lead targets, traced back to evidence
    `nerve_crosstalk_lead_targets.csv` is a hand-curated shortlist. Select an axis to
    see every interaction row supporting it.
    """)
    return


@app.cell
def _lead_target_selector(lead_targets_df, mo):
    axis_select = mo.ui.dropdown(
        options=lead_targets_df["axis"].astype(str).tolist(),
        value=lead_targets_df["axis"].astype(str).iloc[0],
        label="Lead axis",
    )
    mo.vstack(
        [
            mo.ui.table(
                lead_targets_df[
                    ["axis", "target_class", "angle", "best_mag", "best_lrscore",
                     "interfaces", "immune", "n_rows"]
                ],
                selection=None,
                page_size=10,
            ),
            axis_select,
        ]
    )
    return (axis_select,)


@app.cell
def _lead_target_evidence(axis_select, interactions_df, lead_targets_df, mo):
    """Trace a curated axis back to the interaction rows that support it."""
    _row = lead_targets_df[lead_targets_df["axis"].astype(str) == axis_select.value]
    mo.stop(_row.empty, mo.md("*No axis selected.*"))
    _r = _row.iloc[0]
    _a, _b = str(_r["mol_a"]), str(_r["mol_b"])

    # The axis is unordered: match either orientation of the ligand/receptor pair.
    _lig = interactions_df["ligand_complex"].astype(str)
    _rec = interactions_df["receptor_complex"].astype(str)
    _hits = interactions_df[((_lig == _a) & (_rec == _b)) | ((_lig == _b) & (_rec == _a))]

    _cols = [
        "compartment_pair", "direction", "nerve_label", "nerve_cell_type",
        "immune_subtype", "ligand_complex", "receptor_complex", "lrscore",
        "magnitude_rank", "batch_qc_pass",
    ]
    _shown = _hits[[c for c in _cols if c in _hits.columns]].sort_values("magnitude_rank")

    _n_fail = int((~_hits["batch_qc_pass"].astype(bool)).sum()) if len(_hits) else 0
    mo.vstack(
        [
            mo.callout(
                mo.md(
                    f"**{_r['axis']}** — {_r['target_class']} · angle: *{_r['angle']}*\n\n"
                    f"{_r['rationale']}\n\n"
                    f"**Existing drugs:** {_r['existing_drugs']}"
                ),
                kind="info",
            ),
            mo.md(
                f"### Supporting rows: **{len(_hits):,}** "
                f"({_n_fail:,} on a batch-QC-failing group) "
                f"— curated sheet recorded `n_rows = {_r['n_rows']}`"
            ),
            mo.ui.table(_shown.reset_index(drop=True), selection=None, page_size=20),
        ]
    )
    return


@app.cell
def _relay_header(mo):
    mo.md("""
    ---
    ## Panel F - Relay circuits, with cell-type labels
    Pick an immune subtype as the relay hub: leg 1 is the strongest **tumor → immune**
    signalling into it, leg 2 the strongest **immune → nerve** signalling out of it.
    Read together they are candidate **tumor ↦ immune ↦ nerve** relays.
    """)
    return


@app.cell
def _relay_selector(interactions_df, mo):
    _subtypes = sorted(s for s in interactions_df["immune_subtype"].unique() if s)
    relay_hub = mo.ui.dropdown(
        options=_subtypes,
        value=_subtypes[0] if _subtypes else None,
        label="Immune relay hub (subtype)",
    )
    relay_topn = mo.ui.slider(
        start=5, stop=40, value=15, step=5, label="Top-N LR pairs per leg"
    )
    relay_mag = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="magnitude_rank ≤ (both legs)"
    )
    mo.hstack([relay_hub, relay_topn, relay_mag], gap=2)
    return relay_hub, relay_mag, relay_topn


@app.cell
def _relay_view(interactions_df, mo, relay_hub, relay_mag, relay_topn):
    """Two-leg relay through the chosen immune subtype, labelled by nerve cell type."""
    mo.stop(
        relay_hub.value is None,
        mo.md("*No immune subtypes available in the interaction table.*"),
    )
    _hub = relay_hub.value

    def _leg(direction: str, cols: list[str]):
        _d = interactions_df[
            (interactions_df["direction"].astype(str) == direction)
            & (interactions_df["immune_subtype"] == _hub)
            & (interactions_df["magnitude_rank"] <= relay_mag.value)
        ].sort_values("magnitude_rank").head(relay_topn.value)
        return _d[[c for c in cols if c in _d.columns]].reset_index(drop=True)

    _leg1 = _leg(
        "tumor_to_immune",
        ["ligand_complex", "receptor_complex", "lrscore", "magnitude_rank",
         "batch_qc_pass"],
    )
    _leg2 = _leg(
        "immune_to_nerve",
        ["nerve_label", "nerve_cell_type", "ligand_complex", "receptor_complex",
         "lrscore", "magnitude_rank", "batch_qc_pass"],
    )

    _v1 = (
        mo.ui.table(_leg1, selection=None, page_size=15) if len(_leg1)
        else mo.md(f"*No `tumor → {_hub}` pairs at magnitude_rank ≤ {relay_mag.value:.2f}.*")
    )
    _v2 = (
        mo.ui.table(_leg2, selection=None, page_size=15) if len(_leg2)
        else mo.md(f"*No `{_hub} → nerve` pairs at magnitude_rank ≤ {relay_mag.value:.2f}.*")
    )

    _targets = sorted(_leg2["nerve_cell_type"].unique()) if len(_leg2) else []
    mo.vstack(
        [
            mo.callout(
                mo.md(
                    f"**Relay hub: `{_hub}`** — {len(_leg1)} inbound `tumor → {_hub}` "
                    f"pairs, {len(_leg2)} outbound `{_hub} → nerve` pairs "
                    f"(magnitude_rank ≤ {relay_mag.value:.2f})."
                    + (
                        " Outbound signalling reaches: "
                        + ", ".join(f"**{t}**" for t in _targets if t)
                        if _targets else ""
                    )
                ),
                kind="info",
            ),
            mo.hstack(
                [
                    mo.vstack([mo.md(f"### Leg 1 — `tumor → {_hub}`"), _v1]),
                    mo.vstack([mo.md(f"### Leg 2 — `{_hub} → nerve`"), _v2]),
                ],
                widths=[1, 1],
                gap=2,
            ),
            mo.accordion(
                {
                    "How to interpret a relay": mo.md("""
    A relay is a **hypothesis**, not a measured pathway. The two legs are not
    molecularly linked here: leg 1 shows the tumor expressing ligands whose receptors
    are on the immune hub, leg 2 shows that hub expressing ligands whose receptors are
    on nerve clusters. Nothing establishes that the inbound signal *causes* the
    outbound one. Prioritise hubs that are strong on both legs, then look for
    mechanism (signalling or transcriptional evidence) before treating it as a circuit.
    """)
                }
            ),
        ]
    )
    return


@app.cell
def _clinical_header(mo):
    mo.md("""
    ---
    ## Panel G - Clinical association
    """)
    return


@app.cell
def _clinical(clinical_df, mo):
    """Cluster-level clinical association -- reported honestly, including the null."""
    _n = len(clinical_df)
    _sig = clinical_df[clinical_df["padj"].notna() & (clinical_df["padj"] < 0.05)]
    _covs = sorted(clinical_df["covariate"].unique())

    _kind = "warn" if _sig.empty else "success"
    _msg = (
        f"**{len(_sig)} of {_n} tests survive FDR correction "
        f"(`padj < 0.05`).**\n\n"
    )
    if _sig.empty:
        _msg += (
            "There is **no evidence** of nerve-cluster abundance covarying with any "
            "available clinical variable in this cohort. Treat this panel as a "
            "negative result, not as an untested question.\n\n"
            "The test is also underpowered by the metadata: in "
            "`data/external/gdc_clinical.tsv` only `tissue_type`, `gender`, `race` and "
            "`treatment_outcome` vary across the 17 samples. `age_at_index`, "
            "`tumor_grade`, `prior_malignancy` and `project_id` are blank or constant, "
            "so `primary_diagnosis`, `tumor_grade`, `prior_malignancy` and the age "
            "correlation each collapse to a single degenerate row. See "
            "`markdowns/stratifying_open_issues.md` — the fix is a new GDC `cases` API "
            "rule, not a reanalysis."
        )
    _msg += f"\n\nCovariates tested: {', '.join(f'`{c}`' for c in _covs)}."

    mo.vstack(
        [
            mo.callout(mo.md(_msg), kind=_kind),
            mo.ui.table(
                clinical_df.sort_values("padj").reset_index(drop=True),
                selection=None,
                page_size=15,
            ),
        ]
    )
    return


@app.cell
def _export_button(mo):
    export_button = mo.ui.run_button(
        label="Export current filtered view → results/tables/"
    )
    mo.center(export_button)
    return (export_button,)


@app.cell
def _export(
    ambiguous_only,
    celltype_select,
    datetime,
    direction_select,
    export_button,
    filtered_df,
    hashlib,
    hide_qc_fail,
    ligand_search,
    magnitude_slider,
    mo,
    pair_select,
    paths,
    pval_slider,
    receptor_search,
    source_toggle,
    tables_dir,
):
    """Persist the filtered view with FAIR provenance -- hashes every input."""
    mo.stop(not export_button.value, mo.md("*Click the button above to export.*"))

    _ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    _out_path = tables_dir / f"tme_nerve_immune_exploration_{_ts}.csv"
    _prov_path = tables_dir / f"tme_nerve_immune_exploration_{_ts}.provenance.txt"

    filtered_df.to_csv(_out_path, index=False)

    _hashes = {
        _k: hashlib.sha256(_p.read_bytes()).hexdigest()[:12] for _k, _p in paths.items()
    }
    _prov = (
        f"timestamp: {_ts}\n"
        f"notebook: notebooks/04_tme_nerve_immune_explorer.py\n"
        f"cohort: reference (17-sample TCGA-GBM, v1.3.0)\n"
        f"n_rows_exported: {len(filtered_df)}\n"
        "inputs_sha256_12:\n"
        + "".join(
            f"  {_k}: {_h}  ({paths[_k].name})\n" for _k, _h in sorted(_hashes.items())
        )
        + "filters:\n"
        f"  source_toggle_top_pairs: {source_toggle.value}\n"
        f"  compartment_pair: {pair_select.value}\n"
        f"  directions: {direction_select.value}\n"
        f"  nerve_cell_types: {celltype_select.value}\n"
        f"  ambiguous_clusters_only: {ambiguous_only.value}\n"
        f"  hide_qc_failing: {hide_qc_fail.value}\n"
        f"  magnitude_rank_le: {magnitude_slider.value}\n"
        f"  cellphone_pvals_le: {pval_slider.value}\n"
        f"  ligand_contains: {ligand_search.value!r}\n"
        f"  receptor_contains: {receptor_search.value!r}\n"
    )
    _prov_path.write_text(_prov)

    mo.callout(
        mo.md(
            f"**Exported** `{_out_path.name}` ({len(filtered_df):,} rows)\n\n"
            f"**Provenance** `{_prov_path.name}` — records the sha256[:12] of all "
            f"{len(_hashes)} inputs plus every filter value. Keep it next to the CSV; "
            f"it is what makes the export reproducible."
        ),
        kind="success",
    )
    return


@app.cell
def _footer(mo):
    mo.md("""
    ---
    **Scope.** Reference cohort (17-sample TCGA-GBM, v1.3.0) only. The CELLxGENE
    Census replication cohort is excluded — see
    `markdowns/blocker_census_liana_raw_counts.md`.

    **Known limits, all surfaced in the panels above.**
    1. The vasculature (`endothelial`) is annotated but never modelled — no
       perivascular signalling appears anywhere in this project.
    2. `opc` is dropped from the nerve compartment by a config/label naming mismatch.
    3. The immune compartment was re-run after the LR scores were computed; subtype
       cell counts in Panel B do not match the sample sizes behind the LR scores.
    4. Nerve per-patient composition is only available as a purity summary.
    5. No clinical association survives FDR, and the clinical metadata is thin.

    **FAIR.** Inputs are produced by `nerve_tumor_immune_interaction`,
    `annotate_cluster_qc`, `nerve_cluster_annotations`, `immune_cluster_annotations`,
    `immune_cell_subset`, `nerve_batch_qc`, `scrna_annotate` and
    `nerve_clinical_association`. Exports carry a `.provenance.txt` sidecar hashing
    every input.
    """)
    return


if __name__ == "__main__":
    app.run()
