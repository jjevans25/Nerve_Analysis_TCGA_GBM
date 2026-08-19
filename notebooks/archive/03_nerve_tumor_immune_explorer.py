"""
Marimo reactive notebook: Interactive exploration of LIANA ligand–receptor
interactions across the nerve, tumor, and immune compartments in TCGA GBM.
Addresses: Which LR pairs drive nerve↔tumor↔immune crosstalk, across which
compartment interface, in which direction, and for which nerve cluster / immune
subtype — including candidate tumor→immune→nerve relay circuits?
Source rule: workflow/rules/immune.smk → nerve_tumor_immune_interaction.py
(+ annotate_cluster_qc for the _with_qc variants).

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

__generated_with = "0.23.5"
app = marimo.App(
    width="wide",
    app_title="GBM — Nerve↔Tumor↔Immune LR Explorer",
)


@app.cell
def _imports():
    import hashlib
    from datetime import datetime
    from pathlib import Path

    import duckdb
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import yaml

    return Path, datetime, duckdb, hashlib, mo, np, pd, plt, sns, yaml


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
def _hint_helper(mo):
    """Reusable click-to-expand inline definition — uses native <details>."""

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
    """Resolve table paths from config.yaml — no hardcoded absolute paths."""
    project_root = Path(__file__).parent.parent
    _config_path = project_root / "config" / "config.yaml"
    with open(_config_path) as _f:
        config = yaml.safe_load(_f)

    tables_dir = project_root / config["dirs"]["tables"]
    interactions_path = tables_dir / "nerve_tumor_immune_interactions_with_qc.csv"
    top_pairs_path = tables_dir / "nerve_tumor_immune_top_pairs_with_qc.csv"

    _missing = [p for p in (interactions_path, top_pairs_path) if not p.exists()]
    if _missing:
        mo.stop(
            True,
            mo.callout(
                mo.md(
                    "Required artifacts not found:\n\n"
                    + "\n".join(f"- `{p}`" for p in _missing)
                    + "\n\nRun the pipeline first: `snakemake annotate_cluster_qc` "
                    "(which depends on `nerve_tumor_immune_interaction`, which in turn "
                    "depends on `immune_cell_subset` and `immune_cluster_annotations`)."
                ),
                kind="danger",
            ),
        )
    return interactions_path, tables_dir, top_pairs_path


@app.cell
def _header(mo):
    mo.md("""
    # GBM — Nerve↔Tumor↔Immune Ligand–Receptor Explorer
    **Source rule:** `workflow/rules/immune.smk` → `nerve_tumor_immune_interaction.py`
    **Inputs:** `nerve_tumor_immune_interactions_with_qc.csv` (full LIANA table)
    and `nerve_tumor_immune_top_pairs_with_qc.csv` (top-N per pairing).

    Three cross-compartment interfaces are inferred jointly:
    **`nerve-tumor`**, **`immune-nerve`**, **`immune-tumor`**. Each row is one
    ligand–receptor pair for one directional `source → target` group pairing.
    Lower `magnitude_rank` / `specificity_rank` ⇒ stronger evidence.

    > The **Relay circuits** panel near the bottom is the three-way payoff:
    > it stitches `tumor → immune` and `immune → nerve` signals through a chosen
    > immune subtype to surface candidate **tumor ↦ immune ↦ nerve** relays.
    """)
    return


@app.cell
def _glossary(mo):
    """Column definitions — collapsed by default."""
    mo.accordion(
        {
            "Glossary — columns and how to read them": mo.md("""
    | Term | What it is | How to read it |
    |---|---|---|
    | `compartment_pair` | Which two compartments the row spans: `nerve-tumor`, `immune-nerve`, or `immune-tumor`. | Pick the interface you care about in the filter panel. |
    | `direction` | `{source}_to_{target}` at the compartment level, e.g. `tumor_to_nerve`, `immune_to_nerve`. | "Who signals to whom." The ligand is expressed by `source`. |
    | `source`, `target` | Sender / receiver group labels: `malignant`, `nerve_c{N}`, or `immune_{subtype}`. | The directional pair. |
    | `nerve_cluster` | Nerve Leiden cluster id (`nerve_c{N}`) if the row involves nerve, else blank. | Map to biology via `nerve_cluster_markers.csv`. |
    | `immune_subtype` | Immune subtype (microglia / tam / t_cell / nk_cell / dendritic) if the row involves immune, else blank. | From `immune_cluster_annotations.csv`. |
    | `ligand_complex`, `receptor_complex` | LR pair identifiers (LIANA consensus resource). | The molecular contact. |
    | `lrscore` | LIANA aggregate strength in [0,1]. | **Higher = stronger** (dotplot color). |
    | `magnitude_rank` | LIANA aggregate magnitude rank. | **Lower = stronger evidence.** `< 0.05` ≈ top 5%. |
    | `specificity_rank` | LIANA aggregate specificity rank. | **Lower = more group-specific.** |
    | `cellphone_pvals` | CellPhoneDB-style permutation p-value (~1000 permutations). | **Lower = more group-specific.** |
    | `batch_qc_pass` | Combined batch-QC verdict: AND of the nerve-cluster and immune-subtype purity verdicts that apply to the row. | `False` ⇒ at least one side is patient-driven; marked `*`. |
    | `nerve_batch_qc_pass` / `immune_batch_qc_pass` | Per-side batch-QC verdict (blank where that compartment is not involved). | Diagnose which side drives a QC failure. |

    *A `0.05` threshold is convention. Tighten to `0.01` for shortlists, relax to `0.10` for small groups (immune subtypes can be small in this cold tumor).*
    """)
        }
    )
    return


@app.cell
def _load_tables(interactions_path, pd, top_pairs_path):
    """Read both CSVs; coerce categorical/label columns; fill blank compartment slots."""
    interactions_df = pd.read_csv(interactions_path)
    top_pairs_df = pd.read_csv(top_pairs_path)

    for _df in (interactions_df, top_pairs_df):
        _df["nerve_cluster"] = _df["nerve_cluster"].fillna("").astype(str)
        _df["immune_subtype"] = _df["immune_subtype"].fillna("").astype(str)
        _df["compartment_pair"] = _df["compartment_pair"].astype("category")
        _df["direction"] = _df["direction"].astype("category")
        # Coerce QC verdict columns to real booleans. The per-side columns are
        # read as object dtype (blank where that compartment is not involved),
        # so values may arrive as the strings "True"/"False"; normalise them.
        _bool_map = {True: True, "True": True, False: False, "False": False}
        for _c in ("batch_qc_pass", "nerve_batch_qc_pass", "immune_batch_qc_pass"):
            if _c in _df.columns:
                _df[_c] = _df[_c].map(_bool_map)
        if "batch_qc_pass" in _df.columns:
            # Combined verdict applies to every row — fill any parse gaps as pass.
            _df["batch_qc_pass"] = _df["batch_qc_pass"].fillna(True).astype(bool)
        # A readable "interacting unit" label for plotting, per interface.
        _df["pair_unit"] = [
            nc if cp == "nerve-tumor"
            else (isub if cp == "immune-tumor" else f"{nc} × {isub}")
            for cp, nc, isub in zip(
                _df["compartment_pair"].astype(str),
                _df["nerve_cluster"],
                _df["immune_subtype"],
            )
        ]
    return interactions_df, top_pairs_df


@app.cell
def _qc_banner(interactions_df, mo):
    """Banner summarising batch-QC failures by compartment side."""
    _fail = interactions_df[~interactions_df["batch_qc_pass"]]
    if _fail.empty:
        _view = mo.callout(
            mo.md("All groups in this table pass batch-correction QC."),
            kind="success",
        )
    else:
        # A nerve cluster / immune subtype is "QC-failing" if its own per-side
        # verdict is False. Coerce to real booleans first (the columns are read
        # as object dtype because not-involved rows are blank/NaN).
        _nerve_fail_mask = interactions_df["nerve_batch_qc_pass"].map(
            {True: False, False: True}
        ).fillna(False).astype(bool)
        _bad_nerve = sorted(
            {c for c in interactions_df.loc[_nerve_fail_mask, "nerve_cluster"] if c}
        )
        if "immune_batch_qc_pass" in interactions_df.columns:
            _immune_fail_mask = interactions_df["immune_batch_qc_pass"].map(
                {True: False, False: True}
            ).fillna(False).astype(bool)
            _bad_immune = sorted(
                {s for s in interactions_df.loc[_immune_fail_mask, "immune_subtype"] if s}
            )
        else:
            _bad_immune = []
        _view = mo.callout(
            mo.md(
                "**Batch-QC caveat.** Rows whose nerve cluster and/or immune "
                "subtype failed batch-correction purity are marked with `*` in "
                "the plots below and can be dropped with the **Hide QC-failing "
                "rows** checkbox.\n\n"
                f"- QC-failing **nerve clusters**: {', '.join(_bad_nerve) or '—'}\n"
                f"- QC-failing **immune subtypes**: {', '.join(_bad_immune) or '—'}"
            ),
            kind="warn",
        )
    _view
    return


@app.cell
def _overview_stats(duckdb, interactions_df, mo):
    """Headline counts per compartment interface via DuckDB."""
    _conn = duckdb.connect()
    _conn.register("inter", interactions_df)
    overview_df = _conn.execute(
        """
        SELECT
            compartment_pair,
            COUNT(*)                                                       AS n_rows,
            COUNT(DISTINCT ligand_complex)                                 AS n_ligands,
            COUNT(DISTINCT receptor_complex)                              AS n_receptors,
            SUM(CASE WHEN magnitude_rank < 0.05 THEN 1 ELSE 0 END)         AS n_mag_rank_lt_05,
            ROUND(100.0 * SUM(CASE WHEN cellphone_pvals < 0.05 THEN 1 ELSE 0 END)
                  / COUNT(*), 1)                                           AS pct_pval_lt_05
        FROM inter
        GROUP BY compartment_pair
        ORDER BY compartment_pair
        """
    ).df()
    _conn.close()
    mo.md("## Overview by compartment interface")
    return (overview_df,)


@app.cell
def _show_overview(mo, overview_df):
    mo.vstack([mo.ui.table(overview_df, selection=None)])
    return


@app.cell
def _filters(interactions_df, mo):
    """Reactive filter widgets — drive every plot/table downstream."""
    _pairs = sorted(interactions_df["compartment_pair"].cat.categories.tolist())
    _directions = sorted(interactions_df["direction"].cat.categories.tolist())

    source_toggle = mo.ui.switch(
        value=True, label="Use top_pairs (off = full interactions)"
    )
    pair_select = mo.ui.radio(
        options=_pairs, value=_pairs[0], label="Compartment interface"
    )
    direction_select = mo.ui.multiselect(
        options=_directions, value=_directions, label="Direction(s)"
    )
    magnitude_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="magnitude_rank ≤"
    )
    pval_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="cellphone_pvals ≤"
    )
    ligand_search = mo.ui.text(label="Ligand contains", placeholder="e.g. SPP1")
    receptor_search = mo.ui.text(label="Receptor contains", placeholder="e.g. NRXN")
    hide_qc_fail = mo.ui.checkbox(value=False, label="Hide QC-failing rows")
    return (
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
            mo.md("## Filters"),
            mo.md(
                "*Pick a **compartment interface** first — every view below is "
                "scoped to it. Defaults of `0.05` are conventional cutoffs.*"
            ),
            mo.hstack([source_toggle, hide_qc_fail], gap=2),
            pair_select,
            direction_select,
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
    """Apply filters reactively to the chosen source dataframe."""
    _src = top_pairs_df if source_toggle.value else interactions_df
    _df = _src[_src["compartment_pair"].astype(str) == pair_select.value].copy()

    _selected_dirs = direction_select.value or list(
        _df["direction"].astype(str).unique()
    )
    _df = _df[_df["direction"].astype(str).isin(_selected_dirs)]

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
            _df["receptor_complex"].str.contains(receptor_search.value, case=False, na=False)
        ]

    filtered_df = _df.sort_values("magnitude_rank").reset_index(drop=True)
    return (filtered_df,)


@app.cell
def _show_filtered(filtered_df, mo, pair_select):
    mo.vstack(
        [
            mo.md(
                f"### `{pair_select.value}` — filtered rows: **{len(filtered_df):,}**"
            ),
            mo.ui.table(filtered_df, page_size=25, selection=None),
        ]
    )
    return


@app.cell
def _significance_heatmap(filtered_df, magnitude_slider, mo, plt, sns):
    """Counts of LR pairs per interacting unit × direction at the current threshold."""
    mo.stop(
        filtered_df.empty,
        mo.md("*No rows survive current filters — heatmap suppressed.*"),
    )

    # Units that contain any QC-failing row get a `*` tick.
    _fail_units = set(
        filtered_df.loc[~filtered_df["batch_qc_pass"].astype(bool), "pair_unit"]
    )

    def _tick(u):
        return f"{u} *" if u in _fail_units else str(u)

    _counts = (
        filtered_df.groupby(["pair_unit", "direction"], observed=True)
        .size()
        .unstack("direction", fill_value=0)
    )

    _fig, _ax = plt.subplots(
        figsize=(
            max(4, 0.9 * _counts.shape[1] + 3),
            max(3, 0.32 * _counts.shape[0] + 2),
        )
    )
    sns.heatmap(
        _counts,
        annot=True,
        fmt="d",
        cmap="magma_r",
        cbar_kws={"label": f"# LR pairs (magnitude_rank ≤ {magnitude_slider.value:.2f})"},
        ax=_ax,
        yticklabels=[_tick(u) for u in _counts.index],
    )
    _ax.set_title("Significant LR pairs per interacting unit × direction")
    _ax.set_xlabel("Direction (source → target)")
    _ax.set_ylabel("Interacting unit (`*` = batch-QC fail)")
    _fig.tight_layout()
    mo.vstack(
        [
            mo.center(_fig),
            mo.callout(
                mo.md(
                    "**How to read:** each cell counts LR pairs surviving the current "
                    "`magnitude_rank` and `cellphone_pvals` thresholds. For the "
                    "`immune-nerve` interface each row is a `nerve_cluster × immune_subtype` "
                    "unit; for `nerve-tumor` it is the nerve cluster; for `immune-tumor` "
                    "the immune subtype."
                ),
                kind="info",
            ),
        ]
    )
    return


@app.cell
def _topk_slider(mo):
    topk_slider = mo.ui.slider(
        start=5, stop=100, value=25, step=5, label="Top-K LR pairs in dotplot"
    )
    mo.hstack([topk_slider])
    return (topk_slider,)


@app.cell
def _top_lr_dotplot(filtered_df, mo, np, plt, topk_slider):
    """Dotplot of top-K LR pairs from the filtered view (x = interacting unit)."""
    mo.stop(filtered_df.empty, mo.md("*No rows to plot.*"))

    _top = filtered_df.head(topk_slider.value).copy()
    _top["lr_label"] = (
        _top["ligand_complex"].astype(str) + " → " + _top["receptor_complex"].astype(str)
    )
    _top["nlog_mag"] = -np.log10(_top["magnitude_rank"].clip(lower=1e-6))
    _fail_units = set(_top.loc[~_top["batch_qc_pass"].astype(bool), "pair_unit"])

    _x_levels = sorted(_top["pair_unit"].astype(str).unique())
    _y_levels = list(dict.fromkeys(_top["lr_label"].tolist()))
    _x_idx = {v: i for i, v in enumerate(_x_levels)}
    _y_idx = {v: i for i, v in enumerate(_y_levels)}

    _fig, _ax = plt.subplots(
        figsize=(max(5, 0.55 * len(_x_levels) + 4), max(3, 0.3 * len(_y_levels) + 2))
    )
    _sc = _ax.scatter(
        x=[_x_idx[v] for v in _top["pair_unit"].astype(str)],
        y=[_y_idx[v] for v in _top["lr_label"]],
        s=20 + 12 * _top["nlog_mag"].to_numpy(),
        c=_top["lrscore"].to_numpy(),
        cmap="viridis",
        edgecolors="black",
        linewidths=0.3,
    )
    _ax.set_xticks(range(len(_x_levels)))
    _ax.set_xticklabels(
        [f"{v} *" if v in _fail_units else v for v in _x_levels],
        rotation=45,
        ha="right",
        fontsize=8,
    )
    _ax.set_yticks(range(len(_y_levels)))
    _ax.set_yticklabels(_y_levels, fontsize=8)
    _ax.set_xlabel("Interacting unit (`*` = batch-QC fail)")
    _ax.set_ylabel("Ligand → Receptor")
    _ax.set_title(f"Top-{topk_slider.value} LR pairs (size = -log10 magnitude_rank)")
    _fig.colorbar(_sc, ax=_ax, label="lrscore", shrink=0.6)
    _fig.tight_layout()
    mo.center(_fig)
    return


@app.cell
def _relay_header(mo):
    mo.md("""
    ---
    ## Relay circuits — tumor ↦ immune ↦ nerve
    The three-way payoff. Pick an **immune subtype** as the relay hub: the left
    table lists the strongest **`tumor → immune`** LR pairs (tumor signalling
    *into* that immune subtype), the right table the strongest
    **`immune → nerve`** LR pairs (that immune subtype signalling *out to* nerve
    clusters). Read together they are candidate **tumor ↦ immune ↦ nerve** relay
    circuits — a route by which the tumor could reshape nerve cells indirectly
    through immune intermediaries.
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
    """Two-leg relay through the chosen immune subtype."""
    mo.stop(
        relay_hub.value is None,
        mo.md("*No immune subtypes available in the interaction table.*"),
    )
    _hub = relay_hub.value
    _cols = [
        "source",
        "target",
        "nerve_cluster",
        "ligand_complex",
        "receptor_complex",
        "lrscore",
        "magnitude_rank",
        "specificity_rank",
        "batch_qc_pass",
    ]

    # Leg 1: tumor → immune(hub). direction == tumor_to_immune, immune side is hub.
    _leg1 = interactions_df[
        (interactions_df["direction"].astype(str) == "tumor_to_immune")
        & (interactions_df["immune_subtype"] == _hub)
        & (interactions_df["magnitude_rank"] <= relay_mag.value)
    ].sort_values("magnitude_rank").head(relay_topn.value)

    # Leg 2: immune(hub) → nerve. direction == immune_to_nerve, immune side is hub.
    _leg2 = interactions_df[
        (interactions_df["direction"].astype(str) == "immune_to_nerve")
        & (interactions_df["immune_subtype"] == _hub)
        & (interactions_df["magnitude_rank"] <= relay_mag.value)
    ].sort_values("magnitude_rank").head(relay_topn.value)

    _n1, _n2 = len(_leg1), len(_leg2)
    _leg1_view = (
        mo.ui.table(_leg1[_cols].reset_index(drop=True), selection=None, page_size=15)
        if _n1
        else mo.md(f"*No `tumor → {_hub}` pairs at magnitude_rank ≤ {relay_mag.value:.2f}.*")
    )
    _leg2_view = (
        mo.ui.table(_leg2[_cols].reset_index(drop=True), selection=None, page_size=15)
        if _n2
        else mo.md(f"*No `{_hub} → nerve` pairs at magnitude_rank ≤ {relay_mag.value:.2f}.*")
    )

    mo.vstack(
        [
            mo.callout(
                mo.md(
                    f"**Relay hub: `{_hub}`** — "
                    f"{_n1} inbound `tumor → {_hub}` pairs, "
                    f"{_n2} outbound `{_hub} → nerve` pairs "
                    f"(magnitude_rank ≤ {relay_mag.value:.2f})."
                ),
                kind="info",
            ),
            mo.hstack(
                [
                    mo.vstack([mo.md(f"### Leg 1 — `tumor → {_hub}`"), _leg1_view]),
                    mo.vstack([mo.md(f"### Leg 2 — `{_hub} → nerve`"), _leg2_view]),
                ],
                widths=[1, 1],
                gap=2,
            ),
            mo.accordion(
                {
                    "How to interpret a relay": mo.md("""
    A relay is a **hypothesis**, not a measured pathway. Leg 1 shows the tumor
    expressing ligands whose receptors are on the immune hub; Leg 2 shows that same
    immune hub expressing ligands whose receptors are on nerve clusters. When a hub
    receives a strong tumor signal *and* emits a strong nerve signal, it is a
    candidate intermediary by which tumor state could propagate to nerve cells. To
    prioritise, look for hubs that are strong on **both** legs, and cross-check the
    `immune_to_nerve` `nerve_cluster` column against the nerve-cluster biology in
    `nerve_cluster_markers.csv`. Legs are not molecularly linked here — connecting a
    specific inbound ligand to a specific outbound ligand requires downstream
    mechanism (e.g. signalling or transcriptional evidence).
    """)
                }
            ),
        ]
    )
    return


@app.cell
def _export_button(mo):
    export_button = mo.ui.run_button(label="Export current filtered view → results/tables/")
    mo.center(export_button)
    return (export_button,)


@app.cell
def _export(
    datetime,
    direction_select,
    export_button,
    filtered_df,
    hashlib,
    interactions_path,
    ligand_search,
    magnitude_slider,
    mo,
    pair_select,
    pval_slider,
    receptor_search,
    source_toggle,
    tables_dir,
    top_pairs_path,
):
    """Persist the filtered view with FAIR provenance (input hash + filter params)."""
    mo.stop(not export_button.value, mo.md("*Click the button above to export.*"))

    _src_path = top_pairs_path if source_toggle.value else interactions_path
    _src_hash = hashlib.sha256(_src_path.read_bytes()).hexdigest()[:12]
    _ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    _out_path = tables_dir / f"nerve_tumor_immune_exploration_{_ts}.csv"
    _prov_path = tables_dir / f"nerve_tumor_immune_exploration_{_ts}.provenance.txt"

    filtered_df.to_csv(_out_path, index=False)
    _prov = (
        f"timestamp: {_ts}\n"
        f"source: {_src_path.name}\n"
        f"source_sha256_12: {_src_hash}\n"
        f"n_rows_exported: {len(filtered_df)}\n"
        f"filters:\n"
        f"  source_toggle_top_pairs: {source_toggle.value}\n"
        f"  compartment_pair: {pair_select.value}\n"
        f"  directions: {direction_select.value}\n"
        f"  magnitude_rank_le: {magnitude_slider.value}\n"
        f"  cellphone_pvals_le: {pval_slider.value}\n"
        f"  ligand_contains: {ligand_search.value!r}\n"
        f"  receptor_contains: {receptor_search.value!r}\n"
    )
    _prov_path.write_text(_prov)

    mo.callout(
        mo.md(
            f"**Exported** `{_out_path.name}` ({len(filtered_df):,} rows)\n\n"
            f"**Provenance** `{_prov_path.name}` (input sha256[:12] = `{_src_hash}`)"
        ),
        kind="success",
    )
    return


@app.cell
def _footer(mo):
    mo.md("""
    ---
    **FAIR:** input artifacts produced by `workflow/rules/immune.smk`
    (`nerve_tumor_immune_interaction`) and `annotate_cluster_qc` (the `_with_qc`
    batch-QC join). Exports include a `.provenance.txt` sidecar with the input
    hash and every filter value.
    """)
    return


if __name__ == "__main__":
    app.run()
