"""
Marimo reactive notebook: Interactive exploration of LIANA ligand–receptor
interactions between malignant cells and nerve clusters.
Addresses: Which LR pairs drive nerve↔tumor crosstalk in TCGA GBM, in which
direction, and in which nerve cluster?
"""

import marimo

__generated_with = "0.23.5"
app = marimo.App(width="wide", app_title="GBM — Nerve↔Tumor LR Explorer")


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
    interactions_path = tables_dir / "nerve_tumor_interactions.csv"
    top_pairs_path = tables_dir / "nerve_tumor_top_pairs.csv"

    _missing = [p for p in (interactions_path, top_pairs_path) if not p.exists()]
    if _missing:
        mo.stop(
            True,
            mo.callout(
                mo.md(
                    "Required artifacts not found:\n\n"
                    + "\n".join(f"- `{p}`" for p in _missing)
                    + "\n\nRun the pipeline first: `snakemake nerve_tumor_interaction`."
                ),
                kind="danger",
            ),
        )
    return interactions_path, tables_dir, top_pairs_path


@app.cell
def _header(mo):
    mo.md("""
    # GBM — Nerve↔Tumor Ligand–Receptor Explorer
    **Source rule:** `workflow/rules/nerve_cells.smk` → `nerve_tumor_interaction.py`
    **Inputs:** `nerve_tumor_interactions.csv` (full LIANA table) and
    `nerve_tumor_top_pairs.csv` (top-N per cluster × direction)
    **Convention:** lower `magnitude_rank` / `specificity_rank` ⇒ stronger
    evidence. Direction is read as `source → target`.

    > See the **Glossary** accordion below for column definitions and how each is computed.
    """)
    return


@app.cell
def _glossary(mo):
    """Top-level glossary — collapsed by default; expand for definitions."""
    mo.accordion(
        {
            "Glossary — column definitions and how they are computed": mo.md("""
| Term | What it is | How it is computed | How to read it |
|---|---|---|---|
| `cellphone_pvals` | CellPhoneDB-style permutation p-value for the LR pair in `(source, target)`. | Cell-cluster labels are randomly permuted (LIANA default ≈ 1000 permutations); the p-value is the fraction of permutations whose mean(ligand)·mean(receptor) ≥ the observed value. | **Lower ⇒ more cluster-specific.** `0.0` means the observed mean exceeded every permutation. |
| `lr_means` | Test statistic: mean ligand expression in source × mean receptor expression in target. | Computed on the normalized AnnData expression matrix; aggregated by LIANA from the CellPhoneDB method. | Higher = stronger raw signal. Sensitive to library size. |
| `expr_prod` | Raw (un-permuted) per-cluster expression product (ligand × receptor). | Source-cluster ligand mean × target-cluster receptor mean, no normalization. | Same direction as `lr_means`; sanity check. |
| `lrscore` | LIANA aggregate **strength** score in [0, 1]. | Rank-aggregation (RobustRankAggregate) across multiple LR-scoring methods, then normalized. | **Higher = stronger.** Used as the dotplot color and the direction-comparison axes. |
| `magnitude_rank` | LIANA aggregate **magnitude** rank across methods. | RRA rank fraction (0–1) over magnitude-style scores. | **Lower = stronger evidence.** `< 0.05` ≈ top 5% by magnitude. |
| `specificity_rank` | LIANA aggregate **specificity** rank. | RRA over specificity-style scores. | **Lower = more cluster-specific.** |
| `direction` | `malignant_to_nerve` or `nerve_to_malignant` (read as `source → target`). | Set by which side is `source` vs `target` in the LIANA call. | Use to ask "who is signaling to whom?" |
| `nerve_cluster` | Leiden cluster id from `nerve_cells.h5ad` (36 clusters). | Re-clustering of the non-malignant nerve-cell subspace. | Categorical id; map to biology via `nerve_cluster_markers.csv`. |
| `source`, `target` | Sender and receiver cell-type labels. | From the LIANA call (`source_labels`, `target_labels`). | `(source, target)` defines the directional pair. |

*A `0.05` threshold is convention, not a hard rule. Tighten to `0.01` for confirmatory shortlists, relax to `0.10` for small clusters.*
""")
        }
    )
    return


@app.cell
def _load_tables(interactions_path, pd, top_pairs_path):
    """Read both CSVs; coerce categorical columns for memory + plotting."""
    interactions_df = pd.read_csv(interactions_path)
    top_pairs_df = pd.read_csv(top_pairs_path)

    for _df in (interactions_df, top_pairs_df):
        _df["nerve_cluster"] = _df["nerve_cluster"].astype("category")
        _df["direction"] = _df["direction"].astype("category")
    return interactions_df, top_pairs_df


@app.cell
def _overview_stats(duckdb, interactions_df, mo, top_pairs_df):
    """Headline counts — uses DuckDB on the in-memory dataframes."""
    _conn = duckdb.connect()
    _conn.register("inter", interactions_df)
    _conn.register("topp", top_pairs_df)
    overview_df = _conn.execute(
        """
        SELECT
            'full_interactions' AS table,
            COUNT(*)                                          AS n_rows,
            COUNT(DISTINCT nerve_cluster)                     AS n_clusters,
            COUNT(DISTINCT ligand_complex)                    AS n_ligands,
            COUNT(DISTINCT receptor_complex)                  AS n_receptors,
            ROUND(100.0 * SUM(CASE WHEN cellphone_pvals < 0.05 THEN 1 ELSE 0 END)
                  / COUNT(*), 1)                              AS pct_pval_lt_05,
            SUM(CASE WHEN magnitude_rank < 0.05 THEN 1 ELSE 0 END)
                                                              AS n_mag_rank_lt_05
        FROM inter
        UNION ALL
        SELECT
            'top_pairs',
            COUNT(*),
            COUNT(DISTINCT nerve_cluster),
            COUNT(DISTINCT ligand_complex),
            COUNT(DISTINCT receptor_complex),
            ROUND(100.0 * SUM(CASE WHEN cellphone_pvals < 0.05 THEN 1 ELSE 0 END)
                  / COUNT(*), 1),
            SUM(CASE WHEN magnitude_rank < 0.05 THEN 1 ELSE 0 END)
        FROM topp
        """
    ).df()
    _conn.close()

    mo.md("## Overview")
    return (overview_df,)


@app.cell
def _show_overview(mo, overview_df):
    mo.vstack(
        [
            mo.accordion(
                {
                    "What do these summary columns mean?": mo.md("""
- **`pct_pval_lt_05`** — % of LR rows with `cellphone_pvals < 0.05`
  (CellPhoneDB-style permutation p-value; lower = more cluster-specific).
- **`n_mag_rank_lt_05`** — count of LR rows with `magnitude_rank < 0.05`
  (top 5% of LR pairs by aggregate magnitude across LIANA methods).
""")
                }
            ),
            mo.ui.table(overview_df, selection=None),
        ]
    )
    return


@app.cell
def _filters(interactions_df, mo):
    """Reactive filter widgets — drive every plot/table downstream."""
    _clusters = sorted(interactions_df["nerve_cluster"].cat.categories.tolist())

    source_toggle = mo.ui.switch(
        value=True, label="Use top_pairs (off = full interactions)"
    )
    cluster_select = mo.ui.multiselect(
        options=_clusters, value=_clusters, label="Nerve cluster(s)"
    )
    direction_radio = mo.ui.radio(
        options=["both", "malignant_to_nerve", "nerve_to_malignant"],
        value="both",
        label="Direction",
    )
    magnitude_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="magnitude_rank ≤"
    )
    pval_slider = mo.ui.slider(
        start=0.0, stop=1.0, value=0.05, step=0.01, label="cellphone_pvals ≤"
    )
    ligand_search = mo.ui.text(label="Ligand contains", placeholder="e.g. NLGN1")
    receptor_search = mo.ui.text(label="Receptor contains", placeholder="e.g. NRXN")
    return (
        cluster_select,
        direction_radio,
        ligand_search,
        magnitude_slider,
        pval_slider,
        receptor_search,
        source_toggle,
    )


@app.cell
def _show_filters(
    cluster_select,
    direction_radio,
    hint,
    ligand_search,
    magnitude_slider,
    mo,
    pval_slider,
    receptor_search,
    source_toggle,
):
    mo.vstack(
        [
            mo.md("## Filters"),
            mo.md(
                "*Defaults of `0.05` are conventional cutoffs — tighten to `0.01` "
                "for confirmatory shortlists, or relax to `0.10` for small clusters.*"
            ),
            mo.hstack([source_toggle, direction_radio], gap=2),
            mo.md(
                "*`top_pairs` is the pre-filtered top-N per `(cluster × direction)` "
                "shortlist; toggle off to see the unfiltered LIANA table.*"
            ),
            cluster_select,
            mo.hstack([magnitude_slider, pval_slider], gap=2),
            mo.hstack(
                [
                    hint(
                        "magnitude_rank",
                        "LIANA aggregate magnitude rank across methods (RRA). Lower = stronger evidence.",
                    ),
                    hint(
                        "cellphone_pvals",
                        "CellPhoneDB permutation p-value (~1000 cell-label permutations). Lower = more cluster-specific.",
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
    cluster_select,
    direction_radio,
    interactions_df,
    ligand_search,
    magnitude_slider,
    pval_slider,
    receptor_search,
    source_toggle,
    top_pairs_df,
):
    """Apply filters reactively to the chosen source dataframe."""
    _src = top_pairs_df if source_toggle.value else interactions_df
    _df = _src.copy()

    _selected_clusters = cluster_select.value or list(
        _df["nerve_cluster"].cat.categories
    )
    _df = _df[_df["nerve_cluster"].isin(_selected_clusters)]

    if direction_radio.value != "both":
        _df = _df[_df["direction"] == direction_radio.value]

    _df = _df[
        (_df["magnitude_rank"] <= magnitude_slider.value)
        & (_df["cellphone_pvals"] <= pval_slider.value)
    ]

    if ligand_search.value:
        _df = _df[
            _df["ligand_complex"].str.contains(
                ligand_search.value, case=False, na=False
            )
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
def _show_filtered(filtered_df, mo):
    mo.vstack(
        [
            mo.md(f"### Filtered rows: **{len(filtered_df):,}**"),
            mo.ui.table(filtered_df, page_size=25, selection=None),
        ]
    )
    return


@app.cell
def _significance_heatmap(filtered_df, magnitude_slider, mo, plt, pval_slider, sns):
    """Counts of LR pairs per (nerve_cluster, direction) at current threshold."""
    mo.stop(
        filtered_df.empty,
        mo.md("*No rows survive current filters — heatmap suppressed.*"),
    )

    _counts = (
        filtered_df.groupby(["nerve_cluster", "direction"], observed=True)
        .size()
        .unstack("direction", fill_value=0)
    )

    _fig, _ax = plt.subplots(
        figsize=(max(4, 0.6 * _counts.shape[1] + 3), max(3, 0.3 * _counts.shape[0] + 2))
    )
    sns.heatmap(
        _counts,
        annot=True,
        fmt="d",
        cmap="magma_r",
        cbar_kws={"label": f"# LR pairs (magnitude_rank ≤ {magnitude_slider.value:.2f})"},
        ax=_ax,
    )
    _ax.set_title("Significant LR pairs per nerve cluster × direction")
    _ax.set_xlabel("Direction")
    _ax.set_ylabel("Nerve cluster")
    _fig.tight_layout()
    mo.vstack(
        [
            mo.center(_fig),
            mo.callout(
                mo.md(
                    f"**How to read:** each cell counts LR pairs that survive *both* "
                    f"`magnitude_rank ≤ {magnitude_slider.value:.2f}` *and* "
                    f"`cellphone_pvals ≤ {pval_slider.value:.2f}`. "
                    "Darker = more LR pairs supported in that cluster × direction."
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
    return (topk_slider,)


@app.cell
def _show_topk(mo, topk_slider):
    mo.hstack([topk_slider])
    return


@app.cell
def _top_lr_dotplot(filtered_df, mo, np, plt, topk_slider):
    """Dotplot of top-K LR pairs from the filtered view."""
    mo.stop(filtered_df.empty, mo.md("*No rows to plot.*"))

    _top = filtered_df.head(topk_slider.value).copy()
    _top["lr_label"] = (
        _top["ligand_complex"].astype(str) + " → " + _top["receptor_complex"].astype(str)
    )
    _top["cluster_str"] = _top["nerve_cluster"].astype(str)
    _top["nlog_mag"] = -np.log10(_top["magnitude_rank"].clip(lower=1e-6))

    _x_levels = sorted(_top["cluster_str"].unique())
    _y_levels = list(dict.fromkeys(_top["lr_label"].tolist()))
    _x_idx = {v: i for i, v in enumerate(_x_levels)}
    _y_idx = {v: i for i, v in enumerate(_y_levels)}

    _fig, _ax = plt.subplots(
        figsize=(max(5, 0.5 * len(_x_levels) + 4), max(3, 0.3 * len(_y_levels) + 2))
    )
    _sc = _ax.scatter(
        x=[_x_idx[v] for v in _top["cluster_str"]],
        y=[_y_idx[v] for v in _top["lr_label"]],
        s=20 + 12 * _top["nlog_mag"].to_numpy(),
        c=_top["lrscore"].to_numpy(),
        cmap="viridis",
        edgecolors="black",
        linewidths=0.3,
    )
    _ax.set_xticks(range(len(_x_levels)))
    _ax.set_xticklabels(_x_levels, rotation=45, ha="right")
    _ax.set_yticks(range(len(_y_levels)))
    _ax.set_yticklabels(_y_levels, fontsize=8)
    _ax.set_xlabel("Nerve cluster")
    _ax.set_ylabel("Ligand → Receptor")
    _ax.set_title(f"Top-{topk_slider.value} LR pairs (size = -log10 magnitude_rank)")
    _fig.colorbar(_sc, ax=_ax, label="lrscore", shrink=0.6)
    _fig.tight_layout()
    mo.vstack(
        [
            mo.center(_fig),
            mo.callout(
                mo.md("""
**How to read this dotplot**

- **x-axis** — nerve cluster id (Leiden).
- **y-axis** — `ligand_complex → receptor_complex`.
- **Dot size** — `−log10(magnitude_rank)`. Larger ⇒ stronger aggregate magnitude.
- **Dot color** — `lrscore` (LIANA aggregate strength, 0–1). Brighter ⇒ stronger.
- **Worked example:** a large, bright dot at cluster `12` means this LR pair is
  among the top-ranked by magnitude *and* has a high aggregate strength in
  cluster `12`.
"""),
                kind="info",
            ),
        ]
    )
    return


@app.cell
def _direction_compare_selector(interactions_df, mo):
    _clusters = sorted(interactions_df["nerve_cluster"].cat.categories.tolist())
    compare_cluster = mo.ui.dropdown(
        options=_clusters, value=_clusters[0], label="Cluster for direction comparison"
    )
    annotate_n = mo.ui.slider(
        start=0, stop=20, value=8, step=1, label="Annotate top-N by combined lrscore"
    )
    return annotate_n, compare_cluster


@app.cell
def _show_compare_selector(annotate_n, compare_cluster, mo):
    mo.hstack([compare_cluster, annotate_n], gap=2)
    return


@app.cell
def _direction_comparison(
    annotate_n,
    compare_cluster,
    interactions_df,
    mo,
    plt,
):
    """Scatter malignant→nerve vs nerve→malignant lrscore for matching LR pairs."""
    _df = interactions_df[
        interactions_df["nerve_cluster"] == compare_cluster.value
    ].copy()
    mo.stop(_df.empty, mo.md("*No rows for selected cluster.*"))

    _m2n = _df[_df["direction"] == "malignant_to_nerve"][
        ["ligand_complex", "receptor_complex", "lrscore"]
    ].rename(columns={"lrscore": "lrscore_m2n"})
    _n2m = _df[_df["direction"] == "nerve_to_malignant"][
        ["ligand_complex", "receptor_complex", "lrscore"]
    ].rename(
        # swap so a malignant→nerve LR pair lines up with its reciprocal
        # (in nerve→malignant the ligand of one direction is the receptor of the other)
        columns={
            "lrscore": "lrscore_n2m",
            "ligand_complex": "receptor_complex",
            "receptor_complex": "ligand_complex",
        }
    )
    _joined = _m2n.merge(
        _n2m, on=["ligand_complex", "receptor_complex"], how="outer"
    ).fillna(0.0)

    _joined["combined"] = _joined["lrscore_m2n"] + _joined["lrscore_n2m"]
    _joined = _joined.sort_values("combined", ascending=False).reset_index(drop=True)

    _fig, _ax = plt.subplots(figsize=(7, 6))
    _ax.scatter(
        _joined["lrscore_m2n"],
        _joined["lrscore_n2m"],
        s=18,
        alpha=0.6,
        edgecolors="black",
        linewidths=0.2,
    )
    _lim = max(_joined[["lrscore_m2n", "lrscore_n2m"]].to_numpy().max(), 0.01) * 1.05
    _ax.plot([0, _lim], [0, _lim], "--", color="grey", linewidth=0.8)
    _ax.set_xlim(0, _lim)
    _ax.set_ylim(0, _lim)
    _ax.set_xlabel("lrscore (malignant → nerve)")
    _ax.set_ylabel("lrscore (nerve → malignant)")
    _ax.set_title(f"Bidirectional LR signaling — {compare_cluster.value}")

    for _, _row in _joined.head(annotate_n.value).iterrows():
        _ax.annotate(
            f"{_row['ligand_complex']}→{_row['receptor_complex']}",
            xy=(_row["lrscore_m2n"], _row["lrscore_n2m"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )
    _fig.tight_layout()
    mo.vstack(
        [
            mo.accordion(
                {
                    "Why ligand and receptor are swapped, and why many points land on the diagonal": mo.md("""
**Why ligand and receptor are swapped for one direction.** A given LR pair `A → B`
is stored once per direction in the LIANA output. For `malignant → nerve` it appears
as `(ligand=A, receptor=B)`; for `nerve → malignant` the LIANA resource often stores
the same molecular contact under the swapped roles `(ligand=B, receptor=A)`, because
pairs like neurexin/neuroligin are bidirectional adhesion molecules. To compare the
*same* contact across directions, the `nerve_to_malignant` rows have their
`ligand_complex` and `receptor_complex` columns swapped before joining on
`(ligand_complex, receptor_complex)`.

**Why many points land exactly on the diagonal.** This is **expected**, not a bug.
The dominant LIANA component score is `mean(ligand in source) × mean(receptor in target)`.
Multiplication is commutative, so for a symmetric resource pair the swap yields
`mean(A in malignant) × mean(B in nerve)` for both directions — the two scores are
*identical by construction*. Diagonal points therefore mean "this contact is
bidirectional and the chart cannot distinguish a direction for it." Look **off the
diagonal** for direction-specific signal.
""")
                }
            ),
            mo.center(mo.mpl.interactive(_fig)),
            mo.callout(
                mo.md("""
**How to read this chart — three regions, three meanings.**

- **On the diagonal (`y = x`).** Symmetric pairs (e.g., NLGN/NRXN). The LIANA score
  is mathematically identical in both directions because the underlying mean-based
  statistic is commutative. The pair is real and often strong, but the chart cannot
  resolve a direction for it.
- **On the x-axis (`y = 0`).** The reciprocal `nerve → malignant` row was missing
  from the LIANA output (filtered below threshold or absent in the resource for that
  direction); the zero comes from `fillna(0)`, not a measured null.
- **On the y-axis (`x = 0`).** Same as above, in the opposite direction.
- **Off the diagonal, away from the axes.** The interesting signal — these pairs
  have a real direction-specific difference. Points **above** the diagonal favor
  `nerve → malignant`; points **below** favor `malignant → nerve`.

When interpreting any cluster, focus on the off-diagonal mass first; treat the
diagonal pile-up as confirmatory evidence of bidirectional contacts rather than as
direction-specific findings.
"""),
                kind="info",
            ),
        ]
    )
    return


@app.cell
def _export_button(mo):
    export_button = mo.ui.run_button(label="Export current view → results/tables/")
    mo.center(export_button)
    return (export_button,)


@app.cell
def _export(
    cluster_select,
    datetime,
    direction_radio,
    export_button,
    filtered_df,
    hashlib,
    interactions_path,
    ligand_search,
    magnitude_slider,
    mo,
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
    _out_path = tables_dir / f"nerve_tumor_exploration_{_ts}.csv"
    _prov_path = tables_dir / f"nerve_tumor_exploration_{_ts}.provenance.txt"

    filtered_df.to_csv(_out_path, index=False)

    _prov = (
        f"timestamp: {_ts}\n"
        f"source: {_src_path.name}\n"
        f"source_sha256_12: {_src_hash}\n"
        f"n_rows_exported: {len(filtered_df)}\n"
        f"filters:\n"
        f"  source_toggle_top_pairs: {source_toggle.value}\n"
        f"  direction: {direction_radio.value}\n"
        f"  nerve_clusters: {cluster_select.value}\n"
        f"  magnitude_rank_le: {magnitude_slider.value}\n"
        f"  cellphone_pvals_le: {pval_slider.value}\n"
        f"  ligand_contains: {ligand_search.value!r}\n"
        f"  receptor_contains: {receptor_search.value!r}\n"
    )
    _prov_path.write_text(_prov)

    mo.callout(
        mo.md(
            f"**Exported** `{_out_path.name}` ({len(filtered_df):,} rows)\n\n"
            f"**Provenance** `{_prov_path.name}` (input sha256[:12] = `{_src_hash}`)\n\n"
            f"*The `.provenance.txt` records the source CSV hash and every filter value. "
            f"Keep it next to the export — it is what makes the file reproducible.*"
        ),
        kind="success",
    )
    return


@app.cell
def _footer(mo):
    mo.md("""
    ---
    ### Column legend (with how each value is computed)
    | column | meaning | how it is computed |
    |---|---|---|
    | `source`, `target` | sender → receiver cell types | LIANA `source_labels` / `target_labels` |
    | `ligand_complex`, `receptor_complex` | LR pair identifiers | LIANA resource (CellPhoneDB-derived) |
    | `lr_means` | mean expression product across source/target | mean(ligand in source) × mean(receptor in target) on normalized AnnData |
    | `cellphone_pvals` | CellPhoneDB-style permutation p-value | fraction of ~1000 cell-label permutations with statistic ≥ observed |
    | `expr_prod` | raw expression product (ligand × receptor) | same as `lr_means` but un-normalized — sanity check |
    | `lrscore` | LIANA aggregate strength score (0–1, higher = stronger) | RRA rank-aggregation across LR methods, normalized |
    | `magnitude_rank` | LIANA aggregate magnitude rank (lower = better) | RRA rank fraction (0–1) over magnitude-style scores |
    | `specificity_rank` | LIANA aggregate specificity rank (lower = better) | RRA rank fraction over specificity-style scores |
    | `direction` | `malignant_to_nerve` or `nerve_to_malignant` | set by which side is `source` vs `target` |
    | `nerve_cluster` | recipient/sender nerve cluster id | Leiden cluster id from `nerve_cells.h5ad` (36 clusters) |

    **FAIR:** input artifacts produced by `workflow/rules/nerve_cells.smk`
    (rule `nerve_tumor_interaction`). Exports include a `.provenance.txt`
    sidecar with input hash and filter parameters.
    """)
    return


if __name__ == "__main__":
    app.run()
