"""
Marimo reactive notebook: Nerve-cell GSEA enrichment explorer.
Addresses: Which biological themes characterize each nerve-cell cluster, and
which clusters share enrichment profiles?

Data source: results/tables/nerve_enrichment_with_qc.csv (the
`annotate_cluster_qc` Snakemake rule's annotated copy of
`nerve_enrichment.csv`, with `batch_qc_pass` and purity columns joined in).
"""

import marimo

__generated_with = "0.23.1"
app = marimo.App(width="wide", app_title="Nerve Enrichment Explorer")


@app.cell
def _imports():
    import re
    import uuid
    from datetime import datetime, timezone
    from pathlib import Path

    import duckdb
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import yaml
    from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
    from scipy.spatial.distance import pdist, squareform

    return (
        Path,
        dendrogram,
        datetime,
        duckdb,
        leaves_list,
        linkage,
        mo,
        np,
        pd,
        pdist,
        plt,
        re,
        sns,
        squareform,
        timezone,
        uuid,
        yaml,
    )


@app.cell
def _load_config(Path, yaml):
    """Load project config — all paths sourced here, no hardcoded strings."""
    _cfg_path = Path(__file__).parent.parent / "config" / "config.yaml"
    with open(_cfg_path) as _fh:
        config = yaml.safe_load(_fh)
    enrichment_path = Path(config["dirs"]["tables"]) / "nerve_enrichment_with_qc.csv"
    markers_path = Path(config["dirs"]["tables"]) / "nerve_cluster_markers_with_qc.csv"
    msigdb_release = config.get("msigdb", {}).get("release", "unknown")
    return config, enrichment_path, markers_path, msigdb_release


@app.cell
def _header(mo, msigdb_release):
    mo.md(
        f"""
# Nerve-Cell GSEA Enrichment Explorer

**Source:** `results/tables/nerve_enrichment.csv` — GO BP + MF gene sets,
MSigDB release `{msigdb_release}`, computed via `gseapy.prerank` against
cluster Wilcoxon z-scores from the `nerve_cell_heterogeneity` rule.

Four views: (1) term-by-cluster heatmap, (2) interactive top-K per cluster,
(3) cluster similarity from enrichment profiles, (4) thematic term-category
roll-up. Filtering uses adjusted-FDR thresholds (the `Adjusted P-value`
column equals FDR q-val from prerank).
"""
    )
    return


@app.cell
def _check_data(enrichment_path, mo):
    _has_data = (
        enrichment_path.exists()
        and sum(1 for line in enrichment_path.read_text().splitlines() if line.strip())
        > 1
    )
    mo.stop(
        not _has_data,
        mo.callout(
            mo.md(
                f"`{enrichment_path}` is empty or missing. Run the "
                "`nerve_cell_heterogeneity` rule first."
            ),
            kind="warn",
        ),
    )
    return


@app.cell
def _load_enrichment(duckdb, enrichment_path, np, pd):
    """Load CSV via DuckDB; parse Overlap into integer columns; add -log10 padj.

    The annotated `_with_qc.csv` adds `batch_qc_pass` and three purity columns
    so the notebook can flag clusters that failed the batch-correction QC.
    """
    _conn = duckdb.connect()
    _raw = _conn.execute(
        f"""
        SELECT cluster,
               gene_set_library,
               Term,
               CAST("Adjusted P-value" AS DOUBLE) AS padj,
               "Overlap" AS overlap_str,
               batch_qc_pass,
               dominant_sample_fraction,
               n_contributing_samples,
               normalised_entropy
        FROM read_csv('{enrichment_path}', header=true)
        """
    ).df()
    _conn.close()

    _parts = _raw["overlap_str"].str.split("/", n=1, expand=True)
    _raw["lead_size"] = pd.to_numeric(_parts[0], errors="coerce").fillna(0).astype(int)
    _raw["set_size"] = pd.to_numeric(_parts[1], errors="coerce").fillna(0).astype(int)
    _raw["cluster"] = pd.to_numeric(_raw["cluster"], errors="coerce").astype("Int64")
    _raw["neg_log10_padj"] = -np.log10(_raw["padj"].clip(lower=1e-300))

    enrichment_df = _raw.sort_values(["cluster", "padj"]).reset_index(drop=True)
    cluster_ids = sorted(enrichment_df["cluster"].dropna().unique().tolist())
    libraries = sorted(enrichment_df["gene_set_library"].unique().tolist())

    _qc_lookup = (
        enrichment_df.dropna(subset=["cluster"])
        .drop_duplicates(subset=["cluster"])
        .set_index("cluster")["batch_qc_pass"]
        .astype(bool)
    )
    failing_clusters = sorted(int(c) for c in _qc_lookup.index[~_qc_lookup])
    return cluster_ids, enrichment_df, failing_clusters, libraries


@app.cell
def _qc_banner(enrichment_df, failing_clusters, mo):
    """Banner: list clusters that failed batch-correction QC with purity context."""
    if not failing_clusters:
        _qc_banner_view = mo.callout(
            mo.md(
                "All clusters pass batch-correction QC "
                "(`batch_qc_pass == True` for every cluster in this table)."
            ),
            kind="success",
        )
    else:
        _qc_table = (
            enrichment_df[enrichment_df["cluster"].isin(failing_clusters)]
            .drop_duplicates(subset=["cluster"])
            .sort_values("cluster")[
                [
                    "cluster",
                    "dominant_sample_fraction",
                    "n_contributing_samples",
                    "normalised_entropy",
                ]
            ]
            .reset_index(drop=True)
        )
        _ids = ", ".join(str(c) for c in failing_clusters)
        _qc_banner_view = mo.vstack(
            [
                mo.callout(
                    mo.md(
                        f"""
**Batch-QC caveat** — clusters {_ids} fail the batch-correction QC in
`nerve_cluster_sample_purity.csv` and are marked with `*` in the heatmaps
and `[batch-QC fail]` in the cluster picker below. Per-cluster signals from
these clusters may reflect a single patient rather than a cohort program;
treat them as exploratory only.
"""
                    ),
                    kind="warn",
                ),
                mo.ui.table(_qc_table, selection=None),
            ]
        )
    _qc_banner_view
    return


@app.cell
def _qc_label_helpers(failing_clusters):
    """Tick-label / dropdown-label helpers that mark failing clusters."""
    _failing = set(int(c) for c in failing_clusters)

    def cluster_tick(value):
        try:
            _i = int(value)
        except (TypeError, ValueError):
            return str(value)
        return f"{_i} *" if _i in _failing else str(_i)

    def cluster_dropdown_label(value):
        try:
            _i = int(value)
        except (TypeError, ValueError):
            return str(value)
        return f"{_i} [batch-QC fail]" if _i in _failing else str(_i)

    return cluster_dropdown_label, cluster_tick


@app.cell
def _global_summary(enrichment_df, mo, pd):
    _summary = (
        enrichment_df.groupby("gene_set_library")
        .agg(
            n_rows=("Term", "size"),
            n_unique_terms=("Term", "nunique"),
            n_clusters=("cluster", "nunique"),
            n_significant=("padj", lambda s: int((s < 0.05).sum())),
            median_padj=("padj", "median"),
            min_padj=("padj", "min"),
        )
        .reset_index()
    )
    _total = pd.DataFrame(
        [
            {
                "gene_set_library": "ALL",
                "n_rows": int(len(enrichment_df)),
                "n_unique_terms": int(enrichment_df["Term"].nunique()),
                "n_clusters": int(enrichment_df["cluster"].nunique()),
                "n_significant": int((enrichment_df["padj"] < 0.05).sum()),
                "median_padj": float(enrichment_df["padj"].median()),
                "min_padj": float(enrichment_df["padj"].min()),
            }
        ]
    )
    summary_df = pd.concat([_summary, _total], ignore_index=True)
    mo.vstack(
        [
            mo.md("### Global summary"),
            mo.ui.table(summary_df, selection=None),
        ]
    )
    return (summary_df,)


@app.cell
def _section_heatmap(mo):
    mo.md(
        """
---
## 1. Term-by-cluster heatmap

Pivot the top-N enriched terms per cluster (within the chosen library) into a
`Term × cluster` matrix of `-log10(FDR)`. Optional row clustering reorders
terms so co-enriched groups appear together.
"""
    )
    return


@app.cell
def _heatmap_controls(libraries, mo):
    library_pick = mo.ui.dropdown(
        options=libraries,
        value=libraries[0],
        label="Gene-set library",
    )
    top_n = mo.ui.slider(start=1, stop=10, step=1, value=5, label="Top-N terms / cluster")
    padj_cut = mo.ui.slider(
        start=0.001, stop=0.5, step=0.001, value=0.05, label="FDR threshold"
    )
    cluster_terms_switch = mo.ui.switch(value=True, label="Cluster terms by similarity")
    mo.hstack([library_pick, top_n, padj_cut, cluster_terms_switch], gap=1)
    return cluster_terms_switch, library_pick, padj_cut, top_n


@app.cell
def _heatmap_view(
    cluster_terms_switch,
    cluster_tick,
    enrichment_df,
    leaves_list,
    library_pick,
    linkage,
    padj_cut,
    plt,
    sns,
    top_n,
):
    _sub = enrichment_df[
        (enrichment_df["gene_set_library"] == library_pick.value)
        & (enrichment_df["padj"] < padj_cut.value)
    ]
    _keep_terms = (
        _sub.sort_values(["cluster", "padj"])
        .groupby("cluster")
        .head(top_n.value)["Term"]
        .unique()
        .tolist()
    )
    _matrix = (
        enrichment_df[
            (enrichment_df["gene_set_library"] == library_pick.value)
            & (enrichment_df["Term"].isin(_keep_terms))
        ]
        .pivot_table(
            index="Term",
            columns="cluster",
            values="neg_log10_padj",
            aggfunc="max",
            fill_value=0.0,
        )
        .astype(float)
    )

    if cluster_terms_switch.value and _matrix.shape[0] >= 2:
        _link = linkage(_matrix.values, method="average", metric="correlation")
        _matrix = _matrix.iloc[leaves_list(_link), :]

    _h = max(4.0, 0.22 * _matrix.shape[0])
    _w = max(8.0, 0.35 * _matrix.shape[1])
    heatmap_fig, _ax = plt.subplots(figsize=(_w, _h))
    if _matrix.empty:
        _ax.text(
            0.5, 0.5, "No terms pass the chosen FDR threshold", ha="center", va="center"
        )
        _ax.set_axis_off()
    else:
        sns.heatmap(
            _matrix,
            ax=_ax,
            cmap="viridis",
            cbar_kws={"label": "-log10 FDR"},
            linewidths=0.0,
        )
        _ax.set_xticklabels(
            [cluster_tick(c) for c in _matrix.columns],
            rotation=_ax.get_xticklabels()[0].get_rotation() if _ax.get_xticklabels() else 0,
        )
        _ax.set_xlabel("Cluster (leiden) — `*` marks batch-QC-failing clusters")
        _ax.set_ylabel(library_pick.value)
        _ax.set_title(
            f"Top-{top_n.value} terms per cluster (FDR < {padj_cut.value:g}) — "
            f"{_matrix.shape[0]} unique terms × {_matrix.shape[1]} clusters"
        )
    heatmap_fig.tight_layout()
    heatmap_fig
    return (heatmap_fig,)


@app.cell
def _section_per_cluster(mo):
    mo.md(
        """
---
## 2. Top-K terms per cluster (interactive)

Pick a cluster; see its strongest BP and MF enrichments side by side as
horizontal bar charts of `-log10(FDR)`. The companion table shows the
leading-edge / set-size overlap.
"""
    )
    return


@app.cell
def _per_cluster_controls(cluster_dropdown_label, cluster_ids, mo):
    _options = {cluster_dropdown_label(c): str(c) for c in cluster_ids}
    cluster_pick = mo.ui.dropdown(
        options=_options,
        value=cluster_dropdown_label(cluster_ids[0]),
        label="Cluster",
    )
    top_k = mo.ui.slider(start=1, stop=10, step=1, value=10, label="Top K per library")
    include_mf = mo.ui.checkbox(value=True, label="Include GO MF panel")
    mo.hstack([cluster_pick, top_k, include_mf], gap=1)
    return cluster_pick, include_mf, top_k


@app.cell
def _per_cluster_view(
    cluster_pick,
    enrichment_df,
    include_mf,
    libraries,
    mo,
    np,
    pd,
    plt,
    top_k,
):
    _chosen = int(cluster_pick.value)
    _libs = libraries if include_mf.value else libraries[:1]
    _panel_data = {
        _lib: (
            enrichment_df[
                (enrichment_df["cluster"] == _chosen)
                & (enrichment_df["gene_set_library"] == _lib)
            ]
            .sort_values("padj")
            .head(top_k.value)
        )
        for _lib in _libs
    }

    _n_panels = len(_libs)
    per_cluster_fig, _axes = plt.subplots(
        _n_panels, 1, figsize=(11, max(3.0, 0.45 * top_k.value) * _n_panels)
    )
    if _n_panels == 1:
        _axes = [_axes]
    _palette = {
        "GO_Biological_Process_2024": "#4C72B0",
        "GO_Molecular_Function_2024": "#DD8452",
    }
    _fdr_line = float(-np.log10(0.05))
    for _ax, _lib in zip(_axes, _libs):
        _df = _panel_data[_lib]
        if _df.empty:
            _ax.text(
                0.5,
                0.5,
                f"No {_lib} hits for cluster {_chosen}",
                ha="center",
                va="center",
            )
            _ax.set_axis_off()
            continue
        _prefix = _lib.split("_")[0] + "_"
        _terms = [t.replace(_prefix, "", 1) for t in _df["Term"].tolist()]
        _positions = list(range(len(_df)))[::-1]
        _ax.barh(
            _positions,
            _df["neg_log10_padj"].values,
            color=_palette.get(_lib, "#888"),
            edgecolor="white",
        )
        _ax.set_yticks(_positions)
        _ax.set_yticklabels(_terms, fontsize=8)
        _ax.set_xlabel("-log10 FDR")
        _ax.set_title(f"{_lib} — cluster {_chosen}", fontsize=10)
        _ax.axvline(_fdr_line, linestyle="--", color="grey", linewidth=0.7)
    per_cluster_fig.tight_layout()

    per_cluster_table = (
        pd.concat(list(_panel_data.values()), ignore_index=True)
        if _panel_data
        else pd.DataFrame(columns=["gene_set_library", "Term", "padj", "overlap_str"])
    )
    mo.vstack(
        [
            per_cluster_fig,
            mo.ui.table(
                per_cluster_table[["gene_set_library", "Term", "padj", "overlap_str"]],
                selection=None,
            ),
        ]
    )
    return per_cluster_fig, per_cluster_table


@app.cell
def _section_similarity(mo):
    mo.md(
        """
---
## 3. Cluster similarity from enrichment profiles

Build a `cluster × term` matrix from significant enrichments (FDR < 0.05 by
default), then compute pairwise similarity. **Cosine** weights by `-log10(FDR)`
(graded); **Jaccard** binarises (presence/absence). The dendrogram on the
right shows how clusters group.
"""
    )
    return


@app.cell
def _similarity_controls(mo):
    metric_pick = mo.ui.dropdown(
        options=["cosine", "jaccard"], value="cosine", label="Distance metric"
    )
    min_score = mo.ui.slider(
        start=0.0, stop=5.0, step=0.1, value=1.3, label="Min -log10 FDR for inclusion"
    )
    mo.hstack([metric_pick, min_score], gap=1)
    return metric_pick, min_score


@app.cell
def _similarity_view(
    cluster_tick,
    dendrogram,
    enrichment_df,
    linkage,
    metric_pick,
    min_score,
    pd,
    pdist,
    plt,
    sns,
    squareform,
):
    _sig = enrichment_df[enrichment_df["neg_log10_padj"] >= min_score.value]
    _matrix = _sig.pivot_table(
        index="cluster",
        columns="Term",
        values="neg_log10_padj",
        aggfunc="max",
        fill_value=0.0,
    )

    similarity_fig = plt.figure(figsize=(13, 6))
    if _matrix.shape[0] < 2 or _matrix.shape[1] == 0:
        _ax = similarity_fig.add_subplot(111)
        _ax.text(
            0.5,
            0.5,
            "Not enough significant enrichments at the chosen threshold "
            "to compute cluster similarity.",
            ha="center",
            va="center",
            wrap=True,
        )
        _ax.set_axis_off()
    else:
        if metric_pick.value == "jaccard":
            _data = (_matrix.values > 0).astype(float)
            _dist = pdist(_data, metric="jaccard")
        else:
            _dist = pdist(_matrix.values, metric="cosine")
        _link = linkage(_dist, method="average")
        _order = dendrogram(_link, no_plot=True)["leaves"]
        _ordered = [_matrix.index[i] for i in _order]
        _sim_df = pd.DataFrame(
            1.0 - squareform(_dist),
            index=_matrix.index,
            columns=_matrix.index,
        ).loc[_ordered, _ordered]

        _gs = similarity_fig.add_gridspec(1, 2, width_ratios=[3, 2], wspace=0.3)
        _ax_heat = similarity_fig.add_subplot(_gs[0, 0])
        sns.heatmap(
            _sim_df,
            ax=_ax_heat,
            cmap="rocket_r",
            vmin=0,
            vmax=1,
            cbar_kws={"label": f"1 − {metric_pick.value} distance"},
            square=True,
            linewidths=0.0,
            xticklabels=[cluster_tick(c) for c in _sim_df.columns],
            yticklabels=[cluster_tick(c) for c in _sim_df.index],
        )
        _ax_heat.set_title(
            f"Cluster–cluster similarity ({metric_pick.value}, "
            f"{_matrix.shape[1]} terms ≥ -log10 FDR {min_score.value:g})"
        )
        _ax_heat.set_xlabel("cluster (`*` = batch-QC fail)")
        _ax_heat.set_ylabel("cluster (`*` = batch-QC fail)")

        _ax_dend = similarity_fig.add_subplot(_gs[0, 1])
        dendrogram(
            _link,
            labels=[cluster_tick(c) for c in _matrix.index],
            ax=_ax_dend,
            color_threshold=0,
            above_threshold_color="#444",
        )
        _ax_dend.set_title("Hierarchical linkage (average)")
        _ax_dend.set_xlabel("cluster")
        _ax_dend.set_ylabel("distance")
    similarity_fig.tight_layout()
    similarity_fig
    return (similarity_fig,)


@app.cell
def _section_themes(mo):
    mo.md(
        """
---
## 4. Term-category (theme) roll-up

Each MSigDB term is matched against a curated regex set per theme
(axon/dendrite, synapse, myelin/glia, immune, metabolic, cell-cycle,
signaling, transcription, vasculature). Per-cluster theme score = mean
`-log10(FDR)` across significant matching terms. Cell annotation shows the
count of contributing terms.
"""
    )
    return


@app.cell
def _theme_definitions():
    THEMES = {
        "axon_dendrite": [r"AXON", r"DENDRIT", r"NEURITE", r"GROWTH_CONE"],
        "synapse": [r"SYNAP", r"NEUROTRANSMIT", r"POSTSYNAP", r"PRESYNAP"],
        "myelin_glia": [r"MYELIN", r"GLIAL", r"OLIGODENDR", r"ASTROCYT", r"SCHWANN"],
        "ependymal": [r"CILIUM", r"CILIARY", r"CILIA", r"CILIOGENESIS", r"AXONEMAL", r"DYNEIN", r"EPENDYM"],
        "neuron_dev": [r"NEURON_DIFFERENTIATION", r"NEUROGENESIS", r"NEURONAL_DEV"],
        "immune": [
            r"IMMUNE",
            r"INTERFERON",
            r"CYTOKINE",
            r"COMPLEMENT",
            r"INFLAMMA",
            r"INTERLEUKIN",
            r"ANTIGEN",
            r"T_CELL",
            r"B_CELL",
        ],
        "metabolic": [
            r"METABOL",
            r"GLYCOLY",
            r"OXIDATIVE_PHOSPHO",
            r"LIPID",
            r"CHOLESTEROL",
            r"FATTY_ACID",
            r"AMINO_ACID",
        ],
        "cell_cycle": [r"MITOTIC", r"CELL_CYCLE", r"DNA_REPLIC", r"CHROMOSOM"],
        "signaling": [r"SIGNALING", r"KINASE_ACTIVITY", r"GPCR", r"RECEPTOR_BINDING"],
        "transcription": [r"TRANSCRIPTION", r"DNA_BINDING", r"CHROMATIN", r"HISTONE"],
        "vasculature": [r"ANGIOGEN", r"VASCULAR", r"ENDOTHEL", r"BLOOD_VESSEL"],
    }
    return (THEMES,)


@app.cell
def _theme_score_view(THEMES, enrichment_df, mo, plt, re, sns):
    _sig = enrichment_df[enrichment_df["padj"] < 0.05].copy()

    def _assign_themes(term):
        _hits = [
            theme
            for theme, patterns in THEMES.items()
            if any(re.search(p, term) for p in patterns)
        ]
        return _hits or ["other"]

    _sig["themes"] = _sig["Term"].map(_assign_themes)
    _exploded = _sig.explode("themes").rename(columns={"themes": "theme"})

    _score_pivot = _exploded.pivot_table(
        index="cluster",
        columns="theme",
        values="neg_log10_padj",
        aggfunc="mean",
        fill_value=0.0,
    )
    _count_pivot = _exploded.pivot_table(
        index="cluster",
        columns="theme",
        values="Term",
        aggfunc="count",
        fill_value=0,
    )
    _theme_order = list(THEMES.keys()) + (
        ["other"] if "other" in _score_pivot.columns else []
    )
    _theme_order = [t for t in _theme_order if t in _score_pivot.columns]
    _score_pivot = _score_pivot[_theme_order]
    _count_pivot = _count_pivot[_theme_order]

    theme_fig, _ax = plt.subplots(
        figsize=(
            max(8, 0.7 * _score_pivot.shape[1]),
            max(6, 0.3 * _score_pivot.shape[0]),
        )
    )
    sns.heatmap(
        _score_pivot,
        ax=_ax,
        cmap="mako_r",
        annot=_count_pivot.values,
        fmt="d",
        annot_kws={"size": 8},
        cbar_kws={"label": "Mean -log10 FDR (significant matching terms)"},
        linewidths=0.3,
        linecolor="white",
    )
    _ax.set_xlabel("theme")
    _ax.set_ylabel("cluster (leiden)")
    _ax.set_title("Per-cluster theme score (annotation = number of contributing terms)")
    theme_fig.tight_layout()

    dominant_theme_df = _score_pivot.idxmax(axis=1).rename("top_theme").reset_index()
    dominant_theme_df["top_theme_score"] = _score_pivot.max(axis=1).round(3).values
    dominant_theme_df["n_significant_terms"] = (
        _sig.groupby("cluster")["Term"]
        .count()
        .reindex(dominant_theme_df["cluster"])
        .fillna(0)
        .astype(int)
        .values
    )

    mo.vstack(
        [
            theme_fig,
            mo.md("**Dominant theme per cluster** (highest mean -log10 FDR)"),
            mo.ui.table(dominant_theme_df, selection=None),
        ]
    )
    return dominant_theme_df, theme_fig


@app.cell
def _provenance(datetime, enrichment_path, mo, msigdb_release, timezone, uuid):
    run_id = str(uuid.uuid4())
    timestamp_utc = datetime.now(timezone.utc).isoformat()
    mo.callout(
        mo.md(
            f"""
**FAIR provenance**

- `run_id`: `{run_id}`
- `timestamp_utc`: `{timestamp_utc}`
- `source`: `{enrichment_path}`
- `msigdb_release`: `{msigdb_release}`
- All paths sourced from `config/config.yaml`; no hardcoded paths.
"""
        ),
        kind="success",
    )
    return run_id, timestamp_utc


if __name__ == "__main__":
    app.run()
