"""
Marimo reactive notebook: TME x nerve x immune ligand-receptor exploration for the
CELLxGENE Census GBM replication cohort (gbm_cellxgene_56c4912d).
Addresses: Which nerve cell types interface with which immune subtypes in an
independent 170-donor cohort, and which signalling axes reproduce across BOTH
Census arms? Cross-arm agreement is the replication evidence here — the pinned
v1.3.0 reference is no longer used as a comparator (see the Panel F removal note).
Source rules: workflow/rules/datasets.smk -> ds_nerve_tumor_immune_interaction,
ds_annotate_cluster_qc, ds_nerve_cluster_annotations; and
workflow/rules/leads.smk -> nerve_immune_lead_axes.

Part of the Census notebook set: 01_census_cohort_qc (cohort composition and QC),
02_census_nerve_enrichment (what the nerve clusters are enriched for),
06_census_compartment_audit (the Test Oracle — where this notebook's compartment
claims are demonstrated), and this one (the three-way LR analysis).

Its reference-cohort predecessor, 04_tme_nerve_immune_explorer, was archived on
2026-08-19 to notebooks/archive/ together with the rest of the reference set. Their
nerve compartment was built by the logic the 2026-08-05/06 fix removed, and it is
not correctable: the reference is pinned and unreproducible, and its v1.1.0
frozen-barcode insulator bypasses the fixed mask entirely. They are a record of what
was believed, not a comparator. The one capability lost with them is clinical
association — this cohort's gdc_clinical.tsv is a generated stub.

REQUIRES the 2026-07-26 normalization fix. Before it, this cohort's LIANA tables
were computed on raw UMI counts and every row had an empty specificity_rank.
See markdowns/blocker_census_liana_raw_counts.md.

GOVERNING DEPENDENCY: the 2026-08-05/06 compartment-integrity fix. Every table
read here was rebuilt from ds_scrna_annotate down and both Census arms now pass
10/10 compartment gates (nerve 94.8-95.4% neural, was ~11%; tumor 93% malignant,
was ~48%; immune purity 99.6%). Nothing computed before 2026-08-05 is
trustworthy, and the compartment definitions this notebook reports are read from
config at runtime rather than hardcoded, because they changed three times during
that work. See markdowns/post_compartment_fix_next_steps.md.
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
    import sys
    from datetime import datetime
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import yaml

    # Nerve LIANA group ids are NOT uniformly `nerve_c{N}`: neurons are carried as
    # one pooled group, `nerve_neuron`. Slicing the prefix off by hand leaves that
    # id intact, it then matches no annotation row, and the group vanishes from
    # every rollup in silence. fair_utils owns the parsing for exactly this reason
    # (see its docstring, and the same guard in annotate_cluster_qc.py) — import it
    # rather than reimplementing. fair_utils pulls in only stdlib + numpy + pandas.
    sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))
    from fair_utils import nerve_group_key, nerve_group_sort_key

    return (
        Path, datetime, hashlib, json, mo, nerve_group_key, nerve_group_sort_key,
        np, os, pd, plt, sns, yaml,
    )


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
        "lr_provenance": ds_prov / "nerve_tumor_immune_interaction_provenance.json",
        # Root-level, NOT ds_tables: this shortlist carries both Census arms in one
        # table (rank_full/full_best_mag and rank_capped/capped_best_mag), so it
        # belongs to neither arm's namespace. Panel E reads it. Produced by
        # `rule nerve_immune_lead_axes` (workflow/rules/leads.smk).
        "lead_axes_postfix": tables_dir / "nerve_immune_lead_axes_postfix.csv",
    }

    _missing = {k: p for k, p in paths.items() if not p.exists()}
    if _missing:
        _msg = (
            f"Required artifacts not found for cohort `{dataset}`:\n\n"
            + "\n".join(f"- `{p}`" for p in _missing.values())
        )
        # Targets go FIRST: --allowed-rules/--forcerun/--quiet all take nargs='+'
        # and swallow anything placed after them. And always via run_snakemake.sh —
        # a bare `snakemake` reverts to the venv's packages and the workflow/envs
        # pins go unenforced.
        _msg += (
            "\n\nBuild them with:\n\n```\nscripts/run_snakemake.sh \\\n"
            f"  results/tables/{dataset}/nerve_cluster_annotations.csv \\\n"
            f"  results/tables/{dataset}/annotation_summary.csv \\\n"
            "  results/tables/nerve_immune_lead_axes_postfix.csv \\\n"
            "  --use-conda --cores all --rerun-triggers mtime\n```"
        )
        mo.stop(True, mo.callout(mo.md(_msg), kind="danger"))
    return config, dataset, ds_tables, paths


@app.cell
def _header(dataset, mo):
    mo.md(f"""
    # GBM Census — Nerve × Immune Replication Explorer

    **Cohort: `{dataset}`** — the CELLxGENE Census GBM 10x replication cohort
    (170 donors).

    **Read alongside:** `01_census_cohort_qc` (what the cohort is made of, and how
    unevenly), `02_census_nerve_enrichment` (what the nerve clusters are enriched
    for), and **`06_census_compartment_audit`** — which is where the compartment
    claims quoted below are actually demonstrated rather than asserted.

    Its reference-cohort predecessor, `04_tme_nerve_immune_explorer`, was **archived on
    2026-08-19** along with the rest of the reference notebook set
    (`notebooks/archive/`). Its nerve compartment was built by the logic the
    compartment fix removed and cannot be corrected, so it is a historical record, not
    a comparator. This cohort still has **no usable clinical metadata** (its
    `gdc_clinical.tsv` is a generated stub), which is the one capability that was lost
    with it.

    > **Every table below was rebuilt by the 2026-08-05/06 compartment-integrity
    > fix.** Before it, the "nerve" compartment of this cohort was 59% malignant
    > and 27% myeloid, and the `immune` compartment was 96% myeloid because
    > `t_cell` was never in `immune_cells.source_labels`. Both arms now pass 10/10
    > compartment gates — **rendered in `06_census_compartment_audit`, which also
    > shows that 5 of 23 nerve clusters still sit below the 0.80 neural bar despite
    > the aggregate passing at 95.4%.** Two consequences to carry while reading:
    > **(i)** any earlier statement of the form *"immune cells signal to X"* meant
    > *"myeloid cells signal to X"* and has to be re-tested; **(ii)** nothing
    > computed before 2026-08-05 — including the curated 40-axis shortlist — is
    > usable. See `markdowns/post_compartment_fix_next_steps.md`.

    > **This notebook also depends on the earlier 2026-07-26 normalization fix.**
    > Before it, the LIANA tables here were computed on raw UMI counts —
    > `X.max() = 53027` where log1p data peaks near 9 — and every one of the 71,189
    > rows had an empty `specificity_rank`, with 13,492 infinite log-fold-changes.
    > Those tables were void. The fix normalizes to counts-per-10k + log1p before
    > the LIANA call and adds a hard guard that refuses an unnormalized matrix.
    """)
    return


@app.cell
def _load_tables(json, nerve_group_key, paths, pd):
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
    lead_axes_df = pd.read_csv(paths["lead_axes_postfix"])
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

    # The pooled neuron group has no annotation row and never will: neurons are
    # carried as ONE LIANA group (`nerve_neuron`, purity key `neuron`) rather than
    # per Leiden cluster, because ~4.3k neurons across 170 donors is a group, not a
    # cluster set — see config nerve_cells.neuron_labels. nerve_cluster_annotations
    # only ever iterates Leiden clusters, so synthesize the row here rather than
    # letting an unlabelled group drop silently out of every rollup below.
    if "neuron" not in set(nerve_ann_df["cluster_key"]):
        nerve_ann_df = pd.concat([nerve_ann_df, pd.DataFrame([{
            "cluster": "neuron",
            "cluster_key": "neuron",
            "label": "neuron | pooled excitatory + inhibitory",
            "nerve_cell_type": "neuron",
            "nerve_score_type": "neuron",
            "nerve_type_agrees": True,
            "top_markers": "",
            "interpretation": (
                "Pooled subtype-confirmed neurons (excitatory + inhibitory), carried "
                "as a single LIANA group outside the Leiden clustering. Donor-dominated "
                "— see Panel B before reading any neuron-side axis."
            ),
        }])], ignore_index=True)

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
        # Via the shared helper, not str.replace("nerve_c", ""): that leaves the
        # pooled `nerve_neuron` id intact, it matches no annotation row, and 2,679
        # neuron rows lose their cell type without raising. Same helper the pipeline
        # uses in annotate_cluster_qc.py.
        out["cluster_key"] = out["nerve_cluster"].map(nerve_group_key).str.strip()
        out = out.merge(nerve_ann_df[_ann_cols], on="cluster_key", how="left")

        # Notebook-side twin of the [FAIR-ALERT] guard in annotate_cluster_qc.py:
        # a nerve group with no label means the group set moved again. Fail loudly
        # rather than render an empty cell type.
        _orphans = sorted(set(
            out.loc[(out["nerve_cluster"] != "") & out["label"].isna(), "cluster_key"]
        ))
        if _orphans:
            raise RuntimeError(
                f"[FAIR-ALERT] nerve groups with no annotation row: {_orphans}. "
                f"The nerve group set has changed; update nerve_cluster_annotations "
                f"or the pooled-group synthesis above before trusting any panel."
            )

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
        immune_ann_df,
        immune_cluster_purity_df,
        immune_purity_df,
        interactions_df,
        lead_axes_df,
        lr_prov,
        nerve_ann_df,
        nerve_purity_df,
        top_pairs_df,
    )


@app.cell
def _integrity_banner(dataset, interactions_df, lr_prov, mo):
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
                f"signature of the raw-counts defect is absent.\n\n"
                f"*This checks the normalization only.* Compartment integrity — "
                f"whether the cells in each group are what the group is named after "
                f"— is a separate question, evidenced by "
                f"`results/tables/{dataset}/compartment_audit_gates.csv` "
                f"(10/10 PASS on both Census arms, audit in enforcing mode)."
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
                f"```\nscripts/run_snakemake.sh \\\n"
                f"  results/tables/{dataset}/nerve_tumor_immune_interactions_with_qc.csv \\\n"
                f"  --use-conda --cores all --rerun-triggers mtime \\\n"
                f"  --forcerun ds_nerve_tumor_immune_interaction\n```"
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
def _compartment_census(annotation_df, config, interactions_df, lr_prov, mo, pd, plt):
    """Cohort annotation census vs the compartments the LR analysis actually saw."""
    _p = lr_prov["parameters"]
    _groups = set(interactions_df["source"].astype(str)) | set(
        interactions_df["target"].astype(str)
    )
    _immune_groups = {g.removeprefix("immune_") for g in _groups if g.startswith("immune_")}

    # READ FROM CONFIG, never hardcode. Both compartment definitions moved on
    # 2026-08-05/06 and a literal here silently misreports the cohort: the nerve
    # panel lost astrocyte / opc / generic neuron / ependymal, and the immune panel
    # went from the single string "microglia" to eight labels. A hardcoded copy of
    # the old definitions overstated nerve 5x and understated immune 3x.
    _nerve_labels = set(config["nerve_cells"]["cell_types"])
    _immune_cfg = config["immune_cells"]
    _immune_labels = set(
        _immune_cfg.get("source_labels") or [_immune_cfg["source_label"]]
    )

    def _compartment(ct: str) -> str:
        if ct in _nerve_labels:
            return "nerve"
        if ct in _immune_labels:
            return "immune"
        return "neither"

    _c = annotation_df.copy()
    _c["compartment"] = _c["cell_type_predicted"].map(_compartment)
    _c["feeds_a_compartment"] = _c["compartment"] != "neither"
    _c["pct_of_cohort"] = (100 * _c["n_cells"] / _c["n_cells"].sum()).round(2)
    census_df = _c[["cell_type_predicted", "n_cells", "pct_of_cohort",
                    "compartment", "feeds_a_compartment"]].sort_values(
        "n_cells", ascending=False).reset_index(drop=True)

    _total = int(census_df["n_cells"].sum())
    _unmod = census_df.loc[~census_df["feeds_a_compartment"]]
    _n_unmod = int(_unmod["n_cells"].sum())
    _n_nerve_lbl = int(census_df.loc[census_df["compartment"] == "nerve", "n_cells"].sum())
    _n_imm_lbl = int(census_df.loc[census_df["compartment"] == "immune", "n_cells"].sum())

    _fig, _ax = plt.subplots(figsize=(8, 3.4))
    _colors = {"nerve": "#4C72B0", "immune": "#DD8452", "neither": "#BBBBBB"}
    _ax.barh(census_df["cell_type_predicted"][::-1], census_df["n_cells"][::-1],
             color=[_colors[c] for c in census_df["compartment"][::-1]])
    _ax.set_xlabel("cells (marker-argmax annotation)")
    _ax.set_title("Annotation census — grey = label feeds neither nerve nor immune")
    _fig.tight_layout()

    # Three compartments, three columns — no "— total modelled" summary row. A total
    # has no group count, so that row left `groups_in_LR_table` blank, and a blank
    # cell in a table of counts reads as "zero groups" or "not measured" rather than
    # "not applicable". Same reasoning as Panel E's refusal to render an all-null
    # column. Keeping it out also holds `groups_in_LR_table` to a true integer dtype
    # instead of object, so it sorts numerically. The total moves to the label below.
    _compartment_tbl = pd.DataFrame({
        "compartment": ["tumor (CNV-malignant)", "nerve", "immune"],
        "cells_in_LR_run": [
            int(_p["n_tumor_cells"]), int(_p["n_nerve_cells"]),
            int(_p["n_immune_cells"]),
        ],
        "groups_in_LR_table": [
            1, int(_p["n_nerve_clusters"]), len(_immune_groups),
        ],
    })
    _n_modelled = int(_compartment_tbl["cells_in_LR_run"].sum())

    # Rendered as markdown, NOT mo.ui.table. A mo.ui.table stretches to the full page
    # width while its columns size to their content, so three short columns leave a
    # wide empty band on the right that reads as a fourth, blank column; marimo's CSS
    # shrink-wraps a markdown table to its content instead. With three fixed rows the
    # search / explore / export toolbar earns nothing either. Built by hand rather
    # than DataFrame.to_markdown, which needs `tabulate` — installed here but declared
    # in no environment file, so importing it would add an unpinned dependency.
    _tbl_md = "\n".join(
        ["| compartment | cells_in_LR_run | groups_in_LR_table |", "|:---|---:|---:|"]
        + [f"| {_r.compartment} | {_r.cells_in_LR_run:,} | {_r.groups_in_LR_table:,} |"
           for _r in _compartment_tbl.itertuples(index=False)]
    )

    mo.vstack([
        mo.center(_fig),
        # Full width and stacked — NOT mo.hstack(widths=[3, 2]). marimo puts
        # `min-width: 0` on the column wrapper it generates but not on the vstack
        # nested inside it, so a vstack holding a mo.ui.table cannot shrink below
        # the table's min-content width (toolbar + type badges + summary charts).
        # With the row at `flex-wrap: nowrap` both columns overflow their tracks and
        # render on top of each other. `page_size` is pinned to the row count because
        # the census carries 18-19 labels and the default 10 hides a third of them
        # behind a pager that the bar chart above already renders in full.
        mo.md("**Annotation census**"),
        mo.ui.table(census_df, selection=None, page_size=len(census_df)),
        mo.md(f"**Compartments the LR analysis saw** — {_n_modelled:,} cells "
              f"modelled in total\n\n{_tbl_md}"),
        mo.callout(
            mo.md(
                "**Compartment membership is read from `config.yaml` at run time**, not "
                "hardcoded here — nerve = "
                + ", ".join(f"`{lbl}`" for lbl in sorted(_nerve_labels))
                + f"; immune = {len(_immune_labels)} labels "
                + ", ".join(f"`{lbl}`" for lbl in sorted(_immune_labels))
                + ".\n\n"
                f"By label: **{_n_nerve_lbl:,}** cells feed nerve, **{_n_imm_lbl:,}** feed "
                f"immune, **{_n_unmod:,} of {_total:,} ({100 * _n_unmod / _total:.1f}%)** feed "
                f"neither.\n\n"
                "**These label counts are upper bounds, not compartment sizes.** The "
                "compartments additionally drop CNV-malignant cells, which is why nerve "
                f"lands at {int(_p['n_nerve_cells']):,} against {_n_nerve_lbl:,} labelled. "
                "And the **tumor compartment is CNV-derived, so it cuts across every label "
                "in this chart** — a grey bar is not evidence that those cells are absent "
                "from the analysis, only that their *label* feeds neither of the two "
                "label-defined compartments.\n\n"
                "Two things to carry before reading any result below:\n\n"
                "- **Four neural labels are excluded from the nerve compartment by "
                "decision, not by biology.** Measured against the CELLxGENE author "
                "annotation on the full arm, the non-malignant cells carrying each label "
                "are `inhibitory_neuron` 99.5% / `oligodendrocyte` 95.8% / "
                "`excitatory_neuron` 83.9% truly neural, against generic `neuron` 40.6%, "
                "`ependymal` 31.6%, `opc` 18.9%, `astrocyte` 10.0%. The last four are "
                "masked out. They remain annotated in every artifact — **their absence "
                "from the interaction tables is a masking decision and must be reported "
                "as such.** OPCs in particular (relevant to the neuron–glioma interface) "
                "are the largest scientific cost; see §2.4 and §3.2 of "
                "`markdowns/post_compartment_fix_next_steps.md`.\n"
                "- **The immune compartment is no longer myeloid-only.** It was a single "
                "`microglia` label until 2026-08-05, making it 96.4% myeloid / 3.1% "
                "lymphoid; it now takes eight labels and roughly doubled in size. Any "
                'earlier statement of the form *"immune cells signal to X"* was a '
                "**myeloid** finding and has to be re-tested, not merely rescaled."
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
def _patient_purity(
    immune_cluster_purity_df, mo, nerve_group_sort_key, nerve_purity_df, np, plt,
):
    """Per-cluster patient-dominance for both compartments."""
    def _panel(ax, df, title, *, sort_key=None):
        # The nerve purity table mixes numeric Leiden ids with the named `neuron`
        # group, so a plain sort is lexical (0, 1, 10, ..., 9, neuron).
        d = (df.sort_values("cluster", key=lambda s: s.map(sort_key))
             if sort_key else df.sort_values("cluster")).copy()
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
           "Nerve groups — patient dominance (red = fails batch QC)",
           sort_key=nerve_group_sort_key)
    _panel(_axes[1], immune_cluster_purity_df,
           "Immune Leiden clusters — patient dominance")
    _fig.tight_layout()

    _n_fail = int((~nerve_purity_df["pass_overall"].astype(bool)).sum())
    _neu = nerve_purity_df.loc[nerve_purity_df["cluster"].astype(str) == "neuron"]
    _neu_txt = ""
    if not _neu.empty:
        _r = _neu.iloc[0]
        _neu_txt = (
            f"\n\n**The pooled `neuron` group is donor-dominated.** "
            f"{int(_r['n_cells']):,} cells, dominant-sample fraction "
            f"**{float(_r['dominant_sample_fraction']):.4f}** across only "
            f"**{int(_r['n_contributing_samples'])}** donors — on the full arm it is "
            f"0.5032 across 7, i.e. one donor supplies roughly half of all neurons. "
            f"**Any neuron-side ligand–receptor axis is substantially one patient's "
            f"biology.** Report it with that stated, or restrict to axes that survive a "
            f"leave-that-donor-out check. The cohort holds ~3.4k neurons at ~20 per "
            f"donor, a ceiling no pipeline correction can lift — §2.2 and §3.6 of "
            f"`markdowns/post_compartment_fix_next_steps.md`."
        )
    mo.vstack([
        mo.center(_fig),
        mo.callout(
            mo.md(
                f"**{_n_fail} of {len(nerve_purity_df)} nerve groups fail batch QC** "
                f"(dominant-sample fraction ≥ 0.5 or too few contributing donors). "
                f"Their rows are marked `*` and can be dropped in Panel C."
                + _neu_txt
                + "\n\n*Why this is a bar chart and not the `cluster × patient` heatmap "
                "used in the archived notebook 04: this cohort has 170 donors, so that matrix is "
                "unreadable. The purity summary carries the same verdict legibly.*"
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

    # Build the matrices first so the figure can be sized from them: each panel's
    # height tracks its own row count, so y tick labels never collide, and the
    # panels stack vertically so no colorbar sits against a neighbour's labels.
    _mats = [
        (
            _d,
            _df[_df["direction"].astype(str) == _d]
            .groupby(["nerve_cell_type", "immune_subtype"], observed=True)
            .size().unstack("immune_subtype", fill_value=0),
        )
        for _d in sorted(_df["direction"].astype(str).unique())
    ]
    _heights = [max(2.0, 1.2 + 0.42 * _m.shape[0]) for _, _m in _mats]
    _width = max(6.5, 1.8 + 0.62 * max(_m.shape[1] for _, _m in _mats))
    _fig, _axes = plt.subplots(
        len(_mats), 1,
        figsize=(_width, sum(_heights)),
        squeeze=False,
        gridspec_kw={"height_ratios": _heights},
    )
    for _i, (_d, _counts) in enumerate(_mats):
        _ax = _axes[_i][0]
        sns.heatmap(_counts, annot=True, fmt="d", cmap="mako_r", ax=_ax,
                    cbar_kws={"label": "# LR pairs"})
        _ax.set_title(_d)
        _ax.set_ylabel("nerve cell type")
        # Only the bottom panel gets the axis label; tick labels stay on every
        # panel because the immune subtypes present can differ per direction.
        _ax.set_xlabel("immune subtype" if _i == len(_mats) - 1 else "")
        _ax.tick_params(axis="y", rotation=0)
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
def _lead_axes_header(mo):
    mo.md("""
    ---
    ## Panel E — Post-fix lead axes, verified against live
    The §2.1 re-derivation. Each axis carries the numbers from
    `nerve_immune_lead_axes_postfix.csv` — built by `rule nerve_immune_lead_axes` —
    **and** the same quantities recomputed here independently from this cohort's LR
    table. They are expected to agree on every row; a disagreement is a regression.
    """)
    return


@app.cell
def _lead_axes_compute(dataset, interactions_df, lead_axes_df, mo, np, pd):
    """Recompute each curated axis against the live table for the ACTIVE arm."""
    # Explicit map, never a `dataset.endswith("_full")` test — a suffix heuristic
    # would silently mislabel any future arm as "capped".
    # (active_rank, active_mag, active_n_rows, other_rank, other_arm_label)
    _ARM_COLUMNS = {
        "gbm_cellxgene_56c4912d_full": (
            "rank_full", "full_best_mag", "full_n_rows", "rank_capped", "capped"),
        "gbm_cellxgene_56c4912d": (
            "rank_capped", "capped_best_mag", "capped_n_rows", "rank_full", "full"),
    }
    mo.stop(
        dataset not in _ARM_COLUMNS,
        mo.callout(
            mo.md(
                f"**This shortlist does not cover `{dataset}`.** "
                f"`nerve_immune_lead_axes_postfix.csv` carries ranks for "
                f"{', '.join(f'`{k}`' for k in _ARM_COLUMNS)} only. Regenerate it to "
                f"include this cohort before reading a lead axis here."
            ),
            kind="danger",
        ),
    )
    # No leading underscore on names that leave the cell: marimo treats `_name` as
    # cell-local and will not expose it to downstream cells.
    _rank_col, _mag_col, _n_rows_col, _other_rank_col, other_arm = _ARM_COLUMNS[dataset]

    # The key is (axis, nerve_side), NOT axis: 144 distinct axes over 183 rows, 39 of
    # them ranked separately on both nerve sides (APP|CD74 is rank 1 on each), and the
    # ranks restart per side. Keying on axis alone silently collapses those pairs.
    _dupes = int(lead_axes_df.duplicated(["axis", "nerve_side"]).sum())
    if _dupes:
        raise RuntimeError(
            f"[FAIR-ALERT] {_dupes} duplicated (axis, nerve_side) rows in "
            f"nerve_immune_lead_axes_postfix.csv — the panel's key assumption is broken."
        )

    # Independent recomputation, kept as a REGRESSION CHECK rather than a caveat.
    #
    # `rule nerve_immune_lead_axes` now produces this file, and both magnitude and
    # row-count columns are arm-labelled and fully populated, so the active arm's
    # curated numbers and the numbers recomputed here should agree on every row.
    # A disagreement is no longer "the generator used a different rule" — it means
    # the rule and this notebook have drifted apart, and it is a bug to chase.
    #
    # The cost of keeping this: the selection rule below duplicates the one in
    # workflow/scripts/nerve_immune_lead_axes.py, so the two must be edited together.
    # That duplication is the point — an oracle that shares an implementation with the
    # thing it checks cannot catch anything (CLAUDE.md, Test Oracle).
    #
    # The rule: cellphone_pvals <= 0.05 and magnitude_rank <= 0.05, restricted to the
    # axis's own interfaces, plus a per-compartment-side clause (see below).
    live_mag_max = 0.05
    live_pval_max = 0.05
    neuron_group = "nerve_neuron"
    _lig = interactions_df["ligand_complex"].astype(str).to_numpy()
    _rec = interactions_df["receptor_complex"].astype(str).to_numpy()
    _mag = interactions_df["magnitude_rank"].to_numpy()
    _pv = interactions_df["cellphone_pvals"].to_numpy()
    _qc = interactions_df["batch_qc_pass"].astype(bool).to_numpy()
    _nerve_cluster = interactions_df["nerve_cluster"].fillna("").astype(str).to_numpy()

    # The two compartment sides are aggregated under DIFFERENT rules, so the oracle
    # has to branch the same way the rule does. Applying the oligodendrocyte rule to a
    # neuron row would "disagree" on every one of the 55 — the neuron half is built
    # WITHOUT the batch-QC filter (the pooled group fails donor QC in the full arm)
    # and restricted to the pooled `nerve_neuron` group.
    _NEURON_SIDE = "neuron (pooled; donor-caveated)"

    _rows = []
    for _r in lead_axes_df.itertuples(index=False):
        _a, _, _b = str(_r.axis).partition("|")
        _a, _b = _a.strip(), _b.strip()
        # Either orientation: the file records an unordered axis, and measured across
        # the 183 rows 109 are ligand-first and 70 receptor-first.
        _m = ((_lig == _a) & (_rec == _b)) | ((_lig == _b) & (_rec == _a))
        _m &= (_pv <= live_pval_max) & (_mag <= live_mag_max)
        if str(_r.compartment_side) == _NEURON_SIDE:
            _m &= _nerve_cluster == neuron_group
        else:
            _m &= _qc
        # Deliberately NO restriction to the row's own `interfaces`. The rule derives
        # that column FROM these rows rather than filtering on it, so imposing it here
        # would not mirror the rule — and it silently undercounts, because each half's
        # `interfaces` is computed on one arm while this recomputation runs on both
        # (it cost 5/183 on the capped arm and 1/183 on the full one).
        _rows.append({
            "axis": _r.axis,
            "nerve_side": _r.nerve_side,
            "tier_v2": _r.tier_v2,
            "rank_active": getattr(_r, _rank_col),
            f"rank_{other_arm}": getattr(_r, _other_rank_col),
            "curated_best_mag": getattr(_r, _mag_col),
            "live_best_mag": float(_mag[_m].min()) if _m.any() else np.nan,
            "curated_n_rows": int(getattr(_r, _n_rows_col)),
            "live_n_rows": int(_m.sum()),
            "interfaces": _r.interfaces,
            "immune": _r.immune,
            "compartment_side": _r.compartment_side,
            "agents_flagged": _r.agents_flagged,
            "glioma_trials": _r.glioma_trials,
            "withdrawn": _r.withdrawn,
        })
    lead_live_df = pd.DataFrame(_rows)

    # Both magnitude columns are now populated on every row, so every row is
    # comparable in both arms — there is no "other-arm only" class left to flag.
    lead_live_df["n_rows_agrees"] = (
        lead_live_df["curated_n_rows"] == lead_live_df["live_n_rows"])
    lead_live_df["mag_agrees"] = np.isclose(
        lead_live_df["curated_best_mag"].astype(float),
        lead_live_df["live_best_mag"].astype(float),
        rtol=1e-6, equal_nan=False,
    )
    lead_live_df = lead_live_df.sort_values(
        ["rank_active", "nerve_side"]).reset_index(drop=True)
    return lead_live_df, live_mag_max, other_arm


@app.cell
def _lead_axes_filters(lead_live_df, mo):
    """Widgets local to Panel E — deliberately separate from the Panel C/D block."""
    _tiers = sorted(lead_live_df["tier_v2"].astype(str).unique())
    _sides = sorted(lead_live_df["nerve_side"].astype(str).unique())
    tier_select = mo.ui.multiselect(options=_tiers, value=_tiers, label="Druggability tier")
    side_select = mo.ui.multiselect(options=_sides, value=_sides, label="Nerve side")
    disagree_only = mo.ui.checkbox(
        value=False, label="Only axes where curated ≠ live")
    return disagree_only, side_select, tier_select


@app.cell
def _lead_axes_view(
    dataset, disagree_only, lead_axes_df, lead_live_df, live_mag_max, mo,
    other_arm, side_select, tier_select,
):
    """Render the shortlist and assert it against an independent recomputation."""
    _df = lead_live_df[
        lead_live_df["tier_v2"].astype(str).isin(tier_select.value or [])
        & lead_live_df["nerve_side"].astype(str).isin(side_select.value or [])
    ]
    if disagree_only.value:
        _bad = ~lead_live_df["mag_agrees"] | ~lead_live_df["n_rows_agrees"]
        _df = _df[_bad.reindex(_df.index, fill_value=False)]

    _n = len(lead_live_df)
    _n_rows_ok = int(lead_live_df["n_rows_agrees"].sum())
    _mag_ok = int(lead_live_df["mag_agrees"].sum())
    _regression = _mag_ok < _n or _n_rows_ok < _n
    _n_neuron = int(lead_live_df["compartment_side"].astype(str)
                    .str.contains("donor-caveated").sum())
    _n_cxcr4 = int(lead_axes_df["axis"].str.contains("CXCR4", na=False).sum())
    _n_lrp1 = int(lead_axes_df["axis"].str.contains("LRP1", na=False).sum())
    _n_s1pr1 = int(lead_axes_df["axis"].str.contains("S1PR1", na=False).sum())

    _cols = [
        "rank_active", f"rank_{other_arm}", "axis", "nerve_side", "tier_v2",
        "curated_best_mag", "live_best_mag", "mag_agrees",
        "curated_n_rows", "live_n_rows", "n_rows_agrees",
        "interfaces", "immune", "agents_flagged", "glioma_trials", "withdrawn",
    ]

    mo.vstack([
        mo.hstack([tier_select, side_select, disagree_only], gap=2),
        mo.md(f"### Lead axes for `{dataset}` — showing **{len(_df):,}** of {_n}"),
        mo.ui.table(_df[[c for c in _cols if c in _df.columns]],
                    selection=None, page_size=20),
        mo.callout(
            mo.md(
                f"**This shortlist is usable where its predecessor was not.** All "
                f"{_n} axes are present in **both** Census arms and **none** of their "
                f"magnitudes match the pinned v1.3.0 reference — the exact inverse of the "
                f"2026-07-15 `nerve_crosstalk_lead_targets.csv`, which was withdrawn from "
                f"this notebook on 2026-08-07 because it was reference-derived (34 of its "
                f"40 magnitudes reproduced against that unaudited cohort, 0 against this "
                f"one).\n\n"
                f"**The live columns are a regression check, not a caveat.** "
                f"`rule nerve_immune_lead_axes` produces this table, and the recomputation "
                f"here re-derives it independently from `{dataset}`'s LR table: "
                f"`cellphone_pvals ≤ {live_mag_max}` and `magnitude_rank ≤ {live_mag_max}`, "
                f"batch-QC-passing on the oligodendrocyte side, restricted to the pooled "
                f"`nerve_neuron` group and **without** the QC filter on the neuron side — "
                f"because those are the two rules the shortlist is actually built under. "
                f"Magnitude agrees on **{_mag_ok}/{_n}** and row counts on "
                f"**{_n_rows_ok}/{_n}**. Both should read {_n}/{_n}: a `False` in "
                f"`mag_agrees` or `n_rows_agrees` now means **the rule and this notebook "
                f"have drifted apart**, and is a bug to chase rather than a caveat to read "
                f"past.\n\n"
                + (
                    f"**⚠️ REGRESSION: {_n - _mag_ok} magnitude and {_n - _n_rows_ok} "
                    f"row-count mismatches.** Tick *Only axes where curated ≠ live* to "
                    f"isolate them, then reconcile "
                    f"`workflow/scripts/nerve_immune_lead_axes.py` against this cell.\n\n"
                    if _regression else ""
                )
                + f"**The two compartment sides are not comparable to each other.** "
                f"`full_best_mag` / `capped_best_mag` now say which arm they came from, "
                f"but the oligodendrocyte side is aggregated over QC-passing rows from all "
                f"nerve clusters and the neuron side over unfiltered rows from the pooled "
                f"neuron group alone. The arm is explicit; the selection rule still differs "
                f"by `compartment_side`, and no column can paper over that.\n\n"
                f"**{_n_neuron} axes sit on the pooled neuron group**, marked "
                f"`neuron (pooled; donor-caveated)`. One donor supplies roughly half of "
                f"all ~4.3k neurons, so those axes are substantially one patient's "
                f"biology — see Panel B and §2.2.\n\n"
                f"**On the three axes that started this investigation:** CXCR4 appears on "
                f"**{_n_cxcr4}** axes here and LRP1 on **{_n_lrp1}**, but S1PR1 on "
                f"**{_n_s1pr1}**. The corrected tables contain zero S1PR1 rows at all — it "
                f"did not move compartments, it failed to reach significance anywhere once "
                f"the compartments were clean.\n\n"
                f"*`min_pval` is deliberately not rendered.* It is `0.0` in every row (the "
                f"permutation floor at `n_perms=1000`), so it cannot discriminate, and a "
                f"column that is entirely constant invites a reader to infer meaning from "
                f"it. `withdrawn` **is** now rendered — it was empty in all {_n} rows until "
                f"the 2026-08-19 annotation snapshot recovered it from the `[WITHDRAWN]` "
                f"tags inside `agents_flagged`."
            ),
            kind="danger" if _regression else "info",
        ),
    ])
    return


# Panel F — whole-table overlap with the pinned v1.3.0 reference — was removed on
# 2026-08-19. It rendered a Jaccard/Spearman comparison against that reference, under
# a callout explaining at length that the reference is not a valid comparator: its
# nerve compartment was built by the exact logic this pipeline removed, it carries no
# author annotation so has never been audited, it was scored against
# differently-normalized data (defect D8), and it is structurally unreproducible.
#
# A panel whose own callout tells the reader not to use its numbers is an invitation
# to use them anyway. Cross-arm agreement — every axis in Panel E clears the bar in
# BOTH Census arms — is the replication evidence this cohort actually has, and it does
# not need the reference to stand up. Policy is unchanged (§2.5 option (a)): v1.3.0
# stays pinned and is simply no longer compared against.
#
# `ds_cohort_concordance` still runs as a pipeline rule and still writes its summary
# and shared-pairs table; this notebook just no longer reads them.


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
def _footer(dataset, mo, nerve_purity_df):
    _n_fail = int((~nerve_purity_df["pass_overall"].astype(bool)).sum())
    mo.md(f"""
    ---
    **Cohort.** `{dataset}` — CELLxGENE Census GBM 10x, 170 donors. The reference
    cohort's notebook is archived at
    `notebooks/archive/04_tme_nerve_immune_explorer.py`. That reference was *not*
    rebuilt by the compartment fix and is **not** a valid comparator for anything here
    — replication evidence on this page is cross-arm agreement, not agreement with it.

    **Sibling notebooks.** `01_census_cohort_qc`, `02_census_nerve_enrichment`,
    `06_census_compartment_audit`.

    **Limits carried by this cohort, all surfaced above.**
    1. No clinical metadata — `gdc_clinical.tsv` here is a generated stub, so the
       clinical-association panel from the archived notebook 04 has no counterpart and
       cannot be rebuilt for this cohort.
    2. The lead axes in Panel E come from `nerve_immune_lead_axes_postfix.csv`, built
       by `rule nerve_immune_lead_axes` and Census-derived across both arms (the
       withdrawn 2026-07-15 shortlist was not). Its four drug columns are a **pinned
       snapshot re-derived from ChEMBL_37 and the ClinicalTrials.gov v2 API** by
       `rule refresh_drug_annotation`, adopted 2026-08-19 — every value traces to a
       recorded database release, but it is a snapshot rather than a live query and
       glioma-trial coverage is bounded by a hand-curated 20-agent list. See
       `reference/drug_annotation/MANIFEST.json`.
    3. `astrocyte`, `opc`, generic `neuron` and `ependymal` are masked out of the
       nerve compartment by decision (Panel A). Their absence from the interaction
       tables is **not** evidence that they do not participate in crosstalk.
    4. {_n_fail} of {len(nerve_purity_df)} nerve groups fail batch QC, and the pooled
       `neuron` group is donor-dominated on both arms (Panel B).
    5. The two compartment sides in Panel E are aggregated under different selection
       rules (oligodendrocyte: QC-passing, all clusters; neuron: unfiltered, pooled
       group only). Magnitudes are comparable *within* a side, not across the two.

    **FAIR.** Inputs from `ds_nerve_tumor_immune_interaction`, `ds_annotate_cluster_qc`,
    `ds_nerve_cluster_annotations` and `nerve_immune_lead_axes`. Exports carry a
    `.provenance.txt` sidecar hashing every input.
    """)
    return


if __name__ == "__main__":
    app.run()
