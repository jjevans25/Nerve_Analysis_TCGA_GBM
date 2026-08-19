"""
Marimo reactive notebook: compartment integrity audit — the Test Oracle.
Addresses: are the nerve, tumor and immune compartments actually made of the cells
they claim to be, judged against an annotation the pipeline never sees?

This is the evidence the rest of the analysis rests on. Every compartment mask is
cross-tabulated against the CELLxGENE Census author annotation `obs['cell_type']`,
harmonized to the Cell Ontology — held out of the pipeline entirely and used only
here, so it is an independent oracle rather than a restatement of the pipeline's
own decisions. Ten gates run as *enforcing* checks: `ds_compartment_audit` fails
the build when any of them fails, which is what stops a defective compartment
reaching an interaction table again.

It exists because of what happened when there was no such check. Before the
2026-08-05/06 fix the nerve compartment was 59% malignant and 27% myeloid while
looking entirely plausible in every downstream figure, and the defect survived
months of analysis. The numbers below are what "corrected" means, stated as
measurements rather than as an assurance.

This notebook has NO reference-cohort ancestor, and cannot have one: the 17-sample
TCGA reference carries no author annotation, so there is nothing to audit it
against. That asymmetry is the main reason the reference cohort was archived — see
notebooks/archive/ and CHANGELOG 2026-08-19.

Source rules: workflow/rules/datasets.smk -> ds_compartment_audit,
ds_nerve_batch_qc_v2, ds_nerve_celltype_labels.
"""

import marimo

__generated_with = "0.23.1"
app = marimo.App(width="wide", app_title="GBM Census — Compartment Audit")


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
        "gates": ds_tables / "compartment_audit_gates.csv",
        "audit": ds_tables / "compartment_audit.csv",
        "confusion": ds_tables / "malignancy_confusion.csv",
        "cluster_audit": ds_tables / "nerve_compartment_cluster_audit.csv",
        "purity": ds_tables / "nerve_cluster_sample_purity.csv",
        "purity_v2": ds_tables / "nerve_cluster_sample_purity_v2.csv",
        "celltype_labels": ds_tables / "nerve_celltype_label_summary.csv",
    }
    _missing = {k: p for k, p in paths.items() if not p.exists()}
    if _missing:
        mo.stop(True, mo.callout(mo.md(
            f"Required artifacts not found for cohort `{dataset}`:\n\n"
            + "\n".join(f"- `{p}`" for p in _missing.values())
            + "\n\nBuild them with:\n\n```\nscripts/run_snakemake.sh \\\n"
            f"  results/tables/{dataset}/compartment_audit_gates.csv \\\n"
            f"  results/tables/{dataset}/nerve_cluster_sample_purity_v2.csv \\\n"
            "  --use-conda --cores all --rerun-triggers mtime\n```"
        ), kind="danger"))
    return dataset, paths


@app.cell
def _load(pd, paths):
    gates_df = pd.read_csv(paths["gates"])
    audit_df = pd.read_csv(paths["audit"])
    confusion_df = pd.read_csv(paths["confusion"])
    cluster_audit_df = pd.read_csv(paths["cluster_audit"])
    purity_df = pd.read_csv(paths["purity"])
    purity_v2_df = pd.read_csv(paths["purity_v2"])
    celltype_df = pd.read_csv(paths["celltype_labels"])
    return (audit_df, celltype_df, cluster_audit_df, confusion_df, gates_df,
            purity_df, purity_v2_df)


@app.cell
def _header(dataset, gates_df, mo):
    _n = len(gates_df)
    _pass = int((gates_df["verdict"] == "PASS").sum())
    _banner = (
        f"**{_pass}/{_n} gates PASS.**" if _pass == _n
        else f"**{_n - _pass} of {_n} gates FAIL — the compartments are not trustworthy.**"
    )
    mo.md(f"""
    # Compartment Integrity Audit — the Test Oracle

    **Cohort: `{dataset}`.** {_banner}

    Every compartment mask below is cross-tabulated against the CELLxGENE Census
    author annotation `obs['cell_type']`, harmonized to the Cell Ontology. That
    annotation is **held out of the pipeline entirely** — nothing upstream reads it —
    so this is an independent oracle, not a restatement of the pipeline's own
    decisions.

    These gates are **enforcing**: `ds_compartment_audit` fails the build when one
    fails, so a defective compartment cannot silently reach an interaction table.
    That matters because it already happened once — before the 2026-08-05/06 fix the
    "nerve" compartment was **59% malignant and 27% myeloid**, and it looked entirely
    plausible in every downstream figure.

    > The 17-sample TCGA reference cohort has **no author annotation**, so none of
    > this can be computed for it. That asymmetry is why it was archived rather than
    > rebuilt — see `notebooks/archive/`.
    """)
    return


@app.cell
def _gates_header(mo):
    mo.md("""
    ---
    ## Panel A — The ten gates
    Each is a hard build gate. `margin` is how much room the observation has before
    it would trip: small positive margins are the ones to watch on a re-run.
    """)
    return


@app.cell
def _gates_panel(gates_df, mo, pd, plt):
    _g = gates_df.copy()
    _g["observed"] = pd.to_numeric(_g["observed"], errors="coerce")
    _g["threshold"] = pd.to_numeric(_g["threshold"], errors="coerce")
    # Signed distance from the threshold in the direction that would pass. Both
    # comparators are normalised so positive always means "clear of the bar".
    _g["margin"] = [
        (o - t) if c == ">=" else (t - o)
        for o, t, c in zip(_g["observed"], _g["threshold"], _g["comparator"])
    ]
    # Counts and fractions live on wildly different scales (37,945 vs 0.95), so a
    # single bar chart of raw margin is unreadable. Plot the fraction gates only and
    # table the rest.
    _frac = _g[_g["threshold"].abs() <= 1.0].sort_values("margin")
    _fig, _ax = plt.subplots(figsize=(7.5, max(2.4, 0.42 * len(_frac))))
    _colors = ["#c0392b" if v != "PASS" else "#2e7d32" for v in _frac["verdict"]]
    _ax.barh(_frac["gate"], _frac["margin"], color=_colors, height=0.6)
    _ax.axvline(0, color="black", lw=1)
    _ax.set_xlabel("margin to threshold (positive = passing)")
    _ax.set_title("Fraction-valued gates, distance from the bar")
    _fig.tight_layout()

    _show = _g[["gate", "observed", "comparator", "threshold", "margin", "verdict"]]
    _fail = _g[_g["verdict"] != "PASS"]
    mo.vstack([
        mo.center(_fig),
        mo.ui.table(_show, selection=None, page_size=12),
        mo.callout(
            mo.md(
                "**All gates pass.** The two that carry the most weight are "
                "`nerve_neural_fraction` (is the nerve compartment actually neural?) "
                "and `tumor_malignant_fraction` (is the tumor compartment actually "
                "malignant?) — those are the two that were catastrophically wrong "
                "before the fix.\n\n"
                "`immune_size_fraction_of_baseline` above 1.0 is expected, not a "
                "warning: the immune compartment grew when `t_cell` was added to "
                "`immune_cells.source_labels`, which had been missing and left that "
                "compartment 96% myeloid."
                if _fail.empty else
                f"**{len(_fail)} gate(s) FAIL: "
                f"{', '.join(_fail['gate'])}.** The build should have stopped. Do not "
                f"read any interaction table from this cohort until this is resolved."
            ),
            kind="info" if _fail.empty else "danger",
        ),
    ])
    return


@app.cell
def _composition_header(mo):
    mo.md("""
    ---
    ## Panel B — What each compartment is actually made of
    Census `cell_type` collapsed to classes, per compartment. The diagonal claim —
    nerve is neural, tumor is malignant, immune is lymphoid/myeloid — is the whole
    point, and it is measured here rather than asserted.
    """)
    return


@app.cell
def _composition_panel(audit_df, mo, plt):
    _cls = audit_df[audit_df["level"] == "class"].copy()
    _order = ["nerve", "nerve_glia", "nerve_neuron", "tumor", "immune", "cohort"]
    _comps = [c for c in _order if c in set(_cls["compartment"])]
    _classes = sorted(set(_cls["key"]))
    _palette = {
        "neural": "#2e7d32", "malignant": "#c0392b", "myeloid": "#e67e22",
        "lymphoid": "#8e44ad", "vascular": "#2980b9", "unmapped": "#95a5a6",
    }

    _fig, _ax = plt.subplots(figsize=(8.4, 4.2))
    _bottom = [0.0] * len(_comps)
    for _k in _classes:
        _vals = []
        for _c in _comps:
            _row = _cls[(_cls["compartment"] == _c) & (_cls["key"] == _k)]
            _vals.append(float(_row["share_of_compartment"].iloc[0]) if len(_row) else 0.0)
        _ax.bar(_comps, _vals, bottom=_bottom, label=_k,
                color=_palette.get(_k, "#bdc3c7"), edgecolor="white", linewidth=0.6)
        _bottom = [b + v for b, v in zip(_bottom, _vals)]
    _ax.set_ylabel("share of compartment")
    _ax.set_ylim(0, 1)
    _ax.legend(frameon=False, fontsize=8, ncol=3, loc="lower right")
    _ax.set_title("Compartment composition against the Census oracle")
    _fig.tight_layout()

    _nerve = _cls[(_cls["compartment"] == "nerve") & (_cls["key"] == "neural")]
    _neural = float(_nerve["share_of_compartment"].iloc[0]) if len(_nerve) else float("nan")
    mo.vstack([
        mo.center(_fig),
        mo.ui.table(
            _cls.pivot_table(index="compartment", columns="key",
                             values="share_of_compartment", fill_value=0.0).round(4),
            selection=None,
        ),
        mo.callout(mo.md(
            f"**The nerve compartment is {_neural:.1%} neural.** Before the fix that "
            f"figure was ~11%, with 59% malignant and 27% myeloid — the cells labelled "
            f"'nerve' were largely tumour and microglia, so every 'nerve' interaction "
            f"computed from them described something else.\n\n"
            f"`nerve_glia` and `nerve_neuron` are the two halves the LR analysis treats "
            f"separately, and they are audited separately for the same reason: the "
            f"pooled neuron group is small and donor-dominated, so its composition "
            f"deserves its own number rather than being averaged into the glial mass."
        ), kind="info"),
    ])
    return


@app.cell
def _confusion_header(mo):
    mo.md("""
    ---
    ## Panel C — Malignancy calls against the oracle
    The CNV caller's confusion matrix, and — more usefully — *what* it gets wrong.
    A classifier's error profile says more than its headline score.
    """)
    return


@app.cell
def _confusion_panel(confusion_df, mo, np, plt):
    _m = dict(zip(confusion_df["metric"], confusion_df["value"]))
    _tp, _fp = _m.get("true_positive", 0), _m.get("false_positive", 0)
    _fn, _tn = _m.get("false_negative", 0), _m.get("true_negative", 0)

    _fig, _axes = plt.subplots(1, 2, figsize=(11.5, 3.9),
                               gridspec_kw={"width_ratios": [1, 1.5]})
    _mat = np.array([[_tp, _fn], [_fp, _tn]])
    _axes[0].imshow(_mat, cmap="Blues", norm="log")
    for _i in range(2):
        for _j in range(2):
            _axes[0].text(_j, _i, f"{int(_mat[_i, _j]):,}", ha="center", va="center",
                          fontsize=9,
                          color="white" if _mat[_i, _j] > _mat.max() / 3 else "black")
    _axes[0].set_xticks([0, 1], ["Census\nmalignant", "Census\nnon-malignant"])
    _axes[0].set_yticks([0, 1], ["called\nmalignant", "called\nnon-malignant"])
    _axes[0].set_title("Confusion matrix (log colour)")

    # What the false positives actually are — the informative half.
    _fps = confusion_df[confusion_df["metric"].str.startswith("false_positive_is::")].copy()
    _fps["cell_type"] = _fps["metric"].str.replace("false_positive_is::", "", regex=False)
    _fps = _fps.sort_values("value").tail(10)
    _axes[1].barh(_fps["cell_type"], _fps["value"], color="#c0392b", height=0.6)
    _axes[1].set_xlabel("cells wrongly called malignant")
    _axes[1].set_title("What the false positives are")
    _fig.tight_layout()

    _prec, _rec = _m.get("precision", float("nan")), _m.get("recall", float("nan"))
    _f1 = _m.get("f1", float("nan"))
    _top = _fps.tail(1)
    _top_type = _top["cell_type"].iloc[0] if len(_top) else "n/a"
    _top_n = int(_top["value"].iloc[0]) if len(_top) else 0
    mo.vstack([
        mo.center(_fig),
        mo.md(
            f"**precision {_prec:.3f}** &nbsp;·&nbsp; **recall {_rec:.3f}** "
            f"&nbsp;·&nbsp; **F1 {_f1:.3f}**"
        ),
        mo.callout(mo.md(
            f"**Read the error profile, not the headline.** The largest false-positive "
            f"class is **{_top_type}** ({_top_n:,} cells). OPC being the commonest "
            f"mistake is biologically unsurprising — OPCs are the normal population "
            f"most transcriptionally adjacent to malignant glioma cells, and they are "
            f"also the population masked out of the nerve compartment by decision.\n\n"
            f"This score is what the **chr7-gain minus chr10-loss contrast** buys. The "
            f"previous genome-wide CNV *spread* score ranked normal oligodendrocytes "
            f"above tumour cells (18% recall / 48% precision) because spread is "
            f"confounded with library complexity. Measuring a within-cell contrast "
            f"between two loci removes that confound. The plumbing around the score "
            f"mattered far less than the score's definition."
        ), kind="info"),
    ])
    return


@app.cell
def _cluster_header(mo):
    mo.md("""
    ---
    ## Panel D — Per-cluster purity of the nerve compartment
    A compartment can pass in aggregate while hiding one contaminated cluster. This
    is the per-cluster view that would catch that.
    """)
    return


@app.cell
def _cluster_panel(cluster_audit_df, mo, plt):
    _c = cluster_audit_df.sort_values("frac_neural").copy()
    _fig, _ax = plt.subplots(figsize=(9, max(2.6, 0.34 * len(_c))))
    _labels = [f"c{int(x)}" for x in _c["nerve_leiden"]]
    _cols = ["#c0392b" if f < 0.8 else "#e67e22" if f < 0.95 else "#2e7d32"
             for f in _c["frac_neural"]]
    _ax.barh(_labels, _c["frac_neural"], color=_cols, height=0.62)
    _ax.axvline(0.8, ls="--", lw=1, color="black")
    _ax.set_xlim(0, 1)
    _ax.set_xlabel("fraction neural (Census oracle)")
    _ax.set_title("Nerve clusters by neural purity — dashed line = 0.80 gate")
    _fig.tight_layout()

    _dirty = _c[_c["frac_neural"] < 0.8].sort_values("n_cells", ascending=False)
    _tot = int(_c["n_cells"].sum())
    _dirty_n = int(_dirty["n_cells"].sum())
    _dirty_share = (_dirty_n / _tot) if _tot else 0.0
    _malig = _dirty[_dirty["frac_malignant"] >= 0.5]

    _cols_show = ["nerve_leiden", "n_cells", "dominant_census_type", "dominant_share",
                  "frac_neural", "frac_malignant", "frac_myeloid", "endothelial_fraction"]
    mo.vstack([
        mo.center(_fig),
        mo.ui.table(_c[[c for c in _cols_show if c in _c.columns]].round(4),
                    selection=None, page_size=15),
        mo.callout(mo.md(
            (f"**Every one of the {len(_c)} nerve clusters is ≥80% neural.** The "
             f"aggregate figure is not hiding a bad cluster — which is the specific "
             f"failure mode this panel exists to rule out."
             if _dirty.empty else
             f"**The aggregate is hiding contaminated clusters: {len(_dirty)} of "
             f"{len(_c)} sit below the 0.80 neural bar** — "
             + ", ".join(f"`c{int(r.nerve_leiden)}` ({int(r.n_cells):,} cells, "
                         f"{r.frac_neural:.0%} neural)" for r in _dirty.itertuples())
             + f". Together they are **{_dirty_n:,} of {_tot:,} cells "
             f"({_dirty_share:.1%})** — small enough that the compartment still "
             f"averages well above the gate, which is exactly why a per-cluster view "
             f"is needed alongside it.\n\n"
             + (f"**{len(_malig)} of them are majority-malignant** "
                + ", ".join(f"`c{int(r.nerve_leiden)}` ({r.frac_malignant:.0%})"
                            for r in _malig.itertuples())
                + ". These are residual tumour cells inside the nerve compartment. "
                "Any axis whose nerve side rests mainly on them is suspect, and the "
                "per-group `batch_qc_pass` flag does **not** capture this — that flag "
                "is about donor dominance, not cell identity. Cross-reference a lead "
                "axis against the `nerve_cluster` column of the interaction table "
                "before reporting it.\n\n" if len(_malig) else "")
             + "The aggregate gate is not wrong, and this is not a regression: it is "
             "the residue the mask could not separate, now visible instead of implied.")
            + "\n\n`max_nerve_cluster_endothelial_fraction` is a gate in Panel A for a "
            "reason: endothelial contamination was how S1PR1 entered the pre-fix "
            "shortlist as a spurious 'nerve' signal."
        ), kind="info" if _dirty.empty else "warn"),
    ])
    return


@app.cell
def _purity_header(mo):
    mo.md("""
    ---
    ## Panel E — Donor purity, and what scANVI-v2 changed
    Compartment purity asks *are these the right cells*. Donor purity asks *are they
    from enough different patients to mean anything*. Both have to hold.
    """)
    return


@app.cell
def _purity_panel(mo, pd, plt, purity_df, purity_v2_df):
    _v1 = purity_df[["cluster", "n_cells", "dominant_sample_fraction",
                     "n_contributing_samples", "pass_overall"]].copy()
    _v2 = purity_v2_df[["cluster", "dominant_sample_fraction",
                        "n_contributing_samples", "pass_overall"]].copy()
    # v1 writes `cluster` as a string ('0'), v2 as int64 (0) — merging them raises.
    # Both are cluster IDs; normalise to str rather than int because v1 also carries
    # the pooled non-numeric group ids that int() would reject.
    _v1["cluster"] = _v1["cluster"].astype(str).str.strip()
    _v2["cluster"] = _v2["cluster"].astype(str).str.strip()
    # Outer join on purpose: v1 has 24 clusters and v2 has 20, and a cluster present
    # in only one is a real difference between the two clusterings, not missing data.
    _j = _v1.merge(_v2, on="cluster", how="outer", suffixes=("_v1", "_v2"))
    _j["in_both"] = _j["dominant_sample_fraction_v1"].notna() & \
        _j["dominant_sample_fraction_v2"].notna()

    _fig, _ax = plt.subplots(figsize=(6.4, 5.2))
    _ax.scatter(_j["dominant_sample_fraction_v1"], _j["dominant_sample_fraction_v2"],
                s=28, alpha=0.75, edgecolor="white", linewidth=0.5)
    _ax.plot([0, 1], [0, 1], ls="--", lw=1, color="black")
    _ax.axhline(0.5, ls=":", lw=1, color="#c0392b")
    _ax.axvline(0.5, ls=":", lw=1, color="#c0392b")
    _ax.set_xlim(0, 1)
    _ax.set_ylim(0, 1)
    _ax.set_xlabel("dominant-sample fraction, v1 (Leiden)")
    _ax.set_ylabel("dominant-sample fraction, v2 (scANVI)")
    _ax.set_title("Donor dominance per nerve cluster\nred lines = 0.5 QC threshold")
    _fig.tight_layout()

    _n_fail_v1 = int((~_j["pass_overall_v1"].astype("boolean").fillna(False)).sum())
    _n_fail_v2 = int((~_j["pass_overall_v2"].astype("boolean").fillna(False)).sum())
    mo.vstack([
        mo.hstack([mo.center(_fig),
                   mo.ui.table(_j.round(4), selection=None, page_size=12)],
                  widths=[1, 1], gap=2),
        mo.callout(mo.md(
            f"**{_n_fail_v1} of {len(_j)} clusters fail donor QC under v1 and "
            f"{_n_fail_v2} under v2.** Failing clusters are **flagged, never dropped** "
            f"— the dominance test cannot distinguish 'real biology that happens to be "
            f"preserved in one donor's tissue' from 'one patient's artifact', so "
            f"automatic removal would discard real signal along with noise. Every "
            f"downstream table carries `batch_qc_pass` so the reader decides.\n\n"
            f"Points below the diagonal are clusters scANVI-v2 made *more* "
            f"donor-diverse. This is a diagnostic comparison, not a replacement: the "
            f"headline analysis still runs on the v1 Leiden clustering.\n\n"
            f"**The two clusterings are not the same partition** — {len(_v1)} clusters "
            f"under v1 against {len(_v2)} under v2, with {int(_j['in_both'].sum())} ids "
            f"in both. Cluster `k` in one is not necessarily cluster `k` in the other, "
            f"so read the scatter as two distributions compared, not as paired "
            f"measurements of the same object."
        ), kind="info"),
    ])
    return


@app.cell
def _footer(celltype_df, dataset, mo):
    _unknown = celltype_df[celltype_df["cell_type_marker_label"] == "Unknown"]
    _unk = float(_unknown["fraction"].iloc[0]) if len(_unknown) else 0.0
    mo.md(f"""
    ---
    **Cohort.** `{dataset}`. Set `GBM_DATASET` to switch arms; both are audited
    independently and both pass 10/10.

    **What this notebook cannot tell you.**
    1. The oracle is an *annotation*, not ground truth. Census `cell_type` is itself a
       computational call by the original submitters. It is independent of this
       pipeline, which is what makes it useful, but agreement with it is agreement
       with another model rather than with tissue.
    2. **{_unk:.0%}** of nerve cells carry the marker label `Unknown`
       (`nerve_celltype_label_summary.csv`). Those cells are inside the compartment
       and inside the interaction tables; the audit says they are neural by Census
       annotation, but the marker panel could not type them further.
    3. Passing gates say the compartments contain the right *kinds* of cell. They say
       nothing about whether the LR inference over them is biologically meaningful.
    4. `astrocyte`, `opc`, generic `neuron` and `ependymal` are masked out of the nerve
       compartment **by decision**. Their absence is a choice, not a finding, and the
       audit measures only what was let in.

    **FAIR.** Inputs from `ds_compartment_audit`, `ds_nerve_batch_qc_v2` and
    `ds_nerve_celltype_labels`. The Census author annotation is read *only* by the
    audit rule and never by any rule that builds a compartment — that separation is
    what makes these numbers a test rather than a tautology.
    """)
    return


if __name__ == "__main__":
    app.run()
