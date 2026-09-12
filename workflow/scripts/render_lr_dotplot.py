"""Render a LIANA ligand-receptor dotplot from an existing interactions table.

Biological question addressed: which specific ligand-receptor pairs carry the
strongest cross-compartment signal, and between which cell groups — the per-pair
detail behind the compartment-level counts in the significance heatmap.

Split out of `nerve_tumor_immune_interaction.py` on 2026-09-11. The figure used to
be an output of the LIANA rule itself, which meant a purely cosmetic fix to the plot
could only be applied by re-running a seeded 1000-permutation test over ~950k cells
— rewriting the headline interaction tables and cascading through every downstream
rule to change a figure. Reading the table back from disk decouples the two: this
renders in seconds and cannot perturb the data it draws.

**The bug this fixes.** `li.pl.dotplot` facets by `source`, with `target` on the
x-axis inside each facet. The old call sized the figure from the number of cell
labels in the *cohort* (`0.5 * combined.obs["cell_label"].nunique()`), which has
nothing to do with the plotted grid: the top-N subset typically spans ~16 sources
and ~4 targets, i.e. ~64 x-positions crammed into ~13 inches. The facet strip titles
then overlapped into unreadable runs — `erve_c0erve_c1erve_c1` — because each strip
had ~0.8 in for an 8-character label. Width is now derived from the grid actually
being drawn, with a floor per facet wide enough for its own strip text.
"""

import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

import liana as li  # noqa: E402

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

# Layout constants, in inches. Derived from the failure mode rather than guessed:
# a facet must be at least as wide as its own strip label, and each x tick inside a
# facet needs room for a rotated category name.
IN_PER_TARGET = 0.28        # width per x-axis category within one facet
# 0.28 not 0.42: plotnine rotates the x tick labels to vertical here, so a
# category needs roughly a text height, not a text width. At 0.42 the two-way
# plot hit the MAX_W cap and rendered as a 6630x1072 px strip.
IN_PER_STRIP_CHAR = 0.085   # width per character of the facet strip title
IN_PER_ROW = 0.55           # height per DISTINCT interaction on the y-axis
MARGIN_W = 4.0              # y-axis labels + legends
MARGIN_H = 2.0
MIN_H = 5.0                 # the legend stack needs this much regardless of rows
MAX_W = 44.0                # keep the PNG under a sane pixel count at dpi=150


def figure_size(plot_df: pd.DataFrame) -> tuple[float, float]:
    """Size the canvas from the facet grid actually drawn, not the cohort's size."""
    sources = plot_df["source"].astype(str).unique()
    n_targets = max(1, plot_df["target"].astype(str).nunique())

    # Each facet must fit BOTH its x categories and its own strip title.
    longest_label = max((len(s) for s in sources), default=8)
    per_facet = max(n_targets * IN_PER_TARGET, longest_label * IN_PER_STRIP_CHAR)
    width = min(MAX_W, MARGIN_W + len(sources) * per_facet)

    # Height is driven by DISTINCT interactions, not by the number of rows selected.
    # The y-axis is one row per ligand->receptor pair, and the top-N rows collapse
    # hard: the full arm's top 25 rows are only 3 distinct axes spread across 14
    # sources. Sizing on n_rows left ~75% of the canvas empty.
    n_interactions = max(
        1, len(plot_df[["ligand_complex", "receptor_complex"]].drop_duplicates()))
    height = max(MIN_H, MARGIN_H + n_interactions * IN_PER_ROW)
    return float(width), float(height)


def main() -> None:
    params = snakemake.params  # noqa: F821
    top_n = int(getattr(params, "top_n", 25))
    log = snakemake.log[0]  # noqa: F821
    rule = snakemake.rule  # noqa: F821

    lr = pd.read_csv(snakemake.input.lr_table, low_memory=False)  # noqa: F821
    if lr.empty:
        raise ValueError(
            f"[FAIR-ALERT] {snakemake.input.lr_table} is empty — "  # noqa: F821
            "nothing to plot, and an empty figure would read as 'no interactions'."
        )

    n_pairs = min(top_n, len(lr))
    # The on-disk CSV is ordered for readability, not by effect. The plot wants one
    # cohort-wide top list, so rank globally here (matching the original behaviour).
    plot_df = lr.sort_values("magnitude_rank").head(n_pairs).copy()

    width, height = figure_size(plot_df)
    n_facets = plot_df["source"].nunique()
    n_interactions = len(
        plot_df[["ligand_complex", "receptor_complex"]].drop_duplicates())
    log_transformation(
        log, rule,
        f"Plotting top-{n_pairs} of {len(lr)} rows "
        f"({n_interactions} distinct interactions): {n_facets} source facets x "
        f"{plot_df['target'].nunique()} targets -> figure {width:.1f}x{height:.1f} in",
    )

    out = Path(snakemake.output.dotplot)  # noqa: F821
    out.parent.mkdir(parents=True, exist_ok=True)

    try:
        fig = li.pl.dotplot(
            liana_res=plot_df,
            colour="magnitude_rank",
            size="specificity_rank",
            inverse_colour=True,
            inverse_size=True,
            top_n=n_pairs,
            orderby="magnitude_rank",
            orderby_ascending=True,
            figure_size=(width, height),
        )
        # Shrink the strip text as a second line of defence: the width calculation
        # above assumes ~8.5 pt, and plotnine's default is larger.
        try:
            from plotnine import element_text, theme
            fig = fig + theme(
                strip_text_x=element_text(size=8),
                axis_text_x=element_text(size=7),
                axis_text_y=element_text(size=8),
            )
        except ImportError:
            log_transformation(
                log, rule,
                "plotnine not importable for theme override; "
                "relying on the width calculation alone",
                status="WARNING",
            )
        # `limitsize=False` is required, not optional: plotnine hard-refuses any
        # canvas over 25 inches, and a legible facet grid here is ~27 in wide (14
        # source facets x 4 targets). Without it the save raises and this rule
        # silently degrades to the fallback bar chart — which is exactly what the
        # first run of this script did.
        fig.save(str(out), dpi=150, bbox_inches="tight", verbose=False,
                 limitsize=False)
    except Exception as exc:  # plotnine API differences across versions
        log_transformation(
            log, rule,
            f"WARNING: LIANA dotplot failed ({exc}); writing fallback bar chart",
            status="WARNING",
        )
        fig, ax = plt.subplots(figsize=(10, MARGIN_H + n_pairs * 0.35))
        rev = plot_df.iloc[::-1]
        label = (
            rev["source"].astype(str) + " → " + rev["target"].astype(str)
            + " | " + rev["ligand_complex"].astype(str)
            + "→" + rev["receptor_complex"].astype(str)
        )
        ax.barh(label.values,
                -np.log10(rev["magnitude_rank"].clip(lower=1e-6).values))
        ax.set_xlabel("-log10(magnitude_rank)")
        ax.set_title(f"Top {n_pairs} LR pairs (LIANA+ consensus)")
        fig.tight_layout()
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)

    verify_artifact(str(out), min_size_bytes=4096)

    prov = stamp_artifact(
        output_path=str(out),
        rule_name=rule,
        input_paths=[snakemake.input.lr_table],  # noqa: F821
        tool_versions={"liana": li.__version__, "pandas": pd.__version__},
        parameters={
            "top_n": n_pairs,
            "n_source_facets": int(n_facets),
            "n_targets": int(plot_df["target"].nunique()),
            "n_distinct_interactions": int(n_interactions),
            "figure_size_in": [round(width, 2), round(height, 2)],
            "ranked_by": "magnitude_rank",
        },
        description=("LIANA ligand-receptor dotplot rendered from an existing "
                     "interactions table; no inference is performed here"),
        ontology_operation="operation:0337",  # EDAM: Visualisation
    )
    write_provenance(prov, snakemake.output.provenance)  # noqa: F821

    log_transformation(log, rule, "Complete", status="SUCCESS",
                       artifact_paths=[str(out)])


main()
