# workflow/rules/common.smk
# Shared helpers imported by the main Snakefile before any other rules.

import os
from pathlib import Path


def p(*parts) -> str:
    """Join path parts — avoids string literals starting with '/' in rule directives."""
    return str(Path(*parts))


# ---------------------------------------------------------------------------
# Baseline pin (see CHANGELOG 2026-07-21 and 2026-07-26)
#
# `data/processed/nerve_cells.h5ad` was deleted by a failed `nerve_cell_subset`
# job and is NOT reproducible: the retrained scVI latent (39,779 genes)
# re-clusters differently from published v1.3.0 (20,420), so the frozen
# `cl15_split_v1_3_0.csv` barcodes no longer map to one source cluster and the
# integrity guard hard-fails. Researcher decision: pin the surviving v1.3.0
# tables and never regenerate them.
#
# When `baseline.pinned` is true, the rules that would rebuild those artifacts
# are NOT DEFINED (see nerve_cells.smk / immune.smk). A rule that does not exist
# cannot be scheduled by any invocation, so the pin holds without relying on
# `--allowed-rules` being remembered on the command line. Snakemake then treats
# the surviving artifacts as plain source files.
#
# NOTE: `ancient()` does NOT solve this. It only suppresses timestamp-driven
# re-runs; the missing-file path (dag.py -> Job.missing_output -> f.exists()) is
# ancient-blind, so a deleted file still schedules its producer.
# ---------------------------------------------------------------------------
BASELINE_PINNED = bool(config.get("baseline", {}).get("pinned", False))


# ---------------------------------------------------------------------------
# Annotation marker panels (see the `annotation_markers:` config block).
#
# Kept out of `nerve_cells.markers` on purpose: that key is also a param of
# scrna_qc, which runs once per sample (170 donors x 2 arms), so touching it
# invalidates every QC artifact and cascades into the ~11 h scVI train. These
# panels are merged over the neural ones inside scrna_annotate, and win on
# conflict.
# ---------------------------------------------------------------------------
ANNOTATION_MARKERS = config.get("annotation_markers", {}).get("panels", {})
ANNOTATE_LEIDEN_RESOLUTION = config.get("annotation_markers", {}).get(
    "leiden_resolution", config["nerve_cells"]["leiden_resolution"]
)
ANNOTATE_AMBIGUOUS_MARGIN = config.get("annotation_markers", {}).get("ambiguous_margin", 0.0)
ANNOTATE_PANEL_COMPARTMENT = config.get("annotation_markers", {}).get("panel_compartment", {})

# The immune compartment is a union of lineages, not one label. `source_label`
# (str) stays supported so an older config keeps parsing.
IMMUNE_SOURCE_LABELS = config["immune_cells"].get(
    "source_labels", [config["immune_cells"].get("source_label", "microglia")]
)


def pinned_target(*parts) -> list:
    """A frozen v1.3.0 artifact — requested only if it still exists on disk.

    When pinned, the producing rules are not defined, so `rule all` must not
    demand an artifact that is absent: Snakemake would raise "No rule to
    produce". Returns [] for the artifacts lost with `nerve_cells.h5ad`.
    """
    path = p(*parts)
    return [path] if (not BASELINE_PINNED or Path(path).exists()) else []
