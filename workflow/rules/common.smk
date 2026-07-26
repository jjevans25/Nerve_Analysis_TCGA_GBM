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


def pinned_target(*parts) -> list:
    """A frozen v1.3.0 artifact — requested only if it still exists on disk.

    When pinned, the producing rules are not defined, so `rule all` must not
    demand an artifact that is absent: Snakemake would raise "No rule to
    produce". Returns [] for the artifacts lost with `nerve_cells.h5ad`.
    """
    path = p(*parts)
    return [path] if (not BASELINE_PINNED or Path(path).exists()) else []
