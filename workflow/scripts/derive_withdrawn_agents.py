"""Recover the `withdrawn` column of the nerve-immune drug annotation snapshot.

Answers: which of an axis's approved agents have been withdrawn from market? That
matters for target triage — an axis whose only approved agents are withdrawn is a
dead repurposing lead, and the shortlist tiers it `1b_approved_withdrawn_only`
rather than `1_approved_available`.

The 2026-08-07 snapshot carries that column EMPTY in all 144 axes. It is a bug in
the out-of-band generator, not a property of the data: `axis_drug`
(claude_science/code_export/02_chembl_druggability_pipeline.py:182) sourced it from
`DRUG.withdrawn_agents`, which was populated from ChEMBL's bulk chembl_id lookup —
and that lookup, per the export's own README, "does not carry withdrawal status".
The corrected set came from a later by-name `compound_search` sweep and was never
written back into that column.

The information was never lost, only misplaced: the same sweep tagged the agent
strings inline, so `agents_flagged` still reads e.g.

    CALM1: BENZIODARONE (approved) [WITHDRAWN]; PRENYLAMINE (approved) [WITHDRAWN]

This script parses those tags back into the column the generator intended, in the
format `axis_drug` would have produced (`GENE: AGENT; AGENT` joined by ' | ').
It is a data correction over a committed snapshot, not a re-query — nothing here
touches ChEMBL, so the result is deterministic and offline.

Writes a NEW dated snapshot. The 2026-08-07 file is left untouched: rewriting a
dated record would falsify what that artifact contained on that date.
"""

# NOTE: no `from __future__ import annotations` — Snakemake's `script:` directive
# prepends a preamble, so it would stop being the first statement. See
# download_gene_positions.py.

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

RULE = "derive_withdrawn_agents"

log = snakemake.log[0]  # type: ignore[name-defined]

WITHDRAWN_TAG = "[WITHDRAWN]"
# Strips the parenthetical phase/approval-year that annot_agents left on the name,
# e.g. "PROBUCOL (appr 1977)" -> "PROBUCOL".
_PAREN = re.compile(r"\s*\(.*?\)")


def _withdrawn_for_axis(agents_flagged: str) -> str:
    """Per-gene withdrawn agent names, in the generator's intended output format."""
    if not isinstance(agents_flagged, str) or not agents_flagged.strip():
        return ""
    out = []
    for chunk in agents_flagged.split(" | "):
        gene, _, rest = chunk.partition(": ")
        names = [
            _PAREN.sub("", name).replace(WITHDRAWN_TAG, "").strip()
            for name in rest.split(";")
            if WITHDRAWN_TAG in name
        ]
        if names:
            out.append(f"{gene}: {'; '.join(names)}")
    return " | ".join(out)


source = pd.read_csv(snakemake.input.snapshot)  # type: ignore[name-defined]
n_source = len(source)

if source["withdrawn"].notna().any():
    raise ValueError(
        f"[FAIR-ALERT] {snakemake.input.snapshot} already has "  # type: ignore[name-defined]
        f"{int(source['withdrawn'].notna().sum())} non-empty `withdrawn` values. This "
        f"rule exists to recover a column that is empty in the 2026-08-07 snapshot; "
        f"running it over an already-corrected file would silently overwrite real data."
    )

corrected = source.copy()
corrected["withdrawn"] = corrected["agents_flagged"].map(_withdrawn_for_axis)
# Empty string is the generator's own "no withdrawn agents" value, but writing it for
# 131 of 144 rows would turn a genuinely-absent value into a present-but-blank one.
# Keep those as NA so the column reads the same way every other optional column does.
corrected.loc[corrected["withdrawn"] == "", "withdrawn"] = pd.NA

# --- Guards: this must be a pure column recovery, nothing else may move ----------
untouched = ["axis", "tier_v2", "agents_flagged", "glioma_trials"]
for col in untouched:
    if not source[col].fillna("\x00").equals(corrected[col].fillna("\x00")):
        raise ValueError(
            f"[FAIR-ALERT] column `{col}` changed while deriving `withdrawn`. This rule "
            f"must only populate `withdrawn`; every other column is passed through."
        )
if len(corrected) != n_source:
    raise ValueError(f"[FAIR-ALERT] row count changed: {n_source} -> {len(corrected)}")

n_filled = int(corrected["withdrawn"].notna().sum())
tagged = corrected["agents_flagged"].fillna("").str.contains(re.escape(WITHDRAWN_TAG))
if n_filled != int(tagged.sum()):
    raise ValueError(
        f"[FAIR-ALERT] {n_filled} axes got a `withdrawn` value but "
        f"{int(tagged.sum())} carry a {WITHDRAWN_TAG} tag — the parser dropped or "
        f"invented rows."
    )
# Every 1b axis is withdrawn-only by definition, so it must carry a value.
missing_1b = corrected[(corrected["tier_v2"] == "1b_approved_withdrawn_only")
                       & corrected["withdrawn"].isna()]
if len(missing_1b):
    raise ValueError(
        f"[FAIR-ALERT] {len(missing_1b)} axes are tiered 1b_approved_withdrawn_only but "
        f"got no `withdrawn` value: {missing_1b['axis'].tolist()[:10]}. The tier and the "
        f"agent tags disagree, which means one of them is wrong."
    )

Path(snakemake.output.snapshot).parent.mkdir(parents=True, exist_ok=True)  # type: ignore[name-defined]
corrected.to_csv(snakemake.output.snapshot, index=False)  # type: ignore[name-defined]
verify_artifact(snakemake.output.snapshot, min_size_bytes=1024)  # type: ignore[name-defined]

agents = sorted({
    name.strip()
    for value in corrected["withdrawn"].dropna()
    for chunk in value.split(" | ")
    for name in chunk.partition(": ")[2].split(";")
})
summary = {
    "n_axes": n_source,
    "n_withdrawn_populated": n_filled,
    "distinct_withdrawn_agents": agents,
    "source_snapshot": str(snakemake.input.snapshot),  # type: ignore[name-defined]
    "derivation": (
        "Parsed [WITHDRAWN] tags out of agents_flagged; phase parentheticals stripped. "
        "No network access, no ChEMBL re-query — a correction over committed data."
    ),
}
log_transformation(log, RULE,
                   f"{n_filled}/{n_source} axes populated; agents={agents}")

prov = stamp_artifact(
    output_path=snakemake.output.snapshot,  # type: ignore[name-defined]
    rule_name=RULE,
    input_paths=[snakemake.input.snapshot],  # type: ignore[name-defined]
    tool_versions={"pandas": pd.__version__},
    parameters=summary,
    description=("Nerve-immune axis drug annotation with the `withdrawn` column "
                 "recovered from the inline [WITHDRAWN] agent tags"),
    ontology_operation="operation:3096",  # EDAM: Editing / data correction
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]
log_transformation(log, RULE, "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.snapshot])  # type: ignore[name-defined]
