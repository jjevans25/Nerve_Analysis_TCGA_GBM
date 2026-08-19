"""Mint a fresh nerve-immune drug annotation snapshot from ChEMBL and ClinicalTrials.gov.

Answers: for each cross-arm reproducible signalling axis, is there a drug that hits
either end of it, and has anything in that class been tried in glioma? That is the
`tier_v2` / `agents_flagged` / `glioma_trials` / `withdrawn` block of the lead-axes
shortlist.

The pinned 2026-08-07 snapshot cannot be recomputed: it came from a Claude Science
MCP connector whose response caches did not survive, and one of its inputs (a
prior-session artifact addressed only by UUID) is permanently lost. This rule is the
replacement path — it re-derives the same four columns from the two public REST APIs,
with no key and no MCP.

IT WILL NOT REPRODUCE THE PINNED SNAPSHOT, and is not supposed to. Both databases have
moved since 2026-08-07, the original ChEMBL release was never recorded, and the lost
fallback artifact cannot be reconstructed. Output therefore goes to a NEW dated file;
adopting it means pointing `config.lead_axes.drug_annotation` at that file deliberately,
after a diff. Nothing here mutates the pinned snapshot.

Failure semantics matter more than usual here. A partially-resolved snapshot is worse
than no snapshot: a gene whose lookup failed silently drops to `4_no_chembl_target`,
which reads downstream as "this target has no chemistry" — a scientific claim, produced
by a network error. So every gene and every agent must resolve or return a definite
empty result, and any hard failure aborts the rule with no output file written.

Curated inputs (`AGENTS`, `AG2GENE`) come from
`reference/drug_annotation/curated_agent_map_2026-08-07.json`. They bound trial coverage:
an axis shows no trial if its agent is not on that list, which is not the same as no
trial existing.
"""

# NOTE: no `from __future__ import annotations` — Snakemake's `script:` directive
# prepends a preamble, so it would stop being the first statement.

import json
import time
from datetime import datetime, timezone
from pathlib import Path
import sys

import pandas as pd
import requests

sys.path.insert(0, "workflow/scripts")
from fair_utils import (  # noqa: E402
    file_sha256,
    log_transformation,
    stamp_artifact,
    verify_artifact,
    write_provenance,
)

RULE = "refresh_drug_annotation"

log = snakemake.log[0]  # type: ignore[name-defined]
P = snakemake.params  # type: ignore[name-defined]

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
CTGOV = "https://clinicaltrials.gov/api/v2/studies"
MAX_RETRIES = int(getattr(P, "max_retries", 4))
TIMEOUT = int(getattr(P, "timeout_s", 60))
PAUSE = float(getattr(P, "pause_s", 0.34))  # be polite to EBI

# Tier ordering. Keyed on the first two characters, which is why '1_' and '1b' are
# spelled out — see axis_tier2 in the original generator.
TIER_ORDER = {"1_": 1, "1b": 1.5, "2_": 2, "3_": 3, "4_": 4}


def _get(url: str, params: dict) -> dict:
    """GET returning parsed JSON, with exponential backoff. Raises on exhaustion."""
    last = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.get(url, params=params, timeout=TIMEOUT,
                                headers={"Accept": "application/json"})
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException, ValueError) as exc:
            last = exc
            if attempt < MAX_RETRIES:
                wait = 2 ** attempt
                log_transformation(log, RULE,
                                   f"attempt {attempt}/{MAX_RETRIES} for {url} failed "
                                   f"({str(exc)[:120]}); retrying in {wait}s",
                                   status="RETRY")
                time.sleep(wait)
    raise RuntimeError(
        f"[FAIR-ALERT] {url} failed after {MAX_RETRIES} attempts: {last}. Refusing to "
        f"write a partial snapshot — an unresolved gene silently becomes "
        f"'4_no_chembl_target', which reads as a scientific finding rather than a "
        f"network error."
    )


def _axis_genes(axis: str) -> list:
    """Split an undirected axis into genes, expanding LIANA subunit complexes."""
    return [g for end in str(axis).split("|") for g in end.split("_") if g]


# Exact gene-symbol match against a target's component synonyms. NOT `/target/search`:
# that endpoint is full-text, rejects `target_type`/`organism` filters with a 400, and
# 400s outright on short symbols (`q=C3` fails even bare). This filter is both more
# robust and more precise — it resolves each gene to exactly one human target rather
# than a fuzzy ranked list, so a near-miss can never be mistaken for a hit.
_SYNONYM_FILTER = "target_components__target_component_synonyms__component_synonym__iexact"


def _chembl_targets(gene: str) -> list:
    """Human ChEMBL targets for a gene symbol, single-protein first then complexes."""
    for target_type in ("SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN FAMILY"):
        payload = _get(f"{CHEMBL}/target", {
            _SYNONYM_FILTER: gene,
            "organism__iexact": "Homo sapiens",
            "target_type": target_type,
            "format": "json", "limit": 5,
        })
        time.sleep(PAUSE)
        hits = payload.get("targets", [])
        if hits:
            return [{"id": t["target_chembl_id"], "name": t.get("pref_name")}
                    for t in hits]
    return []


def _chembl_mechanism_molecules(target_id: str) -> list:
    payload = _get(f"{CHEMBL}/mechanism", {
        "target_chembl_id": target_id, "format": "json", "limit": 60})
    time.sleep(PAUSE)
    return [m.get("molecule_chembl_id") for m in payload.get("mechanisms", [])
            if m.get("molecule_chembl_id")]


def _phase(raw) -> float:
    """max_phase as a number.

    ChEMBL_37 serialises it as the STRING '4.0', so a naive `phase == 4` is always
    False and every approved drug silently falls through to `3_chembl_target_no_agent`.
    That produced a snapshot claiming zero approved agents across all 156 genes —
    plausible-looking, entirely wrong. Coerce, and treat unparseable as "no phase".
    """
    try:
        return float(raw)
    except (TypeError, ValueError):
        return float("nan")


def _chembl_molecule(molecule_id: str) -> dict:
    """Name, phase and safety flags. Queried per molecule so withdrawal status is real.

    The original generator's bug lived here: its bulk chembl_id lookup did not carry
    `withdrawn_flag`, so the `withdrawn` column came back empty for every axis. The
    per-molecule endpoint does carry it (verified: PROBUCOL -> withdrawn_flag True).
    """
    payload = _get(f"{CHEMBL}/molecule/{molecule_id}", {"format": "json"})
    time.sleep(PAUSE)
    return {
        "name": payload.get("pref_name"),
        "phase": _phase(payload.get("max_phase")),
        "withdrawn": bool(payload.get("withdrawn_flag")),
        "bbw": bool(payload.get("black_box_warning")),
    }


def _glioma_trials(agent: str) -> dict:
    """Glioma trials naming this agent as an intervention.

    Parameter is `query.intr`, not `query.intv`. Verified against the generator's own
    reported result: plerixafor returns 4 trials including NCT01977677, NCT01339039
    and NCT00669669.
    """
    payload = _get(CTGOV, {
        "query.cond": "glioma OR glioblastoma", "query.intr": agent,
        "pageSize": 30, "countTotal": "true",
        "fields": "NCTId,BriefTitle,OverallStatus,Phase",
    })
    time.sleep(PAUSE)
    ncts = [
        s.get("protocolSection", {}).get("identificationModule", {}).get("nctId")
        for s in payload.get("studies", [])
    ]
    ncts = [n for n in ncts if n]
    return {"total": int(payload.get("totalCount", len(ncts))), "ncts": ncts}


# --- Inputs -------------------------------------------------------------------
curated = json.loads(Path(snakemake.input.curated).read_text())  # type: ignore[name-defined]
AGENTS = curated["agents"]
AG2GENE = curated["agent_to_gene"]

table = pd.read_csv(snakemake.input.table)  # type: ignore[name-defined]
axes = sorted(table["axis"].dropna().unique())
genes = sorted({g for a in axes for g in _axis_genes(a)})
log_transformation(log, RULE, f"{len(axes)} axes, {len(genes)} genes, {len(AGENTS)} agents")

# --- ChEMBL release, so this snapshot is pinnable in a way the original was not ---
status = _get(f"{CHEMBL}/status", {"format": "json"})
chembl_release = status.get("chembl_db_version") or status.get("chembl_release")
if not chembl_release:
    raise RuntimeError(
        "[FAIR-ALERT] ChEMBL /status returned no release version. Recording it is the "
        "whole point of this rule — the missing pin is why the 2026-08-07 snapshot is "
        "unreproducible. Refusing to mint another unpinned one."
    )
log_transformation(log, RULE, f"ChEMBL release {chembl_release}")

# --- Gene -> targets -> mechanisms -> molecules --------------------------------
gene_info = {}
for i, gene in enumerate(genes, 1):
    targets = _chembl_targets(gene)
    molecules = []
    for target in targets[:2]:
        molecules.extend(_chembl_mechanism_molecules(target["id"]))
    seen, details = set(), []
    for molecule_id in molecules:
        if molecule_id in seen:
            continue
        seen.add(molecule_id)
        details.append(_chembl_molecule(molecule_id))
    approved = sorted({d["name"] for d in details
                       if d["phase"] == 4 and d["name"] and not d["withdrawn"]})
    withdrawn = sorted({d["name"] for d in details
                        if d["phase"] == 4 and d["name"] and d["withdrawn"]})
    clinical = sorted({d["name"] for d in details
                       if 1 <= d["phase"] <= 3 and d["name"]})
    blackbox = {d["name"] for d in details if d["bbw"] and d["name"]}
    if approved:
        tier = "1_approved_available"
    elif withdrawn:
        tier = "1b_approved_withdrawn_only"
    elif clinical:
        tier = "2_clinical"
    elif targets:
        tier = "3_chembl_target_no_agent"
    else:
        tier = "4_no_chembl_target"
    gene_info[gene] = {"tier": tier, "approved": approved, "withdrawn": withdrawn,
                       "clinical": clinical, "blackbox": blackbox}
    if i % 25 == 0:
        log_transformation(log, RULE, f"resolved {i}/{len(genes)} genes")

# --- Agent -> glioma trials ----------------------------------------------------
gene_trials = {}
for agent in AGENTS:
    result = _glioma_trials(agent)
    if not result["total"]:
        continue
    cell = f"{agent} ({result['total']}: {', '.join(result['ncts'][:3])})"
    for gene in AG2GENE.get(agent, "").replace("(probe)", "").split("/"):
        if gene.strip():
            gene_trials.setdefault(gene.strip(), []).append(cell)

# --- Roll gene facts up to axes ------------------------------------------------
rows = []
for axis in axes:
    members = [g for g in _axis_genes(axis) if g in gene_info]
    tier = (min((gene_info[g]["tier"] for g in members),
                key=lambda t: TIER_ORDER[t[:2]]) if members else "4_no_chembl_target")
    flagged, withdrawn_cells = [], []
    for gene in _axis_genes(axis):
        info = gene_info.get(gene)
        if not info:
            continue
        if info["approved"] or info["withdrawn"]:
            names = [f"{n} [WITHDRAWN]" for n in info["withdrawn"]] + [
                f"{n} [black-box]" if n in info["blackbox"] else n
                for n in info["approved"]]
            flagged.append(f"{gene}: {'; '.join(names)}")
        elif info["clinical"]:
            flagged.append(f"{gene}: [clinical] {'; '.join(info['clinical'])}")
        if info["withdrawn"]:
            withdrawn_cells.append(f"{gene}: {'; '.join(info['withdrawn'])}")
    trials = sorted({c for g in _axis_genes(axis) for c in gene_trials.get(g, [])})
    rows.append({
        "axis": axis,
        "tier_v2": tier,
        "agents_flagged": " | ".join(flagged) or pd.NA,
        "glioma_trials": " | ".join(trials) or pd.NA,
        "withdrawn": " | ".join(withdrawn_cells) or pd.NA,
    })

snapshot = pd.DataFrame(rows).sort_values("axis").reset_index(drop=True)
if len(snapshot) != len(axes):
    raise RuntimeError(
        f"[FAIR-ALERT] built {len(snapshot)} rows for {len(axes)} axes.")

# Plausibility guard. A snapshot with no druggable axis at all is not a finding about
# GBM biology — it is the signature of a parsing failure that resolved targets fine but
# lost every agent, which is exactly what a string-vs-number `max_phase` did on the
# first run: 133 axes tiered 3_chembl_target_no_agent, 0 approved, no error raised.
# Drift from the pinned snapshot is expected; total collapse of the top tiers is not.
_druggable = int(snapshot.tier_v2.isin(
    ["1_approved_available", "1b_approved_withdrawn_only", "2_clinical"]).sum())
if not _druggable:
    raise RuntimeError(
        f"[FAIR-ALERT] 0 of {len(snapshot)} axes resolved to an approved or clinical "
        f"agent. The pinned snapshot has 88. Targets resolved but agents did not, which "
        f"means the mechanism/molecule parse is broken — not that these targets are "
        f"undruggable. Refusing to write a snapshot that would read as a scientific "
        f"finding. Tiers built: {snapshot.tier_v2.value_counts().to_dict()}"
    )

out = Path(snakemake.output.snapshot)  # type: ignore[name-defined]
out.parent.mkdir(parents=True, exist_ok=True)
snapshot.to_csv(out, index=False)
verify_artifact(str(out), min_size_bytes=1024)

manifest = {
    "artifact": out.name,
    "sha256": file_sha256(str(out)),
    "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
    "chembl_release": chembl_release,
    "chembl_api": CHEMBL,
    "clinicaltrials_api": CTGOV,
    "n_axes": int(len(snapshot)),
    "n_genes": len(genes),
    "curated_agent_map": str(snakemake.input.curated),  # type: ignore[name-defined]
    "tier_v2_counts": {str(k): int(v)
                       for k, v in snapshot.tier_v2.value_counts().items()},
    "expect_drift_from_pinned_snapshot": True,
    "notes": [
        "This snapshot is NOT expected to match the 2026-08-07 pinned one. Both "
        "databases have moved, the original ChEMBL release was never recorded, and the "
        "prior-session fallback artifact it used is permanently lost.",
        "Adopting it means pointing config.lead_axes.drug_annotation here deliberately, "
        "after a diff. This rule never mutates the pinned snapshot.",
        "Glioma-trial coverage is bounded by the curated agent list: an empty cell means "
        "no named trial for a LISTED agent, not that no trial exists.",
    ],
}
Path(snakemake.output.manifest).write_text(json.dumps(manifest, indent=2) + "\n")  # type: ignore[name-defined]

log_transformation(log, RULE,
                   f"{len(snapshot)} axes; ChEMBL {chembl_release}; "
                   f"tiers={manifest['tier_v2_counts']}")
prov = stamp_artifact(
    output_path=str(out),
    rule_name=RULE,
    input_paths=[snakemake.input.table, snakemake.input.curated],  # type: ignore[name-defined]
    tool_versions={"pandas": pd.__version__, "requests": requests.__version__},
    parameters=manifest,
    description="Refreshed nerve-immune axis drug annotation from ChEMBL + ClinicalTrials.gov",
    ontology_operation="operation:2422",  # EDAM: Data retrieval
)
write_provenance(prov, snakemake.output.provenance)  # type: ignore[name-defined]
log_transformation(log, RULE, "Complete", status="SUCCESS", artifact_paths=[str(out)])
