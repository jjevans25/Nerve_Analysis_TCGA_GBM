# =============================================================
# Lead axes — post-compartment-fix nerve–immune shortlist
#
# Un-wildcarded on purpose: the output carries BOTH census arms in one table
# (rank_full/best_mag and rank_capped/capped_best_mag), so it belongs to neither
# arm's `{dataset}` namespace and sits at the results/tables/ root.
#
# Closes the provenance gap recorded in CHANGELOG 2026-08-07: this table was
# generated out of band, had no producing rule, and `results/` is gitignored —
# so it could neither be rebuilt nor restored from git.
# =============================================================

_LEAD = config["lead_axes"]
_DRUG_DIR = "reference/drug_annotation"


def _arm_table(arm_key, filename):
    """Path to a per-arm cohort-namespaced table."""
    return os.path.join(config["dirs"]["tables"], _LEAD[arm_key], filename)


rule derive_withdrawn_agents:
    """Recover the `withdrawn` column, empty in the 2026-08-07 snapshot due to a
    generator bug (it was read from ChEMBL's bulk chembl_id lookup, which does not
    carry withdrawal status).

    Pure data correction over committed data — parses the inline [WITHDRAWN] tags
    already present in `agents_flagged`. No network access. Writes a NEW dated
    snapshot; the 2026-08-07 file is never modified.
    """
    input:
        snapshot = os.path.join(_DRUG_DIR, "nerve_immune_axis_drug_annotation_2026-08-07.csv"),
    output:
        snapshot   = os.path.join(_DRUG_DIR, "nerve_immune_axis_drug_annotation_2026-08-19.csv"),
        provenance = os.path.join(config["dirs"]["provenance"], "derive_withdrawn_agents_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "derive_withdrawn_agents.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    script:
        "../scripts/derive_withdrawn_agents.py"


rule refresh_drug_annotation:
    """Mint a NEW dated drug-annotation snapshot from public ChEMBL + ClinicalTrials.gov.

    OPT-IN — deliberately not in `rule all`. Requires network access and takes a while
    (one ChEMBL round-trip per gene, per target and per molecule, politely paced).

    Its output WILL differ from the pinned snapshot and is not supposed to match: both
    databases have moved since 2026-08-07, the original ChEMBL release was never
    recorded, and one of the original's inputs is permanently lost. Adopting a refreshed
    snapshot means pointing `config.lead_axes.drug_annotation` at it deliberately, after
    a diff — this rule never mutates the pinned file.

    Fails loudly rather than emitting a partial snapshot: an unresolved gene silently
    becomes `4_no_chembl_target`, which reads downstream as "this target has no
    chemistry" — a scientific claim manufactured by a network error.

    The `_refreshed_` infix is load-bearing: it keeps this rule's output namespace
    disjoint from the pinned/derived snapshots, so a date wildcard can never resolve to
    a committed file and shadow it (2026-08-19 would otherwise match both this rule and
    `derive_withdrawn_agents`). It also makes re-queried snapshots visually distinct
    from derived ones in a directory listing.

    Invoke explicitly, e.g.:
      scripts/run_snakemake.sh \\
        reference/drug_annotation/nerve_immune_axis_drug_annotation_refreshed_$(date +%F).csv \\
        --use-conda --cores 1 --allowed-rules refresh_drug_annotation
    """
    input:
        table   = os.path.join(config["dirs"]["tables"], "nerve_immune_lead_axes_postfix.csv"),
        curated = os.path.join(_DRUG_DIR, "curated_agent_map_2026-08-07.json"),
    output:
        snapshot   = os.path.join(_DRUG_DIR, "nerve_immune_axis_drug_annotation_refreshed_{refresh_date}.csv"),
        manifest   = os.path.join(_DRUG_DIR, "refresh_manifest_{refresh_date}.json"),
        provenance = os.path.join(config["dirs"]["provenance"], "refresh_drug_annotation_{refresh_date}.json"),
    wildcard_constraints:
        refresh_date = r"\d{4}-\d{2}-\d{2}",
    log:
        os.path.join(config["dirs"]["logs"], "refresh_drug_annotation_{refresh_date}.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 2000,
        threads = 1,
    params:
        max_retries = 4,
        timeout_s   = 60,
        pause_s     = 0.34,
    script:
        "../scripts/refresh_drug_annotation.py"


rule nerve_immune_lead_axes:
    """Cross-arm reproducible nerve–immune signalling axes, tiered by druggability.

    An axis is a lead only if it is significant in both census arms. The output is
    two independently-built halves concatenated — the oligodendrocyte side from the
    full arm (batch-QC-passing) and the pooled-neuron side from the capped arm (no
    QC filter, because that group fails donor QC in the full arm). The selection
    rules genuinely differ between halves; see the script docstring.

    Drug columns come from a pinned snapshot, not from a live query — they are not
    recomputable from anything in this repo.
    """
    input:
        full_lr    = _arm_table("full_arm",   "nerve_tumor_immune_interactions_with_qc.csv"),
        capped_lr  = _arm_table("capped_arm", "nerve_tumor_immune_interactions_with_qc.csv"),
        full_ann   = _arm_table("full_arm",   "nerve_cluster_annotations.csv"),
        capped_ann = _arm_table("capped_arm", "nerve_cluster_annotations.csv"),
        drug       = _LEAD["drug_annotation"],
    output:
        table      = os.path.join(config["dirs"]["tables"], "nerve_immune_lead_axes_postfix.csv"),
        provenance = os.path.join(config["dirs"]["provenance"], "nerve_immune_lead_axes_provenance.json"),
    log:
        os.path.join(config["dirs"]["logs"], "nerve_immune_lead_axes.log"),
    conda:
        "../envs/scrna.yaml",
    resources:
        mem_mb  = 8000,
        threads = 1,
    params:
        pval_max                = _LEAD["pval_max"],
        magnitude_rank_max      = _LEAD["magnitude_rank_max"],
        neuron_group            = _LEAD["neuron_group"],
        compartment_side_oligo  = _LEAD["compartment_side_oligo"],
        compartment_side_neuron = _LEAD["compartment_side_neuron"],
    script:
        "../scripts/nerve_immune_lead_axes.py"
