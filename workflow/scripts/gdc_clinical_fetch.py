"""Fetch clinical metadata from the GDC public REST API for all sample file UUIDs."""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, write_provenance

log     = snakemake.log[0]
samples = snakemake.params.sample_ids
base    = snakemake.params.api_base.rstrip("/")
fields  = snakemake.params.api_fields

log_transformation(log, "gdc_clinical_fetch",
    f"Querying GDC API for {len(samples)} file UUIDs at {base}")

# ---------------------------------------------------------------------------
# Step 1: resolve file UUIDs → case UUIDs via /files endpoint
# ---------------------------------------------------------------------------
files_url = f"{base}/files"
files_payload = {
    "filters": {
        "op": "in",
        "content": {"field": "file_id", "value": samples},
    },
    "fields": "file_id,cases.case_id",
    "size": len(samples),
    "format": "JSON",
}

resp = requests.post(files_url, json=files_payload, timeout=60)
resp.raise_for_status()
files_hits = resp.json()["data"]["hits"]

file_to_case: dict[str, str] = {}
for hit in files_hits:
    fid = hit["file_id"]
    cases = hit.get("cases", [])
    if cases:
        file_to_case[fid] = cases[0]["case_id"]

log_transformation(log, "gdc_clinical_fetch",
    f"Resolved {len(file_to_case)}/{len(samples)} file UUIDs to case UUIDs")

case_ids = list(file_to_case.values())

# ---------------------------------------------------------------------------
# Step 2: fetch clinical fields for each case via /cases endpoint
# ---------------------------------------------------------------------------
cases_url = f"{base}/cases"
cases_payload = {
    "filters": {
        "op": "in",
        "content": {"field": "case_id", "value": case_ids},
    },
    "fields": ",".join(fields),
    "expand": "diagnoses,diagnoses.treatments,samples,demographic",
    "size": len(case_ids),
    "format": "JSON",
}

resp2 = requests.post(cases_url, json=cases_payload, timeout=60)
resp2.raise_for_status()
cases_hits = resp2.json()["data"]["hits"]

log_transformation(log, "gdc_clinical_fetch",
    f"Received clinical data for {len(cases_hits)} cases")

# ---------------------------------------------------------------------------
# Step 3: flatten nested JSON into one row per sample
# ---------------------------------------------------------------------------
case_id_to_info: dict[str, dict] = {}
for hit in cases_hits:
    cid = hit.get("case_id", "")
    diag = hit.get("diagnoses", [{}])[0] if hit.get("diagnoses") else {}
    demo = hit.get("demographic", {})
    samp = hit.get("samples", [{}])[0] if hit.get("samples") else {}
    tx   = diag.get("treatments", [{}])[0] if diag.get("treatments") else {}

    case_id_to_info[cid] = {
        "case_id":           cid,
        "project_id":        hit.get("project", {}).get("project_id", pd.NA),
        "primary_diagnosis": diag.get("primary_diagnosis", pd.NA),
        "tumor_grade":       diag.get("tumor_grade", pd.NA),
        "prior_malignancy":  diag.get("prior_malignancy", pd.NA),
        "tissue_type":       samp.get("tissue_type", pd.NA),
        "treatment_outcome": tx.get("treatment_outcome", pd.NA),
        "gender":            demo.get("gender", pd.NA),
        "race":              demo.get("race", pd.NA),
        "age_at_index":      demo.get("age_at_index", pd.NA),
    }

rows = []
for file_uuid in samples:
    case_id = file_to_case.get(file_uuid, pd.NA)
    info    = case_id_to_info.get(case_id, {}) if case_id is not pd.NA else {}
    row     = {"file_uuid": file_uuid, "case_id": case_id}
    row.update(info)
    rows.append(row)

    if not info:
        log_transformation(log, "gdc_clinical_fetch",
            f"WARNING: no clinical data for file UUID {file_uuid}", status="WARNING")

clinical_df = pd.DataFrame(rows)
clinical_df.to_csv(snakemake.output.clinical, sep="\t", index=False)

n_missing = clinical_df["primary_diagnosis"].isna().sum()
log_transformation(log, "gdc_clinical_fetch",
    f"Wrote {len(clinical_df)} rows; {n_missing} samples lack primary_diagnosis")

# ---------------------------------------------------------------------------
# FAIR provenance
# ---------------------------------------------------------------------------
prov = {
    "rule":              "gdc_clinical_fetch",
    "api_endpoint":      cases_url,
    "query_timestamp":   datetime.now(timezone.utc).isoformat(),
    "n_samples":         len(samples),
    "n_cases_resolved":  len(file_to_case),
    "n_cases_returned":  len(cases_hits),
    "fields_requested":  fields,
    "output":            str(Path(snakemake.output.clinical).resolve()),
}
Path(snakemake.output.provenance).write_text(json.dumps(prov, indent=2))

log_transformation(log, "gdc_clinical_fetch", "Complete", status="SUCCESS",
                   artifact_paths=[snakemake.output.clinical])
