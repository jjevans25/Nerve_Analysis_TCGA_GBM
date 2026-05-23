"""Fetch clinical metadata from the GDC public REST API for all sample file UUIDs."""

import json
import sys
import time
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
timeout_s   = int(snakemake.params.request_timeout)
max_retries = int(snakemake.params.max_retries)


def _post_with_retry(url: str, payload: dict) -> requests.Response:
    """POST with exponential backoff on transient (5xx, network, timeout) failures.

    Mirrors ``download_msigdb_gmt.py:_download_with_retry``. 4xx is *not* retried
    (client error — retrying won't help). Sleep grows 1s, 2s, 4s, 8s, 16s, ...
    """
    last_exc: BaseException | None = None
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.post(url, json=payload, timeout=timeout_s)
            if 400 <= resp.status_code < 500:
                resp.raise_for_status()  # client error -> raise immediately, no retry
            if resp.status_code >= 500:
                log_transformation(
                    log, "gdc_clinical_fetch",
                    f"HTTP {resp.status_code} on POST {url} "
                    f"(attempt {attempt}/{max_retries})",
                    status="WARNING",
                )
                resp.raise_for_status()
            return resp
        except (requests.RequestException, OSError) as exc:
            # Distinguish client (don't retry) from server / network (retry).
            status = getattr(getattr(exc, "response", None), "status_code", None)
            if status is not None and 400 <= status < 500:
                raise
            last_exc = exc
            if attempt < max_retries:
                wait = 2 ** (attempt - 1)
                log_transformation(
                    log, "gdc_clinical_fetch",
                    f"  retry {attempt}/{max_retries} after {wait}s: {exc}",
                    status="WARNING",
                )
                time.sleep(wait)
            else:
                log_transformation(
                    log, "gdc_clinical_fetch",
                    f"All {max_retries} attempts failed for {url}",
                    status="ERROR",
                )
    raise RuntimeError(
        f"[FAIR-ALERT] gdc_clinical_fetch: POST {url} failed after "
        f"{max_retries} attempts; last error: {last_exc}"
    )


log_transformation(log, "gdc_clinical_fetch",
    f"Querying GDC API for {len(samples)} file UUIDs at {base} "
    f"(timeout={timeout_s}s, max_retries={max_retries})")

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

resp = _post_with_retry(files_url, files_payload)
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

resp2 = _post_with_retry(cases_url, cases_payload)
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
