"""Print all provenance records currently registered in the provenance directory."""

import json
from pathlib import Path

provenance_dir = Path(snakemake.params.provenance_dir)
records = sorted(provenance_dir.glob("*.json"))

if not records:
    print("[INFO] No provenance records found yet.")

for r in records:
    data = json.loads(r.read_text())
    print(f"  {data.get('snakemake_rule','?'):30s}  {data.get('created_at','?')}  {r.name}")
