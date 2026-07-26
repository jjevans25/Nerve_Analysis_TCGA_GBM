"""Export a Marimo notebook to HTML and write a FAIR provenance record."""

import json
import os
import subprocess
import uuid
from datetime import datetime
from pathlib import Path


def main() -> None:
    notebook = snakemake.input.notebook  # noqa: F821
    html_out = snakemake.output.html  # noqa: F821
    prov_out = snakemake.output.provenance  # noqa: F821
    log_path = snakemake.log[0]  # noqa: F821

    # Optional `params.env` lets a cohort-namespaced rule tell the notebook which
    # dataset to read (marimo has no argv passthrough). Additive: rules that do
    # not set it are unaffected.
    env = dict(os.environ)
    env.update({k: str(v) for k, v in getattr(snakemake.params, "env", {}).items()})  # noqa: F821

    result = subprocess.run(
        ["marimo", "export", "html", notebook, "-o", html_out],
        capture_output=True,
        text=True,
        env=env,
    )

    Path(log_path).write_text(result.stdout + result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"marimo export failed — see {log_path}")

    input_hash = subprocess.run(
        ["md5", "-q", notebook], capture_output=True, text=True
    ).stdout.strip()

    prov = {
        "artifact":   html_out,
        "run_id":     str(uuid.uuid4()),
        "timestamp":  datetime.utcnow().isoformat() + "Z",
        "tool":       "marimo==0.23.1",
        "rule":       snakemake.rule,  # noqa: F821
        "input_hash": input_hash,
    }
    Path(prov_out).write_text(json.dumps(prov, indent=2))


main()
