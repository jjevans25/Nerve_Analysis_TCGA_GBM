"""Assert that `--use-conda` actually enforces the declared environment.

Not a scientific rule — an environment guard. Every artifact in the Census arms
carries a standing [FAIR-ALERT] because the `claude_science` venv shadowed conda
activation, so `workflow/envs/scrna.yaml` pins were declared but never enforced
and each rule silently ran against whatever the venv happened to contain (see
markdowns/task_conda_env_enforcement.md, Defect 1).

Reproducibility claims are only as good as the environment behind them, so this
runs as a rule under the same mechanism as every other rule and fails loudly if
the interpreter is not the one the workflow declared.
"""

import json
import os
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

sys.path.insert(0, "workflow/scripts")
from fair_utils import log_transformation, stamp_artifact, verify_artifact, write_provenance

log = snakemake.log[0]
expected: dict[str, str] = dict(snakemake.params.expected_versions)

checks: list[dict[str, object]] = []


def record(name: str, observed: object, ok: bool, detail: str = "") -> None:
    """Append one pass/fail row; the rule fails if any row is False."""
    checks.append({"check": name, "observed": str(observed), "pass": bool(ok), "detail": detail})


# --- 1. The interpreter must live inside a conda env, not the venv -----------
exe = Path(sys.executable).resolve()
in_conda = "/.snakemake/conda/" in str(exe) or bool(os.environ.get("CONDA_PREFIX"))
in_venv = "claude_science" in str(exe)
record("interpreter_is_conda_env", exe, in_conda and not in_venv,
       "sys.executable must be $CONDA_ENV/bin/python, not the claude_science venv")

# --- 2. VIRTUAL_ENV must not be leaking into the job shell -------------------
record("virtual_env_unset", os.environ.get("VIRTUAL_ENV", "<unset>"),
       not os.environ.get("VIRTUAL_ENV"),
       "an exported VIRTUAL_ENV re-shadows conda activation for every child")

# --- 3. site-packages must resolve inside the env ----------------------------
env_root = str(exe.parent.parent)
sp_in_env = [p for p in sys.path if p.startswith(env_root) and "site-packages" in p]
record("site_packages_in_env", sp_in_env[:1] or "<none>", bool(sp_in_env))

# --- 4. Declared pins must be the versions actually imported -----------------
# igraph and torch are the canaries from the task doc: the venv carries its own
# igraph, and anaconda base a third copy, so agreement here means the whole env
# is enforced rather than just the packages an earlier hotfix happened to add.
for pkg, want in expected.items():
    try:
        got = version(pkg)
    except PackageNotFoundError:
        record(f"version::{pkg}", "<not installed>", False, f"declared pin {want}")
        continue
    record(f"version::{pkg}", got, got == want, f"declared pin {want}")

# --- 5. The packages the workflow imports must resolve inside the env --------
for mod in ("igraph", "scanpy", "anndata", "torch", "liana"):
    try:
        imported = __import__(mod)
        loc = getattr(imported, "__file__", "") or ""
        record(f"import_location::{mod}", loc, loc.startswith(env_root),
               "must load from the conda env, not the venv or anaconda base")
    except ImportError as exc:
        record(f"import_location::{mod}", f"<ImportError: {exc}>", False)

failed = [c for c in checks if not c["pass"]]
report = {
    "checked_at_utc": __import__("datetime").datetime.now(
        __import__("datetime").timezone.utc).isoformat(),
    "sys_executable": str(exe),
    "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
    "virtual_env": os.environ.get("VIRTUAL_ENV", ""),
    "pass": not failed,
    "n_checks": len(checks),
    "n_failed": len(failed),
    "checks": checks,
}

out = Path(snakemake.output.report)
out.parent.mkdir(parents=True, exist_ok=True)
with open(out, "w") as fh:
    json.dump(report, fh, indent=2)
verify_artifact(out, min_size_bytes=32)

for c in checks:
    log_transformation(log, "conda_env_smoke_test",
                       f"{'PASS' if c['pass'] else 'FAIL'}  {c['check']} = {c['observed']}",
                       status="SUCCESS" if c["pass"] else "ERROR")

prov = stamp_artifact(
    output_path=out,
    rule_name="conda_env_smoke_test",
    input_paths=[],
    tool_versions={"python": sys.version.split()[0]},
    parameters={"expected_versions": expected, "n_failed": len(failed)},
    description="Environment enforcement guard: --use-conda resolves to the declared env",
    ontology_operation="operation:2428",  # EDAM: Validation
)
write_provenance(prov, snakemake.output.provenance)

if failed:
    raise RuntimeError(
        "[FAIR-ALERT] the declared conda environment is NOT being enforced:\n"
        + "\n".join(f"  - {c['check']} = {c['observed']} ({c['detail']})" for c in failed)
        + "\n\nRun the pipeline via scripts/run_snakemake.sh, which unsets VIRTUAL_ENV "
          "and strips the venv from PATH. See markdowns/task_conda_env_enforcement.md."
    )

log_transformation(log, "conda_env_smoke_test",
                   f"All {len(checks)} environment checks passed", status="SUCCESS",
                   artifact_paths=[str(out)])
