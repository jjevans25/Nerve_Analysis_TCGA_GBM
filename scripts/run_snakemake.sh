#!/usr/bin/env bash
# Invoke Snakemake so that `--use-conda` actually enforces workflow/envs/*.yaml.
#
# THE DEFECT (markdowns/task_conda_env_enforcement.md, Defect 1)
# --------------------------------------------------------------
# `claude_science` is a venv whose `bin/python` is a symlink into /opt/anaconda3.
# When it is ACTIVE in the calling shell, `VIRTUAL_ENV` is exported and
# `claude_science/bin` sits first on PATH. Snakemake builds each `script:` job as
#
#     source /opt/anaconda3/bin/activate <env>; set -euo pipefail; python <job>.py
#
# and that child shell inherits the caller's environment. `conda activate`
# prepends the env's bin, but the still-exported VIRTUAL_ENV keeps resolving
# `python` back to the venv — so every rule ran against the venv's packages while
# the run *reported* the declared conda env. The pins in `workflow/envs/*.yaml`
# were never enforced, which is the standing [FAIR-ALERT] on every artifact in
# the Census arms.
#
# THE FIX
# -------
# Unset VIRTUAL_ENV and strip the venv from PATH before starting Snakemake, then
# call the venv's snakemake by absolute path (its shebang still points at the
# venv interpreter, so Snakemake itself runs unchanged). Job subshells inherit
# the cleaned environment, conda activation resolves normally, and the declared
# env wins.
#
# Verified 2026-08-05 against the acceptance criteria in the task doc:
#   sys.executable -> $CONDA_ENV/bin/python
#   igraph 0.11.8, torch 2.12.0, infercnvpy 0.4.3   (all from the conda env)
#   .snakemake/conda/ hashes unchanged — no environment rebuild
#
# Neither candidate approach in the task doc was needed: Snakemake was not
# reinstalled and the venv was not rebuilt.
#
# Usage:  scripts/run_snakemake.sh <targets...> [snakemake flags...]
#         Targets must come FIRST — several Snakemake flags take nargs='+' and
#         will otherwise swallow them (--allowed-rules, --quiet, ...).

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_BIN="${PROJECT_ROOT}/claude_science/bin"
SNAKEMAKE_BIN="${VENV_BIN}/snakemake"

if [[ ! -x "${SNAKEMAKE_BIN}" ]]; then
    echo "[FAIR-ALERT] snakemake not found at ${SNAKEMAKE_BIN}" >&2
    exit 1
fi

# Strip the venv from PATH so no child resolves `python` to it.
CLEAN_PATH="$(printf '%s' "${PATH}" | tr ':' '\n' | grep -v -F "${VENV_BIN}" | paste -sd: -)"

# numba's OpenMP pool SIGSEGVs against torch's libomp once both are loaded; this
# must be set before any import that pulls in numba (scanpy -> umap-learn).
export NUMBA_THREADING_LAYER="${NUMBA_THREADING_LAYER:-workqueue}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"
export PYTORCH_MPS_HIGH_WATERMARK_RATIO="${PYTORCH_MPS_HIGH_WATERMARK_RATIO:-0.0}"

# Only VIRTUAL_ENV and the venv's PATH entry are removed. CONDA_PREFIX and
# CONDA_DEFAULT_ENV must be left alone: `conda activate` reads CONDA_PREFIX to
# find the *current* env's deactivate scripts, and unsetting it makes
# activate.py raise inside posixpath.join before the job ever starts.
cd "${PROJECT_ROOT}"
exec env -u VIRTUAL_ENV PATH="${CLEAN_PATH}" "${SNAKEMAKE_BIN}" --use-conda "$@"
