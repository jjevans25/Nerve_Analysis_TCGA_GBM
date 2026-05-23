"""Out-of-band launcher for the scANVI re-train (Step 3 of the v2 plan).

The Snakemake DAG is currently blocked by an unrelated upstream env issue in
`scrna_integration` (torchvision circular-import in its conda env). This
launcher invokes `workflow/scripts/nerve_scanvi_retrain.py` directly with a
Snakemake-like namespace, so the v2 work isn't gated on fixing the upstream
env. Same parameters as the Snakemake rule.

Launch from project root:

    .snakemake/conda/ef772b6aee7a5ae3b33d107abc89938f_/bin/python \\
        scripts/run_scanvi_retrain.py

Expected runtime: 1.5–3 hours on M4 Max MPS. Outputs:
- results/models/nerve_scanvi_model/
- data/processed/nerve_cells_v2.h5ad
- results/figures/nerve_scanvi_training_curves.png
- provenance/nerve_scanvi_retrain_provenance.json
"""

from __future__ import annotations

import os
import runpy
import types
from pathlib import Path

import yaml

PROJ = Path(__file__).resolve().parent.parent
os.chdir(str(PROJ))

with open(PROJ / "config/config.yaml") as f:
    config = yaml.safe_load(f)

snakemake = types.SimpleNamespace()
snakemake.input = types.SimpleNamespace(
    counts_h5ad="data/processed/nerve_cells_counts_labeled.h5ad",
)
snakemake.output = types.SimpleNamespace(
    model_dir="results/models/nerve_scanvi_model",
    latent_h5ad="data/processed/nerve_cells_v2.h5ad",
    curves="results/figures/nerve_scanvi_training_curves.png",
    provenance="provenance/nerve_scanvi_retrain_provenance.json",
)
snakemake.log = ["logs/nerve_scanvi_retrain.log"]
snakemake.params = types.SimpleNamespace(
    device=config["hardware"]["device"],
    precision=config["hardware"]["precision"],
    n_latent=config["scrna"]["n_latent"],
    n_layers=config["scrna"]["n_layers"],
    random_seed=config["scrna"]["random_seed"],
    batch_key=config["nerve_scanvi"]["batch_key"],
    labels_key=config["nerve_scanvi"]["labels_key"],
    unlabeled_category=config["nerve_scanvi"]["unlabeled_category"],
    scvi_max_epochs=config["nerve_scanvi"]["scvi_max_epochs"],
    scanvi_max_epochs=config["nerve_scanvi"]["scanvi_max_epochs"],
    n_samples_per_label=config["nerve_scanvi"]["n_samples_per_label"],
)

# Ensure output dirs exist
for path in [
    "results/models",
    "results/figures",
    "data/processed",
    "provenance",
    "logs",
]:
    Path(path).mkdir(parents=True, exist_ok=True)

runpy.run_path(
    "workflow/scripts/nerve_scanvi_retrain.py",
    init_globals={"snakemake": snakemake},
)
