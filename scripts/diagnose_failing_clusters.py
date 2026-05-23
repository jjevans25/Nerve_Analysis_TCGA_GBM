"""Diagnostic profile for the 6 batch-QC-failing nerve clusters.

Outputs structured per-cluster evidence (QC, markers, dominant-patient
clinical) to stdout for hand-synthesis into markdowns/failing_cluster_diagnosis.md.

Not a Snakemake rule — one-shot remediation-decision step.
"""

import json
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

FAILING = ["13", "15", "19", "21", "22", "23"]
# size-matched passing controls (242-2306 cells in failing → pick passing
# clusters in the same range): cl17 (1168), cl18 (1159), cl20 (868), cl12 (2480)
CONTROLS = ["12", "17", "18", "20"]

# Gene-set classifiers (uppercased symbols / glob-style prefixes).
GENESETS = {
    "stress_dissociation": {"FOS", "JUN", "JUNB", "JUND", "EGR1", "EGR2",
                            "ATF3", "DUSP1", "IER2", "IER3", "BTG2", "ZFP36"},
    "heatshock":           {"HSPA1A", "HSPA1B", "HSP90AA1", "HSPB1", "DNAJB1"},
    "mito_ribo":           None,  # prefix-matched: MT-, MTRNR, RPL, RPS
    "ambient_blood":       {"HBB", "HBA1", "HBA2", "ALB", "TTR", "HBD"},
    "cellcycle":           {"MKI67", "TOP2A", "CDK1", "AURKB", "CCNB1", "CCNB2",
                            "PCNA", "MCM2", "MCM4", "STMN1", "BIRC5"},
    "excitatory_neuron":   {"SLC17A7", "SLC17A6", "NRGN", "RBFOX3", "SYT1",
                            "SNAP25", "CAMK2A", "GRIN2B", "GRIA1"},
    "inhibitory_neuron":   {"GAD1", "GAD2", "PVALB", "SST", "VIP", "RELN",
                            "LHX6", "DLX1", "DLX5", "ADARB2"},
    "cortical_layer":      {"RORB", "FOXP2", "BCL11B", "CUX2", "SATB2",
                            "FEZF2", "TBR1", "POU3F2"},
    "synaptic_adhesion":   {"CNTNAP2", "CNTNAP4", "CNTNAP5", "NRXN1", "NRXN3",
                            "NLGN1", "NLGN3", "DSCAM", "TENM2", "TENM3"},
    "ependymal_ciliated":  {"FOXJ1", "RSPH1", "PIFO", "TMEM231"},
    "ciliated_dynein":     None,  # prefix-matched: DNAH, CFAP, DNAI
    "reactive_astro_MES":  {"CHI3L1", "CD44", "SERPINA3", "VIM", "GFAP", "S100B"},
    "reactive_astro_A1":   {"LCN2", "GBP2", "C3"},
    "reactive_astro_A2":   {"GLUL", "S100A6", "EMP1", "S100A10", "PTX3"},
    "opc_olig_lineage":    {"PDGFRA", "CSPG4", "SOX10", "OLIG1", "OLIG2",
                            "MBP", "MOG", "PLP1", "MAG"},
}


def classify_gene(g: str) -> list[str]:
    """Return all gene-set labels matched by gene symbol g."""
    g = g.upper()
    hits = []
    # prefix-matched
    if g.startswith(("MT-", "MTRNR")) or re.match(r"^RPL\d", g) or re.match(r"^RPS\d", g):
        hits.append("mito_ribo")
    if g.startswith(("DNAH", "CFAP", "DNAI")):
        hits.append("ciliated_dynein")
    # set-matched
    for label, members in GENESETS.items():
        if members is not None and g in members:
            hits.append(label)
    return hits


# ---------- load -------------------------------------------------------------
print("=" * 76)
print("Loading inputs")
print("=" * 76)

ad_path = ROOT / "data/processed/nerve_cells.h5ad"
adata = ad.read_h5ad(ad_path, backed="r")

# Materialise the needed obs columns into a small DataFrame so we don't drag
# the whole AnnData into RAM.
OBS_COLS = ["nerve_leiden", "sample_id", "pct_counts_mt", "n_genes_by_counts",
            "total_counts", "n_genes", "cnv_score", "cell_type_confidence",
            "cell_type_predicted"]
obs = adata.obs[OBS_COLS].copy()
obs["nerve_leiden"] = obs["nerve_leiden"].astype(str)
obs["sample_id"]    = obs["sample_id"].astype(str)
del adata  # close backed handle

markers = pd.read_csv(ROOT / "results/tables/nerve_cluster_markers_with_qc.csv")
markers["cluster"] = markers["cluster"].astype(str)

annot = pd.read_csv(ROOT / "results/tables/nerve_cluster_annotations.csv")
annot["cluster"] = annot["cluster"].astype(str)

purity = pd.read_csv(ROOT / "results/tables/nerve_cluster_sample_purity.csv")
purity["cluster"] = purity["cluster"].astype(str)

clinical = pd.read_csv(ROOT / "data/external/gdc_clinical.tsv", sep="\t")
print(f"  obs        {obs.shape}")
print(f"  markers    {markers.shape}")
print(f"  annot      {annot.shape}")
print(f"  purity     {purity.shape}")
print(f"  clinical   {clinical.shape} cols={list(clinical.columns)}")
print()

# Map file_uuid (clinical) → case_id; build sample_id → row lookup.
# Per-cell sample_id values are the file_uuid (matching the data/raw filenames).
clinical_lookup = clinical.set_index("file_uuid")

# ---------- per-cluster QC profile -------------------------------------------
QC_COLS = ["pct_counts_mt", "n_genes_by_counts", "total_counts", "cnv_score",
           "cell_type_confidence"]


def cluster_qc(cl: str) -> dict:
    sub = obs[obs["nerve_leiden"] == cl]
    summary = {"n_cells": int(len(sub))}
    for c in QC_COLS:
        v = sub[c].to_numpy()
        summary[c + "_median"] = float(np.nanmedian(v))
        summary[c + "_p90"]    = float(np.nanpercentile(v, 90))
    ct_counts = sub["cell_type_predicted"].astype(str).value_counts(normalize=True).round(3)
    summary["cell_type"] = ct_counts.head(3).to_dict()
    return summary


print("=" * 76)
print("Per-cluster QC profiles (median / p90)")
print("=" * 76)
print(f"{'cluster':>7}  {'n':>6}  {'pct_mt':>9}  {'n_genes':>9}  {'total_cnt':>10}  "
      f"{'cnv':>7}  {'ct_conf':>8}")
print("-" * 76)


def fmt_qc_row(cl: str, kind: str, qc: dict) -> None:
    print(f"{cl:>4} {kind:>2}  {qc['n_cells']:>6}  "
          f"{qc['pct_counts_mt_median']:>5.2f}/{qc['pct_counts_mt_p90']:>5.2f}  "
          f"{qc['n_genes_by_counts_median']:>5.0f}/{qc['n_genes_by_counts_p90']:>5.0f}  "
          f"{qc['total_counts_median']:>5.0f}/{qc['total_counts_p90']:>5.0f}  "
          f"{qc['cnv_score_median']:>6.3f}  "
          f"{qc['cell_type_confidence_median']:>6.3f}")


qc_results = {}
for cl in FAILING:
    qc = cluster_qc(cl)
    qc_results[cl] = qc
    fmt_qc_row(cl, "F", qc)
for cl in CONTROLS:
    qc = cluster_qc(cl)
    qc_results[cl] = qc
    fmt_qc_row(cl, "C", qc)

# Cohort-wide median for context
cohort = {}
for c in QC_COLS:
    v = obs[c].to_numpy()
    cohort[c] = float(np.nanmedian(v))
print(f"\ncohort median: pct_mt={cohort['pct_counts_mt']:.2f}  "
      f"n_genes={cohort['n_genes_by_counts']:.0f}  "
      f"total_counts={cohort['total_counts']:.0f}  "
      f"cnv={cohort['cnv_score']:.3f}  ct_conf={cohort['cell_type_confidence']:.3f}")
print()

# ---------- marker biological coherence --------------------------------------
print("=" * 76)
print("Top-10 markers per failing cluster + gene-set classification")
print("=" * 76)
marker_classifications = {}
for cl in FAILING:
    sub = markers[markers["cluster"] == cl].copy()
    # Take top 10 by adjusted p-value (genes are typically ranked by score; defensive sort)
    if "scores" in sub.columns:
        sub = sub.sort_values("scores", ascending=False)
    top = sub.head(10)
    gene_col = "gene_symbol" if "gene_symbol" in top.columns else "names"
    top_genes = top[gene_col].tolist()
    labels_per_gene = {g: classify_gene(str(g)) for g in top_genes}
    # aggregate
    label_counts = {}
    for hits in labels_per_gene.values():
        for h in hits:
            label_counts[h] = label_counts.get(h, 0) + 1
    marker_classifications[cl] = {
        "top_genes": top_genes,
        "per_gene_labels": labels_per_gene,
        "label_counts": label_counts,
    }
    print(f"\ncluster {cl} top markers:")
    for g in top_genes:
        labs = labels_per_gene[g]
        print(f"   {g:<14}  -> {labs if labs else '(unclassified)'}")
    print(f"   label totals: {label_counts}")

# ---------- dominant-patient clinical cross-reference ------------------------
print("\n" + "=" * 76)
print("Dominant-patient clinical metadata for failing clusters")
print("=" * 76)
dom_summary = {}
for cl in FAILING:
    row = purity[purity["cluster"] == cl].iloc[0]
    dom = row["dominant_sample"]
    frac = float(row["dominant_sample_fraction"])
    n_contrib = int(row["n_contributing_samples"])
    clin = clinical_lookup.loc[dom] if dom in clinical_lookup.index else None
    dom_summary[cl] = {
        "dominant_sample": dom,
        "dominant_fraction": frac,
        "n_contributing_samples": n_contrib,
        "clinical": dict(clin) if clin is not None else None,
    }
    print(f"\ncluster {cl}: dom={dom[:8]} ({frac:.2f}) n_contrib={n_contrib}")
    if clin is not None:
        for k in ["case_id", "primary_diagnosis", "tumor_grade", "gender",
                  "age_at_index", "tissue_type"]:
            if k in clin.index:
                print(f"   {k}: {clin[k]}")
    else:
        print("   (not found in gdc_clinical.tsv)")

# Cohort distribution of relevant clinical columns
print("\n--- cohort clinical distribution ---")
for k in ["primary_diagnosis", "tumor_grade", "gender", "tissue_type"]:
    if k in clinical.columns:
        print(f"{k}: {clinical[k].value_counts(dropna=False).to_dict()}")
print(f"age_at_index: median={clinical['age_at_index'].median():.0f} "
      f"min={clinical['age_at_index'].min()} max={clinical['age_at_index'].max()}")

# ---------- canonical scores (from existing annot) ---------------------------
print("\n" + "=" * 76)
print("Canonical cell-type marker scores (from nerve_cluster_annotations.csv)")
print("=" * 76)
score_cols = [c for c in annot.columns if c.startswith("score_")]
print("cluster | " + " | ".join(s.replace("score_", "")[:5] for s in score_cols))
for cl in FAILING:
    row = annot[annot["cluster"] == cl].iloc[0]
    vals = [f"{row[c]:+.2f}" for c in score_cols]
    print(f"   {cl:>3}  | " + " | ".join(vals))
for cl in CONTROLS:
    row = annot[annot["cluster"] == cl].iloc[0]
    vals = [f"{row[c]:+.2f}" for c in score_cols]
    print(f"   {cl:>3}c | " + " | ".join(vals))

# ---------- save structured JSON for the report writer -----------------------
out = {
    "failing": FAILING,
    "controls": CONTROLS,
    "qc": qc_results,
    "markers": marker_classifications,
    "dominant": dom_summary,
    "cohort_medians": cohort,
}
out_path = ROOT / "results/tables/failing_cluster_diagnosis.json"
out_path.write_text(json.dumps(out, indent=2, default=str))
print(f"\nWrote structured evidence: {out_path}")
