#!/usr/bin/env python3
"""
scVI / scANVI Annotation — Cisplatin_GC_Vehicle (4 conditions)

Author: Xin Wang
Date:   2026-04-14

Runs scVI/scANVI label transfer onto the 4-condition query dataset
(Vehicle, Vehicle-GC, Cisplatin, Cisplatin-GC) using one reference atlas:

  - Lake 2025 snRNA-seq: mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2.h5ad

Pipeline:
  1. Load the Lake reference h5ad; harmonise gene names and cell-type labels
  2. Load query h5ad (Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.h5ad)
  3. Subset to shared genes with query
  4. Train scVI on the reference
  5. Merge reference + query; train scANVI (semi-supervised)
  6. Predict cell types on query cells
  7. Save CSV and diagnostic plots

Prerequisites: Run Cisplatin_GC_Vehicle_Annotation.Rmd first to generate
               the query h5ad.
"""

import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib
matplotlib.use("Agg")          # non-interactive backend for HPC
import matplotlib.pyplot as plt
import torch

# ── Reproducibility ─────────────────────────────────────────────────────────
scvi.settings.seed = 0
sc.settings.verbosity = 1

print("scvi-tools version:", scvi.__version__)
print("GPU available:", torch.cuda.is_available())

# ── Paths ────────────────────────────────────────────────────────────────────
ref_dir   = "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Datasets/Reference/"
ann_dir   = "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/"

lake_path  = os.path.join(ref_dir, "mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2.h5ad")
query_path = os.path.join(ann_dir, "Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.h5ad")
out_dir    = ann_dir
os.makedirs(out_dir, exist_ok=True)

log_dir = os.path.join(ann_dir, "logs")
os.makedirs(log_dir, exist_ok=True)

for p in [lake_path, query_path]:
    print(f"{os.path.basename(p):60s}  exists={os.path.exists(p)}")

# ── Load reference dataset ───────────────────────────────────────────────────
print("\n=== Loading Lake reference ===")
adata_lake = sc.read_h5ad(lake_path)
print(adata_lake)
print(adata_lake.obs.columns.tolist())

# ── Cell-type label column ───────────────────────────────────────────────────
lake_label_col = "SubclassLevel1"     # update if needed

print("\nLake cell types:")
print(adata_lake.obs[lake_label_col].value_counts())

# ── Harmonise gene names (Ensembl → symbol if needed) ───────────────────────
def ensembl_to_symbol(adata):
    """Convert Ensembl IDs in var_names to gene symbols via adata.var['feature_name']."""
    if not adata.var_names[0].startswith("ENSMUS"):
        print("  var_names already look like gene symbols — skipping conversion.")
        return adata
    if "feature_name" not in adata.var.columns:
        raise ValueError("No 'feature_name' column found in adata.var.")
    symbols = adata.var["feature_name"].astype(str).values.copy()

    # Strip trailing _ENSMUSG... suffix
    symbols = np.array([re.sub(r'_ENSMUSG\d+$', '', s) for s in symbols])

    fallback_mask = (symbols == "nan") | (symbols == "") | (symbols == "None")
    n_fallback = fallback_mask.sum()
    if n_fallback:
        print(f"  WARNING: {n_fallback} genes have no feature_name — keeping Ensembl ID.")
        symbols[fallback_mask] = adata.var_names[fallback_mask]

    adata.var_names = symbols
    adata.var_names_make_unique()
    print(f"  Converted {adata.n_vars} var_names to gene symbols ({n_fallback} kept as Ensembl IDs).")
    return adata

print("\nHarmonising gene names...")
adata_lake = ensembl_to_symbol(adata_lake)
print("Lake var sample:", adata_lake.var_names[:5].tolist())

# Sanity check: Lake labels populated
n_lake_labeled = adata_lake.obs[lake_label_col].notna().sum()
print(f"\nLake obs with '{lake_label_col}' populated: {n_lake_labeled} / {adata_lake.n_obs}")
if n_lake_labeled == 0:
    raise ValueError(
        f"Column '{lake_label_col}' is entirely empty in adata_lake. "
        "Check the correct column name."
    )

# ── Load query dataset ───────────────────────────────────────────────────────
print("\n=== Loading query dataset ===")
adata_query = sc.read_h5ad(query_path)
print(adata_query)
print("Conditions:", adata_query.obs["DataSet"].value_counts().to_dict())

# ── scVI / scANVI annotation using Lake reference ────────────────────────────
LABEL_KEY = "cell_type"
ref_name  = "Lake"

print(f"\n{'='*60}")
print(f"  {ref_name} reference")
print(f"{'='*60}")

# Unified label column
adata_lake.obs[LABEL_KEY] = adata_lake.obs[lake_label_col]

# Shared genes
shared = adata_lake.var_names.intersection(adata_query.var_names)
print(f"  Shared genes (Lake ∩ query): {len(shared)}")

ref_sub   = adata_lake[:, shared].copy()
query_sub = adata_query[:, shared].copy()

# ── scVI on reference ────────────────────────────────────────────────────────
scvi.model.SCVI.setup_anndata(ref_sub, labels_key=LABEL_KEY)
scvi_model = scvi.model.SCVI(ref_sub)
scvi_model.train()
scvi_model.save(os.path.join(out_dir, f"scvi_model_{ref_name}"), overwrite=True)
print(f"  scVI model saved.")

# Reference UMAP
ref_sub.obsm["X_scVI"] = scvi_model.get_latent_representation()
sc.pp.neighbors(ref_sub, use_rep="X_scVI")
sc.tl.umap(ref_sub)
fig, ax = plt.subplots(figsize=(9, 6))
sc.pl.umap(ref_sub, color=LABEL_KEY, title=f"{ref_name} reference — cell type", ax=ax, show=False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"{ref_name}_scVI_UMAP.pdf"), bbox_inches="tight")
plt.close()

# ── Prepare query for scANVI ─────────────────────────────────────────────────
scvi.model.SCVI.prepare_query_anndata(query_sub, scvi_model)
query_sub.obs[LABEL_KEY] = pd.Categorical(["Unknown"] * query_sub.n_obs)

# ── Merge ref + query ────────────────────────────────────────────────────────
combined = ref_sub.concatenate(
    query_sub,
    batch_key="batch",
    batch_categories=["ref", "query"],
)

# ── scANVI ───────────────────────────────────────────────────────────────────
scvi.model.SCANVI.setup_anndata(
    combined,
    labels_key=LABEL_KEY,
    batch_key="batch",
    unlabeled_category="Unknown",
)
scanvi_model = scvi.model.SCANVI(combined)
scanvi_model.train(max_epochs=50)
scanvi_model.save(os.path.join(out_dir, f"scanvi_model_{ref_name}"), overwrite=True)
print(f"  scANVI model saved.")

# ── Predict on query ─────────────────────────────────────────────────────────
query_mask  = combined.obs["batch"] == "query"
soft_preds  = scanvi_model.predict(combined[query_mask], soft=True)
label_names = scanvi_model.adata_manager.get_state_registry("labels")["categorical_mapping"]

prob_df     = pd.DataFrame(soft_preds, index=combined.obs_names[query_mask], columns=label_names)
pred_labels = prob_df.idxmax(axis=1)
pred_conf   = prob_df.max(axis=1)

combined.obs.loc[query_mask, f"scanvi_{ref_name}_label"]      = pred_labels.values
combined.obs.loc[query_mask, f"scanvi_{ref_name}_confidence"] = pred_conf.values

# Confidence histogram
pred_conf.hist(bins=50, figsize=(7, 4))
plt.title(f"scANVI confidence — {ref_name} reference")
plt.xlabel("Confidence")
plt.ylabel("Cell count")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"scanvi_{ref_name}_confidence.pdf"))
plt.close()

# Combined UMAP
combined.obsm["X_scANVI"] = scanvi_model.get_latent_representation()
sc.pp.neighbors(combined, use_rep="X_scANVI")
sc.tl.umap(combined)
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
sc.pl.umap(combined, color=LABEL_KEY, title=f"{ref_name} — cell type",   ax=axes[0], show=False)
sc.pl.umap(combined, color="batch",   title=f"{ref_name} — ref / query", ax=axes[1], show=False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"{ref_name}_combined_scANVI_UMAP.pdf"), bbox_inches="tight")
plt.close()

# ── Save predictions ─────────────────────────────────────────────────────────
preds = combined.obs.loc[
    query_mask,
    [f"scanvi_{ref_name}_label", f"scanvi_{ref_name}_confidence"],
].copy()

extra_cols = [
    c for c in adata_query.obs.columns
    if c in ("predicted.id", "prediction.score.max", "Raw_cell_type", "Raw_cell_type_Confidence")
]
if extra_cols:
    preds = preds.join(adata_query.obs[extra_cols], how="left")

csv_path = os.path.join(out_dir, "Cisplatin_GC_Vehicle_scanvi_Lake_annotations.csv")
preds.to_csv(csv_path)
print("Predictions saved to:", csv_path)
print(preds.head())
print(f"\n{ref_name} done.")
print("\nDone.")
