#!/usr/bin/env python3
"""
scVI / scANVI Annotation — PT Subclusters (Lake PT reference)

Author: Xin Wang
Date:   2026-06-30

Runs scVI/scANVI label transfer onto the PT subcluster query dataset
using PT cells from the Lake 2025 snRNA-seq atlas as reference:

  - Lake 2025 snRNA-seq: mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2.h5ad
    (subsetted to SubclassLevel1 == "PT"; labels = SubclassLevel2)

Pipeline:
  1. Load the Lake reference h5ad; subset to PT cells; harmonise gene names
  2. Load query h5ad (PT_Subclustered_Annotated.h5ad)
  3. Subset to shared genes
  4. Train scVI on the Lake PT reference
  5. Merge reference + query; train scANVI (semi-supervised)
  6. Predict PT subtypes on query cells
  7. Save CSV and diagnostic plots

Prerequisites: Run PT_Export_h5ad.R first to generate the query h5ad.
"""

import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import torch

# ── Reproducibility ─────────────────────────────────────────────────────────
scvi.settings.seed = 0
sc.settings.verbosity = 1

print("scvi-tools version:", scvi.__version__)
print("GPU available:", torch.cuda.is_available())

# ── Paths ────────────────────────────────────────────────────────────────────
ref_dir   = "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Datasets/Reference/"
out_dir   = "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/"

lake_path  = os.path.join(ref_dir, "mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2.h5ad")
query_path = os.path.join(out_dir, "PT_Subclustered_Annotated.h5ad")

os.makedirs(out_dir, exist_ok=True)
log_dir = os.path.join(out_dir, "logs")
os.makedirs(log_dir, exist_ok=True)

for p in [lake_path, query_path]:
    print(f"{os.path.basename(p):60s}  exists={os.path.exists(p)}")

# ── Load Lake reference and subset to PT cells ───────────────────────────────
print("\n=== Loading Lake reference ===")
adata_lake = sc.read_h5ad(lake_path)
print(adata_lake)
print("Lake obs columns:", adata_lake.obs.columns.tolist())

lake_l1_col = "SubclassLevel1"
lake_l2_col = "SubclassLevel2"

print(f"\nLake {lake_l1_col} distribution:")
print(adata_lake.obs[lake_l1_col].value_counts())

print("\n=== Subsetting Lake reference to PT cells ===")
adata_lake_pt = adata_lake[adata_lake.obs[lake_l1_col] == "PT"].copy()
print(f"Lake PT cells: {adata_lake_pt.n_obs}")

print(f"\nLake PT {lake_l2_col} distribution:")
print(adata_lake_pt.obs[lake_l2_col].value_counts())

del adata_lake

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
adata_lake_pt = ensembl_to_symbol(adata_lake_pt)
print("Lake PT var sample:", adata_lake_pt.var_names[:5].tolist())

# Sanity check: Lake PT labels populated
n_lake_labeled = adata_lake_pt.obs[lake_l2_col].notna().sum()
print(f"\nLake PT obs with '{lake_l2_col}' populated: {n_lake_labeled} / {adata_lake_pt.n_obs}")
if n_lake_labeled == 0:
    raise ValueError(
        f"Column '{lake_l2_col}' is entirely empty in adata_lake_pt. "
        "Check the correct column name."
    )

# ── Load query dataset ───────────────────────────────────────────────────────
print("\n=== Loading PT query dataset ===")
adata_query = sc.read_h5ad(query_path)
print(adata_query)
print("Query obs columns:", adata_query.obs.columns.tolist())

if "DataSet" in adata_query.obs.columns:
    print("DataSet distribution:", adata_query.obs["DataSet"].value_counts().to_dict())

# ── Use raw counts for scVI if available ────────────────────────────────────
# zellkonverter exports Seurat counts → h5ad layers["counts"]
if "counts" in adata_query.layers:
    print("Using 'counts' layer from query for scVI.")
    adata_query.X = adata_query.layers["counts"].copy()
else:
    print("No 'counts' layer found in query — using X as-is.")

# ── scVI / scANVI annotation using Lake PT reference ─────────────────────────
LABEL_KEY = "cell_type"
ref_name  = "Lake_PT"

print(f"\n{'='*60}")
print(f"  {ref_name} reference (SubclassLevel2 labels)")
print(f"{'='*60}")

# Unified label column
adata_lake_pt.obs[LABEL_KEY] = adata_lake_pt.obs[lake_l2_col].astype(str)

# Shared genes
shared = adata_lake_pt.var_names.intersection(adata_query.var_names)
print(f"  Shared genes (Lake PT ∩ query): {len(shared)}")

ref_sub   = adata_lake_pt[:, shared].copy()
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
sc.pl.umap(ref_sub, color=LABEL_KEY, title=f"{ref_name} reference — PT subtype", ax=ax, show=False)
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
sc.pl.umap(combined, color=LABEL_KEY, title=f"{ref_name} — PT subtype",   ax=axes[0], show=False)
sc.pl.umap(combined, color="batch",   title=f"{ref_name} — ref / query",  ax=axes[1], show=False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"{ref_name}_combined_scANVI_UMAP.pdf"), bbox_inches="tight")
plt.close()

# Query-only UMAP coloured by prediction and existing metadata
query_cells = combined[query_mask].copy()
query_cells.obsm["X_scANVI"] = scanvi_model.get_latent_representation()[query_mask]
sc.pp.neighbors(query_cells, use_rep="X_scANVI")
sc.tl.umap(query_cells)

color_cols = [f"scanvi_{ref_name}_label"]
for col in ["DataSet", "seurat_clusters", "STP_LakeL2", "PT_subtype_Manual"]:
    if col in query_cells.obs.columns:
        color_cols.append(col)

n_plots = len(color_cols)
fig, axes = plt.subplots(1, n_plots, figsize=(8 * n_plots, 6))
if n_plots == 1:
    axes = [axes]
for ax, col in zip(axes, color_cols):
    sc.pl.umap(query_cells, color=col, title=col, ax=ax, show=False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, f"{ref_name}_query_scANVI_UMAP.pdf"), bbox_inches="tight")
plt.close()

# ── Save predictions ─────────────────────────────────────────────────────────
preds = combined.obs.loc[
    query_mask,
    [f"scanvi_{ref_name}_label", f"scanvi_{ref_name}_confidence"],
].copy()

# Also save full probability matrix
prob_df.index = combined.obs_names[query_mask]
prob_csv = os.path.join(out_dir, f"PT_scanvi_{ref_name}_probabilities.csv")
prob_df.to_csv(prob_csv)
print("Probability matrix saved to:", prob_csv)

# Attach existing metadata columns for context
extra_cols = [
    c for c in adata_query.obs.columns
    if c in ("DataSet", "seurat_clusters", "STP_LakeL2", "STP_LakeL2.max",
             "PT_subtype_Manual", "FinalAnnotation_HC")
]
if extra_cols:
    preds = preds.join(adata_query.obs[extra_cols], how="left")

csv_path = os.path.join(out_dir, "PT_scanvi_Lake_annotations.csv")
preds.to_csv(csv_path)
print("Predictions saved to:", csv_path)
print(preds.head(10))

print(f"\nPrediction distribution (scANVI {ref_name}):")
print(preds[f"scanvi_{ref_name}_label"].value_counts())

print(f"\n{ref_name} done.")
print("\nDone.")
