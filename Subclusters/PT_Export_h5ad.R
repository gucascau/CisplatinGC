#!/usr/bin/env Rscript
# =============================================================================
# PT_Export_h5ad.R
# Author: Xin Wang
# Date:   2026-06-30
# Description:
#   Export PT_Subclustered_Annotated.RDS → PT_Subclustered_Annotated.h5ad
#   for downstream scVI / scANVI annotation.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(zellkonverter)
})

OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/"

cat("=== Loading PT_Subclustered_Annotated.RDS ===\n")
pt <- readRDS(paste0(OutDir, "PT_Subclustered_Annotated.RDS"))

cat("PT object summary:\n")
print(pt)
cat("\nMetadata columns:\n")
print(colnames(pt@meta.data))
cat("\nDataSet distribution:\n")
print(table(pt@meta.data$DataSet))
cat("\nSTP_LakeL2 distribution:\n")
print(table(pt@meta.data$STP_LakeL2, useNA = "always"))
cat("\nPT_subtype_Manual distribution:\n")
print(table(pt@meta.data$PT_subtype_Manual, useNA = "always"))

cat("\n=== Joining layers and exporting to h5ad ===\n")
pt <- JoinLayers(pt)
sce <- as.SingleCellExperiment(pt)

h5ad_path <- paste0(OutDir, "PT_Subclustered_Annotated.h5ad")
writeH5AD(sce, file = h5ad_path)
cat("h5ad written to:", h5ad_path, "\n")
