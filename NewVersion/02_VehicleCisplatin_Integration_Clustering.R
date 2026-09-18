#!/usr/bin/env Rscript
## ============================================================
## 02_VehicleCisplatin_Integration_Clustering.R
## Author: Xin Wang
## Email:  xin.wang@nationwidechildrens.org
##
## Description:
##   Stage 2 of the Vehicle + Cisplatin pipeline. Normalizes, regresses out
##   percent.mt, and integrates the two conditions with Harmony, then
##   clusters and runs UMAP.
##
##   Parameter provenance: the reference repo actually has TWO Harmony
##   clustering passes for this project —
##     (a) the first-pass integration in Integration/Cisplatin_Vehicle_GC_Integration.r
##         (NormalizeData/ScaleData with no vars.to.regress, PCA npcs=30,
##          RunHarmony(group.by.vars="DataSet"), dims=1:30, FindClusters
##          resolution=0.4, algorithm=1), and
##     (b) the later re-clustering pass in
##         Annotation/Cisplatin_GC_Vehicle_Annotation_Integration.Rmd, run
##         directly on the Vehicle/Cisplatin(/Cisplatin-GC) subset right
##         before annotation, which regresses out percent.mt, uses PCA
##         npcs=30, RunHarmony(group.by.vars="DataSet", reduction="pca"),
##         FindNeighbors/FindClusters on dims=1:20 at resolution=1
##         (algorithm=1), and RunUMAP(dims=1:20, min.dist=0.1, spread=1.5,
##          n.neighbors=20). That Rmd explicitly calls this "the same
##          parameters as the main paper".
##   Since this script's position in a from-scratch Vehicle+Cisplatin
##   pipeline is exactly analogous to pass (b) — one clustering pass that
##   directly feeds annotation — this script follows pass (b)'s parameters
##   (including the percent.mt regression). ASSUMPTION flagged: pass (b) is
##   treated as the canonical/current strategy per its own in-file comment;
##   if pass (a)'s parameters were intended instead, only the four lines
##   marked "# pass (b) parameter" below would need to change.
##
##   NOTE ON "SCTransform": the task brief for this pipeline mentioned
##   SCTransform / vars.to.regress. A repo-wide search found SCTransform()
##   is never actually called anywhere in this codebase (only
##   `library(sctransform)` is loaded, unused). The real, current
##   normalization strategy in the reference pipeline is standard
##   NormalizeData() + ScaleData(vars.to.regress = "percent.mt"), which is
##   what this script replicates.
##
## Input : /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/QC/
##           VehicleCisplatin_scrna_QC_Doublet_Singlet.RDS   (from script 01)
## Output: /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Clustering/
##           VehicleCisplatin_Harmony_Clustered.RDS
##           ElbowPlot, DimPlots (combined + split by DataSet)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})

options(future.globals.maxSize = 1e15)

## ---- Directories -------------------------------------------------------
InDir  <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/QC/"
OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Clustering/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
setwd(OutDir)

RDSFile <- paste0(InDir, "VehicleCisplatin_scrna_QC_Doublet_Singlet.RDS")

cat("Reading:", RDSFile, "\n")
scrna <- readRDS(RDSFile)
cat("Cells:", ncol(scrna), "  Genes:", nrow(scrna), "\n")

DefaultAssay(scrna) <- "RNA"

# Defensive: drop any stale factor levels before this object heads into
# clustering / downstream table()/ggplot calls.
scrna@meta.data$DataSet <- droplevels(factor(scrna@meta.data$DataSet))
cat("DataSet levels:", paste(levels(scrna@meta.data$DataSet), collapse = ", "), "\n")

## ---- Normalize + regress out percent.mt --------------------------------
cat("\n=== Normalization (percent.mt regressed out, matching reference re-clustering pass) ===\n")

mt_genes     <- grep("^mt-", rownames(scrna), value = TRUE, ignore.case = TRUE)
all_genes    <- rownames(scrna)
non_mt_genes <- setdiff(all_genes, mt_genes)

# percent.mt should already exist from script 01; recompute defensively
if (!"percent.mt" %in% colnames(scrna@meta.data)) {
  cat("percent.mt missing — recomputing...\n")
  scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-")
}

scrna <- scrna %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst")

VariableFeatures(scrna) <- setdiff(VariableFeatures(scrna), mt_genes)

scrna <- ScaleData(scrna, features = non_mt_genes, vars.to.regress = "percent.mt")  # pass (b) parameter

## ---- PCA + Harmony -------------------------------------------------------
cat("\n=== PCA + Harmony integration ===\n")

scrna <- RunPCA(scrna, verbose = TRUE, npcs = 30)  # pass (b) parameter

elbow_p <- ElbowPlot(scrna, ndims = 50)
ggsave(paste0(OutDir, "VehicleCisplatin_ElbowPlot_PCA.pdf"), plot = elbow_p, height = 4, width = 6)

scrna <- RunHarmony(
  scrna,
  group.by.vars   = "DataSet",   # pass (b) parameter — matches reference group.by.vars
  reduction       = "pca",
  reduction.save  = "harmony",
  verbose         = TRUE
)

# Seurat v5 keeps per-DataSet split layers after merge; join them now
# (matches the reference's placement of JoinLayers() right after
# RunHarmony, before any downstream feature-level operation).
scrna <- JoinLayers(scrna)

## ---- Clustering + UMAP ----------------------------------------------------
cat("\n=== Clustering + UMAP ===\n")

scrna <- FindNeighbors(scrna, reduction = "harmony", dims = 1:20)   # pass (b) parameter
scrna <- FindClusters(scrna, resolution = 1, algorithm = 1)         # pass (b) parameter

scrna <- RunUMAP(
  scrna,
  reduction   = "harmony",
  dims        = 1:20,
  min.dist    = 0.1,
  spread      = 1.5,
  n.neighbors = 20
)

cat("Cluster sizes:\n")
print(table(scrna$seurat_clusters))

## ---- Save + basic DimPlots ------------------------------------------------
cat("\nSaving clustered object...\n")
saveRDS(scrna, file = paste0(OutDir, "VehicleCisplatin_Harmony_Clustered.RDS"))

DimplotScrna <- DimPlot(scrna, reduction = "umap", label = TRUE) +
  ggtitle("Vehicle + Cisplatin — Harmony clusters")
ggsave(filename = paste0(OutDir, "VehicleCisplatin_DimPlot.pdf"),
       plot = DimplotScrna, width = 6, height = 5)

DimplotScrnaSplit <- DimPlot(scrna, reduction = "umap", split.by = "DataSet", label = TRUE) +
  ggtitle("Vehicle + Cisplatin — Harmony clusters, split by DataSet")
ggsave(filename = paste0(OutDir, "VehicleCisplatin_DimPlot_SplitByDataSet.pdf"),
       plot = DimplotScrnaSplit, width = 12, height = 5)

cat("\nDone. Output:", paste0(OutDir, "VehicleCisplatin_Harmony_Clustered.RDS"), "\n")
