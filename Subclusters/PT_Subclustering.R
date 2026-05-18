#!/usr/bin/env Rscript
# =============================================================================
# PT_Subclustering.R
# Author: Xin Wang
# Date:   2026-05-18
# Description:
#   Select Proximal Tubule (PT) cells from the final high-confidence annotated
#   dataset and perform subclustering to resolve PT subtypes (S1, S2, S3,
#   injured/aPT, cycling/cycPT, damaged/dPT).
#
#   Strategy mirrors the whole-dataset annotation pipeline:
#     1. Subset PT cells from FinalAnnotation_HC == "PT"
#     2. Re-run Seurat pipeline with Harmony batch correction
#     3. Seurat label transfer from Lake snRNA reference (SubclassLevel2)
#     4. FeaturePlots / DotPlots for PT subtype markers
#     5. Manual annotation framework (fill in cluster-to-subtype mapping)
#     6. Save annotated Seurat object + h5ad for downstream scVI
#
#   Input : IntegrateL1/Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS
#   Output: Results/Subclusters/PT/
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(sctransform)
  library(ggthemes)
  library(Matrix)
  library(limma)
  library(tidyverse)
  library(data.table)
  library(cowplot)
  library(ggsci)
  library(ggplot2)
  library(gridExtra)
  library(harmony)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(here)
  library(future)
  library(zellkonverter)
})

options(future.globals.maxSize = 1e15)
set.seed(10000)

# =============================================================================
# Directories
# =============================================================================
AnnDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/IntegrateL1/"
RefDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Datasets/Reference/"
OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
setwd(OutDir)

# =============================================================================
# Helper: convert Ensembl IDs -> gene symbols in a Seurat object
# =============================================================================
ensembl_to_symbol_seurat <- function(obj) {
  default_assay <- DefaultAssay(obj)
  avail_layers  <- tryCatch(Layers(obj[[default_assay]]), error = function(e) character(0))
  message("Available layers in '", default_assay, "': ", paste(avail_layers, collapse = ", "))

  preferred    <- c("counts", "data", "scale.data")
  layer_to_use <- preferred[preferred %in% avail_layers]
  if (length(layer_to_use) == 0) layer_to_use <- avail_layers
  layer_to_use <- layer_to_use[1]

  counts <- tryCatch(
    GetAssayData(obj, layer = layer_to_use),
    error = function(e) GetAssayData(obj, slot = layer_to_use)
  )

  if (is.null(counts) || nrow(counts) == 0)
    stop("Could not extract a non-empty matrix from assay '", default_assay,
         "' (layers: ", paste(avail_layers, collapse = ", "), ")")

  genes <- rownames(counts)

  if (!grepl("^ENSMUS", genes[1])) {
    message("Gene names already look like symbols — no conversion needed.")
    return(obj)
  }
  message("Detected Ensembl IDs — converting to gene symbols via org.Mm.eg.db ...")

  genes_clean <- sub("\\.[0-9]+$", "", genes)

  sym_map <- mapIds(
    org.Mm.eg.db,
    keys      = genes_clean,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  sym_map[is.na(sym_map)] <- genes[is.na(sym_map)]
  sym_map <- make.unique(as.character(sym_map))

  stopifnot("sym_map length != nrow(counts)" = length(sym_map) == nrow(counts))
  rownames(counts) <- sym_map

  new_obj <- CreateSeuratObject(
    counts    = counts,
    meta.data = obj@meta.data,
    project   = Project(obj)
  )
  DefaultAssay(new_obj) <- "RNA"
  message(sprintf("  Renamed %d genes (sample: %s ...)", length(sym_map),
                  paste(head(sym_map, 5), collapse = ", ")))
  return(new_obj)
}

# =============================================================================
# Step 1: Load the full annotated object and subset PT cells
# =============================================================================
cat("=== Loading full annotated object ===\n")
scrna_full <- readRDS(
  paste0(AnnDir, "Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS")
)

cat("FinalAnnotation_HC distribution:\n")
print(table(scrna_full@meta.data$FinalAnnotation_HC, useNA = "always"))

# Enforce dataset order
scrna_full@meta.data$DataSet <- ordered(
  factor(scrna_full@meta.data$DataSet),
  levels = c("Vehicle", "Vehicle-GC", "Cisplatin", "Cisplatin-GC")
)

# Subset PT cells, excluding Vehicle-GC
pt <- subset(scrna_full, subset = FinalAnnotation_HC == "PT" & DataSet != "Vehicle-GC")
cat("\nPT cells per dataset:\n")
print(table(pt@meta.data$DataSet))

rm(scrna_full)
gc()

# =============================================================================
# Step 2: Re-run Seurat pipeline on PT subset
# =============================================================================
cat("\n=== Re-running Seurat pipeline on PT subset ===\n")

DefaultAssay(pt) <- "RNA"
pt <- JoinLayers(pt)

pt <- NormalizeData(pt)
pt <- FindVariableFeatures(pt, selection.method = "vst", nfeatures = 3000)
pt <- ScaleData(pt)
pt <- RunPCA(pt, verbose = TRUE, npcs = 30)

# Elbow plot to assess PC dimensionality
ggsave(
  filename = "PT_ElbowPlot.pdf",
  plot     = ElbowPlot(pt, ndims = 30),
  height   = 4, width   = 5
)

# Harmony batch correction by DataSet
pt <- RunHarmony(pt, group.by.vars = "DataSet", reduction.save = "harmony")
pt <- RunUMAP(pt,      dims = 1:20, reduction = "harmony")
pt <- FindNeighbors(pt, dims = 1:20, reduction = "harmony")
pt <- FindClusters(pt, resolution = 0.4, algorithm = 1)

cat("PT clusters at resolution 0.4:\n")
print(table(pt@meta.data$seurat_clusters))

# Initial UMAP
ggsave(
  filename = "PT_DimPlot_clusters_res0.4.pdf",
  plot     = DimPlot(pt, reduction = "umap", label = TRUE) + ggtitle("PT subclusters (res 0.4)"),
  height   = 5, width = 6
)
ggsave(
  filename = "PT_DimPlot_clusters_res0.4_SplitByDataSet.pdf",
  plot     = DimPlot(pt, reduction = "umap", split.by = "DataSet", label = TRUE, ncol = 3),
  height   = 5, width = 16
)

saveRDS(pt, file = paste0(OutDir, "PT_Subclustered_res0.4.RDS"))

# =============================================================================
# Step 3: Seurat label transfer — Lake snRNA reference (PT cells, SubclassLevel2)
# =============================================================================
cat("\n=== Seurat label transfer from Lake snRNA reference (PT cells only, SubclassLevel2) ===\n")

reference_lake <- readRDS(paste0(RefDir, "LakesnRNA_seurat.rds"))
reference_lake <- UpdateSeuratObject(reference_lake)

reference_lake <- ensembl_to_symbol_seurat(reference_lake)

# Subset Lake reference to PT-annotated cells only for a focused transfer
cat("Lake SubclassLevel1 distribution before subsetting:\n")
print(table(reference_lake@meta.data$SubclassLevel1))

reference_lake_PT <- subset(reference_lake, subset = SubclassLevel1 == "PT")
cat("Lake PT cells for reference:", ncol(reference_lake_PT), "\n")
cat("Lake PT SubclassLevel2 subtypes:\n")
print(table(reference_lake_PT@meta.data$SubclassLevel2))

rm(reference_lake)
gc()

reference_lake_PT <- NormalizeData(reference_lake_PT)
reference_lake_PT <- FindVariableFeatures(reference_lake_PT, selection.method = "vst", nfeatures = 4000)
reference_lake_PT <- ScaleData(reference_lake_PT)

features_lake <- intersect(VariableFeatures(reference_lake_PT), VariableFeatures(pt))
cat("Lake PT shared variable features:", length(features_lake), "\n")

# Transfer SubclassLevel2 (PT-S1, PT-S1/2, PT-S2, PT-S3, aPT, dPT, cycPT, ...)
anchors_lake <- FindTransferAnchors(
  reference = reference_lake_PT,
  query     = pt,
  features  = features_lake,
  dims      = 1:30
)

predictions_lake_L2 <- TransferData(
  anchorset = anchors_lake,
  refdata   = reference_lake_PT@meta.data[["SubclassLevel2"]],
  dims      = 1:30
)

colnames(predictions_lake_L2) <- sub("^predicted\\.id$",           "STP_LakeL2",     colnames(predictions_lake_L2))
colnames(predictions_lake_L2) <- sub("^prediction\\.score\\.max$",  "STP_LakeL2.max", colnames(predictions_lake_L2))
colnames(predictions_lake_L2) <- sub("^prediction\\.score\\.(.+)$", "STP_LakeL2.\\1", colnames(predictions_lake_L2))

pt <- AddMetaData(pt, metadata = predictions_lake_L2)

rm(reference_lake_PT, anchors_lake, predictions_lake_L2)
gc()

cat("\nLake L2 transfer distribution (PT cells):\n")
print(table(pt@meta.data$STP_LakeL2, useNA = "always"))

# Ordered factor levels for L2
LakeL2_PT_Order <- c(
  "PT-S1", "PT-S1/2", "PT-S2", "PT-S3",
  "aPT", "dPT", "cycPT",
  "PapE"
)

pt@meta.data$STP_LakeL2 <- factor(pt@meta.data$STP_LakeL2)

ggsave(
  filename = "PT_DimPlot_STP_LakeL2.pdf",
  plot     = DimPlot(pt, group.by = "STP_LakeL2", label = TRUE, repel = TRUE) +
             ggtitle("Lake L2 label transfer (PT)") + scale_color_hue(),
  height   = 5, width = 8
)
ggsave(
  filename = "PT_DimPlot_STP_LakeL2_SplitByDataSet.pdf",
  plot     = DimPlot(pt, group.by = "STP_LakeL2", split.by = "DataSet",
                     label = TRUE, repel = TRUE, ncol = 3) +
             ggtitle("Lake L2 label transfer by dataset") + scale_color_hue(),
  height   = 5, width = 16
)

saveRDS(pt, file = paste0(OutDir, "PT_Subclustered_LakeTransfer.RDS"))

# =============================================================================
# Step 4: PT subtype marker visualization
# =============================================================================
cat("\n=== Generating PT subtype marker plots ===\n")

# Canonical PT subtype markers
PT_markers <- list(
  PT_S1    = c("Slc5a2", "Slc34a1", "Lrp2", "Cubn", "Slc27a2", "Slc17a3"),
  PT_S2    = c("Slc34a1", "Lrp2", "Slc22a30", "Slc16a9"),
  PT_S3    = c("Slc22a6", "Slc22a8", "Slc22a30", "Slc16a9"),
  aPT      = c("Havcr1", "Vcam1", "Vim", "Cdh6", "Ly6c1"),
  dPT      = c("Sox4", "Tgfb1", "Cldn1", "Fn1", "Hspa1a"),
  cycPT    = c("Mki67", "Top2a", "Pcna", "Cdk1", "Cenpa"),
  General  = c("Lrp2", "Slc34a1", "Slc13a3", "Slc5a2",
               "Slc22a6", "Slc22a8", "Havcr1")
)

# Combined feature vector for DotPlot (deduplicated, ordered)
all_PT_features <- unique(c(
  PT_markers$General,
  PT_markers$PT_S1,
  PT_markers$PT_S2,
  PT_markers$PT_S3,
  PT_markers$aPT,
  PT_markers$dPT,
  PT_markers$cycPT
))
all_PT_features <- intersect(all_PT_features, rownames(pt))

# DotPlot across clusters
ggsave(
  filename = "PT_MarkerDotPlot_ByCluster.pdf",
  plot     = DotPlot(pt, features = all_PT_features, group.by = "seurat_clusters",
                     cols = "RdYlBu") +
             RotatedAxis() + ggtitle("PT subtype markers by cluster"),
  height   = 6, width = 14
)

# DotPlot grouped by Lake L2 prediction
ggsave(
  filename = "PT_MarkerDotPlot_BySTP_LakeL2.pdf",
  plot     = DotPlot(pt, features = all_PT_features, group.by = "STP_LakeL2",
                     cols = "RdYlBu") +
             RotatedAxis() + ggtitle("PT subtype markers by Lake L2 transfer"),
  height   = 6, width = 14
)

# FeaturePlots for each subtype marker group
MarkerDir <- file.path(OutDir, "MarkerFeaturePlots")
dir.create(MarkerDir, showWarnings = FALSE, recursive = TRUE)

for (subtype in names(PT_markers)) {
  markers_present <- intersect(PT_markers[[subtype]], rownames(pt))
  if (length(markers_present) == 0) {
    message("Skipping ", subtype, " — no markers found in object")
    next
  }
  subtype_dir <- file.path(MarkerDir, subtype)
  dir.create(subtype_dir, showWarnings = FALSE, recursive = TRUE)

  ncols <- min(3, length(markers_present))

  fp <- FeaturePlot(pt, features = markers_present, ncol = ncols)
  ggsave(
    filename = file.path(subtype_dir, paste0(subtype, "_FeaturePlot.pdf")),
    plot     = fp,
    width    = 6 * ncols,
    height   = 5 * ceiling(length(markers_present) / ncols)
  )

  fp_split <- FeaturePlot(pt, features = markers_present,
                          split.by = "DataSet", ncol = ncols)
  ggsave(
    filename = file.path(subtype_dir, paste0(subtype, "_FeaturePlot_SplitByDataSet.pdf")),
    plot     = fp_split,
    width    = 6 * ncols,
    height   = 5 * ceiling(length(markers_present) / ncols)
  )

  dp <- DotPlot(pt, features = markers_present) + RotatedAxis() +
        ggtitle(paste(subtype, "markers"))
  ggsave(
    filename = file.path(subtype_dir, paste0(subtype, "_DotPlot.pdf")),
    plot     = dp,
    width    = 3 + length(markers_present) * 0.8,
    height   = 5
  )

  message("Saved plots for: ", subtype)
}

# =============================================================================
# Step 5: Find cluster markers (for manual annotation support)
# =============================================================================
cat("\n=== Finding cluster markers ===\n")

Idents(pt) <- "seurat_clusters"
pt_markers <- FindAllMarkers(
  pt,
  only.pos    = TRUE,
  min.pct     = 0.25,
  logfc.threshold = 0.25
)

write.csv(pt_markers,
          file = paste0(OutDir, "PT_Cluster_Markers.csv"),
          row.names = FALSE)

# Top 5 markers per cluster for quick review
top5 <- pt_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)
write.csv(top5,
          file = paste0(OutDir, "PT_Cluster_Top5Markers.csv"),
          row.names = FALSE)

cat("Top 5 markers per PT cluster:\n")
print(top5)

# =============================================================================
# Step 6: Manual annotation
# =============================================================================
# -----------------------------------------------------------------------------
# INSTRUCTIONS: Inspect the DotPlots, FeaturePlots, and Lake L2 transfer, then
# fill in the CellType mapping below. Current placeholders keep cluster numbers.
# PT subtype labels to use:
#   PT-S1  : high Slc5a2, Slc34a1, Lrp2 (S1 segment)
#   PT-S2  : moderate Slc34a1, Lrp2 (S2 segment)
#   PT-S3  : high Slc22a6, Slc22a8 (S3 segment)
#   aPT    : high Havcr1/Kim-1, Vcam1 (injured/adaptive)
#   dPT    : high Sox4, Tgfb1, Cldn1 (failed repair / damaged)
#   cycPT  : high Mki67, Top2a (cycling)
# -----------------------------------------------------------------------------
cat("\n=== Manual annotation (placeholder — update mappings below) ===\n")

FinalCluster <- as.character(levels(factor(pt@meta.data$seurat_clusters)))
CellType_PT  <- data.frame(cluster = FinalCluster, subtype = FinalCluster)

# ---- UPDATE THESE MAPPINGS AFTER INSPECTING PLOTS ----
# Example (adjust cluster numbers based on your data):
# CellType_PT[CellType_PT$cluster %in% c("0", "1"),   "subtype"] <- "PT-S1"
# CellType_PT[CellType_PT$cluster %in% c("2", "3"),   "subtype"] <- "PT-S2"
# CellType_PT[CellType_PT$cluster %in% c("4"),        "subtype"] <- "PT-S3"
# CellType_PT[CellType_PT$cluster %in% c("5"),        "subtype"] <- "aPT"
# CellType_PT[CellType_PT$cluster %in% c("6"),        "subtype"] <- "dPT"
# CellType_PT[CellType_PT$cluster %in% c("7"),        "subtype"] <- "cycPT"
# -------------------------------------------------------

pt@meta.data$PT_subtype_Manual <- sapply(
  as.character(pt@meta.data$seurat_clusters),
  function(x) CellType_PT[CellType_PT$cluster == x, "subtype"]
)

cat("Manual PT subtype distribution (placeholder):\n")
print(table(pt@meta.data$PT_subtype_Manual))

ggsave(
  filename = "PT_DimPlot_ManualAnnotation.pdf",
  plot     = DimPlot(pt, group.by = "PT_subtype_Manual",
                     label = TRUE, repel = TRUE, split.by = "DataSet", ncol = 3) +
             scale_color_hue() + ggtitle("PT manual subtype annotation"),
  height   = 5, width = 16
)

# =============================================================================
# Step 7: Save annotated PT Seurat object
# =============================================================================
cat("\n=== Saving annotated PT Seurat object ===\n")

saveRDS(pt,
        file = paste0(OutDir, "PT_Subclustered_Annotated.RDS"))
cat("Saved:", paste0(OutDir, "PT_Subclustered_Annotated.RDS"), "\n")

# =============================================================================
# Step 8: Export to h5ad for downstream scVI / scANVI annotation
# =============================================================================
cat("\n=== Exporting to h5ad for scVI ===\n")

pt_export <- JoinLayers(pt)
sce <- as.SingleCellExperiment(pt_export)

writeH5AD(
  sce,
  file = paste0(OutDir, "PT_Subclustered_Annotated.h5ad")
)
cat("h5ad written to:", paste0(OutDir, "PT_Subclustered_Annotated.h5ad"), "\n")
cat("Next step: run PT_scvi_annotation.ipynb for scANVI-based subtype annotation\n")
