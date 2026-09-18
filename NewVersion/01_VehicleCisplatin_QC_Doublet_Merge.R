#!/usr/bin/env Rscript
## ============================================================
## 01_VehicleCisplatin_QC_Doublet_Merge.R
## Author: Xin Wang
## Email:  xin.wang@nationwidechildrens.org
##
## Description:
##   Stage 1 of the Vehicle + Cisplatin (non-GC pair) clustering/annotation
##   pipeline. Mirrors the QC / doublet-detection / merge strategy used in
##   the canonical 4-condition pipeline (Integration/Cisplatin_Vehicle_GC_Integration.r),
##   scoped to just the Vehicle and Cisplatin raw matrices, with ONE
##   deliberate deviation: percent.mt QC thresholds are condition-specific
##   instead of a single global cutoff.
##       Vehicle   cells must pass percent.mt < 20
##       Cisplatin cells must pass percent.mt < 40
##   (The reference script used a single global "percent.mt < 50" filter
##   at the per-sample QC step; everything else here — min.cells/min.features,
##   nFeature/nCount bounds, rDNA bound, per-sample preprocessing before
##   doublet detection, DoubletFinder parameters, and the merge strategy —
##   tracks the reference as closely as possible.)
##
## Raw matrix source (confirmed from PreProcess/PipseekerPipeline_Vehicle.sh
## and PreProcess/PipseekerPipeline_Cisplatin.sh — NOT from the stale macOS
## paths in Integration/Cisplastin_Vehicle_Integration.Rmd):
##   PipSeeker "full" output-path for each sample is
##     .../RawDZscRNAseq/Results/<id>/          (id = Vehicle | Cisplatin)
##   which on THIS cluster resolves (confirmed present on disk) to:
##     /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Preprocess/<id>/raw_matrix/
##       matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
##   This matches exactly the ReadMtx() path convention used in
##   Integration/Cisplatin_Vehicle_GC_Integration.r ( Indir/<i>/raw_matrix/... ).
##
## Input : raw PipSeeker matrices for Vehicle and Cisplatin
##           /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Preprocess/Vehicle/raw_matrix/
##           /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Preprocess/Cisplatin/raw_matrix/
## Output: /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/QC/
##           VehicleCisplatin_scrna_merged_withdoublets.RDS   (merged, doublet_id added, pre-removal)
##           VehicleCisplatin_scrna_QC_Doublet_Singlet.RDS    (final: doublets removed — feeds script 02)
##           QC violin plots, per-sample doublet diagnostic plots
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(DoubletFinder)
})

set.seed(10000)  # matches reference seed

## ---- Directories -------------------------------------------------------
RawDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Preprocess/"
OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/QC/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
setwd(OutDir)

Conditions <- c("Vehicle", "Cisplatin")

## Condition-specific mito QC thresholds — the one deliberate deviation from
## the reference pipeline (which used a single global percent.mt < 50 at
## this stage). Everything else in the per-sample filter below matches the
## reference's per-sample subset() criteria in Cisplatin_Vehicle_GC_Integration.r.
MitoThreshold <- c(Vehicle = 20, Cisplatin = 40)

## ---- Step 1: per-sample read-in, QC, and pre-processing -----------------
cat("=== Step 1: Reading raw matrices and applying per-sample QC ===\n")

scrna.list <- list()

for (i in Conditions) {
  cat("\n--- Sample:", i, "---\n")
  mat_dir <- paste0(RawDir, i, "/raw_matrix/")

  counts <- ReadMtx(
    mtx      = paste0(mat_dir, "matrix.mtx.gz"),
    features = paste0(mat_dir, "features.tsv.gz"),
    cells    = paste0(mat_dir, "barcodes.tsv.gz")
  )
  counts <- as(counts, "dgCMatrix")

  obj <- CreateSeuratObject(
    counts       = counts,
    project      = i,
    min.cells    = 3,     # matches reference
    min.features = 200    # matches reference
  )

  obj <- RenameCells(obj, add.cell.id = i)
  obj$DataSet <- rep(i, length(obj$orig.ident))

  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj$rDNA       <- PercentageFeatureSet(obj, pattern = "^Rp[sl][[:digit:]]")

  mt_cut <- MitoThreshold[[i]]
  cat("  Cells before QC filter:", ncol(obj), "\n")
  cat("  Applying: nFeature_RNA > 150 & nCount_RNA >= 150 & nFeature_RNA < 10000 &",
      "nCount_RNA < 20000 & percent.mt <", mt_cut, "& rDNA < 40\n")

  obj <- subset(
    obj,
    subset = nFeature_RNA > 150 &
      nCount_RNA >= 150 &
      nFeature_RNA < 10000 &
      nCount_RNA < 20000 &
      percent.mt < mt_cut &
      rDNA < 40
  )
  cat("  Cells after QC filter:", ncol(obj), "\n")

  scrna.list[[i]] <- obj
}

## ---- Step 2: merge samples ------------------------------------------------
cat("\n=== Step 2: Merging Vehicle + Cisplatin into one object ===\n")

scrna <- merge(
  x = scrna.list[[1]],
  y = scrna.list[2:length(scrna.list)],
  project = "VehicleCisplatin"
)

scrna@meta.data$DataSet <- ordered(
  factor(scrna@meta.data$DataSet),
  levels = c("Vehicle", "Cisplatin")
)
cat("DataSet counts after merge:\n")
print(table(scrna@meta.data$DataSet))

# Recompute percent.mt / rDNA at the merged-object level (matches reference
# pattern of recomputing after merge, even though per-sample QC already used
# these same metrics before subsetting).
scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-")
scrna[["rDNA"]]       <- PercentageFeatureSet(scrna, pattern = "^Rp[sl][[:digit:]]")

QCReport <- VlnPlot(
  scrna,
  features = c("rDNA", "percent.mt", "nCount_RNA", "nFeature_RNA"),
  split.by = "DataSet",
  pt.size  = 0,
  ncol     = 1
)
ggsave(filename = paste0(OutDir, "VehicleCisplatin_SampleQC_reports.pdf"),
       plot = QCReport, height = 12, width = 8)

## ---- Step 3: doublet detection (DoubletFinder), per sample ---------------
cat("\n=== Step 3: DoubletFinder — per-sample doublet detection ===\n")

DoubletDir <- paste0(OutDir, "DoubletFinder/")
dir.create(DoubletDir, showWarnings = FALSE, recursive = TRUE)
setwd(DoubletDir)

FindDoublets <- function(library_id, seurat_aggregate) {
  seurat_obj <- subset(seurat_aggregate, idents = library_id)

  # Seurat v5 keeps per-DataSet split layers after merge; a per-sample
  # subset still needs its layers joined before any feature-level op
  # (NormalizeData/ScaleData/RunPCA read layer data directly, and stale
  # empty layers from the other sample can otherwise trip these calls).
  seurat_obj <- JoinLayers(seurat_obj)

  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

  ## pK identification (no ground truth) — matches reference exactly
  sweep.res.list_kidney <- paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
  sweep.stats_kidney    <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney          <- find.pK(sweep.stats_kidney)
  pK <- bcmvn_kidney %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))

  seurat_doublets <- doubletFinder(
    seurat_obj, PCs = 1:20, pN = 0.25, pK = pK,
    nExp = round(0.05 * length(seurat_obj@active.ident)),
    reuse.pANN = FALSE, sct = FALSE
  )

  DF.class <- names(seurat_doublets@meta.data) %>% str_subset("DF.classifications")
  pANN     <- names(seurat_doublets@meta.data) %>% str_subset("pANN")

  p1 <- ggplot(bcmvn_kidney, aes(x = pK, y = BCmetric)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("pKmax=", pK)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- DimPlot(seurat_doublets, group.by = DF.class)
  p3 <- FeaturePlot(seurat_doublets, features = pANN)

  ggsave(filename = paste0(library_id, "_pKmax_distribution.pdf"), p1, height = 5, width = 6)
  ggsave(filename = paste0(library_id, "_Dimplotdoublet_distribution.pdf"), p2, height = 5, width = 5)
  ggsave(filename = paste0(library_id, "_Featureplotdoublet_distribution.pdf"), p3, height = 5, width = 5)

  df_doublet_barcodes <- as.data.frame(cbind(
    rownames(seurat_doublets@meta.data),
    seurat_doublets@meta.data[[DF.class]]
  ))
  return(df_doublet_barcodes)
}

orig.ident <- levels(Idents(scrna))
cat("Samples for doublet detection:", paste(orig.ident, collapse = ", "), "\n")

list.doublet.bc <- lapply(orig.ident, function(x) FindDoublets(x, seurat_aggregate = scrna))

doublet_id <- list.doublet.bc %>%
  bind_rows() %>%
  dplyr::rename("doublet_id" = "V2") %>%
  tibble::column_to_rownames(var = "V1")

cat("\nDoublet call distribution:\n")
print(table(doublet_id))

scrna <- AddMetaData(scrna, doublet_id)

cat("\nSaving merged object with doublet calls (pre-removal)...\n")
setwd(OutDir)
saveRDS(scrna, file = paste0(OutDir, "VehicleCisplatin_scrna_merged_withdoublets.RDS"))

## ---- Step 4: remove doublets, save final QC object ------------------------
cat("\n=== Step 4: Removing doublets, saving final QC/merge object ===\n")

scrna <- readRDS(paste0(OutDir, "VehicleCisplatin_scrna_merged_withdoublets.RDS"))
scrna <- subset(scrna, subset = (doublet_id == "Singlet"))

# Drop unused DataSet factor levels (defensive — none expected here since
# only Vehicle/Cisplatin were ever loaded, but this mirrors the convention
# used throughout the rest of this pipeline whenever a Seurat object is
# subsetted).
scrna@meta.data$DataSet <- droplevels(factor(scrna@meta.data$DataSet))

VQCReport <- VlnPlot(
  scrna,
  features = c("rDNA", "percent.mt", "nCount_RNA", "nFeature_RNA"),
  pt.size  = 0,
  ncol     = 1
)
ggsave(filename = paste0(OutDir, "VehicleCisplatin_SampleQC_reports_PostDoublet.pdf"),
       plot = VQCReport, height = 12, width = 4)

cat("Final cell counts (post-QC, post-doublet-removal):\n")
print(table(scrna@meta.data$DataSet))

cat("Saving final QC/doublet/merge object...\n")
saveRDS(scrna, file = paste0(OutDir, "VehicleCisplatin_scrna_QC_Doublet_Singlet.RDS"))

cat("\nDone. Output:", paste0(OutDir, "VehicleCisplatin_scrna_QC_Doublet_Singlet.RDS"), "\n")
