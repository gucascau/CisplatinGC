#!/usr/bin/env Rscript
## ============================================================
## 03_VehicleCisplatin_MultiMethod_Annotation.R
## Author: Xin Wang
## Email:  xin.wang@nationwidechildrens.org
##
## Description:
##   Stage 3 of the Vehicle + Cisplatin pipeline. Reproduces the same
##   multi-method annotation strategy used in
##   Annotation/Cisplatin_GC_Vehicle_Annotation.Rmd / .R:
##     1. Seurat label transfer from the Mouse Kidney Atlas (MKA) reference
##        -> Seaurat_Transfer_Predicted_MKA (+ .max / per-class scores)
##     2. Seurat label transfer from the Lake 2025 snRNA reference
##        (SubclassLevel1 column, already L1-level)
##        -> Seaurat_Transfer_Predicted_Lake (+ .max / per-class scores)
##     3. Manual/marker-guided cluster annotation -> Raw_cell_type
##        (+ Raw_cell_type_Confidence = 1, deterministic)
##     4. scANVI label transfer (MKA reference AND Lake reference,
##        each run separately) -> scanvi_MKA_label / scanvi_MKA_confidence
##        and scanvi_Lake_label / scanvi_Lake_confidence.
##        scANVI itself runs OUTSIDE R (python/scvi-tools) — this script
##        exports the h5ad scANVI needs, and re-imports the resulting CSVs,
##        exactly mirroring the boundary used in the reference pipeline
##        (Cisplatin_GC_Vehicle_Annotation.R lines ~376-483, and the
##        scanvi_annotation_*.py scripts).
##
## Column-naming note: this deliberately keeps the "Seaurat_Transfer..."
## spelling (typo carried over from the reference R script) so that
## downstream scripts (04, and any future harmonisation code copy-pasted
## from the reference) can be reused verbatim without a rename step.
##
## Input : /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Clustering/
##           VehicleCisplatin_Harmony_Clustered.RDS        (from script 02)
##         /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Datasets/Reference/
##           MKA_seurat.rds, LakesnRNA_seurat.rds           (Seurat label transfer)
##           MouseKidneyATLAS_MKA_updated.h5ad               (scANVI, external)
##           mouse_kidney_snRNAseq_Lake2025_bioRxiv_V2.h5ad   (scANVI, external)
##
## Output: /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/
##           VehicleCisplatin_Raw_Predict_MKA_Lake_Annotation.RDS   (after Seurat transfer, before manual/scANVI)
##           VehicleCisplatin_Raw_Predict_MKA_Lake_Annotation.h5ad  (export for external scANVI run)
##           MarkerFeaturePlots/<CellType>/...                      (manual-annotation review plots)
##           VehicleCisplatin_Combined_Predict_MKA_Lake_scanvi_Annotation.RDS   (final, all 5 methods)
##
##   ** EXTERNAL STEP REQUIRED BEFORE THE scANVI READ-BACK SECTION RUNS **
##   After this script writes the h5ad above, submit the two adapted scANVI
##   jobs (outside R, on an HPC node with the cell2loc_env conda env):
##     sbatch NewVersion/03b_VehicleCisplatin_scvi_annotation_MKAOnly.sh
##     sbatch NewVersion/03c_VehicleCisplatin_scvi_annotation_LakeOnly.sh
##   These already point at this script's h5ad output and are restricted to
##   the Vehicle/Cisplatin query (no GC groups involved at all in this
##   pipeline). They write:
##     VehicleCisplatin_scanvi_MKA_annotations.csv
##     VehicleCisplatin_scanvi_Lake_annotations.csv
##   into the Annotation/ OutDir above, each with a first column "X" holding
##   the query cell barcode with a "-query" suffix (AnnData.concatenate()'s
##   default naming) plus scanvi_<ref>_label / scanvi_<ref>_confidence —
##   matching the read-back format below exactly. The read-back section is
##   guarded by file.exists() so this script can be run in two passes: once
##   to produce the h5ad, and again (after the two scANVI jobs finish) to
##   complete the merge.
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(zellkonverter)
})

options(future.globals.maxSize = 1e15)

## ---- Directories -------------------------------------------------------
InDir  <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Clustering/"
RefDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Datasets/Reference/"
OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
setwd(OutDir)

RDSFile <- paste0(InDir, "VehicleCisplatin_Harmony_Clustered.RDS")

cat("Reading:", RDSFile, "\n")
scrna <- readRDS(RDSFile)
DefaultAssay(scrna) <- "RNA"
scrna@meta.data$DataSet <- droplevels(factor(scrna@meta.data$DataSet))
cat("DataSet levels:", paste(levels(scrna@meta.data$DataSet), collapse = ", "), "\n")

## ---- Helper: Ensembl ID -> gene symbol (verbatim from reference) --------
## Detects whether rownames are Ensembl IDs (ENSMUSG...) and, if so,
## converts them to gene symbols using org.Mm.eg.db. Needed because the MKA
## / Lake reference objects can be Ensembl-indexed.
ensembl_to_symbol_seurat <- function(obj) {
  default_assay <- DefaultAssay(obj)
  avail_layers  <- tryCatch(Layers(obj[[default_assay]]), error = function(e) character(0))
  message("Available layers in '", default_assay, "': ", paste(avail_layers, collapse = ", "))

  preferred <- c("counts", "data", "scale.data")
  layer_to_use <- preferred[preferred %in% avail_layers]
  if (length(layer_to_use) == 0) layer_to_use <- avail_layers
  layer_to_use <- layer_to_use[1]

  counts <- tryCatch(
    GetAssayData(obj, layer = layer_to_use),
    error = function(e) GetAssayData(obj, slot = layer_to_use)
  )
  if (is.null(counts) || nrow(counts) == 0)
    stop("Could not extract a non-empty matrix from assay '", default_assay, "'")

  genes <- rownames(counts)
  if (!grepl("^ENSMUS", genes[1])) {
    message("Gene names already look like symbols — no conversion needed.")
    return(obj)
  }
  message("Detected Ensembl IDs — converting to gene symbols via org.Mm.eg.db ...")

  genes_clean <- sub("\\.[0-9]+$", "", genes)
  sym_map <- mapIds(org.Mm.eg.db, keys = genes_clean, column = "SYMBOL",
                     keytype = "ENSEMBL", multiVals = "first")
  sym_map[is.na(sym_map)] <- genes[is.na(sym_map)]
  sym_map <- make.unique(as.character(sym_map))
  stopifnot(length(sym_map) == nrow(counts))
  rownames(counts) <- sym_map

  new_obj <- CreateSeuratObject(counts = counts, meta.data = obj@meta.data, project = Project(obj))
  DefaultAssay(new_obj) <- "RNA"
  return(new_obj)
}

## ==========================================================================
## Step 1: Seurat label transfer — MKA reference
## ==========================================================================
cat("\n=== Step 1a: Seurat label transfer — MKA reference ===\n")

reference_mka <- readRDS(paste0(RefDir, "MKA_seurat.rds"))
reference_mka <- UpdateSeuratObject(reference_mka)
reference_mka <- ensembl_to_symbol_seurat(reference_mka)

reference_mka <- NormalizeData(reference_mka)
reference_mka <- FindVariableFeatures(reference_mka, selection.method = "vst", nfeatures = 4000)
reference_mka <- ScaleData(reference_mka)

scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 4000)

features_mka <- intersect(VariableFeatures(reference_mka), VariableFeatures(scrna))
cat("MKA shared variable features:", length(features_mka), "\n")

anchors_mka <- FindTransferAnchors(
  reference = reference_mka, query = scrna, features = features_mka, dims = 1:30
)
predictions_mka <- TransferData(
  anchorset = anchors_mka, refdata = reference_mka$author_cell_type, dims = 1:30
)
colnames(predictions_mka) <- sub("^predicted.id$",          "Seaurat_Transfer_Predicted_MKA",     colnames(predictions_mka))
colnames(predictions_mka) <- sub("^prediction\\.score\\.max$", "Seaurat_Transfer_Predicted_MKA.max", colnames(predictions_mka))
colnames(predictions_mka) <- sub("^prediction\\.score\\.(.+)$","Seaurat_Transfer_Predicted_MKA.\\1", colnames(predictions_mka))

scrna <- AddMetaData(scrna, metadata = predictions_mka)
cat("MKA transfer distribution:\n")
print(table(scrna$Seaurat_Transfer_Predicted_MKA, useNA = "always"))

## ==========================================================================
## Step 1b: Seurat label transfer — Lake snRNA reference
## ==========================================================================
cat("\n=== Step 1b: Seurat label transfer — Lake reference ===\n")

reference_lake <- readRDS(paste0(RefDir, "LakesnRNA_seurat.rds"))
reference_lake <- UpdateSeuratObject(reference_lake)
reference_lake <- ensembl_to_symbol_seurat(reference_lake)

reference_lake <- NormalizeData(reference_lake)
reference_lake <- FindVariableFeatures(reference_lake, selection.method = "vst", nfeatures = 4000)
reference_lake <- ScaleData(reference_lake)

features_lake <- intersect(VariableFeatures(reference_lake), VariableFeatures(scrna))
cat("Lake shared variable features:", length(features_lake), "\n")

lake_label_col <- "SubclassLevel1"  # L1-level label, matches reference

anchors_lake <- FindTransferAnchors(
  reference = reference_lake, query = scrna, features = features_lake, dims = 1:30
)
predictions_lake <- TransferData(
  anchorset = anchors_lake, refdata = reference_lake@meta.data[[lake_label_col]], dims = 1:30
)
colnames(predictions_lake) <- sub("^predicted\\.id$",          "Seaurat_Transfer_Predicted_Lake",     colnames(predictions_lake))
colnames(predictions_lake) <- sub("^prediction\\.score\\.max$", "Seaurat_Transfer_Predicted_Lake.max", colnames(predictions_lake))
colnames(predictions_lake) <- sub("^prediction\\.score\\.(.+)$","Seaurat_Transfer_Predicted_Lake.\\1", colnames(predictions_lake))

scrna <- AddMetaData(scrna, metadata = predictions_lake)
cat("Lake transfer distribution:\n")
print(table(scrna$Seaurat_Transfer_Predicted_Lake, useNA = "always"))

p_mka  <- DimPlot(scrna, group.by = "Seaurat_Transfer_Predicted_MKA",  label = TRUE, repel = TRUE) + ggtitle("MKA transfer")
p_lake <- DimPlot(scrna, group.by = "Seaurat_Transfer_Predicted_Lake", label = TRUE, repel = TRUE) + ggtitle("Lake transfer")
ggsave(paste0(OutDir, "VehicleCisplatin_Reference_Comparison_DimPlot.pdf"),
       plot = p_mka + p_lake, height = 5, width = 12)

cat("Saving intermediate (post-Seurat-transfer) object...\n")
saveRDS(scrna, file = paste0(OutDir, "VehicleCisplatin_Raw_Predict_MKA_Lake_Annotation.RDS"))
scrna<- readRDS(paste0(OutDir, "VehicleCisplatin_Raw_Predict_MKA_Lake_Annotation.RDS"))
## ==========================================================================
## Step 2: Manual / marker-guided cluster annotation -> Raw_cell_type
## ==========================================================================
cat("\n=== Step 2: Manual marker-guided annotation ===\n")

Finalcelltypemarkers <- c(
  "Lrp2","Slc34a1","Slc13a3",
  "Slc5a2","Slc34a1","Lrp2",
  "Slc22a6","Slc22a8","Havcr1",
  "Slc12a3","Pvalb","Wnk1","Trpm6","Fxyd2",
  "Umod","Slc12a1","Cldn10","Kcnj1","Ptger3",
  "Klk1","Slc8a1","Calb1",
  "Slc14a2","Fst","Bst1",
  "Aqp2","Hsd11b2","Scnn1g",
  "Atp6v1g3","Aqp6","Slc26a7",
  "Wt1","Nphs1","Nphs2",
  "Cdh5","Igfbp3","Pecam1",
  "Col3a1","Fbn1","Lum","Col1a1",
  "Ireb2","Alas2",
  "S100a8","S100a9","Il1b",
  "C1qa","C1qb","Aif1",
  "Cxcr6","Cd247","Cd3e",
  "Igkc","Cd79a","Cd79b",
  "Ccl5","Nkg7","Cd7"
)

CelltypeMarkers <- list(
  PT         = c("Lrp2","Slc34a1","Slc13a3"),
  PT_S1S2    = c("Slc5a2","Slc34a1","Lrp2"),
  PT_S3      = c("Slc5a1","Slc22a6","Slc22a8","Havcr1"),
  DCT        = c("Slc12a3","Pvalb","Wnk1","Trpm6","Fxyd2"),
  TAL        = c("Umod","Slc12a1","Cldn10","Kcnj1","Ptger3"),
  CNT        = c("Klk1","Slc8a1","Calb1"),
  LOH        = c("Slc14a2","Aqp1","Clcnka"),
  PC         = c("Aqp2","Hsd11b2","Scnn1g"),
  IC         = c("Atp6v1g3","Aqp6","Slc26a7"),
  POD        = c("Wt1","Nphs1","Nphs2"),
  END        = c("Cdh5","Pecam1","Kdr","Emcn","Plvap"),
  FIB        = c("Col1a1","Col3a1","Fbn1","Lum","Pdgfra","Dcn"),
  pDC        = c("Siglech","Bst2","Irf7"),
  Neutrophil = c("S100a8","S100a9","Ly6g","Mpo"),
  Macrophage = c("C1qa","C1qb","Aif1","Adgre1","Cd68"),
  T_cells    = c("Cd3e","Cd247","Trac"),
  B_cells    = c("Cd79a","Cd79b","Ms4a1"),
  NKCD8      = c("Nkg7","Ccl5","Gzmb","Cd8a"),
  ICA        = c("Atp6v1g3","Aqp6","Kit","Slc4a1","Slc26a7"),
  ICB        = c("Atp6v1g3","Slc26a4","Spink8","Hmx2"),
  Mesangial  = c("Pdgfrb","Acta2","Tagln"),
  Pericyte   = c("Pdgfrb","Rgs5","Kcnj8"),
  PEC        = c("Krt8","Krt18","Claudin1"),
  DTL        = c("Aqp1","Unc5d","Adgrl3")
)

Idents(scrna) <- "seurat_clusters"

MarkerDotPlot <- DotPlot(scrna, features = intersect(Finalcelltypemarkers, rownames(scrna))) + RotatedAxis()
ggsave(paste0(OutDir, "VehicleCisplatin_ManualAnnotation_DotPlot.pdf"),
       plot = MarkerDotPlot, height = 5, width = 14)

MarkerDir <- file.path(OutDir, "MarkerFeaturePlots")
dir.create(MarkerDir, showWarnings = FALSE, recursive = TRUE)

for (ct in names(CelltypeMarkers)) {
  markers_present <- intersect(CelltypeMarkers[[ct]], rownames(scrna))
  if (length(markers_present) == 0) {
    message("Skipping ", ct, " — no markers found in object")
    next
  }
  ct_dir <- file.path(MarkerDir, ct)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

  fp <- FeaturePlot(scrna, features = markers_present, ncol = min(3, length(markers_present)))
  ggsave(file.path(ct_dir, paste0(ct, "_FeaturePlot.pdf")), plot = fp,
         width = 6 * min(3, length(markers_present)),
         height = 5 * ceiling(length(markers_present) / 3))

  fp_split <- FeaturePlot(scrna, features = markers_present, split.by = "DataSet",
                           ncol = min(3, length(markers_present)))
  ggsave(file.path(ct_dir, paste0(ct, "_FeaturePlot_SplitByCondition.pdf")), plot = fp_split,
         width = 6 * min(3, length(markers_present)),
         height = 5 * ceiling(length(markers_present) / 3))

  dp <- DotPlot(scrna, features = markers_present) + RotatedAxis() + ggtitle(paste(ct, "markers"))
  ggsave(file.path(ct_dir, paste0(ct, "_DotPlot.pdf")), plot = dp,
         width = 3 + length(markers_present) * 0.6, height = 5)

  message("Saved plots for: ", ct)
}

## -----------------------------------------------------------------------
## MANUAL STEP (matches reference pipeline exactly): inspect the DotPlot
## and per-cell-type FeaturePlots/DotPlots saved above, then fill in the
## cluster -> cell type mapping below using the actual seurat_clusters IDs
## produced by script 02. The reference Annotation.Rmd leaves this same
## step for manual fill-in after visual review — we do not invent cluster
## assignments here. Until this table is filled in, Raw_cell_type falls
## back to the raw cluster number as a placeholder so the pipeline still
## runs end-to-end and downstream scripts have a well-formed column.
## -----------------------------------------------------------------------
FinalCluster <- as.character(levels(factor(scrna$seurat_clusters)))
CellType <- data.frame(cluster = FinalCluster, cell = FinalCluster)

# ASSUMPTION / TODO: update the mappings below after inspecting
# VehicleCisplatin_ManualAnnotation_DotPlot.pdf and MarkerFeaturePlots/.
# Example (uncomment and edit once clusters are reviewed):
# CellType[CellType$cluster %in% c("0","2","5"), 2] <- "PT"
# CellType[CellType$cluster %in% c("1"),         2] <- "TAL"
# CellType[CellType$cluster %in% c("3"),         2] <- "DCT"
# ... etc.

scrna@meta.data$Raw_cell_type <- sapply(
  scrna@meta.data$seurat_clusters,
  function(x) { CellType[CellType$cluster %in% x, 2] }
)
scrna@meta.data$Raw_cell_type_Confidence <- 1  # deterministic manual call, matches reference

cat("Raw_cell_type distribution (placeholder = cluster number until manually reviewed):\n")
print(table(scrna@meta.data$Raw_cell_type))

Idents(scrna) <- "Raw_cell_type"
p_manual <- DimPlot(scrna, split.by = "DataSet", label = TRUE, ncol = 2) +
  ggtitle("Manual annotation (placeholder pending cluster review)")
ggsave(paste0(OutDir, "VehicleCisplatin_ManualAnnotation_DimPlot.pdf"), plot = p_manual, height = 5, width = 12)

cat("Saving intermediate (post-manual-annotation) object...\n")
saveRDS(scrna, file = paste0(OutDir, "VehicleCisplatin_Raw_Predict_MKA_Lake_Customed_Annotation.RDS"))

scrna <- readRDS(paste0(OutDir, "VehicleCisplatin_Raw_Predict_MKA_Lake_Customed_Annotation.RDS"))
## ==========================================================================
## Step 3: Export to h5ad for external scANVI annotation
## ==========================================================================
cat("\n=== Step 3: Export h5ad for scANVI (MKA + Lake references) ===\n")

# Seurat v5 keeps per-sample split layers; as.SingleCellExperiment() calls
# GetAssayData() internally, which errors on multiple layers. JoinLayers()
# collapses them back into single counts/data matrices first — same fix
# used in the reference Annotation.Rmd export step.
scrna <- JoinLayers(scrna)
cat("Layers after JoinLayers:", paste(Layers(scrna[["RNA"]]), collapse = ", "), "\n")

sce <- as.SingleCellExperiment(scrna)
h5ad_path <- paste0(OutDir, "VehicleCisplatin_Raw_Predict_MKA_Lake_Annotation.h5ad")
writeH5AD(sce, file = h5ad_path)
cat("h5ad written to:", h5ad_path, "\n")
cat("Next step (EXTERNAL, outside R): submit\n",
    "  sbatch NewVersion/03b_VehicleCisplatin_scvi_annotation_MKAOnly.sh\n",
    "  sbatch NewVersion/03c_VehicleCisplatin_scvi_annotation_LakeOnly.sh\n",
    "then re-run this script to read back their CSV output from:", OutDir, "\n")

## ==========================================================================
## Step 4: Read back scANVI predictions (run this script again once the
## external scANVI CSVs exist)
## ==========================================================================
mka_csv  <- paste0(OutDir, "VehicleCisplatin_scanvi_MKA_annotations.csv")
lake_csv <- paste0(OutDir, "VehicleCisplatin_scanvi_Lake_annotations.csv")

if (file.exists(mka_csv) && file.exists(lake_csv)) {
  cat("\n=== Step 4: Reading back scANVI predictions ===\n")

  scanvi_annotation_mka <- read.csv(mka_csv) %>%
    mutate(CellID = str_replace(X, "-query", ""))
  rownames(scanvi_annotation_mka) <- scanvi_annotation_mka$CellID
  scanvi_annotation_mka <- scanvi_annotation_mka %>% dplyr::select(-X, -CellID)
  scrna <- AddMetaData(scrna, metadata = scanvi_annotation_mka)

  scanvi_annotation_lake <- read.csv(lake_csv) %>%
    mutate(CellID = str_replace(X, "-query", ""))
  rownames(scanvi_annotation_lake) <- scanvi_annotation_lake$CellID
  scanvi_annotation_lake <- scanvi_annotation_lake %>% dplyr::select(-X, -CellID)
  scrna <- AddMetaData(scrna, metadata = scanvi_annotation_lake)

  cat("scanvi_MKA_label distribution:\n");  print(table(scrna$scanvi_MKA_label,  useNA = "always"))
  cat("scanvi_Lake_label distribution:\n"); print(table(scrna$scanvi_Lake_label, useNA = "always"))

  p_scanvi_mka  <- DimPlot(scrna, group.by = "scanvi_MKA_label",  split.by = "DataSet", label = TRUE, repel = TRUE) +
    ggtitle("scanvi MKA transfer") + scale_color_hue()
  p_scanvi_lake <- DimPlot(scrna, group.by = "scanvi_Lake_label", split.by = "DataSet", label = TRUE, repel = TRUE) +
    ggtitle("scanvi Lake transfer") + scale_color_hue()
  ggsave(paste0(OutDir, "VehicleCisplatin_scanvi_MKA_Annotation_DimPlot.pdf"),  plot = p_scanvi_mka,  height = 5, width = 14)
  ggsave(paste0(OutDir, "VehicleCisplatin_scanvi_Lake_Annotation_DimPlot.pdf"), plot = p_scanvi_lake, height = 5, width = 14)

  cat("Saving final combined (5-method) annotation object...\n")
  saveRDS(scrna, file = paste0(OutDir, "VehicleCisplatin_Combined_Predict_MKA_Lake_scanvi_Annotation.RDS"))
  cat("\nDone. Output:", paste0(OutDir, "VehicleCisplatin_Combined_Predict_MKA_Lake_scanvi_Annotation.RDS"), "\n")
} else {
  cat("\n=== Step 4 SKIPPED: scANVI CSVs not found yet ===\n")
  cat("Run the external scANVI step described above, then re-run this script to\n",
      "produce VehicleCisplatin_Combined_Predict_MKA_Lake_scanvi_Annotation.RDS.\n")
  cat("(Intermediate object with MKA/Lake transfer + manual annotation was already\n",
      "saved to VehicleCisplatin_Raw_Predict_MKA_Lake_Customed_Annotation.RDS.)\n")
}
