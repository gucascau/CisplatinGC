#!/usr/bin/env Rscript
# =============================================================================
# Cisplatin_GC_Vehicle_Annotation.R
# Author: Xin Wang
# Date:   2026-04-14
# Description:
#   Cell type annotation for the 4-condition (Vehicle, Vehicle-GC, Cisplatin,
#   Cisplatin-GC) integrated single-cell dataset. Three annotation methods:
#     1. Seurat label transfer from Mouse Kidney Atlas (MKA) reference
#     2. Seurat label transfer from Lake snRNA reference
#     3. Manual annotation guided by canonical kidney cell-type markers
#   The annotated Seurat object and an h5ad file are saved for downstream use.
#   (scVI/scANVI run separately in Cisplatin_GC_Vehicle_scvi_annotation.ipynb)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(devtools)
  library(sctransform)
  library(ggthemes)
  library(Matrix)
  library(limma)
  library(tidyverse)
  library(biomaRt)
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
  library(DoubletFinder)
  library(reticulate)
  library(zellkonverter)
})

options(future.globals.maxSize = 1e15)

# =============================================================================
# Directories
# =============================================================================
IntegratedDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/IntegratedFourConditions/"
RefDir        <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Datasets/Reference/"
OutDir        <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
setwd(OutDir)

# =============================================================================
# Read the integrated Seurat object
# =============================================================================
scrna <- readRDS(paste0(IntegratedDir, "Cisplatin_GC_Vehicle_singlecell_doublet_harmony_v04082026.RDS"))

scrna@meta.data$DataSet <- ordered(
  factor(scrna@meta.data$DataSet),
  levels = c("Vehicle", "Vehicle-GC", "Cisplatin", "Cisplatin-GC")
)
scrna@meta.data$DataSet %>% table()
scrna@meta.data %>% head()

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
# Step 1a: Seurat label transfer — Mouse Kidney Atlas (MKA)
# =============================================================================
reference_mka <- readRDS(paste0(RefDir, "MKA_seurat.rds"))
reference_mka <- UpdateSeuratObject(reference_mka)
cat("MKA default assay before conversion:", DefaultAssay(reference_mka), "\n")

reference_mka <- ensembl_to_symbol_seurat(reference_mka)
cat("MKA default assay after conversion:", DefaultAssay(reference_mka), "\n")
cat("MKA sample gene names:", head(rownames(reference_mka), 5), "\n")

DefaultAssay(scrna) <- "RNA"

reference_mka <- NormalizeData(reference_mka)
reference_mka <- FindVariableFeatures(reference_mka, selection.method = "vst", nfeatures = 4000)
reference_mka <- ScaleData(reference_mka)

scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 4000)

features_mka <- intersect(VariableFeatures(reference_mka), VariableFeatures(scrna))
cat("MKA shared variable features:", length(features_mka), "\n")

reference_mka@meta.data %>% head()

anchors_mka <- FindTransferAnchors(
  reference = reference_mka,
  query     = scrna,
  features  = features_mka,
  dims      = 1:30
)

predictions_mka <- TransferData(
  anchorset = anchors_mka,
  refdata   = reference_mka$author_cell_type,
  dims      = 1:30
)
scrna <- AddMetaData(scrna, metadata = predictions_mka)
scrna@meta.data %>% head()

# =============================================================================
# Step 1b: Seurat label transfer — Lake snRNA reference
# =============================================================================
reference_lake <- readRDS(paste0(RefDir, "LakesnRNA_seurat.rds"))
reference_lake <- UpdateSeuratObject(reference_lake)
cat("Lake default assay before conversion:", DefaultAssay(reference_lake), "\n")

reference_lake <- ensembl_to_symbol_seurat(reference_lake)
cat("Lake default assay after conversion:", DefaultAssay(reference_lake), "\n")
cat("Lake sample gene names:", head(rownames(reference_lake), 5), "\n")

reference_lake <- NormalizeData(reference_lake)
reference_lake <- FindVariableFeatures(reference_lake, selection.method = "vst", nfeatures = 4000)
reference_lake <- ScaleData(reference_lake)

features_lake <- intersect(VariableFeatures(reference_lake), VariableFeatures(scrna))
cat("LakesnRNA shared variable features:", length(features_lake), "\n")

cat("Lake reference metadata columns:\n")
print(colnames(reference_lake@meta.data))

lake_label_col <- "SubclassLevel1"

reference_lake@meta.data$SubclassLevel1 %>% table()

anchors_lake <- FindTransferAnchors(
  reference = reference_lake,
  query     = scrna,
  features  = features_lake,
  dims      = 1:30
)

predictions_lake <- TransferData(
  anchorset = anchors_lake,
  refdata   = reference_lake@meta.data[[lake_label_col]],
  dims      = 1:30
)

colnames(predictions_lake) <- sub("^predicted\\.id$",          "predicted.id.Lake",         colnames(predictions_lake))
colnames(predictions_lake) <- sub("^prediction\\.score\\.max$", "prediction.score.max.Lake", colnames(predictions_lake))
colnames(predictions_lake) <- sub("^prediction\\.score\\.(.+)$","prediction.score.Lake.\\1", colnames(predictions_lake))

scrna <- AddMetaData(scrna, metadata = predictions_lake)
scrna@meta.data %>% head()

# =============================================================================
# Step 1c: Compare MKA vs Lake transfer predictions
# =============================================================================
table(MKA  = scrna$predicted.id,
      Lake = scrna$predicted.id.Lake) %>% head(10)

p_mka  <- DimPlot(scrna, group.by = "predicted.id",      label = TRUE, repel = TRUE) + ggtitle("MKA transfer")
p_lake <- DimPlot(scrna, group.by = "predicted.id.Lake", label = TRUE, repel = TRUE) + ggtitle("Lake transfer")
ggsave("Reference_Comparison_DimPlot.pdf", plot = p_mka + p_lake, height = 5, width = 25)

saveRDS(scrna,
        file = paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.RDS"))

# =============================================================================
# Step 2: Subclustering for manual annotation
# =============================================================================
Idents(scrna) <- "seurat_clusters"

scrna <- FindSubCluster(scrna, cluster = c("6"), graph.name = "RNA_snn",
                        subcluster.name = "Sub1", resolution = 0.3)
Idents(scrna) <- "Sub1"
scrna <- FindSubCluster(scrna, cluster = c("3"), graph.name = "RNA_snn",
                        subcluster.name = "Sub2", resolution = 0.2)
Idents(scrna) <- "Sub2"
ggsave("Subclustering_DimPlot.pdf",
       plot = DimPlot(scrna, split.by = "DataSet", label = TRUE, ncol = 2),
       height = 5, width = 14)

# =============================================================================
# Step 3: Manual cell-type assignment
# =============================================================================
Finalcelltypemarkers <- c(
  "Lrp2","Slc34a1","Slc13a3",          # Proximal tubule (PT)
  "Slc12a3","Pvalb","Wnk1",            # Distal convoluted tubule (DCT)
  "Umod","Slc12a1","Cldn10",           # Thick ascending limb (TAL)
  "Klk1","Slc8a1","Calb1",             # Connecting tubule (CNT)
  "Slc14a2","Fst","Bst1",              # Thin descending limb (DTL/LOH)
  "Aqp2","Hsd11b2","Scnn1g",           # Principal cells (PC)
  "Atp6v1g3","Aqp6","Slc26a7",         # Intercalated cells (IC)
  "Wt1","Nphs1","Nphs2",               # Podocyte (POD)
  "Cdh5","Igfbp3","Pecam1",            # Endothelium (END)
  "Col3a1","Fbn1","Lum","Col1a1",      # Fibroblast (FIB)
  "Ireb2","Alas2",                     # Plasmacytoid dendritic cells
  "S100a8","S100a9","Il1b",            # Neutrophil
  "C1qa","C1qb","Aif1",               # Macrophage
  "Cxcr6","Cd247","Cd3e",             # T cells
  "Igkc","Cd79a","Cd79b",             # B cells
  "Ccl5","Nkg7","Cd7"                 # NK/CD8 cells
)

ggsave("Cisplatin_GC_Vehicle_MarkerDotPlot_PreAnnotation.pdf",
       plot = DotPlot(scrna, features = Finalcelltypemarkers) + RotatedAxis(),
       height = 5, width = 14)

FinalCluster <- as.character(levels(factor(scrna$Sub2)))
CellType <- data.frame(cluster = FinalCluster, cell = FinalCluster)

CellType[CellType$cluster %in% c("6_0","6_1"), 2] <- "TAL"
CellType[CellType$cluster %in% c("6_2"),        2] <- "DCT"
CellType[CellType$cluster %in% c("11"),         2] <- "LOH"
CellType[CellType$cluster %in% c("0","5","7","8","1","2","4","14"), 2] <- "PT"
CellType[CellType$cluster %in% c("12"),         2] <- "IC"
CellType[CellType$cluster %in% c("3_0","3_1"),  2] <- "PC"
CellType[CellType$cluster %in% c("3_2"),        2] <- "CNT"
CellType[CellType$cluster %in% c("15","9","10"),2] <- "END"
CellType[CellType$cluster %in% c("13"),         2] <- "IMM"

scrna@meta.data$Raw_cell_type <- sapply(
  scrna@meta.data$Sub2,
  function(x) { CellType[CellType$cluster %in% x, 2] }
)
scrna@meta.data$Raw_cell_type_Confidence <- 1

levels(factor(scrna@meta.data$Raw_cell_type))
Idents(scrna) <- "Raw_cell_type"
ggsave("Cisplatin_GC_Vehicle_ManualAnnotation_DimPlot.pdf",
       plot = DimPlot(scrna, split.by = "DataSet", label = TRUE, ncol = 2),
       height = 5, width = 14)

ManualDotPlot <- DotPlot(scrna, features = Finalcelltypemarkers) + RotatedAxis()
ggsave("Cisplatin_GC_Vehicle_ManualAnnotation_DotPlot.pdf",
       plot = ManualDotPlot, height = 5, width = 14)

# =============================================================================
# Save annotated Seurat object (intermediate — before scVI)
# =============================================================================
saveRDS(scrna,
        file = paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_Annotation.RDS"))

# =============================================================================
# Export to h5ad for scVI annotation
# =============================================================================
sce <- as.SingleCellExperiment(scrna)

writeH5AD(
  sce,
  file = paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_Annotation.h5ad")
)
cat("h5ad written to:", paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_Annotation.h5ad"), "\n")
cat("Next step: run Cisplatin_GC_Vehicle_scvi_annotation.ipynb\n")
