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

MKA_Order<- c(c("PTS1", "PTS2", "PTS3", "PTS3T2","PEC","DCT","DCT-CNT","CNT","ATL","CTAL","MTAL","DTL-ATL","DTL","LOH","Podo","ICA","ICB","CD-Trans","PC","Per","MC", "Endo","Asc-Vasa-Recta","Desc-Vasa-Recta","Vas-Efferens","Vas-Afferens", "Glom-Endo", "Fib","Neutro","Macro","T lymph","B lymph","NK","DC"))
DefaultAssay(scrna) <- "RNA"

# generate the Dimplot for the MKA and the Markers
DimPlotRefereMKA<- DimPlot(reference_mka, group.by = "author_cell_type", label = TRUE, repel = TRUE) + ggtitle("MKA reference author_cell_type")
ggsave("Cisplatin_GC_Vehicle_MKA_Reference_DimPlot.pdf", plot = DimPlotRefereMKA, height = 5, width = 9)

# generate the DotPlot for the MKA reference of marker genes
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
DefaultAssay(reference_mka) <- "RNA"
# marker genes shared between MKA and our dataset
markers_mka <- intersect(Finalcelltypemarkers, rownames(reference_mka))
DotPlotMarkerGenesMKD <- DotPlot(reference_mka, features = markers_mka) + RotatedAxis() + ggtitle("MKA reference marker genes")
ggsave("Cisplatin_GC_Vehicle_MKA_Reference_MarkerGenes_DotPlot.pdf", plot = DotPlotMarkerGenesMKD, height = 5, width = 9)

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

colnames(predictions_mka) <- sub("^predicted.id$",          "Seaurat_Transfer_Predicted_MKA",         colnames(predictions_mka))
colnames(predictions_mka) <- sub("^prediction\\.score\\.max$", "Seaurat_Transfer_Predicted_MKA.max", colnames(predictions_mka))
colnames(predictions_mka) <- sub("^prediction\\.score\\.(.+)$","Seaurat_Transfer_Predicted_MKA.\\1", colnames(predictions_mka))

scrna <- AddMetaData(scrna, metadata = predictions_mka)
scrna@meta.data %>% head()
unique(scrna@meta.data$Seaurat_Transfer_Predicted_MKA)


scrna@meta.data$Seaurat_Transfer_Predicted_MKA <- ordered(
  factor(scrna@meta.data$Seaurat_Transfer_Predicted_MKA),
  levels = MKA_Order
)


# =============================================================================
# Step 1b: Seurat label transfer — Lake snRNA reference
# =============================================================================
reference_lake <- readRDS(paste0(RefDir, "LakesnRNA_seurat.rds"))
reference_lake <- UpdateSeuratObject(reference_lake)
cat("Lake default assay before conversion:", DefaultAssay(reference_lake), "\n")

LakeL1Cell_Order<- c("PT","PEC","DCT","CNT","DTL","ATL","TAL","IC","PC","POD","EC","FIB","PapE","Ad","VSM/P","Myeloid","NEU","Lymphoid")
length(LakeL1Cell_Order)
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
reference_lake@meta.data$SubclassLevel1  %>% unique()
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

colnames(predictions_lake) <- sub("^predicted\\.id$",          "Seaurat_Transfer_Predicted_Lake",         colnames(predictions_lake))
colnames(predictions_lake) <- sub("^prediction\\.score\\.max$", "Seaurat_Transfer_Predicted_Lake.max", colnames(predictions_lake))
colnames(predictions_lake) <- sub("^prediction\\.score\\.(.+)$","Seaurat_Transfer_Predicted_Lake.\\1", colnames(predictions_lake))

scrna <- AddMetaData(scrna, metadata = predictions_lake)
scrna@meta.data %>% head()

scrna@meta.data$Seaurat_Transfer_Predicted_Lake <- ordered(
  factor(scrna@meta.data$Seaurat_Transfer_Predicted_Lake),
  levels = LakeL1Cell_Order
)
# =============================================================================
# Step 1c: Compare MKA vs Lake transfer predictions
# =============================================================================
table(MKA  = scrna$Seaurat_Transfer_Predicted_Lake,
      Lake = scrna$Seaurat_Transfer_Predicted_MKA) %>% head(10)

p_mka  <- DimPlot(scrna, group.by = "Seaurat_Transfer_Predicted_MKA",      label = TRUE, repel = TRUE) + ggtitle("MKA transfer")
p_lake <- DimPlot(scrna, group.by = "Seaurat_Transfer_Predicted_Lake", label = TRUE, repel = TRUE) + ggtitle("Lake transfer")
ggsave("Reference_Comparison_DimPlot.pdf", plot = p_mka + p_lake, height = 5, width = 25)

saveRDS(scrna,
        file = paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.RDS"))


# read the RDS
scrna <- readRDS(paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.RDS"))


# =============================================================================
# Step 2: Subclustering for manual annotation
# =============================================================================

scrna <- FindClusters(scrna,  resolution = 0.8,algorithm = 1)

#scrna <- FindClusters(scrna,  resolution = 2,algorithm = 1)

Idents(scrna) <- "seurat_clusters"
DimPlot(scrna, group.by = "seurat_clusters", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("Manual annotation")

#scrna <- FindSubCluster(scrna, cluster = c("5"), graph.name = "RNA_snn",
                        subcluster.name = "Sub1", resolution = 0.2)
#Idents(scrna) <- "Sub1"

scrna <- FindSubCluster(scrna, cluster = c("7"), graph.name = "RNA_snn",
                        subcluster.name = "Sub2", resolution = 0.6)
Idents(scrna) <- "Sub2"
# DimPlot(scrna, split.by = "DataSet", label = TRUE, ncol = 2)
ggsave("Subclustering_DimPlot.pdf",
       plot = DimPlot(scrna, split.by = "DataSet", label = TRUE, ncol = 2),
       height = 5, width = 14)
DimPlot(scrna, group.by = "Sub2", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("Manual annotation")
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
DotPlot(scrna, features = Finalcelltypemarkers) + RotatedAxis()
FinalCluster <- as.character(levels(factor(scrna$Sub2)))
CellType <- data.frame(cluster = FinalCluster, cell = FinalCluster)

CellType[CellType$cluster %in% c("7_0","7_2","7_3"), 2] <- "DCT"
CellType[CellType$cluster %in% c("7_1"),        2] <- "CNT"
CellType[CellType$cluster %in% c("4"),         2] <- "TAL"
CellType[CellType$cluster %in% c("0","1","2","3","6","10","17"), 2] <- "PT"
CellType[CellType$cluster %in% c("12"),         2] <- "IC"
CellType[CellType$cluster %in% c("5"),  2] <- "PC"
CellType[CellType$cluster %in% c("9"),        2] <- "Myeloid"
CellType[CellType$cluster %in% c("8","13","16"),2] <- "END"
CellType[CellType$cluster %in% c("15"),         2] <- "FIB"
CellType[CellType$cluster %in% c("11"),         2] <- "Lymphoid"

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
ggsave("Cisplatin_GC_Vehicle_ManualAnnotation_DotPlots.pdf",
       plot = ManualDotPlot, height = 5, width = 14)


scrna <- FindClusters(scrna,  resolution = 2,algorithm = 1)

Idents(scrna) <- "seurat_clusters"

ggsave("Cisplatin_GC_Vehicle_ManualAnnotation_DimPlot_26clusters.pdf",
       plot = DimPlot(scrna, split.by = "DataSet", label = TRUE, ncol = 2),
       height = 9, width = 14)
ManualDotPlot <- DotPlot(scrna, features = Finalcelltypemarkers) + RotatedAxis()
ggsave("Cisplatin_GC_Vehicle_ManualAnnotation_DotPlot_26clusters.pdf",
       plot = ManualDotPlot, height = 10, width = 14)
# =============================================================================
# Save annotated Seurat object (intermediate — before scVI)
# =============================================================================
saveRDS(scrna,
        file = paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.RDS"))


# =============================================================================
# Export to annotation from scVI annotation
# =============================================================================
AnnotationDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/"


# =============================================================================
scrna <- readRDS(paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.RDS"))
# 1. Export to annotation from scVI annotation
# read the scANVI annotation CSV, this trained with the MKA reference and Lake together:
scanvi_annotation<- read.csv(paste0(AnnotationDir, "Cisplatin_GC_Vehicle_scanvi_annotations.csv"))

# add the CellID from X and remove "-query" string
scanvi_annotation <- scanvi_annotation %>%
  mutate(CellID = str_replace(X, "-query", ""))

# add the row name as the CellID and remove the X column
rownames(scanvi_annotation) <- scanvi_annotation$CellID
scanvi_annotation <- scanvi_annotation %>%
  select(-X, -CellID)
# add the metadata to the Seurat object
scrna <- AddMetaData(scrna,   metadata = scanvi_annotation)

# check the new metadata columns
scrna@meta.data %>% tail()
# 
scrna@meta.data$scanvi_MKA_label %>% table()
unique(scrna@meta.data$scanvi_MKA_label ) %>% length()

length(MKA_Order)
# set up the levels order of scanvi_MKA_label
scrna@meta.data$scanvi_MKA_label <- ordered(
  factor(scrna@meta.data$scanvi_MKA_label),
  levels = MKA_Order
)

scrna@meta.data$scanvi_Lake_label <- ordered(
  factor(scrna@meta.data$scanvi_Lake_label),
  levels = LakeL1Cell_Order
)


# DimPlot(scrna, group.by = "scanvi_MKA_label", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("scanvi_MKA transfer")

# DimPlot(scrna, group.by = "scanvi_Lake_label", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("scanvi_Lake transfer")
# =============================================================================

# =============================================================================
# Export to annotation from scVI annotation, we finally used the individual annoation by two datasets, not the combined one, because the combined one has some weird annotation that does not match with the MKA and Lake annotations, and we are not sure if it is a bug in the scVI annotation or if it is a real biological difference. So we will use the individual annotations for MKA and Lake separately to compare with our Seurat label transfer and manual annotation, and to see if there are any discrepancies or additional insights from the scVI annotations. We will also compare the scVI annotations with the MKA and Lake predictions to see if they are consistent or if there are any differences that we need to investigate further.
# =============================================================================

# 2 . Export to annotation from scVI annotation
# read the scANVI annotation CSV, this trained only with Lake reference:
scanvi_annotation_lake<- read.csv(paste0(AnnotationDir, "Cisplatin_GC_Vehicle_scanvi_Lake_annotations.csv"))
# add the CellID from X and remove "-query" string
scanvi_annotation_lake <- scanvi_annotation_lake %>%
  mutate(CellID = str_replace(X, "-query", ""))
scanvi_annotation_lake %>% head()
# add the row name as the CellID and remove the X column
rownames(scanvi_annotation_lake) <- scanvi_annotation_lake$CellID
scanvi_annotation_lake <- scanvi_annotation_lake %>%
  dplyr::select(-X)

head(scanvi_annotation_lake)

# add the metadata to the Seurat object
scrna <- AddMetaData(scrna,   metadata = scanvi_annotation_lake)
# check the new metadata columns
scrna@meta.data %>% tail()
# change the order of scanvi_Lake_label
scrna@meta.data$scanvi_Lake_label <- ordered(
  factor(scrna@meta.data$scanvi_Lake_label),
  levels = LakeL1Cell_Order
)

CisplatinLakeAnnotationDimplot<- DimPlot(scrna, group.by = "scanvi_Lake_label", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("scanvi_Lake transfer") + scale_color_hue()
ggsave("Cisplatin_GC_Vehicle_scanvi_Lake_Annotation_DimPlot.pdf",
       plot = CisplatinLakeAnnotationDimplot, height = 5, width = 14)


# 3 . Export to annotation from scVI annotation
# read the scANVI annotation CSV, this trained only with MKA reference:
scanvi_annotation_MKA<- read.csv(paste0(AnnotationDir, "Cisplatin_GC_Vehicle_scanvi_MKA_annotations.csv"))
# add the CellID from X and remove "-query" string
scanvi_annotation_MKA <- scanvi_annotation_MKA %>%
  mutate(CellID = str_replace(X, "-query", ""))
scanvi_annotation_MKA %>% head()
# add the row name as the CellID and remove the X column
rownames(scanvi_annotation_MKA) <- scanvi_annotation_MKA$CellID
scanvi_annotation_MKA <- scanvi_annotation_MKA %>%
  dplyr::select(-X)
# add the metadata to the Seurat object
scrna <- AddMetaData(scrna,   metadata = scanvi_annotation_MKA)
# change the order of scanvi_Lake_label
scrna@meta.data$scanvi_MKA_label <- ordered(
  factor(scrna@meta.data$scanvi_MKA_label),
  levels = MKA_Order
)

# check the new metadata columns
scrna@meta.data %>% tail()
CisplatinMKAAnnotationDimplot<- DimPlot(scrna, group.by = "scanvi_MKA_label", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("scanvi_MKA transfer") + scale_color_hue()
ggsave("Cisplatin_GC_Vehicle_scanvi_MKA_Annotation_DimPlot.pdf",
       plot = CisplatinMKAAnnotationDimplot, height = 5, width = 14)

# save all the svanvi annotations in the Seurat object
saveRDS(scrna,
        file = paste0(OutDir, "Cisplatin_GC_Vehicle_Combined_Predict_MKA_Lake_scanvi_Annotation.RDS"))


# =============================================================================
#  Comprehensive sub cluster annatation
# =============================================================================

# Method 1: Subset each major cell type and re-run clustering and marker identification to find subtypes
#  then Use known subtype markers to manually annotate subtypes within each major cell type

# Method 2: Use the subcluster labels from the previous step and cross-reference with the MKA and Lake annotations to assign more specific cell type labels to each subcluster. This can be done by looking at the overlap between the subclusters and the predicted labels from MKA and Lake, as well as checking the expression of known marker genes for specific cell types within each subcluster.

# We have the subcluter annotations
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

lake_label_col <- "SubclassLevel2"

reference_lake@meta.data$SubclassLevel2 %>% table()
#DimPlot(reference_lake, group.by = "SubclassLevel2", label = TRUE, repel = TRUE) + ggtitle("Lake reference SubclassLevel2")

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
# we need to change the column names to avoid overwriting the previous Lake predictions
predictions_lake_2 <- predictions_lake

colnames(predictions_lake)
colnames(predictions_lake_2) <- sub("^predicted.id$",          "Seaurat_Transfer_Predicted_LakeL2",         colnames(predictions_lake_2))
colnames(predictions_lake_2) <- sub("^prediction\\.score\\.max$", "Seaurat_Transfer_Predicted_LakeL2.max", colnames(predictions_lake_2))
colnames(predictions_lake_2) <- sub("^prediction\\.score\\.(.+)$","Seaurat_Transfer_Predicted_LakeL2.\\1", colnames(predictions_lake_2))
colnames(predictions_lake_2)
scrna <- AddMetaData(scrna, metadata = predictions_lake_2)
scrna@meta.data %>% head()
LakeL2Cell_Order <- c("PT-S1",  "PT-S1/2" , "PT-S2", "PT-S3","aPT" ,"dPT" ,  "cycPT", "PapE", "DTL1", "DTL2", "DTL3","aDTL", "DCT","aDCT",  "dDCT",  "dCNT", "ATL", "aTAL","C-TAL" ,"C/M-TAL", "dC-TAL",  "M-TAL", "frTAL" , "dM-TAL", "infTAL", "M-PC" , "dPC" , "C-IC-A" , "IC-B", "M-IC-A" , "EC-AA", "EC-AVR", "EC-DVR", "EC-GC","cycEC","infEC-PTC" ,"dEC-PTC"  ,"dEC-AVR" , "FIB", "dM-FIB", "infFIB", "eaFIB" , "M-FIB",  "VSMC/P", "IMCD","MC", "MON" , "moMAC", "moMAC-INF", "DC","T" ,  "B"  )

# change the order of scanvi_Lake_label
scrna@meta.data$Seaurat_Transfer_Predicted_LakeL2 <- ordered(
  factor(scrna@meta.data$Seaurat_Transfer_Predicted_LakeL2),
  levels = LakeL2Cell_Order
)


# =============================================================================
# Step 1c: Compare MKA vs Lake vs Lake L2 transfer predictions
# =============================================================================
table(MKA  = scrna$Seaurat_Transfer_Predicted_MKA,
      LakeL1 = scrna$Seaurat_Transfer_Predicted_Lake,
      LakeL2 = scrna$Seaurat_Transfer_Predicted_LakeL2) %>% head(10)


p_mka  <- DimPlot(scrna, group.by = "Seaurat_Transfer_Predicted_MKA",      label = TRUE, repel = TRUE) + ggtitle("MKA transfer") + scale_color_hue()
p_lake  <- DimPlot(scrna, group.by = "Seaurat_Transfer_Predicted_Lake", label = TRUE, repel = TRUE) + ggtitle("Lake transfer") + scale_color_hue()
p_lakeL2 <- DimPlot(scrna, group.by = "Seaurat_Transfer_Predicted_LakeL2", label = TRUE, repel = TRUE) + ggtitle("Lake L2 transfer") + scale_color_hue()
ggsave("Reference_Comparison_LabelTransfer_DimPlot.pdf", plot = p_mka + p_lake + p_lakeL2, height = 5, width = 30)

saveRDS(scrna,
        file = paste0(OutDir, "Cisplatin_GC_Vehicle_Combined_Predict_MKA_LakeWithL2_scanvi_Annotation.RDS"))


# Method 3: Use the scVI annotations to further refine the cell type labels, especially for ambiguous clusters. Compare the scVI predictions with the MKA and Lake predictions, and use the scVI annotations to resolve any discrepancies or to provide additional granularity in cell type classification.
# we need to run Cisplatin_GC_Vehicle_scvi_annotation_LakeOnly_ComprehensiveAnnotation.py first to get the scVI annotations, then we can add them to the Seurat object and compare with the MKA and Lake predictions to refine our cell type labels. This will help us to identify any ambiguous clusters and assign more accurate cell type labels based on the combined evidence from all three annotation methods.

#scrna<- readRDS ("Cisplatin_GC_Vehicle_Combined_Predict_MKA_LakeWithL2_scanvi_Annotation.RDS")


# Method 3 . Export to annotation from scVI annotation
# read the scANVI annotation CSV, this trained only with Lake reference:
scanvi_annotation_lakeL2<- read.csv(paste0(AnnotationDir, "Cisplatin_GC_Vehicle_scanvi_LakeL2_annotations.csv"))
# add the CellID from X and remove "-query" string
scanvi_annotation_lakeL2 <- scanvi_annotation_lakeL2 %>%
  mutate(CellID = str_replace(X, "-query", ""))
scanvi_annotation_lakeL2 %>% head()
# add the row name as the CellID and remove the X column
rownames(scanvi_annotation_lakeL2) <- scanvi_annotation_lakeL2$CellID
scanvi_annotation_lakeL2 <- scanvi_annotation_lakeL2 %>%
  dplyr::select(-X)

head(scanvi_annotation_lakeL2)

# add the metadata to the Seurat object
scrna <- AddMetaData(scrna,   metadata = scanvi_annotation_lakeL2)
# check the new metadata columns
levels( factor(scrna@meta.data$scanvi_LakeL2_label))
# add the order of the LakeL2Cell_Order
LakeL2Cell_Order <- c("PT-S1",  "PT-S1/2" , "PT-S2", "PT-S3","aPT" ,"dPT" ,  "cycPT", "PapE", "DTL1", "DTL2", "DTL3","aDTL", "DCT","aDCT",  "dDCT",  "dCNT", "ATL", "aTAL","C-TAL" ,"C/M-TAL", "dC-TAL",  "M-TAL", "frTAL" , "dM-TAL", "infTAL", "M-PC" , "dPC" , "C-IC-A" , "IC-B", "M-IC-A" , "EC-AA", "EC-AVR", "EC-DVR", "EC-GC","cycEC","infEC-PTC" ,"dEC-PTC"  ,"dEC-AVR" , "FIB", "dM-FIB", "infFIB", "eaFIB" , "M-FIB",  "VSMC/P", "IMCD","MC", "MON" , "moMAC", "moMAC-INF", "DC","T" ,  "B"  )
length(LakeL2Cell_Order)


scrna@meta.data %>% tail()
# change the order of scanvi_Lake_label
scrna@meta.data$scanvi_LakeL2_label <- ordered(
  factor(scrna@meta.data$scanvi_LakeL2_label),
  levels = LakeL2Cell_Order
)

CisplatinLakeAnnotationL2Dimplot<- DimPlot(scrna, group.by = "scanvi_LakeL2_label", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("scanvi_Lake transfer") + scale_color_hue()
ggsave("Cisplatin_GC_Vehicle_scanvi_Lake_Annotation_scanvi_LakeL2_label_DimPlot.pdf",
       plot = CisplatinLakeAnnotationL2Dimplot, height = 5, width = 14)


saveRDS(scrna,
        file = paste0(OutDir, "Cisplatin_GC_Vehicle_Combined_Predict_MKA_LakeWithL2_scanvi_LakeL2scanvi_Annotation.RDS"))
# scrna <- readRDS(paste0(OutDir, "Cisplatin_GC_Vehicle_Combined_Predict_MKA_LakeWithL2_scanvi_LakeL2scanvi_Annotation.RDS"))

#DimPlot(scrna, group.by = "scanvi_LakeL2_label", split.by = "DataSet",label = TRUE, repel = TRUE) + ggtitle("scanvi_Lake transfer") + scale_color_hue()

# 
scrna@meta.data$Raw_cell_type %>% table()
scrna@meta.data %>% tail()



# =============================================================================
# Export to h5ad for scVI annotation
# =============================================================================

#srcna <- readRDS(paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_MKA_Lake_Annotation.RDS"))
sce <- as.SingleCellExperiment(scrna)

writeH5AD(
  sce,
  file = paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_Annotation.h5ad")
)
cat("h5ad written to:", paste0(OutDir, "Cisplatin_GC_Vehicle_Raw_Predict_Annotation.h5ad"), "\n")
cat("Next step: run Cisplatin_GC_Vehicle_scvi_annotation.ipynb\n")



