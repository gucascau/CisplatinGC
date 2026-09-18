## ============================================================
## Mitochondrial-rate QC — Vehicle vs Cisplatin (clean, GC groups removed)
## Input : Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation_VehicleCisplatin_Clean.RDS
##         (produced by Annotation/Cisplatin_Subset_VehicleCisplatin_Clean.R)
## Output: CSV summary table (percent.mt / nFeature / nCount per cell type x
##         DataSet), Vln/BoxPlots of percent.mt, and percent.mt vs
##         nFeature_RNA / nCount_RNA scatter plots (colored by DataSet,
##         faceted by cell type) so mito-high populations can be judged as
##         technical (low depth/complexity too) vs biological (injury
##         signal, depth preserved) before deciding whether to exclude them.
## ============================================================

library(dplyr)
library(Seurat)
library(ggplot2)

InDir  <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/IntegrateL1/"
OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/IntegrateL1CisplatinVehicle/QC/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)

RDSFile <- paste0(InDir, "Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation_VehicleCisplatin_Clean.RDS")

cat("Reading RDS...\n")
obj <- readRDS(RDSFile)
cat("Done. Cells:", ncol(obj), "\n")

# percent.mt should already be in meta.data from upstream QC; recompute if
# it's missing for some reason
if (!"percent.mt" %in% colnames(obj@meta.data)) {
  cat("percent.mt not found in meta.data — recomputing from '^mt-' genes...\n")
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
}

# Drop unused DataSet factor levels (in case any Vehicle-GC/Cisplatin-GC
# levels are still lingering) so only real Vehicle/Cisplatin groups appear
obj@meta.data$DataSet <- droplevels(factor(obj@meta.data$DataSet))
cat("DataSet levels present:", paste(levels(obj@meta.data$DataSet), collapse = ", "), "\n")

meta <- obj@meta.data %>%
  dplyr::select(FinalAnnotation_HC, DataSet, percent.mt, nFeature_RNA, nCount_RNA) %>%
  mutate(FinalAnnotation_HC = as.character(FinalAnnotation_HC),
         FinalAnnotation_HC = ifelse(is.na(FinalAnnotation_HC) | FinalAnnotation_HC == "",
                                      "Uncertain", FinalAnnotation_HC))

# ── 0. Whole-sample summary (per DataSet only, all cell types pooled) ────
# Checks whether the mito%/depth differences seen per cell type reflect a
# whole-sample technical difference (capture depth, dissociation quality,
# sequencing batch) rather than something specific to individual clusters.
DataSetSummary <- meta %>%
  group_by(DataSet) %>%
  summarise(
    n_cells           = n(),
    mean_percent_mt   = mean(percent.mt),
    median_percent_mt = median(percent.mt),
    mean_nFeature     = mean(nFeature_RNA),
    median_nFeature   = median(nFeature_RNA),
    mean_nCount       = mean(nCount_RNA),
    median_nCount     = median(nCount_RNA),
    .groups = "drop"
  )

cat("\n=== Whole-sample summary (per DataSet, all cell types pooled) ===\n")
print(DataSetSummary, n = Inf)

write.csv(DataSetSummary,
          file = file.path(OutDir, "VehicleCisplatin_DataSet_Overall_Summary.csv"),
          row.names = FALSE)

# ── 1. Summary table: mito% + depth metrics per cell type x DataSet ──────
MitoSummary <- meta %>%
  group_by(FinalAnnotation_HC, DataSet) %>%
  summarise(
    n_cells           = n(),
    mean_percent_mt   = mean(percent.mt),
    median_percent_mt = median(percent.mt),
    mean_nFeature     = mean(nFeature_RNA),
    mean_nCount       = mean(nCount_RNA),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_percent_mt))

cat("\n=== Mito% summary by cell type x DataSet (highest first) ===\n")
print(MitoSummary, n = Inf)

write.csv(MitoSummary,
          file = file.path(OutDir, "VehicleCisplatin_MitoRate_Summary.csv"),
          row.names = FALSE)

# ── 2. VlnPlot: percent.mt per cell type, split by DataSet ───────────────
p_vln <- VlnPlot(obj,
                  features = "percent.mt",
                  group.by = "FinalAnnotation_HC",
                  split.by = "DataSet",
                  pt.size  = 0) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Mitochondrial % by cell type — Vehicle vs Cisplatin",
       x = NULL, y = "percent.mt")

ggsave(
  filename = file.path(OutDir, "VehicleCisplatin_MitoRate_VlnPlot.pdf"),
  plot     = p_vln,
  width    = 12, height = 6
)

# ── 3. BoxPlot alternative (clearer for many groups/outliers) ────────────
p_box <- ggplot(meta, aes(x = FinalAnnotation_HC, y = percent.mt, fill = DataSet)) +
  geom_boxplot(outlier.size = 0.3, position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Mitochondrial % by cell type — Vehicle vs Cisplatin",
       x = NULL, y = "percent.mt")

ggsave(
  filename = file.path(OutDir, "VehicleCisplatin_MitoRate_BoxPlot.pdf"),
  plot     = p_box,
  width    = 12, height = 6
)

# ── 4. Scatter: percent.mt vs nFeature_RNA / nCount_RNA ──────────────────
# Colored by DataSet, faceted by cell type: lets you see, within each cell
# type, whether high-mito cells in Cisplatin still retain depth/complexity
# (biology) or collapse together with low counts/features (technical).
n_types <- length(unique(meta$FinalAnnotation_HC))
facet_ncol <- min(5, n_types)
scatter_w  <- max(10, facet_ncol * 2.6)
scatter_h  <- max(6, ceiling(n_types / facet_ncol) * 2.6)

p_scatter_feature <- ggplot(meta, aes(x = nFeature_RNA, y = percent.mt, color = DataSet)) +
  geom_point(size = 0.3, alpha = 0.3) +
  facet_wrap(~ FinalAnnotation_HC, ncol = facet_ncol) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold", size = 8)) +
  labs(title = "percent.mt vs nFeature_RNA, by cell type",
       x = "nFeature_RNA", y = "percent.mt")

ggsave(
  filename = file.path(OutDir, "VehicleCisplatin_Mito_vs_nFeature_Scatter.pdf"),
  plot     = p_scatter_feature,
  width    = scatter_w, height = scatter_h
)

p_scatter_count <- ggplot(meta, aes(x = nCount_RNA, y = percent.mt, color = DataSet)) +
  geom_point(size = 0.3, alpha = 0.3) +
  facet_wrap(~ FinalAnnotation_HC, ncol = facet_ncol) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold", size = 8)) +
  labs(title = "percent.mt vs nCount_RNA, by cell type",
       x = "nCount_RNA", y = "percent.mt")

ggsave(
  filename = file.path(OutDir, "VehicleCisplatin_Mito_vs_nCount_Scatter.pdf"),
  plot     = p_scatter_count,
  width    = scatter_w, height = scatter_h
)

cat("\nSaved to:", OutDir, "\n",
    " - VehicleCisplatin_MitoRate_Summary.csv\n",
    " - VehicleCisplatin_MitoRate_VlnPlot.pdf\n",
    " - VehicleCisplatin_MitoRate_BoxPlot.pdf\n",
    " - VehicleCisplatin_Mito_vs_nFeature_Scatter.pdf\n",
    " - VehicleCisplatin_Mito_vs_nCount_Scatter.pdf\n")
