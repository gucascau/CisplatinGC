## ============================================================
## Mitochondrial-rate QC summary, by cell type (FinalAnnotation_HC) and DataSet
## Input : Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS
## Output: CSV summary table + VlnPlot/boxplot of percent.mt per cell type,
##         split by DataSet, so mito-high clusters can be reviewed before
##         deciding whether to exclude them.
## ============================================================

library(dplyr)
library(Seurat)
library(ggplot2)

InDir  <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/IntegrateL1/"
OutDir <- InDir

RDSFile <- paste0(InDir, "Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS")

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

meta <- obj@meta.data %>%
  select(FinalAnnotation_HC, DataSet, percent.mt, nFeature_RNA, nCount_RNA) %>%
  mutate(FinalAnnotation_HC = as.character(FinalAnnotation_HC),
         FinalAnnotation_HC = ifelse(is.na(FinalAnnotation_HC) | FinalAnnotation_HC == "",
                                      "Uncertain", FinalAnnotation_HC))

# ── Summary table: mito% + depth metrics per cell type x DataSet ─────────
MitoSummary <- meta %>%
  group_by(FinalAnnotation_HC, DataSet) %>%
  summarise(
    n_cells         = n(),
    mean_percent_mt = mean(percent.mt),
    median_percent_mt = median(percent.mt),
    mean_nFeature   = mean(nFeature_RNA),
    mean_nCount     = mean(nCount_RNA),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_percent_mt))

cat("\n=== Mito% summary by cell type x DataSet (highest first) ===\n")
print(MitoSummary, n = Inf)

write.csv(MitoSummary,
          file = file.path(OutDir, "FinalAnnotation_MitoRate_Summary.csv"),
          row.names = FALSE)

# ── VlnPlot: percent.mt per cell type, split by DataSet ──────────────────
p_vln <- VlnPlot(obj,
                  features = "percent.mt",
                  group.by = "FinalAnnotation_HC",
                  split.by = "DataSet",
                  pt.size  = 0) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Mitochondrial % by cell type, split by DataSet",
       x = NULL, y = "percent.mt")

ggsave(
  filename = file.path(OutDir, "FinalAnnotation_MitoRate_VlnPlot.pdf"),
  plot     = p_vln,
  width    = 12, height = 6
)

# ── Boxplot alternative (clearer for many groups/outliers) ───────────────
p_box <- ggplot(meta, aes(x = FinalAnnotation_HC, y = percent.mt, fill = DataSet)) +
  geom_boxplot(outlier.size = 0.3, position = position_dodge(width = 0.8)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(title = "Mitochondrial % by cell type, split by DataSet",
       x = NULL, y = "percent.mt")

ggsave(
  filename = file.path(OutDir, "FinalAnnotation_MitoRate_BoxPlot.pdf"),
  plot     = p_box,
  width    = 12, height = 6
)

cat("\nSaved:\n",
    " ", file.path(OutDir, "FinalAnnotation_MitoRate_Summary.csv"), "\n",
    " ", file.path(OutDir, "FinalAnnotation_MitoRate_VlnPlot.pdf"), "\n",
    " ", file.path(OutDir, "FinalAnnotation_MitoRate_BoxPlot.pdf"), "\n")
