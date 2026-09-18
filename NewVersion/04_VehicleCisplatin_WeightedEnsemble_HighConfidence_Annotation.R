#!/usr/bin/env Rscript
## ============================================================
## 04_VehicleCisplatin_WeightedEnsemble_HighConfidence_Annotation.R
## Author: Xin Wang
## Email:  xin.wang@nationwidechildrens.org
##
## Description:
##   Stage 4 (final) of the Vehicle + Cisplatin pipeline. Combines the five
##   per-method L1 annotation columns produced by script 03 into a single
##   weighted-ensemble call, then derives a high-confidence final call,
##   following exactly the logic in:
##     Annotation/Cisplatin_GC_Vehicle_Annotation_Integration.Rmd   (L1 harmonisation)
##     Annotation/Cisplatin_GC_Vehicle_WeightedEnsemble_Annotation.Rmd (weighted ensemble)
##     Annotation/Cisplatin_GC_Vehicle_HighConfidence_Annotation.Rmd   (confidence tiers + FinalAnnotation_HC)
##
##   L1 harmonisation:
##     L1_SeuratMKA  = map_MKA_to_L1(Seaurat_Transfer_Predicted_MKA)   [MKA abbreviations -> shared L1 vocabulary]
##     L1_SeuratLake = Seaurat_Transfer_Predicted_Lake                  [already L1 — SubclassLevel1]
##     L1_Raw        = map_raw_to_L1(Raw_cell_type)                     [pass-through / validate against L1_levels]
##     L1_scanviMKA  = map_MKA_to_L1(scanvi_MKA_label)
##     L1_scanviLake = scanvi_Lake_label                                [already L1]
##
##   Weighted ensemble: each of the 5 methods gets equal weight (0.20,
##   "no published benchmark directly comparing all five methods on mouse
##   kidney data is available" — same rationale/weights as reference).
##   Each method's vote is weighted by method_weight x per-cell confidence
##   score (min-max normalised to [0,1] per method); highest total wins.
##
##   Confidence tiers (identical thresholds to reference):
##     High   : >=3 of 5 methods agree with the weighted winner AND
##              cluster coherence >= 0.70
##     Medium : >=2 methods agree OR cluster coherence >= 0.50
##     Low    : everything else -> falls back to cluster-majority weighted
##              label (safer, cluster-level call) for FinalAnnotation_HC
##     (cluster = seurat_clusters from script 02 — the reference used
##      "Sub2", a sub-clustering column produced by an interactive
##      FindSubCluster() review step upstream that has no equivalent in
##      this from-scratch 2-condition pipeline. ASSUMPTION: seurat_clusters
##      is used as the coherence/majority grouping variable instead.)
##
##   Mes/FIB collapse: both map to "FIB/VSM" in FinalAnnotation_HC, matching
##   the reference's stated rationale (these two labels are frequently
##   confused / low-confidence).
##
## Input : /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/
##           VehicleCisplatin_Combined_Predict_MKA_Lake_scanvi_Annotation.RDS   (from script 03)
## Output: /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/IntegrateL1/
##           VehicleCisplatin_L1Harmonised_Weighted_FinalAnnotation_VehicleCisplatin_Clean.RDS
##           CSV summaries + DimPlots
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
})

options(future.globals.maxSize = 1e15)

## ---- Directories -------------------------------------------------------
InDir  <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/"
OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/IntegrateL1/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)
setwd(OutDir)

RDSFile <- paste0(InDir, "VehicleCisplatin_Combined_Predict_MKA_Lake_scanvi_Annotation.RDS")

cat("Reading:", RDSFile, "\n")
obj <- readRDS(RDSFile)
DefaultAssay(obj) <- "RNA"
obj@meta.data$DataSet <- droplevels(factor(obj@meta.data$DataSet))

required_cols <- c("Seaurat_Transfer_Predicted_MKA", "Seaurat_Transfer_Predicted_Lake",
                    "Raw_cell_type", "scanvi_MKA_label", "scanvi_Lake_label")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing columns in metadata (did script 03's scANVI read-back run?): ",
       paste(missing_cols, collapse = ", "))
}

## ==========================================================================
## Step 1: Harmonise all five methods to a shared L1 vocabulary
## ==========================================================================
cat("\n=== Step 1: L1 harmonisation ===\n")

L1_levels <- c("PT","TAL","DCT","CNT","PC","IC","IC/CNT","POD","END","VSM/P","PapE","Myeloid","Lymphoid")

# Verbatim from Annotation/Cisplatin_GC_Vehicle_Annotation_Integration.Rmd
map_MKA_to_L1 <- function(x) {
  dplyr::case_when(
    x %in% c("PTS1","PTS2","PTS3","PTS3T2")                                      ~ "PT",
    x %in% c("CTAL","MTAL")                                                      ~ "TAL",
    x %in% c("DCT","DCT-CNT")                                                    ~ "DCT",
    x %in% c("CNT")                                                               ~ "CNT",
    x %in% c("ATL","DTL-ATL","DTL","LOH")                                        ~ "DTL",
    x %in% c("PC")                                                                ~ "PC",
    x %in% c("ICA","ICB")                                                         ~ "IC",
    x %in% c("CD-Trans")                                                          ~ "PC",
    x %in% c("Podo","PEC")                                                        ~ "POD",
    x %in% c("Endo","Asc-Vasa-Recta","Desc-Vasa-Recta",
              "Vas-Efferens","Vas-Afferens","Glom-Endo","Per")                   ~ "END",
    x %in% c("Fib")                                                               ~ "FIB",
    x %in% c("MC")                                                                ~ "Mes",
    x %in% c("Neutro","Macro")                                                   ~ "Myeloid",
    x %in% c("T lymph","B lymph","NK","DC")                                      ~ "Lymphoid",
    TRUE ~ "Uncertain"
  )
}

map_raw_to_L1 <- function(x) {
  ifelse(x %in% L1_levels, x, "Uncertain")
}

obj@meta.data <- obj@meta.data %>%
  mutate(
    L1_SeuratMKA  = map_MKA_to_L1(Seaurat_Transfer_Predicted_MKA),
    L1_SeuratLake = Seaurat_Transfer_Predicted_Lake,
    L1_Raw        = map_raw_to_L1(Raw_cell_type),
    L1_scanviMKA  = map_MKA_to_L1(scanvi_MKA_label),
    L1_scanviLake = scanvi_Lake_label
  )

for (col in c("L1_SeuratMKA","L1_SeuratLake","L1_Raw","L1_scanviMKA","L1_scanviLake")) {
  n_unc <- sum(obj@meta.data[[col]] == "Uncertain", na.rm = TRUE)
  cat(col, "-> Uncertain:", n_unc, "(", round(100 * n_unc / ncol(obj), 1), "%)\n")
}

saveRDS(obj, file = paste0(OutDir, "VehicleCisplatin_L1Harmonised_Annotation.RDS"))

## ==========================================================================
## Step 2: Weighted ensemble
## ==========================================================================
cat("\n=== Step 2: Weighted ensemble ===\n")

method_weights <- c(
  L1_SeuratMKA  = 0.20,
  L1_SeuratLake = 0.20,
  L1_Raw        = 0.20,
  L1_scanviMKA  = 0.20,
  L1_scanviLake = 0.20
)
stopifnot(abs(sum(method_weights) - 1) < 1e-9)

method_score_cols <- c(
  L1_SeuratMKA  = "Seaurat_Transfer_Predicted_MKA.max",
  L1_SeuratLake = "Seaurat_Transfer_Predicted_Lake.max",
  L1_Raw        = "Raw_cell_type_Confidence",
  L1_scanviMKA  = "scanvi_MKA_confidence",
  L1_scanviLake = "scanvi_Lake_confidence"
)
missing_score_cols <- setdiff(method_score_cols, colnames(obj@meta.data))
if (length(missing_score_cols) > 0) {
  stop("Missing confidence-score columns: ", paste(missing_score_cols, collapse = ", "))
}

normalize_01 <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  r <- range(x, na.rm = TRUE)
  if (diff(r) < 1e-10) return(rep(1.0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

meta <- obj@meta.data

all_labels <- sort(unique(na.omit(unlist(
  lapply(names(method_weights), function(m) {
    v <- meta[[m]]
    v[!is.na(v) & v != "Uncertain"]
  })
))))
cat("L1 labels in ensemble:", paste(all_labels, collapse = ", "), "\n")

n_cells      <- nrow(meta)
score_matrix <- matrix(0.0, nrow = n_cells, ncol = length(all_labels),
                        dimnames = list(rownames(meta), all_labels))

for (method in names(method_weights)) {
  labels <- meta[[method]]
  scores <- normalize_01(meta[[method_score_cols[method]]])
  w <- method_weights[method]
  weighted <- w * scores

  valid <- !is.na(labels) & labels != "Uncertain" & labels %in% all_labels
  rows  <- which(valid)
  cols  <- match(labels[valid], all_labels)
  lin_idx <- rows + (cols - 1L) * n_cells
  score_matrix[lin_idx] <- score_matrix[lin_idx] + weighted[rows]
}

row_max  <- apply(score_matrix, 1, max)
best_col <- max.col(score_matrix, ties.method = "first")
second_max <- apply(score_matrix, 1, function(r) {
  s <- sort(r, decreasing = TRUE)
  if (length(s) >= 2L) s[2L] else 0
})

obj@meta.data$WeightedAnnotation        <- ifelse(row_max == 0, "Uncertain", all_labels[best_col])
obj@meta.data$WeightedAnnotation_score  <- row_max
obj@meta.data$WeightedAnnotation_margin <- row_max - second_max

cat("Weighted ensemble annotation distribution:\n")
print(table(obj@meta.data$WeightedAnnotation, useNA = "always"))

## ==========================================================================
## Step 3: High-confidence tiers + FinalAnnotation_HC
## ==========================================================================
cat("\n=== Step 3: Confidence tiers ===\n")

method_label_cols <- names(method_weights)

agree_matrix <- sapply(method_label_cols, function(m) {
  as.integer(!is.na(obj@meta.data[[m]]) &
             obj@meta.data[[m]] != "Uncertain" &
             obj@meta.data[[m]] == obj@meta.data$WeightedAnnotation)
})
obj@meta.data$WeightedAnnotation_n_agree <- rowSums(agree_matrix)

cat("Method agreement count (0-5):\n")
print(table(obj@meta.data$WeightedAnnotation_n_agree))

# Cluster coherence — grouped by seurat_clusters (see header note: the
# reference used "Sub2", a manually-reviewed sub-clustering column that has
# no equivalent in this from-scratch pipeline).
coherence_tbl <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  mutate(
    ClusterCoherence = {
      valid_labels <- WeightedAnnotation[WeightedAnnotation != "Uncertain"]
      if (length(valid_labels) == 0) 0
      else max(table(valid_labels)) / n()
    }
  ) %>%
  ungroup() %>%
  pull(ClusterCoherence)
obj@meta.data$ClusterCoherence <- coherence_tbl

cat("Cluster coherence summary (per cell):\n")
print(summary(obj@meta.data$ClusterCoherence))

obj@meta.data$AnnotationTier <- dplyr::case_when(
  obj@meta.data$WeightedAnnotation == "Uncertain"                                              ~ "Uncertain",
  obj@meta.data$WeightedAnnotation_n_agree >= 3 & obj@meta.data$ClusterCoherence >= 0.70        ~ "High",
  obj@meta.data$WeightedAnnotation_n_agree >= 2 | obj@meta.data$ClusterCoherence >= 0.50        ~ "Medium",
  TRUE                                                                                           ~ "Low"
)

cat("Confidence tier distribution:\n")
print(table(obj@meta.data$AnnotationTier, useNA = "always"))

# Cluster-majority weighted label (fallback for "Low" tier cells)
cluster_majority <- obj@meta.data %>%
  filter(WeightedAnnotation != "Uncertain") %>%
  group_by(seurat_clusters, WeightedAnnotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  dplyr::rename(ClusterMajority_Weighted = WeightedAnnotation) %>%
  dplyr::select(seurat_clusters, ClusterMajority_Weighted)

obj@meta.data$ClusterMajority_Weighted <- cluster_majority$ClusterMajority_Weighted[
  match(obj@meta.data$seurat_clusters, cluster_majority$seurat_clusters)
]

obj@meta.data$FinalAnnotation_HC <- dplyr::case_when(
  obj@meta.data$AnnotationTier %in% c("High", "Medium") ~ obj@meta.data$WeightedAnnotation,
  obj@meta.data$AnnotationTier == "Low"                  ~ obj@meta.data$ClusterMajority_Weighted,
  TRUE                                                    ~ "Uncertain"
)

# Mes/FIB collapse — matches reference rationale (frequently confused / low confidence)
obj@meta.data$FinalAnnotation_HC <- dplyr::case_when(
  obj@meta.data$FinalAnnotation_HC %in% c("Mes","FIB") ~ "FIB/VSM",
  TRUE ~ obj@meta.data$FinalAnnotation_HC
)

cat("Final high-confidence annotation distribution:\n")
print(table(obj@meta.data$FinalAnnotation_HC, useNA = "always"))

## ==========================================================================
## Step 4: Summaries, plots, save
## ==========================================================================
cat("\n=== Step 4: Summaries + plots ===\n")

hc_cluster_summary <- obj@meta.data %>%
  group_by(seurat_clusters, FinalAnnotation_HC) %>%
  summarise(
    n_cells      = n(),
    pct_High     = round(100 * mean(AnnotationTier == "High"),   1),
    pct_Medium   = round(100 * mean(AnnotationTier == "Medium"), 1),
    pct_Low      = round(100 * mean(AnnotationTier == "Low"),    1),
    mean_score   = round(mean(WeightedAnnotation_score),   3),
    mean_n_agree = round(mean(WeightedAnnotation_n_agree), 2),
    .groups = "drop"
  ) %>%
  group_by(seurat_clusters) %>%
  slice_max(n_cells, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(hc_cluster_summary,
          file = paste0(OutDir, "VehicleCisplatin_Cluster_HighConfidence_Annotation_Summary.csv"),
          row.names = FALSE)

p_weighted <- DimPlot(obj, group.by = "WeightedAnnotation", split.by = "DataSet", label = TRUE) +
  scale_color_hue() + ggtitle("Weighted ensemble annotation") + theme(legend.position = "right")
ggsave(paste0(OutDir, "VehicleCisplatin_WeightedAnnotation_DimPlot.pdf"), plot = p_weighted, height = 5, width = 10)

p_tier <- DimPlot(obj, group.by = "AnnotationTier", split.by = "DataSet", label = FALSE) +
  scale_colour_manual(values = c(High = "#2ca02c", Medium = "#ff7f0e", Low = "#d62728", Uncertain = "grey80")) +
  ggtitle("Annotation confidence tier") + theme(legend.position = "right")
ggsave(paste0(OutDir, "VehicleCisplatin_AnnotationTier_DimPlot.pdf"), plot = p_tier, height = 5, width = 10)

p_hc <- DimPlot(obj, group.by = "FinalAnnotation_HC", split.by = "DataSet", label = TRUE) +
  scale_color_hue() + ggtitle("Final high-confidence annotation") + theme(legend.position = "right")
ggsave(paste0(OutDir, "VehicleCisplatin_FinalAnnotation_HC_DimPlot.pdf"), plot = p_hc, height = 5, width = 10)

OutRDSFile <- paste0(OutDir, "VehicleCisplatin_L1Harmonised_Weighted_FinalAnnotation_VehicleCisplatin_Clean.RDS")
cat("Saving final annotated object to:", OutRDSFile, "\n")
saveRDS(obj, file = OutRDSFile)

write.csv(
  obj@meta.data[, c(
    "DataSet", "seurat_clusters",
    "Raw_cell_type", "Seaurat_Transfer_Predicted_MKA", "Seaurat_Transfer_Predicted_Lake",
    "scanvi_MKA_label", "scanvi_Lake_label",
    "L1_Raw","L1_SeuratMKA","L1_SeuratLake","L1_scanviMKA","L1_scanviLake",
    "WeightedAnnotation","WeightedAnnotation_score","WeightedAnnotation_margin",
    "WeightedAnnotation_n_agree","ClusterCoherence","AnnotationTier","FinalAnnotation_HC"
  )],
  file = paste0(OutDir, "VehicleCisplatin_PerCell_HighConfidence_Annotation_Summary.csv")
)

cat("\nDone. Output:", OutRDSFile, "\n")
