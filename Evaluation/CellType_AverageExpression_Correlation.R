#!/usr/bin/env Rscript
# =============================================================================
# CellType_AverageExpression_Correlation.R
# Author: Xin Wang
# Date:   2026-05-18
# Description:
#   Evaluate final cell-type annotations by computing average expression per
#   cell type and Pearson correlating those profiles between:
#     1. Our dataset (FinalAnnotation_HC) vs MKA (author_cell_type)
#     2. Our dataset (FinalAnnotation_HC) vs Lake (author_cell_type)
#   Outputs heatmaps + CSV tables of the correlation matrices.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# =============================================================================
# Paths
# =============================================================================
RefDir  <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Datasets/Reference/"
AnnDir  <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/IntegrateL1/"
OutDir  <- "/vast0/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Evaluation/CellTypeCorrelation/"
dir.create(OutDir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Helper: Ensembl IDs -> gene symbols (same logic as annotation script)
# =============================================================================
ensembl_to_symbol_seurat <- function(obj) {
  default_assay <- DefaultAssay(obj)
  avail_layers  <- tryCatch(Layers(obj[[default_assay]]), error = function(e) character(0))

  preferred    <- c("counts", "data", "scale.data")
  layer_to_use <- preferred[preferred %in% avail_layers]
  if (length(layer_to_use) == 0) layer_to_use <- avail_layers
  layer_to_use <- layer_to_use[1]

  counts <- tryCatch(
    GetAssayData(obj, layer = layer_to_use),
    error = function(e) GetAssayData(obj, slot = layer_to_use)
  )

  genes <- rownames(counts)
  if (!grepl("^ENSMUS", genes[1])) {
    message("Gene names already look like symbols — skipping conversion.")
    return(obj)
  }
  message("Converting Ensembl IDs -> gene symbols ...")

  genes_clean <- sub("\\.[0-9]+$", "", genes)
  sym_map <- mapIds(org.Mm.eg.db, keys = genes_clean, column = "SYMBOL",
                    keytype = "ENSEMBL", multiVals = "first")
  sym_map[is.na(sym_map)] <- genes[is.na(sym_map)]
  sym_map <- make.unique(as.character(sym_map))
  rownames(counts) <- sym_map

  new_obj <- CreateSeuratObject(counts = counts, meta.data = obj@meta.data,
                                project = Project(obj))
  DefaultAssay(new_obj) <- "RNA"
  message(sprintf("  Renamed %d genes.", length(sym_map)))
  return(new_obj)
}

# =============================================================================
# Helper: compute average log-normalised expression per cell type
# =============================================================================
compute_avg_expr <- function(obj, group_col) {
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- group_col
  avg <- AverageExpression(obj, assays = "RNA", slot = "data", return.seurat = FALSE)
  as.matrix(avg$RNA)   # genes x cell-types
}

# =============================================================================
# Helper: correlation heatmap between two avg-expression matrices
# =============================================================================
plot_correlation_heatmap <- function(mat_query, mat_ref, title, outfile,
                                     query_label = "Our annotation",
                                     ref_label   = "Reference") {
  shared_genes <- intersect(rownames(mat_query), rownames(mat_ref))
  message(sprintf("Shared genes for '%s': %d", title, length(shared_genes)))

  m1 <- mat_query[shared_genes, , drop = FALSE]
  m2 <- mat_ref[shared_genes,   , drop = FALSE]

  # Pearson correlation: each column is a cell-type profile
  cor_mat <- cor(m1, m2, method = "pearson")   # query cell-types x ref cell-types

  # save CSV
  csv_path <- sub("\\.pdf$", ".csv", outfile)
  write.csv(cor_mat, file = csv_path)
  message("Correlation matrix saved: ", csv_path)

  # heatmap
  color_breaks <- seq(-1, 1, length.out = 101)
  color_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

  pheatmap(
    cor_mat,
    color            = color_palette,
    breaks           = color_breaks,
    main             = title,
    angle_col        = 45,
    fontsize_row     = 9,
    fontsize_col     = 9,
    border_color     = NA,
    clustering_method = "ward.D2",
    filename         = outfile,
    width            = max(7, ncol(cor_mat) * 0.35 + 3),
    height           = max(6, nrow(cor_mat) * 0.35 + 3)
  )
  message("Heatmap saved: ", outfile)

  invisible(cor_mat)
}

# =============================================================================
# 1. Load our final annotated dataset
# =============================================================================
message("Loading our final annotated Seurat object ...")
scrna <- readRDS(paste0(AnnDir, "Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS"))
scrna <- UpdateSeuratObject(scrna)

stopifnot("FinalAnnotation_HC" %in% colnames(scrna@meta.data))
message("Our cell types (FinalAnnotation_HC):")
print(table(scrna$FinalAnnotation_HC))

DefaultAssay(scrna) <- "RNA"
scrna <- NormalizeData(scrna, verbose = FALSE)

message("Computing average expression for our dataset ...")
avg_ours <- compute_avg_expr(scrna, "FinalAnnotation_HC")
message(sprintf("  Matrix: %d genes x %d cell types", nrow(avg_ours), ncol(avg_ours)))

# =============================================================================
# 2. Load MKA reference
# =============================================================================
message("Loading MKA reference ...")
ref_mka <- readRDS(paste0(RefDir, "MKA_seurat.rds"))
ref_mka  <- UpdateSeuratObject(ref_mka)
ref_mka  <- ensembl_to_symbol_seurat(ref_mka)

stopifnot("author_cell_type" %in% colnames(ref_mka@meta.data))
message("MKA cell types (author_cell_type):")
print(table(ref_mka$author_cell_type))

DefaultAssay(ref_mka) <- "RNA"
ref_mka <- NormalizeData(ref_mka, verbose = FALSE)

message("Computing average expression for MKA ...")
avg_mka <- compute_avg_expr(ref_mka, "author_cell_type")
message(sprintf("  Matrix: %d genes x %d cell types", nrow(avg_mka), ncol(avg_mka)))

# =============================================================================
# 3. Load Lake snRNA reference
# =============================================================================
message("Loading Lake snRNA reference ...")
ref_lake <- readRDS(paste0(RefDir, "LakesnRNA_seurat.rds"))
ref_lake  <- UpdateSeuratObject(ref_lake)
ref_lake  <- ensembl_to_symbol_seurat(ref_lake)

lake_label_col <- if ("author_cell_type" %in% colnames(ref_lake@meta.data)) {
  "author_cell_type"
} else if ("SubclassLevel1" %in% colnames(ref_lake@meta.data)) {
  "SubclassLevel1"
} else {
  stop("Neither 'author_cell_type' nor 'SubclassLevel1' found in Lake metadata. Columns: ",
       paste(colnames(ref_lake@meta.data), collapse = ", "))
}
message("Using Lake label column: ", lake_label_col)
message("Lake cell types (", lake_label_col, "):")
print(table(ref_lake@meta.data[[lake_label_col]]))

DefaultAssay(ref_lake) <- "RNA"
ref_lake <- NormalizeData(ref_lake, verbose = FALSE)

message("Computing average expression for Lake ...")
avg_lake <- compute_avg_expr(ref_lake, lake_label_col)
message(sprintf("  Matrix: %d genes x %d cell types", nrow(avg_lake), ncol(avg_lake)))

# =============================================================================
# 4. Correlation: Our dataset vs MKA
# =============================================================================
message("\n--- Correlation: Our annotation vs MKA ---")
cor_ours_mka <- plot_correlation_heatmap(
  mat_query   = avg_ours,
  mat_ref     = avg_mka,
  title       = "Pearson correlation: Our annotation vs MKA",
  outfile     = paste0(OutDir, "Correlation_OurAnnotation_vs_MKA.pdf"),
  query_label = "Our annotation (FinalAnnotation_HC)",
  ref_label   = "MKA (author_cell_type)"
)

# =============================================================================
# 5. Correlation: Our dataset vs Lake
# =============================================================================
message("\n--- Correlation: Our annotation vs Lake ---")
cor_ours_lake <- plot_correlation_heatmap(
  mat_query   = avg_ours,
  mat_ref     = avg_lake,
  title       = "Pearson correlation: Our annotation vs Lake",
  outfile     = paste0(OutDir, "Correlation_OurAnnotation_vs_Lake.pdf"),
  query_label = "Our annotation (FinalAnnotation_HC)",
  ref_label   = "Lake (author_cell_type)"
)

# =============================================================================
# 6. Summary: best-matching reference cell type per query cell type
# =============================================================================
summarise_best_match <- function(cor_mat, ref_name) {
  data.frame(
    OurCellType    = rownames(cor_mat),
    BestMatch      = colnames(cor_mat)[apply(cor_mat, 1, which.max)],
    MaxCorrelation = apply(cor_mat, 1, max),
    Reference      = ref_name,
    stringsAsFactors = FALSE
  )
}

summary_mka  <- summarise_best_match(cor_ours_mka,  "MKA")
summary_lake <- summarise_best_match(cor_ours_lake, "Lake")

summary_combined <- bind_rows(summary_mka, summary_lake) %>%
  pivot_wider(names_from = Reference,
              values_from = c(BestMatch, MaxCorrelation),
              names_glue = "{Reference}_{.value}")

write.csv(summary_combined,
          file = paste0(OutDir, "BestMatch_Summary_OurAnnotation_vs_References.csv"),
          row.names = FALSE)
message("Best-match summary saved: ",
        paste0(OutDir, "BestMatch_Summary_OurAnnotation_vs_References.csv"))

# =============================================================================
# 7. Combined dot-style summary plot
# =============================================================================
plot_best_match_bar <- function(df, ref_col_match, ref_col_cor, ref_name, outfile) {
  df <- df %>%
    rename(BestMatch = !!sym(ref_col_match),
           MaxCor    = !!sym(ref_col_cor)) %>%
    arrange(desc(MaxCor)) %>%
    mutate(OurCellType = factor(OurCellType, levels = OurCellType))

  p <- ggplot(df, aes(x = MaxCor, y = OurCellType, fill = MaxCor)) +
    geom_col() +
    geom_text(aes(label = BestMatch), hjust = -0.05, size = 3) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0.5, limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1.35)) +
    labs(title = paste0("Best-match correlation: Our annotation vs ", ref_name),
         x = "Max Pearson correlation", y = "Our cell type",
         fill = "Correlation") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  ggsave(outfile, plot = p,
         width  = 7,
         height = max(4, nrow(df) * 0.4 + 1))
  message("Bar plot saved: ", outfile)
}

plot_best_match_bar(
  df           = summary_combined,
  ref_col_match = "MKA_BestMatch",
  ref_col_cor   = "MKA_MaxCorrelation",
  ref_name      = "MKA",
  outfile       = paste0(OutDir, "BestMatch_Bar_OurAnnotation_vs_MKA.pdf")
)

plot_best_match_bar(
  df           = summary_combined,
  ref_col_match = "Lake_BestMatch",
  ref_col_cor   = "Lake_MaxCorrelation",
  ref_name      = "Lake",
  outfile       = paste0(OutDir, "BestMatch_Bar_OurAnnotation_vs_Lake.pdf")
)

message("\nAll outputs written to: ", OutDir)
message("Done.")
