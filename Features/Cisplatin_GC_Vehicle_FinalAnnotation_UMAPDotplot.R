## ============================================================
## Cisplatin GC Vehicle — Final Annotation UMAP & Dot Plots
## Input : Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS
## Output: UMAPs (per method, combined + split by DataSet)
##         DotPlots (per method + FinalAnnotation_HC)
## ============================================================

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggsci)
library(patchwork)

options(future.globals.maxSize = 1e15)

# ── Directories ──────────────────────────────────────────────
OutDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/IntegrateL1/"
RDSFile <- paste0(OutDir, "Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS")

dir.create(file.path(OutDir, "UMAPs"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OutDir, "DotPlots"), showWarnings = FALSE, recursive = TRUE)

# ── Load object ──────────────────────────────────────────────
cat("Reading RDS...\n")
obj <- readRDS(RDSFile)
cat("Done. Cells:", ncol(obj), "\n")

# ── Shared settings ──────────────────────────────────────────
method_cols <- c("L1_SeuratMKA", "L1_SeuratLake", "L1_Raw",
                 "L1_scanviMKA", "L1_scanviLake")

all_annotation_cols <- c(method_cols, "FinalAnnotation_HC")

# Nephron-axis cell-type order (for DotPlots)
celltype_order <- c(
  "PT", "DTL", "TAL", "DCT", "CNT", "IC/CNT",
  "PC", "IC", "POD", "PapE",
  "END", "FIB", "VSM/P", "Myeloid", "Lymphoid"
)

# Named marker list (grouped by cell type for clear DotPlot labelling)
marker_list <- c(
  # PT
  "Lrp2", "Slc34a1", "Slc13a3", "Slc5a2", "Slc22a6",
  # DTL
  "Aqp1", "Slc14a1",
  # TAL
  "Umod", "Slc12a1", "Cldn10", "Kcnj1", "Ptger3",
  # DCT
  "Slc12a3", "Pvalb", "Wnk1", "Fxyd2",
  # CNT
  "Klk1", "Slc8a1", "Calb1",
  # PC
  "Aqp2", "Hsd11b2", "Scnn1g",
  # IC
  "Atp6v1g3", "Aqp6", "Slc4a1", "Slc26a4",
  # POD
  "Nphs1", "Nphs2", "Wt1", "Synpo",
  # END
  "Pecam1", "Cdh5", "Kdr",
  # FIB
  "Col1a1", "Dcn", "Pdgfra",
  # VSM/P
  "Acta2", "Tagln", "Rgs5", "Pdgfrb",
  # Myeloid
  "C1qa", "C1qb", "Aif1", "Cd68",
  # Lymphoid
  "Ptprc", "Cd3e", "Cd79a", "Nkg7"
)

available_markers <- intersect(marker_list, rownames(obj))
n_markers <- length(available_markers)
cat("Markers available in assay:", n_markers, "of", length(marker_list), "\n")

# ── Helper: factor-order labels for a given column ───────────
# Returns ordered factor levels AND their count (for figure sizing)
order_labels <- function(seurat_obj, col) {
  vals <- as.character(seurat_obj@meta.data[[col]])
  vals[is.na(vals) | vals == ""] <- "Uncertain"
  present <- intersect(celltype_order, unique(vals))
  extra   <- setdiff(unique(vals), c(present, "Uncertain"))
  ordered_levels <- c(present, sort(extra), "Uncertain")
  seurat_obj@meta.data[[col]] <- factor(vals, levels = ordered_levels)
  seurat_obj
}

n_types_for <- function(seurat_obj, col) {
  nlevels(factor(as.character(seurat_obj@meta.data[[col]])))
}

# ── Figure-size helpers ───────────────────────────────────────
# UMAP combined : extra legend width for many cell types
umap_comb_size <- function(n_types) {
  w <- max(9, 7 + ceiling(n_types / 10) * 1.5)
  h <- max(6, 5 + ceiling(n_types / 20) * 0.5)
  c(w, h)
}

# UMAP split : each panel ~4.5 × 4.5 inches + 2.5 for shared legend
umap_split_size <- function(n_datasets, ncols = 3) {
  n_rows <- ceiling(n_datasets / ncols)
  w <- ncols * 4.5 + 2.5
  h <- n_rows  * 4.5 + 0.8   # 0.8 for title
  c(w, h)
}

# DotPlot horizontal : markers on x, cell types on y
dot_horiz_size <- function(n_types, n_markers) {
  w <- max(10, n_markers * 0.32 + 3)
  h <- max(4,  n_types   * 0.40 + 2)
  c(w, h)
}

# DotPlot flipped : cell types on x, markers on y
dot_flip_size <- function(n_types, n_markers) {
  w <- max(7,  n_types   * 0.50 + 3)
  h <- max(8,  n_markers * 0.32 + 2)
  c(w, h)
}

# ── 1. UMAPs ─────────────────────────────────────────────────
cat("\n=== Generating UMAPs ===\n")

n_datasets <- length(unique(obj@meta.data$DataSet))
umap_ncols <- min(n_datasets, 3)

for (col in all_annotation_cols) {
  cat("  UMAP:", col, "\n")
  tmp     <- order_labels(obj, col)
  n_types <- nlevels(tmp@meta.data[[col]])

  # -- 1a. Combined (no split) --
  sz_comb <- umap_comb_size(n_types)

  p_comb <- DimPlot(tmp,
                    group.by   = col,
                    label      = TRUE,
                    label.size = 3,
                    repel      = TRUE,
                    pt.size    = 0.3) +
    ggtitle(paste0(col, " — all cells")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold"))

  ggsave(
    filename = file.path(OutDir, "UMAPs", paste0(col, "_UMAP_Combined.pdf")),
    plot     = p_comb,
    width    = sz_comb[1], height = sz_comb[2]
  )

  # -- 1b. Split by DataSet --
  sz_split <- umap_split_size(n_datasets, umap_ncols)

  p_split <- DimPlot(tmp,
                     group.by   = col,
                     split.by   = "DataSet",
                     label      = TRUE,
                     label.size = 2.5,
                     repel      = TRUE,
                     pt.size    = 0.3,
                     ncol       = umap_ncols) +
    ggtitle(paste0(col, " — split by DataSet")) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right",
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))

  ggsave(
    filename = file.path(OutDir, "UMAPs", paste0(col, "_UMAP_SplitByDataSet.pdf")),
    plot     = p_split,
    width    = sz_split[1], height = sz_split[2]
  )
}

cat("UMAPs saved to:", file.path(OutDir, "UMAPs"), "\n")

# ── 2. Dot Plots ─────────────────────────────────────────────
cat("\n=== Generating DotPlots ===\n")

make_dotplot <- function(seurat_obj, col, tag) {
  tmp     <- order_labels(seurat_obj, col)
  n_types <- nlevels(tmp@meta.data[[col]])

  # Enforce nephron-axis order: set Idents to the ordered factor so
  # DotPlot rows appear in the same sequence as celltype_order
  Idents(tmp) <- tmp@meta.data[[col]]

  sz_h <- dot_horiz_size(n_types, n_markers)
  sz_f <- dot_flip_size(n_types, n_markers)

  # -- Version 1: horizontal, RdYlBu distiller --
  p1 <- DotPlot(tmp,
                features = available_markers,
                group.by = col) +
    scale_colour_distiller(palette = "RdYlBu", direction = -1) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 9),
          legend.position = "right") +
    labs(title = paste0("Marker expression — ", col), x = NULL, y = NULL)

  ggsave(
    filename = file.path(OutDir, "DotPlots", paste0(tag, "_DotPlot_RdYlBu.pdf")),
    plot     = p1,
    width    = sz_h[1], height = sz_h[2]
  )

  # -- Version 2: flipped (markers on y), cell types readable on x --
  p2 <- DotPlot(tmp,
                features = available_markers,
                group.by = col,
                cols     = "RdYlBu") +
    coord_flip() +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 8),
          legend.position = "right") +
    labs(title = paste0("Marker expression (flipped) — ", col), x = NULL, y = NULL)

  ggsave(
    filename = file.path(OutDir, "DotPlots", paste0(tag, "_DotPlot_Flipped.pdf")),
    plot     = p2,
    width    = sz_f[1], height = sz_f[2]
  )

  cat(sprintf("  Saved: %-20s  horiz=%dx%d  flip=%dx%d  (%d types)\n",
              tag, round(sz_h[1]), round(sz_h[2]),
              round(sz_f[1]), round(sz_f[2]), n_types))
}

for (col in all_annotation_cols) {
  make_dotplot(obj, col, tag = col)
}

cat("DotPlots saved to:", file.path(OutDir, "DotPlots"), "\n")
cat("\nAll done.\n")
