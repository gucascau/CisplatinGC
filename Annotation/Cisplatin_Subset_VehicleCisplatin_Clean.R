## ============================================================
## Subset Cisplatin GC/Vehicle object to a clean Vehicle + Cisplatin object
## Input : Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS
## Output: ..._VehicleCisplatin_Clean.RDS
##         Only Vehicle and Cisplatin cells are kept, and unused DataSet
##         factor levels (Vehicle-GC / Cisplatin-GC) are dropped so they
##         cannot reappear as empty categories in downstream table()/ggplot
##         calls.
## ============================================================

library(Seurat)

InDir <- "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/IntegrateL1/"

RDSFile    <- paste0(InDir, "Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation.RDS")
OutRDSFile <- paste0(InDir, "Cisplatin_GC_Vehicle_Reclustered_MarkerAnnot_L1Harmonised_NoVehGC_Weighted_FinalAnnotation_VehicleCisplatin_Clean.RDS")

cat("Reading RDS...\n")
obj <- readRDS(RDSFile)
cat("Done. Cells:", ncol(obj), "\n")
cat("DataSet counts before subsetting:\n")
print(table(obj@meta.data$DataSet))

# Keep only Vehicle and Cisplatin (drop Vehicle-GC / Cisplatin-GC)
obj <- subset(obj, subset = DataSet %in% c("Vehicle", "Cisplatin"))

# Seurat v5 leaves per-dataset split layers (e.g. counts.Vehicle, data.Cisplatin)
# after subsetting. Join them now so downstream FetchData()/DotPlot() calls on
# the saved object can find features without erroring.
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)

# subset() does not drop unused factor levels on its own, so Vehicle-GC /
# Cisplatin-GC would otherwise still show up as empty categories later on
obj@meta.data$DataSet <- droplevels(factor(obj@meta.data$DataSet))

cat("DataSet counts after subsetting:\n")
print(table(obj@meta.data$DataSet))

cat("Saving clean object to:", OutRDSFile, "\n")
saveRDS(obj, file = OutRDSFile)
cat("Done.\n")
