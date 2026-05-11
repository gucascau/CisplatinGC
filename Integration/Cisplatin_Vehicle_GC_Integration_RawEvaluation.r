#!/usr/bin/env Rscript
# ---
# title: "Cisplatin_Vehicle_GC_Integration.r"
# author: "Xin Wang"
# date: "2026-04-06"
# email: xin.wang@nationwidechildrens.org
# output: html_document
# description: this script is to integrate the single cells from Cisplatin, Cisplatin GC, Vehicle and Vehicle GC
#           1. The quality controls of each samples
#           2. The doublet detection of each samples
#           3. The integration of the four samples by harmony
#           4. The candidate markers for each clusters

# loading the packages
# check the session information
sessionInfo()
library(dplyr)
library(Seurat)
library(patchwork)

#install.packages('devtools') #assuming it is not already installed

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
library(devtools)
#install_github('andreacirilloac/updateR')
# install.packages("sctransform")
# install.packages("ggthemes")
# loading the libarary
#library(updateR) # Windows-only, not available on Linux HPC

library("sctransform")
library("ggthemes")
library(Seurat)

# BiocManager::install("limma")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("biomaRt")
# BiocManager::install("methrix")BiocManager::install("org.Hs.eg.db")
# devtools::install_github("saeyslab/nichenetr")
library(limma)
#library(nichenetr)
library(tidyverse)
library(dplyr)

library(biomaRt)
library("data.table")


library(patchwork)
library(cowplot)
#library(umap)     # not installed; Seurat uses uwot internally
#library(installr) # Windows-only, not available on Linux HPC

#install.packages("Matrix", repos = "http://cran.r-project.org")
### we used discre color palettes from ggsci

library("ggsci")
library("ggplot2")
library("gridExtra")


#options(buildtools.check = function(action) TRUE )
#devtools::install_url('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz')
#devtools::install_github('cole-trapnell-lab/monocle3')
#library(monocle3)

library(harmony)
library(Seurat)
#library(SeuratData)
library(tidyverse)

#library(BPCells)
library(ggplot2)
library("devtools")
library("AnnotationDbi")
library("org.Mm.eg.db")  # mouse genome annotation (data is mouse)
#library(tibble)
#library(future)
library(here)
library(patchwork)
library(future)
#library(monocle3) # not installed and not used in this script
#devtools::install_dev("remotes")
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#remotes::install_github(repo ='chris-mcginnis-ucsf/DoubletFinder' )
#install.packages("DoubletFinder")
library(DoubletFinder)
#library(future.callr) # not installed; parallel plan is disabled anyway
library("ggsci")
library("ggplot2")
library("gridExtra")

setwd(
  "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/"
)
Indir <-
  c(
    "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Preprocess/"
  )
Outdir <-
  c(
    "/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/"
  )
## check the files in the current directory
Datafiles <-
  list.files(
    path = Indir,
    recursive = F,
    full.names = F
  )

# setting up the seed 
set.seed(10000)


QCDir<- paste0(Outdir, "QuanlityControls/")
dir.create(QCDir)
setwd(QCDir)
# create an empty list to store the objects
scrna.list = list()

# please evaluate the quality of the data before the integration, and filter the low quality cells, and remove the doublets. The quality control is based on the number of features, number of counts, percentage of mitochondrial genes and percentage of rDNA genes. The doublet detection is based on the DoubletFinder package. The integration is based on the Harmony package. The candidate markers are based on the FindAllMarkers function in Seurat.
for (i in Datafiles) {
  # initiate the seurate object
  scrna.list[[i]] = ReadMtx(
    mtx = paste0(Indir, i, "/raw_matrix/matrix.mtx.gz"),
    features = paste0(Indir, i, "/raw_matrix/features.tsv.gz"),
    cells = paste0(Indir, i, "/raw_matrix/barcodes.tsv.gz")
  )
  
  scrna.list[[i]] <- as(scrna.list[[i]], "dgCMatrix")
  # initiate the seurate object
  #scrna.list[[name2]]<- CreateSeuratObject(counts=rds, project= name2)
  scrna.list[[i]] = CreateSeuratObject(
    counts = scrna.list[[i]],
    project = i,
    min.cells = 3,
    min.features = 200
  )
  # measure the mitochondrial percentage
  scrna.list[[i]][["percent.mt"]] <-
  PercentageFeatureSet(scrna.list[[i]], pattern = "^mt-")
  scrna.list[[i]][["rDNA"]] <-
  PercentageFeatureSet(scrna.list[[i]], pattern = "^Rp[sl][[:digit:]]")
  
  # rename the cell with name id
  scrna.list[[i]] <- RenameCells(scrna.list[[i]], add.cell.id = i)
  ## using the assign the name to the objects
  assign (i, scrna.list[[i]])
  ## create a group that show the name of samples
  scrna.list[[i]]$DataSet <-
    rep(i, length(scrna.list[[i]]$orig.ident))

  #  Please evaluate the quality of the data before the integration, and filter the low quality cells, and remove the doublets. The quality control is based on the number of features, number of counts, percentage of mitochondrial genes and percentage of rDNA genes.
  #  Generate the vlnplot 

    #VlnPlot(scrna, features = c("percent.mt", "nCount_RNA", "nFeature_RNA","rDNA"))
    QCReport <-
        VlnPlot(
     scrna.list[[i]],
     features = c("rDNA", "percent.mt", "nCount_RNA", "nFeature_RNA") ,
     #split.by = "DataSet",
     pt.size=0,
     ncol = 1
    )
  QCReport[[3]] <- QCReport[[3]] + coord_cartesian(ylim = c(0, 20000))
  QCReport[[4]] <- QCReport[[4]] + coord_cartesian(ylim = c(0, 6000))
# save the vlnplot for each sample
  ggsave(filename = paste0(i, "_QC_reports.pdf"), plot = QCReport, height = 12, width = 4)
  

}

scrna.list

# please measure these lists, ignore the Vehicle-GC

scrna <-
  merge(
    x = scrna.list[[1]],
    y = scrna.list[2:3],
    project = "Cisplatin_Vehicle_GC_Integration"
  )

# check the meta data
levels (as.factor(scrna@meta.data$DataSet ))

scrna@meta.data$DataSet <-
  ordered(factor(scrna@meta.data$DataSet),
          levels = c("Vehicle", "Cisplatin","Cisplatin-GC"))

scrna@meta.data %>% tail()
# create a condition details with and without infection

# measure the mitochondrial percentage
scrna[["percent.mt"]] <-
  PercentageFeatureSet(scrna, pattern = "^mt-")
scrna[["rDNA"]] <-
  PercentageFeatureSet(scrna, pattern = "^Rp[sl][[:digit:]]")

#VlnPlot(scrna, features = c("percent.mt", "nCount_RNA", "nFeature_RNA","rDNA"))
 QCReport <-
   VlnPlot(
     scrna,
     features = c("rDNA", "percent.mt", "nCount_RNA", "nFeature_RNA") ,
     split.by = "DataSet",
     pt.size=0,
     ncol = 1
   )
QCReport[[3]] <- QCReport[[3]] + coord_cartesian(ylim = c(0, 15000))
QCReport[[4]] <- QCReport[[4]] + coord_cartesian(ylim = c(0, 5000))

ggsave(filename = "Three_SampleQC_reports.pdf",plot = QCReport, height = 12, width = 4 )

DoubletDir<- paste0(Outdir, "DoubletFinder/")
dir.create(DoubletDir)
setwd(DoubletDir)
# check the doublets
scrna <- readRDS("Cisplatin_GC_Vehicle_scrna_merged_withdoublets.rds")

# summarize the doublet detection results
scrna@meta.data %>% group_by(DataSet) %>% summarise(Doublet = sum(doublet_id == "Doublet"), Singlet = sum(doublet_id == "Singlet"))

scrna@meta.data %>% group_by(DataSet)  %>% head()
