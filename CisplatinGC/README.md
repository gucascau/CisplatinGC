# CisplatinGC — Single-Cell RNA-seq Analysis

Integration and analysis of mouse kidney single-cell RNA-seq data across four treatment conditions:
**Vehicle**, **Vehicle-GC**, **Cisplatin**, and **Cisplatin-GC**.

---

## Project Structure

```
CisplatinGC/
├── PreProcess/                          # FASTQ processing and cell-barcode quantification
│   ├── ExtractSequenceCorrectedBaseOnLength.pl   # Read-length filtering (Perl)
│   ├── PipseekerPipeline_Vehicle.sh              # PipSeeker job for Vehicle
│   ├── PipseekerPipeline_Vehicle-GC.sh           # PipSeeker job for Vehicle-GC
│   ├── PipseekerPipeline_Cisplatin.sh            # PipSeeker job for Cisplatin
│   └── PipseekerPipeline_Cisplatin-GC.sh         # PipSeeker job for Cisplatin-GC
└── Integration/                         # Downstream R analysis
    ├── Cisplatin_Vehicle_GC_Integration.r         # Main integration script (R)
    ├── Cisplatin_Vehicle_GC_Integration.sh        # SLURM job wrapper for R script
    └── Cisplatin_Vehicle_GC_Integration.rmd       # R Markdown version
```

---

## Workflow Overview

### Step 1 — Preprocessing (PreProcess/)

Each sample is processed independently with a two-step approach:

1. **Read-length filtering** (`ExtractSequenceCorrectedBaseOnLength.pl`)  
   Filters paired-end FASTQ reads to a target insert length of **54 bp**, producing `*_Cfiltered_L00{1,2}` files.

2. **Cell-barcode quantification** (PipSeeker v3.3)  
   Runs `pipseeker full` with:
   - Chemistry: `v4`
   - Reference genome: GRCm39 (`pipseeker-gex-reference-GRCm39-2022.04`)
   - Aligner: STAR 2.7.9a

**SLURM resource requirements per job:** `himem` partition, 12 CPUs, 40 h wall time.

---

### Step 2 — Integration and Analysis (Integration/)

Run via SLURM: `sbatch Cisplatin_Vehicle_GC_Integration.sh` (R/4.4.0, 24 CPUs, `himem`).

The main script `Cisplatin_Vehicle_GC_Integration.r` performs the following steps:

#### 2a. Quality Control
- Reads raw matrices from all four samples (`raw_matrix/matrix.mtx.gz`)
- Creates Seurat objects (min 3 cells per gene, min 200 features per cell)
- Calculates mitochondrial (`^mt-`) and ribosomal (`^Rp[sl]`) gene percentages
- Filters cells:
  - `nFeature_RNA`: 150 – 10,000
  - `nCount_RNA`: 150 – 20,000
  - `percent.mt` < 50%
  - `rDNA` < 40%
- Outputs QC violin plots to `QuanlityControls/`

#### 2b. Doublet Detection
- Runs **DoubletFinder** per sample (expected ~5% doublets, pN = 0.25, PCs 1–20)
- Automatic pK selection via BCmetric maximization
- Saves doublet barcodes and UMAP plots to `DoubletFinder/`
- Removes doublets; saves merged object with doublet labels:  
  `Cisplatin_GC_Vehicle_scrna_merged_withdoublets.rds`

#### 2c. Harmony Integration
- Merges all four singlet datasets into one Seurat object
- Normalization → Variable features (VST, 2000 genes) → Scaling → PCA (30 PCs)
- Batch correction with **Harmony** (grouping by `DataSet`)
- UMAP on Harmony embedding (dims 1–30)
- Clustering: `FindNeighbors` + `FindClusters` (resolution = 0.4, Louvain)
- Saves integrated object to `IntegratedFourConditions/`:  
  `Cisplatin_GC_Vehicle_singlecell_doublet_harmony_v<date>.RDS`

#### 2d. Candidate Marker Analysis
- **FeaturePlot** for urothelial and stem-cell markers (e.g., `Krt14`, `Trp63`, `Upk2`)
- **DotPlot** across kidney cell-type marker panel (urothelial, proximal tubule, distal tubule, TAL, CNT, DTL, principal cells, intercalated cells, endothelium, fibroblasts, immune cells)
- Outputs saved to `CandidateMarkers/`

---

## Software Dependencies

| Tool | Version |
|------|---------|
| STAR | 2.7.9a |
| PipSeeker | v3.3.0 |
| R | 4.4.0 |
| Seurat | ≥ 4.x |
| Harmony | latest |
| DoubletFinder | latest (GitHub) |
| monocle3 | latest (GitHub) |

Key R packages: `Seurat`, `harmony`, `DoubletFinder`, `monocle3`, `sctransform`,
`tidyverse`, `ggplot2`, `patchwork`, `cowplot`, `ggsci`, `future`, `biomaRt`

---

## Input / Output

| | Path |
|---|---|
| Raw FASTQ | `Results/RawDZscRNAseq/<SampleID>/` |
| PipSeeker output | `Results/<SampleID>/raw_matrix/` |
| Integration output | `Results/Integration/IntegrationGC/` |

---

## Author

**Xin Wang** — xin.wang@nationwidechildrens.org  
Date: 2026-04-06
