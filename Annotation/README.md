# Annotation Pipeline — Cisplatin / GC / Vehicle scRNA-seq

**Author:** Xin Wang  
**Dataset:** 4-condition mouse kidney single-cell dataset (Vehicle, Vehicle-GC, Cisplatin, Cisplatin-GC)

---

## Overview

Cell type annotation is performed by combining five independent methods into a single consensus label. The pipeline is intentionally redundant: no single method is trusted alone. A cluster-level majority vote determines the final call, and UCell module scores provide marker-gene-based validation independent of any reference.

---

## Scripts (run in order)

| Step | Script | Output |
|------|--------|--------|
| 1 | `Cisplatin_GC_Vehicle_Annotation.Rmd` | `*_Combined_Predict_MKA_LakeWithL2_scanvi_*.RDS` |
| 2 | `Cisplatin_GC_Vehicle_scvi_annotation*.ipynb` | scanvi label CSV files |
| 3 | `Cisplatin_GC_Vehicle_Annotation_Integration.Rmd` | `Cisplatin_GC_Vehicle_Final_Annotated_v04202026.rds` |

---

## Five Annotation Methods

| # | Column | Method | Reference |
|---|--------|--------|-----------|
| 1 | `Seaurat_Transfer_Predicted_MKA` | Seurat label transfer | MKA reference atlas |
| 2 | `Seaurat_Transfer_Predicted_Lake` | Seurat label transfer | Lake 2025 reference |
| 3 | `Raw_cell_type` | Manual annotation | — |
| 4 | `scanvi_MKA_label` | scANVI | MKA reference atlas |
| 5 | `scanvi_Lake_label` | scANVI | Lake 2025 reference |

MKA uses granular subtypes (PTS1–3, CTAL/MTAL, DCT/CNT, etc.) that are mapped to the shared L1 vocabulary before voting. Lake labels are already L1-compatible.

---

## Integration Steps

### 1. Subsetting & Re-clustering
Cells from Vehicle, Cisplatin, and Cisplatin-GC are subset from the full 4-condition object and re-clustered together using standard Seurat workflow: normalisation → HVG selection (mt-genes excluded) → PCA → **Harmony** batch correction (by `DataSet`) → neighbour graph → Leiden clustering (resolution = 1) → UMAP.

### 2. L1 Vocabulary Harmonisation
All five method labels are mapped to a shared 14-class L1 vocabulary:

`PT · TAL · DCT · CNT · ATL · PC · IC · CD · POD · END · FIB · Mes · Myeloid · Lymphoid`

Labels that cannot be mapped become `"Uncertain"` and are excluded from voting.

### 3. Cluster-level Majority Vote
For each Seurat cluster, all per-cell votes from the five methods are pooled. The L1 label with the most votes wins (`ClusterMajority`). The vote fraction (`ClusterMajority_frac`) is retained as a confidence measure.

### 4. UCell Marker-gene Validation
UCell module scores are computed for 23 cell-type signatures (see marker table below). For each cluster the top-scoring signature (`UCell_TopType`) is compared with `ClusterMajority`. Concordance rate is reported.

### 5. Final Annotation
`FinalAnnotation = ClusterMajority` by default. Discordant clusters can be overridden in the `manual_overrides` table inside `Cisplatin_GC_Vehicle_Annotation_Integration.Rmd`. Every override requires a documented reason.

---

## Marker Gene Signatures (UCell & FeaturePlots)

| Cell Type | Key Markers |
|-----------|-------------|
| PT | Lrp2, Slc34a1, Slc13a3 |
| PT_S1S2 | Slc5a2, Slc34a1, Lrp2 |
| PT_S3 | Slc22a6, Slc22a8, Havcr1 |
| TAL | Umod, Slc12a1, Cldn10, Kcnj1 |
| DCT | Slc12a3, Pvalb, Wnk1, Trpm6 |
| CNT | Klk1, Slc8a1, Calb1 |
| LOH | Slc14a2, Aqp1, Clcnka |
| PC | Aqp2, Hsd11b2, Scnn1g |
| IC | Atp6v1g3, Aqp6, Slc26a7 |
| ICA | Atp6v1g3, Aqp6, Kit, Slc4a1, Slc26a7 |
| ICB | Atp6v1g3, Slc26a4, Spink8, Hmx2 |
| POD | Wt1, Nphs1, Nphs2 |
| END | Cdh5, Pecam1, Kdr, Emcn, Plvap |
| FIB | Col1a1, Col3a1, Fbn1, Lum, Pdgfra, Dcn |
| Mesangial | Pdgfrb, Acta2, Tagln |
| Pericyte | Pdgfrb, Rgs5, Kcnj8 |
| PEC | Krt8, Krt18, Claudin1 |
| pDC | Siglech, Bst2, Irf7 |
| Neutrophil | S100a8, S100a9, Ly6g, Mpo |
| Macrophage | C1qa, C1qb, Aif1, Adgre1, Cd68 |
| T_cells | Cd3e, Cd247, Trac |
| B_cells | Cd79a, Cd79b, Ms4a1 |
| NK_CD8 | Nkg7, Ccl5, Gzmb, Cd8a |

---

## Outputs

| File | Description |
|------|-------------|
| `ElbowPlot_PCA.pdf` | PCA elbow plot for PC cutoff selection |
| `*_DimPlot.pdf` | UMAP coloured by each annotation method |
| `Cluster_MajorityVote_Table.csv` | Per-cluster vote counts and majority label |
| `Cluster_UCell_TopType_Table.csv` | Per-cluster top UCell signature |
| `Cluster_MajorityVote_UCell_Concordance.csv` | Concordance between majority vote and UCell |
| `PerCell_Annotation_Summary.csv` | Per-cell table of all five method labels + final annotation |
| `FinalAnnotation_DotPlot.pdf` | Dot plot of marker genes by final cell type |
| `Cluster_UCell_Heatmap.pdf` | UCell score heatmap across clusters |
| `MarkerFeaturePlots/<CellType>/<Gene>_FeaturePlot.pdf` | Per-gene feature plots split by DataSet |
| `Cisplatin_GC_Vehicle_Final_Annotated_v04202026.rds` | Final annotated Seurat object |

---

## Notes

- `Vehicle-GC` is **excluded** from re-clustering to avoid confounding by GC effects, but its per-cell annotation labels are retained in the object.
- Column name `Seaurat_Transfer_Predicted_*` (typo "Seaurat") is intentional for consistency with the upstream annotation script.
- Harmony batch correction is applied over `DataSet` (not `orig.ident`) so that condition-level effects are corrected while biological variation is preserved.
- Genes absent from the object or with all-zero/NA expression are silently skipped in FeaturePlot generation.
