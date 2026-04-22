# Cisplatin GC Vehicle — Cell Type Annotation Pipeline

## Pipeline diagram

```mermaid
flowchart TD

    %% ── Prerequisites ────────────────────────────────────────────────────────
    subgraph PRE ["Prerequisites (run once)"]
        direction LR
        P1["Cisplatin_GC_Vehicle_Annotation.Rmd\n────────────────────────────────────\n• Seurat Transfer — MKA reference\n  → Seaurat_Transfer_Predicted_MKA\n  → Seaurat_Transfer_Predicted_MKA.max\n• Seurat Transfer — Lake 2025 reference\n  → Seaurat_Transfer_Predicted_Lake\n  → Seaurat_Transfer_Predicted_Lake.max\n• Manual annotation\n  → Raw_cell_type\n  → Raw_cell_type_Confidence"]
        P2["scvi_annotation*.ipynb\n────────────────────────────────────\n• scANVI — MKA reference\n  → scanvi_MKA_label\n  → scanvi_MKA_confidence\n• scANVI — Lake 2025 reference\n  → scanvi_Lake_label\n  → scanvi_Lake_confidence"]
    end

    %% ── Step 1 ───────────────────────────────────────────────────────────────
    subgraph S1 ["Step 1 · Annotation_Integration.Rmd"]
        direction TB
        S1A["Subset to Vehicle + Cisplatin + Cisplatin-GC\n(Vehicle-GC excluded from analysis)"]
        S1B["Harmony integration across conditions\n→ Reclustering  res=1, dims 1:20"]
        S1C["Sub-cluster clusters 5 & 14\n→ Sub2 (final cluster column)"]
        S1D["Harmonise all 5 labels → shared L1 vocabulary\n13 categories: PT · TAL · DCT · CNT · PC · IC\nIC/CNT · POD · END · VSM/P · PapE · Myeloid · Lymphoid\n────────────────────────────────────────────────────────\nL1_SeuratMKA   L1_SeuratLake   L1_Raw\nL1_scanviMKA   L1_scanviLake"]
    end

    %% ── Step 2 ───────────────────────────────────────────────────────────────
    subgraph S2 ["Step 2 · FinalAnnotation.Rmd  ← PRIMARY SCRIPT"]
        direction TB
        S2A["PRIMARY — Cluster-level majority vote\n──────────────────────────────────────────────────\nPool all 5-method votes across all cells per Sub2 cluster\nPlurality label  →  ClusterAnnotation\nVote fraction    →  ClusterAnnotation_frac"]
        S2B["VALIDATION — Per-cell equal-weight consensus\n──────────────────────────────────────────────────\nNormalise each confidence score to 0–1 (min-max)\nConsensus score = Σ 0.20 × norm_confidence  (5 methods)\nTop-scoring label  →  WeightedAnnotation"]
        S2C["Concordance check\nDoes WeightedAnnotation = ClusterAnnotation?\n→  Cell_Concordant  TRUE / FALSE"]
        S2D["Confidence tier\n──────────────────────────────────────────────────\nHigh   :  vote fraction ≥ 70%  AND  Concordant\nMedium :  vote fraction ≥ 50%  OR   Concordant\nLow    :  vote fraction < 50%  AND  Discordant\n→  AnnotationConfidence"]
        S2E["FinalAnnotation = ClusterAnnotation\nDiscordant cells flagged, not removed"]
    end

    %% ── Optional ─────────────────────────────────────────────────────────────
    subgraph OPT ["Optional · WeightedEnsemble_Annotation.Rmd"]
        direction TB
        O1["Detailed per-cell consensus analysis\nScore distributions · decision margins\nMethod agreement counts\nSupplementary figures only"]
    end

    %% ── Output ───────────────────────────────────────────────────────────────
    OUT[("Final annotated RDS\n─────────────────────────────────\nFinalAnnotation\nClusterAnnotation_frac\nWeightedAnnotation\nWeightedAnnotation_score\nCell_Concordant\nAnnotationConfidence")]

    %% ── Connections ──────────────────────────────────────────────────────────
    PRE --> |"*_Combined_Predict_MKA_LakeWithL2_scanvi_*.RDS\n+ scanvi CSV files"| S1
    S1  --> |"*_L1Harmonised_NoVehGC_Annotation.RDS"| S2
    S1  -.-> |"optional deeper dive"| OPT
    S2  --> OUT
```

---

## Script run order

| Order | Script | Role | Key output columns |
|---|---|---|---|
| 0a | `Cisplatin_GC_Vehicle_Annotation.Rmd` | Seurat Transfer + manual annotation | `Seaurat_Transfer_Predicted_*`, `Raw_cell_type` |
| 0b | `scvi_annotation*.ipynb` | scANVI predictions | `scanvi_*_label`, `scanvi_*_confidence` |
| 1 | `Cisplatin_GC_Vehicle_Annotation_Integration.Rmd` | Reclustering + L1 harmonisation | `L1_SeuratMKA`, `L1_SeuratLake`, `L1_Raw`, `L1_scanviMKA`, `L1_scanviLake` |
| **2** | **`Cisplatin_GC_Vehicle_FinalAnnotation.Rmd`** | **Cluster majority vote + per-cell consensus + tiers** | `FinalAnnotation`, `ClusterAnnotation_frac`, `WeightedAnnotation`, `Cell_Concordant`, `AnnotationConfidence` |
| (opt) | `Cisplatin_GC_Vehicle_WeightedEnsemble_Annotation.Rmd` | Detailed per-cell consensus analysis for supplementary figures | `WeightedAnnotation_score`, `WeightedAnnotation_margin` |

> **Note:** `FinalAnnotation.Rmd` is self-contained — it computes the per-cell consensus internally and does not depend on the WeightedEnsemble script having been run first.

---

## Design rationale

### Why cluster-level majority vote as the primary label?

Single-cell data is inherently noisy (dropout, ambient RNA). Annotating at the cluster level — pooling votes from all five methods across all cells within each Seurat sub-cluster — averages out per-cell noise. This is the standard approach in published kidney single-cell studies and is the most defensible for a peer-reviewed paper.

### Why keep the per-cell consensus as a validation layer?

The per-cell consensus provides an independent cell-level check. Cells where it disagrees with the cluster label (`Cell_Concordant = FALSE`) may represent transitional states, doublets, or cluster boundary cells. Flagging them rather than silently overriding them preserves this biological information for downstream interpretation.

### Why normalise confidence scores?

Seurat `prediction.score.max`, scANVI posterior probability, and manual confidence are measured on different scales. Min-max normalisation to [0, 1] ensures each method's confidence contributes comparably before summing.

### Why equal weights?

All five methods are weighted equally (0.20 each). No published benchmark directly comparing all five methods on mouse kidney single-cell data is available, making equal weighting the most conservative and reviewer-defensible choice.

| Method column | Confidence column | Weight |
|---|---|---|
| `L1_SeuratMKA` | `Seaurat_Transfer_Predicted_MKA.max` | 0.20 |
| `L1_SeuratLake` | `Seaurat_Transfer_Predicted_Lake.max` | 0.20 |
| `L1_Raw` | `Raw_cell_type_Confidence` | 0.20 |
| `L1_scanviMKA` | `scanvi_MKA_confidence` | 0.20 |
| `L1_scanviLake` | `scanvi_Lake_confidence` | 0.20 |

---

## Key output files

| File | Script | Contents |
|---|---|---|
| `Cluster_MajorityVote_Final.csv` | FinalAnnotation | Cluster vote tallies and winning label |
| `Cluster_Concordance_ClusterMajority_vs_Weighted.csv` | FinalAnnotation | Per-cluster concordance rate |
| `PerCell_FinalAnnotation_Summary.csv` | FinalAnnotation | Per-cell annotation with all layers |
| `FinalAnnotation_DimPlot.pdf` | FinalAnnotation | UMAP coloured by final annotation |
| `AnnotationConfidence_DimPlot.pdf` | FinalAnnotation | UMAP coloured by High / Medium / Low tier |
| `CellConcordance_DimPlot.pdf` | FinalAnnotation | UMAP highlighting discordant cells |
| `Cluster_VoteFraction_Bar.pdf` | FinalAnnotation | Cluster vote strength bar chart |
| `FinalAnnotation_ConfidenceBar.pdf` | FinalAnnotation | Tier breakdown per cell type |
| `WeightedAnnotation_ScoreSummary.csv` | WeightedEnsemble | Per-cell consensus score distributions |
| `WeightedAnnotation_DimPlot.pdf` | WeightedEnsemble | UMAP coloured by per-cell consensus label |

---

## Methods text for paper

### Cell type annotation

Cell types were annotated using five independent methods applied to the
harmonised single-cell dataset: Seurat label transfer with the MKA reference
(Stuart et al., 2019, *Cell*), Seurat label transfer with the Lake 2025
reference (Lake et al., 2025, *Nature*), scANVI with the MKA reference
(Xu et al., 2021, *Molecular Systems Biology*), scANVI with the Lake 2025
reference, and manual annotation guided by established mouse kidney marker
genes. All five predictions were harmonised to a shared L1 cell-type
vocabulary comprising 13 categories (PT, TAL, DCT, CNT, PC, IC, IC/CNT,
POD, END, VSM/P, PapE, Myeloid, Lymphoid).

The primary cell-type label for each Seurat sub-cluster was determined by
plurality vote, pooling predictions from all five methods across all cells
within each cluster. This cluster-level approach is robust to per-cell
transcriptomic noise inherent in single-cell data. The vote fraction —
the proportion of total votes assigned to the winning label — was recorded
as a measure of cluster-level annotation confidence.

To provide independent per-cell validation, a consensus score was computed
for each cell as the sum of each method's normalised prediction confidence,
with all five methods weighted equally (0.20 each). Prediction confidence
scores from each method were first min-max normalised to [0, 1] to ensure
comparability across methods that use different scoring scales (e.g. Seurat
prediction score versus scANVI posterior probability). The label with the
highest total consensus score was designated the per-cell consensus label.

Cells where the per-cell consensus label agreed with the cluster majority
label were classified as concordant. Cluster annotations were assigned a
confidence tier based on vote fraction and concordance: High (vote fraction
≥ 70% and concordant), Medium (vote fraction ≥ 50% or concordant), and Low
(both criteria unmet). Discordant cells were flagged in the metadata but
retained for downstream analyses. All annotation code is available at
[GitHub repository URL].
