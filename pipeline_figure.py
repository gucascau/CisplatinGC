#!/usr/bin/env python3
"""
Methods pipeline figure for cisplatin kidney scRNA-seq study.
Generates a publication-quality workflow diagram.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

# ── Palette ──────────────────────────────────────────────────────────────────
HDR_BG    = '#2B5FA0'   # column header fill
HDR_TXT   = '#FFFFFF'
ROW_BG    = '#3A3A5C'   # left row-label fill
ROW_TXT   = '#FFFFFF'
CELL_BG   = '#F0F5FF'   # content cell fill
CELL_EDGE = '#2B5FA0'
SUB_BLUE  = '#C8DCFA'   # sub-box: reference-based methods
SUB_GRN   = '#C8EDDC'   # sub-box: model-based methods
SUB_ORG   = '#FDE9C8'   # sub-box: validation
SUB_GRY   = '#E8E8EE'   # sub-box: neutral
ARROW_C   = '#444466'
TXT_C     = '#1A1A2E'

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 8.5))
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 22)
ax.set_ylim(0, 8.5)
ax.axis('off')

# ── Helper utilities ──────────────────────────────────────────────────────────
def rect(x, y, w, h, fc, ec, lw=1.2, zorder=1):
    p = FancyBboxPatch((x, y), w, h,
                        boxstyle="round,pad=0.08",
                        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=zorder)
    ax.add_patch(p)

def label(x, y, txt, fs=8.5, ha='center', va='center', color=TXT_C, bold=False):
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, txt, fontsize=fs, ha=ha, va=va, color=color,
            fontweight=weight, multialignment='center', linespacing=1.4)

def horiz_arrow(x1, x2, y):
    ax.annotate('', xy=(x2, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle='->', color=ARROW_C, lw=1.6,
                                mutation_scale=14))

def vert_arrow(x, y1, y2):
    ax.annotate('', xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle='->', color=ARROW_C, lw=1.2,
                                mutation_scale=12))

def sub_box(x, y, w, h, txt, fc=SUB_BLUE, fs=7.8):
    rect(x, y, w, h, fc=fc, ec=CELL_EDGE, lw=0.8, zorder=3)
    label(x + w/2, y + h/2, txt, fs=fs, color=TXT_C)

# ── Column layout  (title, x_left, width) ────────────────────────────────────
COL_H   = 0.62          # header height
CONTENT_Y = 0.35        # content box bottom
HEADER_Y  = 7.55        # header box bottom
CELL_TOP  = HEADER_Y    # = content top

cols = [
    ('Pre-Processing',         1.55,  3.10),
    ('QC & Doublet Removal',   4.90,  2.90),
    ('Integration',            8.05,  2.90),
    ('Cell Type Annotation',  11.20,  4.40),
    ('Visualization\n& Analysis', 15.85, 5.90),
]

# ── Left row label ────────────────────────────────────────────────────────────
rect(0.10, CONTENT_Y, 1.30, CELL_TOP - CONTENT_Y, fc=ROW_BG, ec=ROW_BG, lw=0)
label(0.75, (CONTENT_Y + CELL_TOP)/2,
      'Mouse\nKidney\nscRNA-seq\n\nVehicle\nCisplatin\nCisplatin-GC',
      fs=8.0, color=ROW_TXT, bold=False)

# ── Column headers ────────────────────────────────────────────────────────────
for title, x, w in cols:
    rect(x, HEADER_Y, w, COL_H, fc=HDR_BG, ec=HDR_BG, lw=0, zorder=2)
    label(x + w/2, HEADER_Y + COL_H/2, title, fs=10, color=HDR_TXT, bold=True)

# ── Content cells (background) ────────────────────────────────────────────────
for _, x, w in cols:
    rect(x, CONTENT_Y, w, CELL_TOP - CONTENT_Y, fc=CELL_BG, ec=CELL_EDGE, lw=1.0)

# ── Arrows between columns ────────────────────────────────────────────────────
for i in range(len(cols) - 1):
    _, x1, w1 = cols[i]
    _, x2, _  = cols[i + 1]
    horiz_arrow(x1 + w1 + 0.01, x2 - 0.01, (CONTENT_Y + CELL_TOP) / 2)

# ─────────────────────────────────────────────────────────────────────────────
# Col 1 — Pre-Processing
# ─────────────────────────────────────────────────────────────────────────────
x, w = cols[0][1], cols[0][2]
cx = x + w / 2
sub_box(x+0.12, 5.65, w-0.24, 1.65,
        'Read-Length Filtering\n(Perl)\nTarget insert: 54 bp\nProduces *_Cfiltered files',
        fc=SUB_GRY, fs=7.8)
vert_arrow(cx, 5.62, 5.28)
sub_box(x+0.12, 3.40, w-0.24, 1.80,
        'PipSeeker v3.3\nChemistry: v4\nAligner: STAR 2.7.9a\nReference: GRCm39\n(pipseeker-gex-2022.04)',
        fc=SUB_BLUE, fs=7.8)
vert_arrow(cx, 3.37, 3.07)
label(cx, 2.80, 'Count matrix per sample\n(raw_matrix/matrix.mtx.gz)',
      fs=7.2, color='#444466')

# ─────────────────────────────────────────────────────────────────────────────
# Col 2 — QC & Doublet Removal
# ─────────────────────────────────────────────────────────────────────────────
x, w = cols[1][1], cols[1][2]
cx = x + w / 2
# QC filters sub-box
sub_box(x+0.10, 5.35, w-0.20, 1.90,
        'Cell Quality Control\nnFeature: 150 – 10,000\nnCount: 150 – 20,000\n%MT < 50%\nrDNA < 40%',
        fc=SUB_GRY, fs=7.8)
vert_arrow(cx, 5.32, 4.92)
# DoubletFinder sub-box
sub_box(x+0.10, 3.50, w-0.20, 1.30,
        'DoubletFinder\n~5% expected doublets\npN=0.25, PCs 1–20\nauto pK selection',
        fc=SUB_BLUE, fs=7.8)
vert_arrow(cx, 3.47, 3.07)
label(cx, 2.80, 'Singlet cells retained\nper sample',
      fs=7.2, color='#444466')

# ─────────────────────────────────────────────────────────────────────────────
# Col 3 — Integration
# ─────────────────────────────────────────────────────────────────────────────
x, w = cols[2][1], cols[2][2]
cx = x + w / 2
sub_box(x+0.10, 5.35, w-0.20, 1.90,
        'Seurat\nNormalization (NormalizeData)\nHVG: VST, top 2,000 genes\nPCA (30 PCs)',
        fc=SUB_GRY, fs=7.8)
vert_arrow(cx, 5.32, 4.92)
sub_box(x+0.10, 3.50, w-0.20, 1.30,
        'Harmony\nBatch correction by DataSet\nUMAP (dims 1–30)\nLouvain clustering (res = 0.4)',
        fc=SUB_BLUE, fs=7.8)
vert_arrow(cx, 3.47, 3.07)
label(cx, 2.80, 'Integrated Seurat object\n(4 conditions)',
      fs=7.2, color='#444466')

# ─────────────────────────────────────────────────────────────────────────────
# Col 4 — Cell Type Annotation  (5-method consensus)
# ─────────────────────────────────────────────────────────────────────────────
x, w = cols[3][1], cols[3][2]
cx = x + w / 2
bw = (w - 0.30) / 2   # width of each method sub-box

# Method boxes (2 × 2 grid + 1 centered)
# Row A: Seurat transfer ×2
sub_box(x+0.10,          6.10, bw, 1.20,
        'Seurat Label\nTransfer\nMKA Atlas', fc=SUB_BLUE)
sub_box(x+0.10+bw+0.10,  6.10, bw, 1.20,
        'Seurat Label\nTransfer\nLake 2025', fc=SUB_BLUE)
# Row B: scANVI ×2
sub_box(x+0.10,          4.75, bw, 1.20,
        'scANVI\n(scvi-tools)\nMKA Atlas', fc=SUB_GRN)
sub_box(x+0.10+bw+0.10,  4.75, bw, 1.20,
        'scANVI\n(scvi-tools)\nLake 2025', fc=SUB_GRN)
# Manual annotation — centered between the two pairs
sub_box(x+0.10+bw/2,     3.48, bw, 1.10,
        'Manual\nAnnotation\n(canonical markers)', fc=SUB_ORG)

# Convergence arrow + consensus box
vert_arrow(cx, 3.44, 3.02)
sub_box(x+0.10, 2.08, w-0.20, 0.84,
        'Cluster-level Majority Vote  ·  14-class L1 vocabulary',
        fc='#EEE8FF', fs=7.8)
vert_arrow(cx, 2.04, 1.68)
sub_box(x+0.10, 0.90, w-0.20, 0.68,
        'UCell validation  ·  23 cell-type signatures  ·  concordance check',
        fc=SUB_ORG, fs=7.6)

# ─────────────────────────────────────────────────────────────────────────────
# Col 5 — Visualization & Analysis
# ─────────────────────────────────────────────────────────────────────────────
x, w = cols[4][1], cols[4][2]
cx = x + w / 2
bw2 = (w - 0.35) / 2

# Top row: UMAP + DotPlot
sub_box(x+0.12,          5.90, bw2, 1.40,
        'UMAP / DimPlot\nColored by:\n  cell type\n  condition\n  cluster', fc=SUB_BLUE)
sub_box(x+0.12+bw2+0.11, 5.90, bw2, 1.40,
        'DotPlot\nCanonical kidney\ncell-type markers\nper cluster', fc=SUB_BLUE)

# Middle row: FeaturePlot + Signaling
sub_box(x+0.12,          4.32, bw2, 1.40,
        'FeaturePlot\nPer-gene expression\nsplit by condition\n(23 marker signatures)', fc=SUB_GRN)
sub_box(x+0.12+bw2+0.11, 4.32, bw2, 1.40,
        'Signaling Gene\nFeatures\n(NicheNet)\nligand–receptor', fc=SUB_GRN)

# Bottom: candidate markers
sub_box(x+0.12, 2.80, w-0.24, 1.35,
        'Candidate Markers & Pathways\n'
        'Renal regional · AKI · CK · Inflammation · Fibrosis\n'
        'Apoptosis · AMPs · Cytokines · Chemokines\n'
        'Cisplatin nephrotoxicity signatures',
        fc=SUB_ORG, fs=7.8)

# ─────────────────────────────────────────────────────────────────────────────
# Figure title
# ─────────────────────────────────────────────────────────────────────────────
ax.text(11.0, 8.22,
        'Single-Cell RNA-seq Analysis Pipeline — Cisplatin Nephrotoxicity Study',
        fontsize=12, fontweight='bold', ha='center', va='center', color='#1A1A2E')

# ─────────────────────────────────────────────────────────────────────────────
# Legend
# ─────────────────────────────────────────────────────────────────────────────
legend_items = [
    mpatches.Patch(facecolor=SUB_BLUE, edgecolor=CELL_EDGE, label='Reference-based / established tools'),
    mpatches.Patch(facecolor=SUB_GRN,  edgecolor=CELL_EDGE, label='Deep learning / model-based'),
    mpatches.Patch(facecolor=SUB_ORG,  edgecolor=CELL_EDGE, label='Validation / downstream'),
    mpatches.Patch(facecolor=SUB_GRY,  edgecolor=CELL_EDGE, label='Preprocessing / filtering'),
]
ax.legend(handles=legend_items, loc='lower right',
          bbox_to_anchor=(0.99, 0.00),
          fontsize=7.8, framealpha=0.9, ncol=2,
          handlelength=1.2, handletextpad=0.5, columnspacing=1.0)

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
out_prefix = ('/vast0/home/gdzepedaorozcolab/lab/xxw004/Projects/'
              'RawDZscRNAseq/Scripts/Methods_Pipeline')
fig.savefig(out_prefix + '.pdf', dpi=300, bbox_inches='tight')
fig.savefig(out_prefix + '.png', dpi=300, bbox_inches='tight')
print(f"Saved: {out_prefix}.pdf  /  .png")
