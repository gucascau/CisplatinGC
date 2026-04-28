#!/bin/bash
#SBATCH --job-name=CispGC_UMAPDot
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs/CispGC_UMAPDot_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs/CispGC_UMAPDot_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --partition=himem

# ---------------------------------------------------------------------------
# Load modules
# ---------------------------------------------------------------------------
module purge
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.0

# ---------------------------------------------------------------------------
# Create log directory if it does not exist
# ---------------------------------------------------------------------------
mkdir -p /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs

# ---------------------------------------------------------------------------
# Run UMAP + DotPlot script
# ---------------------------------------------------------------------------
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/Annotation"

Rscript "${SCRIPT_DIR}/Cisplatin_GC_Vehicle_FinalAnnotation_UMAPDotplot.R"
