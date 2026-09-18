#!/bin/bash
#SBATCH --job-name=VehCisp_QC
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/QC/logs/VehCisp_QC_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/QC/logs/VehCisp_QC_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=himem

# ---------------------------------------------------------------------------
# Load R module
# ---------------------------------------------------------------------------
module purge
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.0

# ---------------------------------------------------------------------------
# Create log directory if it does not exist
# ---------------------------------------------------------------------------
mkdir -p /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/QC/logs

# ---------------------------------------------------------------------------
# Run Stage 1: QC (condition-specific mito thresholds) + doublet removal + merge
# ---------------------------------------------------------------------------
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/NewVersion"

Rscript "${SCRIPT_DIR}/01_VehicleCisplatin_QC_Doublet_Merge.R"
