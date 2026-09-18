#!/bin/bash
#SBATCH --job-name=VehCisp_Ensemble
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/IntegrateL1/logs/VehCisp_Ensemble_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/IntegrateL1/logs/VehCisp_Ensemble_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
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
mkdir -p /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/IntegrateL1/logs

# ---------------------------------------------------------------------------
# Run Stage 4: L1 harmonisation + weighted ensemble + high-confidence
# FinalAnnotation_HC
# ---------------------------------------------------------------------------
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/NewVersion"

Rscript "${SCRIPT_DIR}/04_VehicleCisplatin_WeightedEnsemble_HighConfidence_Annotation.R"
