#!/bin/bash
#SBATCH --job-name=VehCisp_Cluster
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Clustering/logs/VehCisp_Cluster_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Clustering/logs/VehCisp_Cluster_%j.err
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
mkdir -p /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Clustering/logs

# ---------------------------------------------------------------------------
# Run Stage 2: normalization + Harmony integration + clustering + UMAP
# ---------------------------------------------------------------------------
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/NewVersion"

Rscript "${SCRIPT_DIR}/02_VehicleCisplatin_Integration_Clustering.R"
