#!/bin/bash
#SBATCH --job-name=VehCisp_Annot
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/logs/VehCisp_Annot_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/logs/VehCisp_Annot_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=himem

# ---------------------------------------------------------------------------
# Load R module
# ---------------------------------------------------------------------------
# module purge
# ml GCC/9.3.0
# ml OpenMPI/4.0.3
# ml R/4.4.0
conda activate scrna_env

# ---------------------------------------------------------------------------
# Create log directory if it does not exist
# ---------------------------------------------------------------------------
mkdir -p /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/VehicleCisplatin/Annotation/logs

# ---------------------------------------------------------------------------
# Run Stage 3: Seurat MKA/Lake label transfer + manual annotation + h5ad
# export for scANVI. NOTE: this script must be run TWICE —
#   1st run: writes the h5ad, then exits (scANVI CSVs don't exist yet)
#   ---- in between, submit the two scANVI jobs and wait for them ----
#     sbatch 03b_VehicleCisplatin_scvi_annotation_MKAOnly.sh
#     sbatch 03c_VehicleCisplatin_scvi_annotation_LakeOnly.sh
#   2nd run: reads back the scANVI CSVs and saves the final 5-method object
# ---------------------------------------------------------------------------
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/NewVersion"

Rscript "${SCRIPT_DIR}/03_VehicleCisplatin_MultiMethod_Annotation.R"
