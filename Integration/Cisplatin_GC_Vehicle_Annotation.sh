#!/bin/bash
#SBATCH --job-name=CispGC_Annot
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs/CispGC_Annot_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs/CispGC_Annot_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=himem

# ---------------------------------------------------------------------------
# Load R module (adjust version/module name to match your HPC environment)
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
# Run the annotation script
# ---------------------------------------------------------------------------
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/Integration"

Rscript "${SCRIPT_DIR}/Cisplatin_GC_Vehicle_Annotation.R"
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/Integration"

Rscript "${SCRIPT_DIR}/Cisplatin_GC_Vehicle_Annotation.R"
