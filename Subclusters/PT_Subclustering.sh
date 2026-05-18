#!/bin/bash
#SBATCH --job-name=PT_Subcluster
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/logs/PT_Subcluster_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/logs/PT_Subcluster_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=himem

module purge
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.0

mkdir -p /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/logs

SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/Subclusters"

Rscript "${SCRIPT_DIR}/PT_Subclustering.R"
