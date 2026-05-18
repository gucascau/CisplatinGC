#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --partition=himem
#SBATCH --cpus-per-task=24
#SBATCH --job-name=CellTypeCorr
#SBATCH --output=slurm_%j.out
set -e

# loading the module
module purge
ml GCC/9.3.0
ml OpenMPI/4.0.3
ml R/4.4.0

CWD=/vast0/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/Evaluation

cd $CWD
Rscript $CWD/CellType_AverageExpression_Correlation.R
