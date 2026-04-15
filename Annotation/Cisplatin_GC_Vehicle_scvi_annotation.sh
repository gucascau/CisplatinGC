#!/bin/bash
#SBATCH --job-name=CispGC_scVI
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs/CispGC_scVI_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs/CispGC_scVI_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=himem
# To use a GPU node instead, comment out the two lines above and uncomment:
# #SBATCH --partition=gpu
# #SBATCH --gres=gpu:1

# ---------------------------------------------------------------------------
# Create log directory if it does not exist
# ---------------------------------------------------------------------------
mkdir -p /home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Integration/IntegrationGC/Annotation/logs

# ---------------------------------------------------------------------------
# Activate the conda environment
# ---------------------------------------------------------------------------
source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

# ---------------------------------------------------------------------------
# Print environment info for reproducibility
# ---------------------------------------------------------------------------
echo "Python:      $(python --version)"
echo "scvi-tools:  $(python -c 'import scvi; print(scvi.__version__)')"
echo "Hostname:    $(hostname)"
echo "Start time:  $(date)"
echo ""

# ---------------------------------------------------------------------------
# Run the annotation script
# ---------------------------------------------------------------------------
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts/Integration"

python "${SCRIPT_DIR}/Cisplatin_GC_Vehicle_scvi_annotation.py"

echo ""
echo "End time: $(date)"
