#!/bin/bash
#SBATCH --job-name=PT_scVI_Lake
#SBATCH --output=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/logs/PT_scVI_Lake_%j.out
#SBATCH --error=/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT/logs/PT_scVI_Lake_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --partition=himem
# To use a GPU node instead, comment out the two lines above and uncomment:
# #SBATCH --partition=gpu
# #SBATCH --gres=gpu:1

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
OUT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Results/Subclusters/PT"
SCRIPT_DIR="/home/gdzepedaorozcolab/lab/xxw004/Projects/RawDZscRNAseq/Scripts"

mkdir -p "${OUT_DIR}/logs"

# ---------------------------------------------------------------------------
# Step 1: Export PT_Subclustered_Annotated.RDS → h5ad (if not already done)
# ---------------------------------------------------------------------------
H5AD_PATH="${OUT_DIR}/PT_Subclustered_Annotated.h5ad"

if [ ! -f "${H5AD_PATH}" ]; then
    echo "=== Step 1: Exporting PT RDS to h5ad ==="
    echo "Start time: $(date)"

    module purge
    ml GCC/9.3.0
    ml OpenMPI/4.0.3
    ml R/4.4.0

    Rscript "${SCRIPT_DIR}/Subclusters/PT_Export_h5ad.R"

    echo "h5ad export done: $(date)"
else
    echo "=== Step 1: h5ad already exists — skipping export ==="
    echo "  ${H5AD_PATH}"
fi

# ---------------------------------------------------------------------------
# Step 2: Run scVI / scANVI annotation
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Running PT scVI / scANVI annotation ==="
echo "Start time: $(date)"

source /home/gdbecknelllab/xxw004/Software/miniconda3/etc/profile.d/conda.sh
conda activate cell2loc_env

echo "Python:      $(python --version)"
echo "scvi-tools:  $(python -c 'import scvi; print(scvi.__version__)')"
echo "Hostname:    $(hostname)"

python "${SCRIPT_DIR}/Annotation/PT_scvi_annotation_LakeOnly.py"

echo ""
echo "End time: $(date)"
