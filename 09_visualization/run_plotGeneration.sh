#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J plot_tcga
#SBATCH -o plot_tcga_%j.out
#SBATCH -t 4:00:00
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

# Change to working directory
cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/tcga_endometrial

# Load singularity module (needed for satc_dump)
ml PDC singularity

# Activate conda environment
source ~/.bashrc
conda activate myenv

# Define paths
DATASET_NAME="TCGA_Endometrial"
OUT_FOLDER="./plot_results"
METADATA_PATH="./metadata.tsv"
SATC_FOLDER="./result_satc"
PV_DF_PATH="./result.after_correction.scores.tsv"
SAMPLE_MAPPING="./sample_name_to_id.mapping.txt"
ANCHOR_FILE="./anchors_to_plot.txt"
# satc_dump must be run via singularity due to GLIBC dependencies
# Using wrapper script that calls satc_dump inside the container
SATC_DUMP="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/tcga_endometrial/satc_dump_wrapper.sh"

# Create output folder if it doesn't exist
mkdir -p ${OUT_FOLDER}

# Check if anchor file exists
# Note: Create ${ANCHOR_FILE} with your anchors of interest (one per line, no header)
if [ ! -f ${ANCHOR_FILE} ]; then
    echo "ERROR: ${ANCHOR_FILE} does not exist. Please create it with your anchors to plot (one per line, no header)."
    exit 1
fi

# Run plotGeneration.py
# Remove --skipSATC flag on first run to dump .satc files
# Add --skipSATC flag if rerunning with same anchors
python plotGeneration.py \
    --dsName ${DATASET_NAME} \
    --outFolder ${OUT_FOLDER} \
    --metadataPath ${METADATA_PATH} \
    --satcFolder ${SATC_FOLDER} \
    --pvDfPath ${PV_DF_PATH} \
    --sampleMappingTxt ${SAMPLE_MAPPING} \
    --anchorFile ${ANCHOR_FILE} \
    --satc_dump_file ${SATC_DUMP}

echo "Plot generation completed!"

