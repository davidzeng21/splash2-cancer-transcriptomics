#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J extendor_analysis
#SBATCH -o extendor_analysis_%j.out
#SBATCH -t 8:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

# Change to working directory
cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched

# Load singularity module (needed for satc_dump)
ml PDC singularity

# Activate conda environment
source ~/.bashrc
conda activate myenv

# Define paths
INPUT_FILE="input_CPTAC.txt"
SCORES_FILE="2025-12-04_test.after_correction.scores.tsv"
SATC_FOLDER="2025-12-04_test_satc"
SAMPLE_MAPPING="sample_name_to_id.mapping.txt"
SATC_DUMP_WRAPPER="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/satc_dump_wrapper.sh"
OUT_FOLDER="extendor_analysis_results"

# Number of top anchors to analyze (reduce for faster testing)
TOP_N=10000

# Run the analysis
python extendor_analysis.py \
    --inputFile ${INPUT_FILE} \
    --scoresFile ${SCORES_FILE} \
    --satcFolder ${SATC_FOLDER} \
    --sampleMapping ${SAMPLE_MAPPING} \
    --satcDumpWrapper ${SATC_DUMP_WRAPPER} \
    --outFolder ${OUT_FOLDER} \
    --topN ${TOP_N}

# If you want to rerun without re-dumping SATC files, add --skipSATCDump flag:
# python extendor_analysis.py \
#     --inputFile ${INPUT_FILE} \
#     --scoresFile ${SCORES_FILE} \
#     --satcFolder ${SATC_FOLDER} \
#     --sampleMapping ${SAMPLE_MAPPING} \
#     --satcDumpWrapper ${SATC_DUMP_WRAPPER} \
#     --outFolder ${OUT_FOLDER} \
#     --topN ${TOP_N} \
#     --skipSATCDump

echo "Extendor analysis completed!"

