#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J extendor_top_effect
#SBATCH -o extendor_top_effect_%j.out
#SBATCH -t 4:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 128

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
OUT_FOLDER="extendor_top_effect_size"

# Extract anchors from top effect size file
echo "Extracting anchors from top_effect_size_bin.tsv..."
mkdir -p ${OUT_FOLDER}
tail -n +2 2025-12-04_test.after_correction.scores.top_effect_size_bin.tsv | awk -F'\t' '{print $1}' | sed 's/^[[:space:]]*//' > ${OUT_FOLDER}/anchors_top_effect_size.txt
echo "Extracted $(wc -l < ${OUT_FOLDER}/anchors_top_effect_size.txt) anchors"

# Remove any existing symlink to avoid using wrong data
rm -f ${OUT_FOLDER}/satc_dumps 2>/dev/null

# Run the analysis - MUST run satc_dump with new anchors
python extendor_analysis.py \
    --inputFile ${INPUT_FILE} \
    --scoresFile ${SCORES_FILE} \
    --satcFolder ${SATC_FOLDER} \
    --sampleMapping ${SAMPLE_MAPPING} \
    --satcDumpWrapper ${SATC_DUMP_WRAPPER} \
    --outFolder ${OUT_FOLDER} \
    --anchorFile ${OUT_FOLDER}/anchors_top_effect_size.txt

echo "Top effect size extendor analysis completed!"

