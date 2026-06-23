#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J extendor_2026_02_20
#SBATCH -o extendor_2026_02_20_%j.out
#SBATCH -t 8:00:00
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 256
#SBATCH --mem=400G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

set -euo pipefail

# Change to working directory
cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched

# Load singularity module (needed for satc_dump)
ml PDC singularity

# Activate conda environment
source ~/.bashrc
conda activate myenv

# ============================================================
# Paths
# ============================================================
INPUT_FILE="input_CPTAC_107.txt"
SCORES_FILE="2026-02-20_min_50_samples/2026-02-20_min_50_samples.after_correction.scores.tsv"
SATC_FOLDER="2026-02-20_min_50_samples/2026-02-20_min_50_samples_satc"
SAMPLE_MAPPING="2026-02-20_min_50_samples/sample_name_to_id.mapping.txt"
SATC_DUMP_WRAPPER="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/satc_dump_wrapper.sh"
OUT_FOLDER="2026-02-20_min_50_samples/extendor_analysis"

mkdir -p "${OUT_FOLDER}"

# ============================================================
# Step 1: Filter anchors by effect_size_bin >= 0.2
# ============================================================
echo "=== Step 1: Filtering anchors (effect_size_bin >= 0.2) ==="
echo "Date: $(date)"

ANCHOR_FILE="${OUT_FOLDER}/anchors_es_ge0.2.txt"

if [ ! -f "${ANCHOR_FILE}" ]; then
    awk -F'\t' 'NR > 1 && $4 >= 0.2 {print $1}' "${SCORES_FILE}" > "${ANCHOR_FILE}"
    echo "Filtered to $(wc -l < ${ANCHOR_FILE}) anchors with effect_size_bin >= 0.2"
else
    echo "Anchor file already exists: $(wc -l < ${ANCHOR_FILE}) anchors"
fi

# ============================================================
# Step 2: Run extendor analysis
# ============================================================
echo ""
echo "=== Step 2: Running extendor analysis ==="
echo "Date: $(date)"

python extendor_analysis.py \
    --inputFile "${INPUT_FILE}" \
    --scoresFile "${SCORES_FILE}" \
    --satcFolder "${SATC_FOLDER}" \
    --sampleMapping "${SAMPLE_MAPPING}" \
    --satcDumpWrapper "${SATC_DUMP_WRAPPER}" \
    --outFolder "${OUT_FOLDER}" \
    --anchorFile "${ANCHOR_FILE}"

# --- Alternatives ---
# Skip satc_dump if dumps already exist in ${OUT_FOLDER}/satc_dumps/:
#   add --skipSATCDump

# Skip everything up to and including proportion calculation
# (re-run plotting + summary stats only, using existing extendor_proportions.tsv):
#   replace the command above with:
# python extendor_analysis.py \
#     --inputFile "${INPUT_FILE}" \
#     --outFolder "${OUT_FOLDER}" \
#     --skipPropCalc

echo ""
echo "=== Analysis complete ==="
echo "Date: $(date)"
echo "Output folder: ${OUT_FOLDER}"
