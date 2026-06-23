#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J threshold_sweep
#SBATCH -o threshold_sweep_%j.out
#SBATCH -t 24:00:00
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 256
#SBATCH --mem=440G

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
FULL_SCORES="2025-12-04_test.after_correction.scores.tsv"
INPUT_FILE="input_CPTAC.txt"
SATC_FOLDER="2025-12-04_test_satc"
SAMPLE_MAPPING="sample_name_to_id.mapping.txt"
SATC_DUMP_WRAPPER="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/scripts/project_specific/cptac_lung/satc_dump_wrapper.sh"
OUT_FOLDER="extendor_threshold_sweep"

mkdir -p ${OUT_FOLDER}

# ============================================================
# STEP 1: Filter superset anchors (samples>=100, es>=0.8)
# ============================================================
echo "=== Step 1: Filtering superset anchors (samples>=100, es>=0.8) ==="
echo "Date: $(date)"

SUPERSET_SCORES="${OUT_FOLDER}/superset_scores.tsv"

if [ ! -f "${SUPERSET_SCORES}" ]; then
    awk -F'\t' 'NR==1 || ($8 >= 100 && $4 >= 0.8)' "${FULL_SCORES}" > "${SUPERSET_SCORES}"
    echo "Filtered to $(( $(wc -l < ${SUPERSET_SCORES}) - 1 )) anchors"
else
    echo "Superset scores already exist: $(( $(wc -l < ${SUPERSET_SCORES}) - 1 )) anchors"
fi

# Extract anchor list
ANCHOR_FILE="${OUT_FOLDER}/anchors_superset.txt"
cut -f1 "${SUPERSET_SCORES}" | tail -n +2 > "${ANCHOR_FILE}"
echo "Anchor list: $(wc -l < ${ANCHOR_FILE}) anchors"

# ============================================================
# STEP 2: Run SATC dump for superset anchors
# ============================================================
echo ""
echo "=== Step 2: Running SATC dump ==="
echo "Date: $(date)"

python /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/scripts/main_analysis/extendor_analysis.py \
    --inputFile ${INPUT_FILE} \
    --scoresFile "${SUPERSET_SCORES}" \
    --satcFolder ${SATC_FOLDER} \
    --sampleMapping ${SAMPLE_MAPPING} \
    --satcDumpWrapper ${SATC_DUMP_WRAPPER} \
    --outFolder ${OUT_FOLDER} \
    --anchorFile "${ANCHOR_FILE}"

echo ""
echo "=== Step 2 complete ==="
echo "Date: $(date)"

# ============================================================
# STEP 3: Run threshold sweep
# ============================================================
echo ""
echo "=== Step 3: Running threshold sweep ==="
echo "Date: $(date)"

python3 /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/threshold_sweep.py \
    --proportions "${OUT_FOLDER}/extendor_proportions.tsv" \
    --scores "${SUPERSET_SCORES}" \
    --inputFile "${INPUT_FILE}" \
    --outFolder "${OUT_FOLDER}"

echo ""
echo "=== All steps complete ==="
echo "Date: $(date)"
