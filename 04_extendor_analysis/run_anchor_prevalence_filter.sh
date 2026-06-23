#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J anchor_prev
#SBATCH -o anchor_prevalence_filter_%j.out
#SBATCH -t 1:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

set -euo pipefail

cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched

source /cfs/klemming/home/j/jlzeng/myenv/bin/activate

SCRIPT="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/ref_free_cancer/04_extendor_analysis/anchor_prevalence_filter.py"

echo "=== Anchor prevalence filter ==="
echo "Date: $(date)"

python "${SCRIPT}" \
    --inputFile input_CPTAC_107.txt \
    --dumpFolder 2026-02-20_min_50_samples/extendor_analysis/satc_dumps \
    --propParquet 2026-02-20_min_50_samples/extendor_analysis/extendor_proportions.parquet \
    --outFolder 2026-02-20_min_50_samples/extendor_analysis/anchor_prevalence \
    --minTumorFrac 0.10 \
    --minTumorFracRounding ceil \
    --fcThresholds 2,3

echo ""
echo "=== Done ==="
echo "Date: $(date)"
