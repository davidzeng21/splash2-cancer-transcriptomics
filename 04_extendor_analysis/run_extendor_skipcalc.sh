#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J extendor_skipcalc
#SBATCH -o extendor_skipcalc_%j.out
#SBATCH -t 1:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

set -euo pipefail

cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched

source /cfs/klemming/home/j/jlzeng/myenv/bin/activate

echo "=== Re-running plots and summary stats (skipPropCalc) ==="
echo "Date: $(date)"

python extendor_analysis.py \
    --inputFile input_CPTAC_107.txt \
    --outFolder 2026-02-20_min_50_samples/extendor_analysis \
    --skipPropCalc

echo ""
echo "=== Done ==="
echo "Date: $(date)"
