#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J refilter_no_human_hits
#SBATCH -o refilter_no_human_hits_%j.out
#SBATCH -t 1:00:00
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G

ml PDC
source /cfs/klemming/home/j/jlzeng/myenv/bin/activate

SCRIPT_DIR="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/2026-02-20_min_50_samples/extendor_analysis/blast_tumor_specific"

echo "=== Re-filter no_human_hits (ge3 + ge2) ==="
echo "Date: $(date)"
echo "Script: ${SCRIPT_DIR}/refilter_no_human_hits.py"
echo ""

python3 ${SCRIPT_DIR}/refilter_no_human_hits.py

echo ""
echo "=== Done ==="
echo "Date: $(date)"
