#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J replot_sweep
#SBATCH -o replot_sweep_%j.out
#SBATCH -t 4:00:00
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --mem=240G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jialinz@kth.se

cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched

source ~/.bashrc
conda activate myenv

OUT_FOLDER="extendor_threshold_sweep"

echo "=== Replotting threshold sweep with updated color scale ==="
echo "Date: $(date)"

python3 threshold_sweep.py \
    --proportions "${OUT_FOLDER}/extendor_proportions.tsv" \
    --scores "${OUT_FOLDER}/superset_scores.tsv" \
    --inputFile "input_CPTAC.txt" \
    --outFolder "${OUT_FOLDER}"

echo ""
echo "=== Done ==="
echo "Date: $(date)"
