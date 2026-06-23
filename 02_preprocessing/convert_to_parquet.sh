#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J tsv2parquet
#SBATCH -o tsv2parquet_%j.out
#SBATCH -t 30:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --begin=now
#SBATCH --mail-type=ALL

source ~/.bashrc
conda activate myenv

TSV="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/2026-02-20_min_50_samples/extendor_analysis/extendor_proportions.tsv"
PARQUET="${TSV%.tsv}.parquet"

echo "Converting $TSV -> $PARQUET"
echo "Start: $(date)"

python3 - <<'EOF'
import polars as pl, os

tsv = "/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/2026-02-20_min_50_samples/extendor_analysis/extendor_proportions.tsv"
parquet = tsv.replace('.tsv', '.parquet')

print(f"Streaming {tsv} -> {parquet} with type downcasting ...")
(
    pl.scan_csv(tsv, separator='\t')
    .with_columns([
        pl.col('count_normal').cast(pl.Int32),
        pl.col('count_tumor').cast(pl.Int32),
        pl.col('prop_tumor').cast(pl.Float32),
        pl.col('prop_normal').cast(pl.Float32),
        pl.col('log10FC').cast(pl.Float32),
    ])
    .sink_parquet(parquet)
)
size_gb = os.path.getsize(parquet) / 1e9
print(f"Done! Parquet size: {size_gb:.2f} GB")
EOF

echo "End: $(date)"
