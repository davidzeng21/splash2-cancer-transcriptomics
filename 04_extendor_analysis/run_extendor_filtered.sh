#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J extendor_filtered
#SBATCH -o extendor_filtered_%j.out
#SBATCH -t 4:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=64G

# Change to working directory
cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched

# Load singularity module (needed for satc_dump)
ml PDC singularity

# Activate conda environment
source ~/.bashrc
conda activate myenv

# Define paths
INPUT_FILE="input_CPTAC.txt"
SCORES_FILE="2025-12-04_test.after_correction.scores.filtered_for_blast.tsv"
SATC_FOLDER="2025-12-04_test_satc"
SAMPLE_MAPPING="sample_name_to_id.mapping.txt"
SATC_DUMP_WRAPPER="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/scripts/project_specific/cptac_lung/satc_dump_wrapper.sh"
OUT_FOLDER="extendor_analysis_filtered_blast"

# Create output folder
mkdir -p ${OUT_FOLDER}

# Step 1: Extract anchor list from filtered scores file
echo "=== Step 1: Extracting anchor list ==="
cut -f1 ${SCORES_FILE} | tail -n +2 > ${OUT_FOLDER}/anchors_to_analyze.txt
echo "Extracted $(wc -l < ${OUT_FOLDER}/anchors_to_analyze.txt) anchors"

# Step 2: Run extendor analysis
echo ""
echo "=== Step 2: Running extendor analysis ==="
python /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/scripts/main_analysis/extendor_analysis.py \
    --inputFile ${INPUT_FILE} \
    --scoresFile ${SCORES_FILE} \
    --satcFolder ${SATC_FOLDER} \
    --sampleMapping ${SAMPLE_MAPPING} \
    --satcDumpWrapper ${SATC_DUMP_WRAPPER} \
    --outFolder ${OUT_FOLDER} \
    --anchorFile ${OUT_FOLDER}/anchors_to_analyze.txt

# Step 3: Create BLAST FASTA files
echo ""
echo "=== Step 3: Creating BLAST FASTA files ==="
mkdir -p ${OUT_FOLDER}/blast_analysis
cd ${OUT_FOLDER}/blast_analysis

python3 << 'EOF'
import pandas as pd
import os

# Read the proportions file
df = pd.read_csv('../extendor_proportions.tsv', sep='\t')

print(f"Total extendors: {len(df)}")

# Top tumor-enriched (zero in normal, high in tumor)
tumor_enriched = df[(df['count_normal'] == 0) & (df['count_tumor'] >= 50)].nlargest(50, 'count_tumor')

# Top normal-enriched (zero in tumor, high in normal)  
normal_enriched = df[(df['count_tumor'] == 0) & (df['count_normal'] >= 50)].nlargest(50, 'count_normal')

# Also get top by log2FC (both directions)
df_both = df[(df['count_tumor'] > 0) & (df['count_normal'] > 0)]
if len(df_both) > 0:
    top_tumor_fc = df_both.nlargest(50, 'log2FC')
    top_normal_fc = df_both.nsmallest(50, 'log2FC')
else:
    top_tumor_fc = pd.DataFrame()
    top_normal_fc = pd.DataFrame()

print(f"Tumor-enriched (exclusive): {len(tumor_enriched)}")
print(f"Normal-enriched (exclusive): {len(normal_enriched)}")
print(f"Top tumor by log2FC: {len(top_tumor_fc)}")
print(f"Top normal by log2FC: {len(top_normal_fc)}")

# Write tumor-enriched FASTA
with open('tumor_enriched.fasta', 'w') as f:
    for i, row in enumerate(tumor_enriched.itertuples(), 1):
        # Write full extendor (anchor + target)
        f.write(f">tumor_exclusive_{i:02d}_counts{int(row.count_tumor)}\n")
        f.write(f"{row.anchor}{row.target}\n")

# Write normal-enriched FASTA
with open('normal_enriched.fasta', 'w') as f:
    for i, row in enumerate(normal_enriched.itertuples(), 1):
        f.write(f">normal_exclusive_{i:02d}_counts{int(row.count_normal)}\n")
        f.write(f"{row.anchor}{row.target}\n")

# Write combined file with all differential extendors
with open('top_extendors_combined.fasta', 'w') as f:
    # Tumor exclusive
    for i, row in enumerate(tumor_enriched.itertuples(), 1):
        f.write(f">tumor_exclusive_{i:02d}_counts{int(row.count_tumor)}\n")
        f.write(f"{row.anchor}{row.target}\n")
    # Normal exclusive
    for i, row in enumerate(normal_enriched.itertuples(), 1):
        f.write(f">normal_exclusive_{i:02d}_counts{int(row.count_normal)}\n")
        f.write(f"{row.anchor}{row.target}\n")
    # Top tumor by FC
    for i, row in enumerate(top_tumor_fc.itertuples(), 1):
        f.write(f">tumor_fc_{i:02d}_log2fc{row.log2FC:.2f}\n")
        f.write(f"{row.anchor}{row.target}\n")
    # Top normal by FC
    for i, row in enumerate(top_normal_fc.itertuples(), 1):
        f.write(f">normal_fc_{i:02d}_log2fc{row.log2FC:.2f}\n")
        f.write(f"{row.anchor}{row.target}\n")

total_seqs = len(tumor_enriched) + len(normal_enriched) + len(top_tumor_fc) + len(top_normal_fc)
print(f"\nCreated FASTA files:")
print(f"  - tumor_enriched.fasta ({len(tumor_enriched)} sequences)")
print(f"  - normal_enriched.fasta ({len(normal_enriched)} sequences)")
print(f"  - top_extendors_combined.fasta ({total_seqs} total sequences)")

EOF

echo ""
echo "=== Extendor analysis completed! ==="
echo "Output folder: ${OUT_FOLDER}"
echo ""
echo "Next step: Submit BLAST job with run_blast_human.sh"

