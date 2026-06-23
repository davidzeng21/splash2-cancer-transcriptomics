#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J extendor_blast
#SBATCH -o extendor_blast_%j.out
#SBATCH -t 4:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=64G

# Load BLAST+
ml PDC blast+

cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/extendor_top_effect_size/blast_analysis

# Database paths
REFSEQ_RNA="/sw/data/blast_databases/refseq_rna"
CORE_NT="/sw/data/blast_databases/core_nt"

echo "=== BLAST Analysis of Top Differential Extendors ==="
echo "Date: $(date)"
echo ""

# BLAST against RefSeq RNA (transcriptome)
echo "--- BLASTing tumor-enriched extendors against RefSeq RNA ---"
blastn -query top_extendors_combined.fasta \
    -db ${REFSEQ_RNA} \
    -out blast_refseq_rna.txt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 5 \
    -num_threads 16 \
    -evalue 1e-5 \
    -word_size 11

echo "RefSeq RNA BLAST complete: $(wc -l < blast_refseq_rna.txt) hits"

# BLAST against core_nt (nucleotide database)
echo ""
echo "--- BLASTing tumor-enriched extendors against core_nt ---"
blastn -query top_extendors_combined.fasta \
    -db ${CORE_NT} \
    -out blast_core_nt.txt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 5 \
    -num_threads 16 \
    -evalue 1e-5 \
    -word_size 11

echo "Core NT BLAST complete: $(wc -l < blast_core_nt.txt) hits"

# Create summary report
echo ""
echo "=== Creating Summary Report ==="

python3 << 'EOF'
import pandas as pd

# Column names for BLAST outfmt 6 with stitle
cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

# Read RefSeq RNA results
try:
    df_rna = pd.read_csv('blast_refseq_rna.txt', sep='\t', names=cols)
    df_rna['database'] = 'refseq_rna'
except:
    df_rna = pd.DataFrame(columns=cols + ['database'])

# Read core_nt results
try:
    df_nt = pd.read_csv('blast_core_nt.txt', sep='\t', names=cols)
    df_nt['database'] = 'core_nt'
except:
    df_nt = pd.DataFrame(columns=cols + ['database'])

# Combine
df = pd.concat([df_rna, df_nt], ignore_index=True)

if len(df) > 0:
    # Get best hit per query
    df_best = df.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')
    
    # Separate tumor and normal
    tumor_hits = df_best[df_best['qseqid'].str.contains('tumor')]
    normal_hits = df_best[df_best['qseqid'].str.contains('normal')]
    
    print("\n" + "="*80)
    print("TUMOR-ENRICHED EXTENDORS - Best BLAST Hits")
    print("="*80)
    for _, row in tumor_hits.iterrows():
        print(f"\n{row['qseqid']}:")
        print(f"  Hit: {row['stitle'][:80]}...")
        print(f"  Identity: {row['pident']:.1f}%, Length: {row['length']}, E-value: {row['evalue']:.2e}")
    
    print("\n" + "="*80)
    print("NORMAL-ENRICHED EXTENDORS - Best BLAST Hits")
    print("="*80)
    for _, row in normal_hits.iterrows():
        print(f"\n{row['qseqid']}:")
        print(f"  Hit: {row['stitle'][:80]}...")
        print(f"  Identity: {row['pident']:.1f}%, Length: {row['length']}, E-value: {row['evalue']:.2e}")
    
    # Save full results
    df_best.to_csv('blast_best_hits.tsv', sep='\t', index=False)
    print(f"\n\nFull results saved to blast_best_hits.tsv")
else:
    print("No BLAST hits found!")

EOF

echo ""
echo "=== BLAST Analysis Complete ==="



