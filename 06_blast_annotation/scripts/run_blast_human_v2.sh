#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J blast_human_v2
#SBATCH -o blast_human_v2_%j.out
#SBATCH -t 24:00:00
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=240G

# Load BLAST+
ml PDC blast+

# Set BLASTDB environment variable to include taxonomy database
export BLASTDB=/sw/data/blast_databases

# Change to blast analysis directory
cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/extendor_analysis_filtered_blast/blast_analysis

# Input file
QUERY_FILE="top_extendors_expanded.fasta"

# Database paths
CORE_NT="/sw/data/blast_databases/core_nt"
REFSEQ_RNA="/sw/data/blast_databases/refseq_rna"

# Human taxonomy ID
HUMAN_TAXID=9606

echo "=== BLAST Analysis of Differential Extendors (Human Only) ==="
echo "Date: $(date)"
echo "Query file: ${QUERY_FILE}"
echo "Sequences: $(grep -c '^>' ${QUERY_FILE})"
echo "BLASTDB: ${BLASTDB}"
echo ""

# Verify taxonomy database is accessible
echo "Checking taxonomy database..."
ls -la ${BLASTDB}/taxonomy4blast.sqlite3 ${BLASTDB}/taxdb.btd ${BLASTDB}/taxdb.bti
echo ""

# BLAST against core_nt (human only)
echo "--- BLASTing against core_nt (human taxid ${HUMAN_TAXID}) ---"
echo "Start time: $(date)"
blastn -query ${QUERY_FILE} \
    -db ${CORE_NT} \
    -taxids ${HUMAN_TAXID} \
    -out blast_core_nt_human_v2.txt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 10 \
    -num_threads 64 \
    -evalue 1e-3 \
    -word_size 11

echo "Core NT BLAST complete: $(wc -l < blast_core_nt_human_v2.txt) hits"
echo "End time: $(date)"

# BLAST against RefSeq RNA (human only)
echo ""
echo "--- BLASTing against RefSeq RNA (human taxid ${HUMAN_TAXID}) ---"
echo "Start time: $(date)"
blastn -query ${QUERY_FILE} \
    -db ${REFSEQ_RNA} \
    -taxids ${HUMAN_TAXID} \
    -out blast_refseq_rna_human_v2.txt \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    -max_target_seqs 10 \
    -num_threads 64 \
    -evalue 1e-3 \
    -word_size 11

echo "RefSeq RNA BLAST complete: $(wc -l < blast_refseq_rna_human_v2.txt) hits"
echo "End time: $(date)"

# Create summary report
echo ""
echo "=== Creating Summary Report ==="

python3 << 'EOF'
import pandas as pd

# Column names for BLAST outfmt 6 with stitle
cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

# Read core_nt results
try:
    df_nt = pd.read_csv('blast_core_nt_human_v2.txt', sep='\t', names=cols)
    df_nt['database'] = 'core_nt'
    print(f"Core NT hits: {len(df_nt)}")
except:
    df_nt = pd.DataFrame(columns=cols + ['database'])
    print("No Core NT hits")

# Read RefSeq RNA results
try:
    df_rna = pd.read_csv('blast_refseq_rna_human_v2.txt', sep='\t', names=cols)
    df_rna['database'] = 'refseq_rna'
    print(f"RefSeq RNA hits: {len(df_rna)}")
except:
    df_rna = pd.DataFrame(columns=cols + ['database'])
    print("No RefSeq RNA hits")

# Combine
df = pd.concat([df_nt, df_rna], ignore_index=True)

if len(df) > 0:
    # Get best hit per query per database
    df_best = df.sort_values('bitscore', ascending=False).drop_duplicates(['qseqid', 'database'])
    
    # Separate tumor and normal
    tumor_hits = df_best[df_best['qseqid'].str.contains('tumor')]
    normal_hits = df_best[df_best['qseqid'].str.contains('normal')]
    
    # Save full results
    df_best.to_csv('blast_best_hits_human_v2.tsv', sep='\t', index=False)
    df.to_csv('blast_all_hits_human_v2.tsv', sep='\t', index=False)
    
    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY (Human Only)")
    print("="*80)
    print(f"Total queries with human hits: {df_best['qseqid'].nunique()}")
    print(f"Tumor extendors with human hits: {tumor_hits['qseqid'].nunique()}")
    print(f"Normal extendors with human hits: {normal_hits['qseqid'].nunique()}")
    
    # Show top 20 hits for each category
    print("\n" + "="*80)
    print("TOP 20 TUMOR EXTENDOR HITS")
    print("="*80)
    for _, row in tumor_hits.head(20).iterrows():
        print(f"\n{row['qseqid']} ({row['database']}):")
        stitle = str(row['stitle'])[:100] + "..." if len(str(row['stitle'])) > 100 else row['stitle']
        print(f"  Hit: {stitle}")
        print(f"  Identity: {row['pident']:.1f}%, Length: {row['length']}, E-value: {row['evalue']:.2e}")
    
    print("\n" + "="*80)
    print("TOP 20 NORMAL EXTENDOR HITS")
    print("="*80)
    for _, row in normal_hits.head(20).iterrows():
        print(f"\n{row['qseqid']} ({row['database']}):")
        stitle = str(row['stitle'])[:100] + "..." if len(str(row['stitle'])) > 100 else row['stitle']
        print(f"  Hit: {stitle}")
        print(f"  Identity: {row['pident']:.1f}%, Length: {row['length']}, E-value: {row['evalue']:.2e}")
    
    print(f"\n\nResults saved to:")
    print(f"  - blast_best_hits_human_v2.tsv (best hit per query per database)")
    print(f"  - blast_all_hits_human_v2.tsv (all hits)")
else:
    print("No BLAST hits found!")

EOF

echo ""
echo "=== BLAST Analysis Complete ==="
echo "Date: $(date)"

