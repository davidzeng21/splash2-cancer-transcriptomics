#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J fc_summary_refseq
#SBATCH -o fc_summary_refseq_%j.out
#SBATCH -t 1:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=64G

cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/extendor_analysis_filtered_blast/blast_analysis

echo "=== FC Extendors Summary (RefSeq RNA only) ==="
echo "Date: $(date)"
echo ""
conda activate myenv
python3 << 'EOF'
import pandas as pd

# Read extendor proportions
print("Loading extendor_proportions.tsv...")
prop_df = pd.read_csv('../extendor_proportions.tsv', sep='\t')
print(f"Loaded: {len(prop_df)} extendors")

# Read RefSeq RNA BLAST hits only
print("\nLoading blast_refseq_rna_human_v2.txt...")
cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']
blast_df = pd.read_csv('blast_refseq_rna_human_v2.txt', sep='\t', names=cols)
blast_df['database'] = 'refseq_rna'
print(f"Loaded: {len(blast_df)} BLAST hits")

# Read the FASTA to get mapping of query IDs to anchor+target
print("\nParsing FASTA file...")
fasta_mapping = {}
with open('top_extendors_expanded.fasta', 'r') as f:
    current_id = None
    for line in f:
        if line.startswith('>'):
            current_id = line.strip()[1:]
        else:
            seq = line.strip()
            anchor = seq[:27]
            target = seq[27:]
            fasta_mapping[current_id] = {'anchor': anchor, 'target': target}

print(f"Parsed: {len(fasta_mapping)} sequences")

# Get tumor_fc and normal_fc extendors
tumor_fc_ids = [k for k in fasta_mapping.keys() if k.startswith('tumor_fc_')]
normal_fc_ids = [k for k in fasta_mapping.keys() if k.startswith('normal_fc_')]

print(f"\nTumor FC extendors: {len(tumor_fc_ids)}")
print(f"Normal FC extendors: {len(normal_fc_ids)}")

# Build summary
print("\nBuilding summary...")
results = []

for qid in sorted(tumor_fc_ids) + sorted(normal_fc_ids):
    info = fasta_mapping[qid]
    anchor = info['anchor']
    target = info['target']
    
    # Get proportion data
    prop_row = prop_df[(prop_df['anchor'] == anchor) & (prop_df['target'] == target)]
    
    if len(prop_row) > 0:
        prop_row = prop_row.iloc[0]
        count_tumor = prop_row['count_tumor']
        count_normal = prop_row['count_normal']
        prop_tumor = prop_row['prop_tumor']
        prop_normal = prop_row['prop_normal']
        log2fc = prop_row['log2FC']
    else:
        count_tumor = count_normal = prop_tumor = prop_normal = log2fc = 'NA'
    
    # Get BLAST hits (RefSeq RNA only)
    hits = blast_df[blast_df['qseqid'] == qid].sort_values('bitscore', ascending=False)
    
    if len(hits) > 0:
        best_hit = hits.iloc[0]
        blast_hit = best_hit['stitle']
        blast_pident = best_hit['pident']
        blast_evalue = best_hit['evalue']
        blast_length = best_hit['length']
    else:
        blast_hit = 'No hit'
        blast_pident = blast_evalue = blast_length = 'NA'
    
    results.append({
        'query_id': qid,
        'anchor': anchor,
        'target': target,
        'count_tumor': count_tumor,
        'count_normal': count_normal,
        'prop_tumor': prop_tumor,
        'prop_normal': prop_normal,
        'log2FC': log2fc,
        'blast_hit': blast_hit,
        'blast_pident': blast_pident,
        'blast_evalue': blast_evalue,
        'blast_length': blast_length
    })

# Create DataFrame and save
result_df = pd.DataFrame(results)
result_df.to_csv('fc_extendors_summary_refseq_only.tsv', sep='\t', index=False)
print(f"\nSaved fc_extendors_summary_refseq_only.tsv with {len(result_df)} extendors")

# Print summary
print("\n" + "="*100)
print("TOP 20 TUMOR FOLD CHANGE EXTENDORS (RefSeq RNA hits only)")
print("="*100)
tumor_df = result_df[result_df['query_id'].str.contains('tumor_fc')].head(20)
for _, row in tumor_df.iterrows():
    print(f"\n{row['query_id']}:")
    print(f"  Counts: Tumor={row['count_tumor']}, Normal={row['count_normal']}, log2FC={row['log2FC']:.2f}" if row['log2FC'] != 'NA' else f"  Counts: NA")
    print(f"  BLAST: {str(row['blast_hit'])[:80]}...")
    print(f"         Identity={row['blast_pident']}%, Length={row['blast_length']}, E-value={row['blast_evalue']}")

print("\n" + "="*100)
print("TOP 20 NORMAL FOLD CHANGE EXTENDORS (RefSeq RNA hits only)")
print("="*100)
normal_df = result_df[result_df['query_id'].str.contains('normal_fc')].head(20)
for _, row in normal_df.iterrows():
    print(f"\n{row['query_id']}:")
    print(f"  Counts: Tumor={row['count_tumor']}, Normal={row['count_normal']}, log2FC={row['log2FC']:.2f}" if row['log2FC'] != 'NA' else f"  Counts: NA")
    print(f"  BLAST: {str(row['blast_hit'])[:80]}...")
    print(f"         Identity={row['blast_pident']}%, Length={row['blast_length']}, E-value={row['blast_evalue']}")

# Summary statistics
print("\n" + "="*100)
print("SUMMARY STATISTICS")
print("="*100)
tumor_with_hits = result_df[(result_df['query_id'].str.contains('tumor_fc')) & (result_df['blast_hit'] != 'No hit')]
normal_with_hits = result_df[(result_df['query_id'].str.contains('normal_fc')) & (result_df['blast_hit'] != 'No hit')]
print(f"Tumor FC extendors with RefSeq RNA hits: {len(tumor_with_hits)}/{len(tumor_fc_ids)}")
print(f"Normal FC extendors with RefSeq RNA hits: {len(normal_with_hits)}/{len(normal_fc_ids)}")

# Hit categories
print("\n=== HIT CATEGORIES ===")
print("\nTumor FC hits:")
tumor_hits = result_df[(result_df['query_id'].str.contains('tumor_fc')) & (result_df['blast_hit'] != 'No hit')]['blast_hit']
for hit, count in tumor_hits.value_counts().head(15).items():
    print(f"  {count:3d} - {hit[:70]}...")

print("\nNormal FC hits:")
normal_hits = result_df[(result_df['query_id'].str.contains('normal_fc')) & (result_df['blast_hit'] != 'No hit')]['blast_hit']
for hit, count in normal_hits.value_counts().head(15).items():
    print(f"  {count:3d} - {hit[:70]}...")

EOF

echo ""
echo "=== Analysis Complete ==="
echo "Date: $(date)"
