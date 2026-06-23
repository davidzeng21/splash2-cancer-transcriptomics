#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J fc_summary_combined
#SBATCH -o fc_summary_combined_%j.out
#SBATCH -t 1:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

cd /cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/extendor_analysis_filtered_blast/blast_analysis

echo "=== FC Extendors Summary (RefSeq RNA + Core NT fallback) ==="
echo "Date: $(date)"
echo ""

conda activate myenv

python3 << 'EOF'
import pandas as pd
from collections import Counter

# Read extendor proportions
print("Loading extendor_proportions.tsv...")
prop_df = pd.read_csv('../extendor_proportions.tsv', sep='\t')
print(f"Loaded: {len(prop_df)} extendors")

# Read BLAST results
cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

print("\nLoading blast_refseq_rna_human_v2.txt...")
refseq_df = pd.read_csv('blast_refseq_rna_human_v2.txt', sep='\t', names=cols)
refseq_df['database'] = 'refseq_rna'
print(f"Loaded: {len(refseq_df)} RefSeq RNA hits")

print("\nLoading blast_core_nt_human_v2.txt...")
core_nt_df = pd.read_csv('blast_core_nt_human_v2.txt', sep='\t', names=cols)
core_nt_df['database'] = 'core_nt'
print(f"Loaded: {len(core_nt_df)} Core NT hits")

# Read the FASTA to get mapping
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

# Get FC extendors
tumor_fc_ids = sorted([k for k in fasta_mapping.keys() if k.startswith('tumor_fc_')])
normal_fc_ids = sorted([k for k in fasta_mapping.keys() if k.startswith('normal_fc_')])

print(f"\nTumor FC extendors: {len(tumor_fc_ids)}")
print(f"Normal FC extendors: {len(normal_fc_ids)}")

# Build summary with RefSeq RNA priority, Core NT fallback
print("\nBuilding summary (RefSeq RNA priority, Core NT fallback)...")
results = []

def classify_hit(stitle):
    """Classify BLAST hit into categories"""
    stitle_lower = stitle.lower()
    
    # Ribosomal RNA
    if 'ribosomal' in stitle_lower or 'rrna' in stitle_lower or '45s pre-ribosomal' in stitle_lower:
        if '45s' in stitle_lower or '18s' in stitle_lower or '28s' in stitle_lower or '5.8s' in stitle_lower:
            return 'rRNA (45S/18S/28S/5.8S)'
        elif '5s' in stitle_lower:
            return 'rRNA (5S)'
        return 'rRNA (other)'
    
    # Small nucleolar RNA
    if 'snord' in stitle_lower or 'small nucleolar' in stitle_lower or 'snorna' in stitle_lower:
        return 'snoRNA'
    
    # Small nuclear RNA
    if 'snrna' in stitle_lower or 'small nuclear' in stitle_lower:
        return 'snRNA'
    
    # Transfer RNA
    if 'trna' in stitle_lower or 'transfer rna' in stitle_lower:
        return 'tRNA'
    
    # Satellite/centromeric DNA
    if 'satellite' in stitle_lower or 'alpha satellite' in stitle_lower or 'centromeric' in stitle_lower:
        return 'Satellite/Centromeric DNA'
    
    # Mitochondrial
    if 'mitochondri' in stitle_lower:
        return 'Mitochondrial'
    
    # Chromosome/genomic
    if 'chromosome' in stitle_lower and 'primary assembly' in stitle_lower:
        return 'Genomic DNA'
    
    # Long non-coding RNA
    if 'lncrna' in stitle_lower or 'long non-coding' in stitle_lower or 'linc' in stitle_lower:
        return 'lncRNA'
    
    # MicroRNA
    if 'microrna' in stitle_lower or 'mirna' in stitle_lower or stitle_lower.startswith('mir'):
        return 'miRNA'
    
    # Uncharacterized/predicted
    if 'uncharacterized' in stitle_lower or 'predicted:' in stitle_lower.lower():
        return 'Uncharacterized/Predicted'
    
    # Known protein-coding genes
    if 'mrna' in stitle_lower or 'transcript variant' in stitle_lower:
        return 'mRNA/Protein-coding'
    
    return 'Other'

for qid in tumor_fc_ids + normal_fc_ids:
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
        count_tumor = count_normal = prop_tumor = prop_normal = log2fc = float('nan')
    
    # Try RefSeq RNA first
    refseq_hits = refseq_df[refseq_df['qseqid'] == qid].sort_values('bitscore', ascending=False)
    
    if len(refseq_hits) > 0:
        best_hit = refseq_hits.iloc[0]
        blast_hit = best_hit['stitle']
        blast_pident = best_hit['pident']
        blast_evalue = best_hit['evalue']
        blast_length = best_hit['length']
        blast_source = 'refseq_rna'
    else:
        # Fallback to Core NT
        core_hits = core_nt_df[core_nt_df['qseqid'] == qid].sort_values('bitscore', ascending=False)
        if len(core_hits) > 0:
            best_hit = core_hits.iloc[0]
            blast_hit = best_hit['stitle']
            blast_pident = best_hit['pident']
            blast_evalue = best_hit['evalue']
            blast_length = best_hit['length']
            blast_source = 'core_nt'
        else:
            blast_hit = 'No hit'
            blast_pident = blast_evalue = blast_length = float('nan')
            blast_source = 'none'
    
    # Classify hit
    hit_category = classify_hit(blast_hit) if blast_hit != 'No hit' else 'No hit'
    
    results.append({
        'query_id': qid,
        'anchor': anchor,
        'target': target,
        'count_tumor': count_tumor,
        'count_normal': count_normal,
        'prop_tumor': prop_tumor,
        'prop_normal': prop_normal,
        'log2FC': log2fc,
        'blast_source': blast_source,
        'blast_hit': blast_hit,
        'hit_category': hit_category,
        'blast_pident': blast_pident,
        'blast_evalue': blast_evalue,
        'blast_length': blast_length
    })

# Create DataFrame
result_df = pd.DataFrame(results)
result_df.to_csv('fc_extendors_summary_combined.tsv', sep='\t', index=False)
print(f"\nSaved fc_extendors_summary_combined.tsv with {len(result_df)} extendors")

# Separate tumor and normal
tumor_df = result_df[result_df['query_id'].str.contains('tumor_fc')]
normal_df = result_df[result_df['query_id'].str.contains('normal_fc')]

# ============================================================
# STATISTICS TABLES
# ============================================================
print("\n" + "="*100)
print("SUMMARY STATISTICS TABLES")
print("="*100)

# Table 1: Hit source distribution
print("\n### TABLE 1: BLAST Database Source Distribution ###")
print("-"*60)
print(f"{'Category':<25} {'Tumor FC':<15} {'Normal FC':<15}")
print("-"*60)
tumor_refseq = len(tumor_df[tumor_df['blast_source'] == 'refseq_rna'])
tumor_corент = len(tumor_df[tumor_df['blast_source'] == 'core_nt'])
tumor_none = len(tumor_df[tumor_df['blast_source'] == 'none'])
normal_refseq = len(normal_df[normal_df['blast_source'] == 'refseq_rna'])
normal_corент = len(normal_df[normal_df['blast_source'] == 'core_nt'])
normal_none = len(normal_df[normal_df['blast_source'] == 'none'])

print(f"{'RefSeq RNA':<25} {tumor_refseq:<15} {normal_refseq:<15}")
print(f"{'Core NT (fallback)':<25} {tumor_corент:<15} {normal_corент:<15}")
print(f"{'No hit':<25} {tumor_none:<15} {normal_none:<15}")
print("-"*60)
print(f"{'Total':<25} {len(tumor_df):<15} {len(normal_df):<15}")

# Table 2: Hit category distribution
print("\n### TABLE 2: Hit Category Distribution ###")
print("-"*80)
print(f"{'Category':<35} {'Tumor FC':<12} {'%':<8} {'Normal FC':<12} {'%':<8}")
print("-"*80)

all_categories = sorted(set(result_df['hit_category'].unique()))
for cat in all_categories:
    t_count = len(tumor_df[tumor_df['hit_category'] == cat])
    n_count = len(normal_df[normal_df['hit_category'] == cat])
    t_pct = 100 * t_count / len(tumor_df) if len(tumor_df) > 0 else 0
    n_pct = 100 * n_count / len(normal_df) if len(normal_df) > 0 else 0
    print(f"{cat:<35} {t_count:<12} {t_pct:<8.1f} {n_count:<12} {n_pct:<8.1f}")

print("-"*80)

# Table 3: Detailed hit counts
print("\n### TABLE 3: Top Specific Hits (Tumor FC) ###")
print("-"*90)
tumor_hits = tumor_df[tumor_df['blast_hit'] != 'No hit']['blast_hit'].value_counts().head(15)
for hit, count in tumor_hits.items():
    source = tumor_df[tumor_df['blast_hit'] == hit]['blast_source'].iloc[0]
    print(f"{count:3d} [{source:10}] {hit[:70]}...")

print("\n### TABLE 4: Top Specific Hits (Normal FC) ###")
print("-"*90)
normal_hits = normal_df[normal_df['blast_hit'] != 'No hit']['blast_hit'].value_counts().head(15)
for hit, count in normal_hits.items():
    source = normal_df[normal_df['blast_hit'] == hit]['blast_source'].iloc[0]
    print(f"{count:3d} [{source:10}] {hit[:70]}...")

# Table 5: Log2FC statistics by category
print("\n### TABLE 5: Log2FC Statistics by Hit Category ###")
print("-"*100)
print(f"{'Category':<35} {'Mean log2FC':<15} {'Min':<12} {'Max':<12} {'N':<8}")
print("-"*100)

# Tumor
print("TUMOR FC:")
for cat in tumor_df['hit_category'].unique():
    subset = tumor_df[tumor_df['hit_category'] == cat]['log2FC'].dropna()
    if len(subset) > 0:
        print(f"  {cat:<33} {subset.mean():<15.2f} {subset.min():<12.2f} {subset.max():<12.2f} {len(subset):<8}")

# Normal
print("\nNORMAL FC:")
for cat in normal_df['hit_category'].unique():
    subset = normal_df[normal_df['hit_category'] == cat]['log2FC'].dropna()
    if len(subset) > 0:
        print(f"  {cat:<33} {subset.mean():<15.2f} {subset.min():<12.2f} {subset.max():<12.2f} {len(subset):<8}")

# ============================================================
# Print sample entries
# ============================================================
print("\n" + "="*100)
print("SAMPLE ENTRIES (Top 10 Tumor FC, Top 10 Normal FC)")
print("="*100)

print("\n### TUMOR FC (Top 10 by log2FC) ###")
for _, row in tumor_df.head(10).iterrows():
    print(f"\n{row['query_id']}:")
    print(f"  Counts: Tumor={row['count_tumor']:.0f}, Normal={row['count_normal']:.0f}, log2FC={row['log2FC']:.2f}")
    print(f"  Category: {row['hit_category']}")
    print(f"  Source: {row['blast_source']}")
    print(f"  Hit: {str(row['blast_hit'])[:80]}...")

print("\n### NORMAL FC (Top 10 by -log2FC) ###")
for _, row in normal_df.head(10).iterrows():
    print(f"\n{row['query_id']}:")
    print(f"  Counts: Tumor={row['count_tumor']:.0f}, Normal={row['count_normal']:.0f}, log2FC={row['log2FC']:.2f}")
    print(f"  Category: {row['hit_category']}")
    print(f"  Source: {row['blast_source']}")
    print(f"  Hit: {str(row['blast_hit'])[:80]}...")

# Save summary table as TSV
summary_data = []
for cat in all_categories:
    t_count = len(tumor_df[tumor_df['hit_category'] == cat])
    n_count = len(normal_df[normal_df['hit_category'] == cat])
    summary_data.append({
        'category': cat,
        'tumor_fc_count': t_count,
        'tumor_fc_pct': round(100 * t_count / len(tumor_df), 1),
        'normal_fc_count': n_count,
        'normal_fc_pct': round(100 * n_count / len(normal_df), 1)
    })

summary_table = pd.DataFrame(summary_data)
summary_table.to_csv('fc_extendors_category_summary.tsv', sep='\t', index=False)
print(f"\nSaved fc_extendors_category_summary.tsv")

EOF

echo ""
echo "=== Analysis Complete ==="
echo "Date: $(date)"
