#!/usr/bin/env python3
"""Create FASTA files for BLAST analysis of top differential extendors."""

import pandas as pd
import os

# Read the proportions file
df = pd.read_csv('../extendor_proportions.tsv', sep='\t')

# Top tumor-enriched (zero in normal, high in tumor)
tumor_enriched = df[(df['count_normal'] == 0) & (df['count_tumor'] >= 50)].nlargest(20, 'count_tumor')

# Top normal-enriched (zero in tumor, high in normal)  
normal_enriched = df[(df['count_tumor'] == 0) & (df['count_normal'] >= 50)].nlargest(20, 'count_normal')

# Write tumor-enriched FASTA
with open('tumor_enriched.fasta', 'w') as f:
    for i, row in enumerate(tumor_enriched.itertuples(), 1):
        # Write anchor
        f.write(f">tumor_{i:02d}_counts{int(row.count_tumor)}_anchor\n")
        f.write(f"{row.anchor}\n")
        # Write target  
        f.write(f">tumor_{i:02d}_counts{int(row.count_tumor)}_target\n")
        f.write(f"{row.target}\n")

# Write normal-enriched FASTA
with open('normal_enriched.fasta', 'w') as f:
    for i, row in enumerate(normal_enriched.itertuples(), 1):
        f.write(f">normal_{i:02d}_counts{int(row.count_normal)}_anchor\n")
        f.write(f"{row.anchor}\n")
        f.write(f">normal_{i:02d}_counts{int(row.count_normal)}_target\n")
        f.write(f"{row.target}\n")

# Write combined file with full extendor (anchor+target concatenated)
with open('top_extendors_combined.fasta', 'w') as f:
    for i, row in enumerate(tumor_enriched.itertuples(), 1):
        f.write(f">tumor_{i:02d}_counts{int(row.count_tumor)}_extendor\n")
        f.write(f"{row.anchor}{row.target}\n")
    for i, row in enumerate(normal_enriched.itertuples(), 1):
        f.write(f">normal_{i:02d}_counts{int(row.count_normal)}_extendor\n")
        f.write(f"{row.anchor}{row.target}\n")

print("Created FASTA files:")
print(f"  - tumor_enriched.fasta ({len(tumor_enriched)*2} sequences)")
print(f"  - normal_enriched.fasta ({len(normal_enriched)*2} sequences)")
print(f"  - top_extendors_combined.fasta ({len(tumor_enriched)+len(normal_enriched)} full extendors)")

print("\n=== Top 5 Tumor-Enriched Extendors ===")
for i, row in enumerate(tumor_enriched.head().itertuples(), 1):
    print(f"{i}. {row.anchor} | counts={int(row.count_tumor)}")
    
print("\n=== Top 5 Normal-Enriched Extendors ===")
for i, row in enumerate(normal_enriched.head().itertuples(), 1):
    print(f"{i}. {row.anchor} | counts={int(row.count_normal)}")



