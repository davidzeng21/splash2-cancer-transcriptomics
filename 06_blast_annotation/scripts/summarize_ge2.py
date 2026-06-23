"""
Memory-efficient summariser for ge2 BLAST results.
Processes the 18GB genome file line-by-line (no pandas for large files).
Uses pandas only for the small RNA files.
"""

import os, sys
from collections import defaultdict
import pandas as pd

ODIR = os.path.dirname(os.path.abspath(__file__))

COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
        'stitle', 'sscinames', 'scomnames']

def parse_header(qseqid):
    parts = qseqid.split('|')
    extendor = parts[0]
    meta = {k: v for k, v in (p.split('=') for p in parts[1:])}
    return extendor, float(meta.get('log10FC', 0)), int(meta.get('tumor', 0)), int(meta.get('normal', 0))

# ── Load all query IDs from FASTA ─────────────────────────────────────────────
fasta_path = os.path.join(ODIR, 'tumor_extendors_log10FC_ge2.fasta')
all_queries = {}   # qseqid -> (extendor_id, log10FC, tumor, normal)
with open(fasta_path) as f:
    for line in f:
        if line.startswith('>'):
            q = line[1:].strip()
            all_queries[q] = parse_header(q)
print(f"Total ge2 queries: {len(all_queries)}")

# ── Stream the 18GB genome file: keep best hit per query ──────────────────────
genome_path = os.path.join(ODIR, 'ge2_human_genome.txt')
print(f"\nStreaming {genome_path} ({os.path.getsize(genome_path)/1e9:.1f} GB)...")

best_genome = {}   # qseqid -> list of col values for the row with max bitscore
n_lines = 0
with open(genome_path) as f:
    for line in f:
        n_lines += 1
        if n_lines % 5_000_000 == 0:
            print(f"  ...processed {n_lines/1e6:.0f}M lines, "
                  f"{len(best_genome)} unique queries seen so far")
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 12:
            continue
        qseqid  = parts[0]
        bitscore = float(parts[11])
        if qseqid not in best_genome or bitscore > float(best_genome[qseqid][11]):
            best_genome[qseqid] = parts

print(f"  Done. Total lines: {n_lines:,}, unique queries with genome hits: {len(best_genome)}")

# Write best genome hits
genome_best_path = os.path.join(ODIR, 'ge2_human_genome_best_hits.tsv')
with open(genome_best_path, 'w') as f:
    f.write('\t'.join(COLS) + '\textendor_id\tlog10FC\tcount_tumor\tcount_normal\n')
    for qseqid, parts in sorted(best_genome.items()):
        extendor, fc, tumor, normal = parse_header(qseqid)
        row = parts[:len(COLS)]
        while len(row) < len(COLS):
            row.append('N/A')
        f.write('\t'.join(row) + f'\t{extendor}\t{fc:.4f}\t{tumor}\t{normal}\n')
print(f"Best genome hits saved: {genome_best_path}")

# ── Load RNA files with pandas (small, safe) ─────────────────────────────────
def load_blast_tsv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=COLS)
    return pd.read_csv(path, sep='\t', names=COLS)

df_select = load_blast_tsv(os.path.join(ODIR, 'ge2_refseq_select_rna_human.txt'))
df_rna    = load_blast_tsv(os.path.join(ODIR, 'ge2_refseq_rna_human.txt'))

print(f"\nRefSeq select RNA hits: {len(df_select)}, unique queries: {df_select['qseqid'].nunique() if len(df_select) else 0}")
print(f"RefSeq RNA hits:        {len(df_rna)}, unique queries: {df_rna['qseqid'].nunique() if len(df_rna) else 0}")

# Best hit per query for RNA databases
def best_per_query(df, db_label):
    if len(df) == 0:
        return pd.DataFrame()
    df = df.copy()
    df['database'] = db_label
    df[['extendor_id', 'log10FC', 'count_tumor', 'count_normal']] = \
        pd.DataFrame([parse_header(q) for q in df['qseqid']], index=df.index)
    return df.sort_values('bitscore', ascending=False).drop_duplicates('qseqid')

rna_select_best = best_per_query(df_select, 'refseq_select_rna')
rna_full_best   = best_per_query(df_rna,    'refseq_rna')

rna_best_path = os.path.join(ODIR, 'ge2_rna_best_hits.tsv')
pd.concat([rna_select_best, rna_full_best], ignore_index=True).to_csv(rna_best_path, sep='\t', index=False)
print(f"RNA best hits saved: {rna_best_path}")

# ── Identify queries with NO human hit in ANY database ────────────────────────
hit_in_genome = set(best_genome.keys())
hit_in_rna    = set(df_select['qseqid'].tolist() + df_rna['qseqid'].tolist())
hit_anywhere  = hit_in_genome | hit_in_rna
no_hit_queries = {q: v for q, v in all_queries.items() if q not in hit_anywhere}

print(f"\n{'='*70}")
print(f"SUMMARY - ge2 extendors (log10FC ≥ 2, count_normal ≥ 10)")
print(f"{'='*70}")
print(f"Total queries:                {len(all_queries)}")
print(f"Queries with genome hits:     {len(hit_in_genome)}")
print(f"Queries with RNA hits:        {len(hit_in_rna)}")
print(f"Queries with ANY human hit:   {len(hit_anywhere)}")
print(f"Queries with NO human hit:    {len(no_hit_queries)}  <-- potential non-human")

no_hit_path = os.path.join(ODIR, 'ge2_no_human_hits.txt')
with open(no_hit_path, 'w') as f:
    f.write('extendor\tlog10FC\tcount_tumor\tcount_normal\n')
    for q, (extendor, fc, tumor, normal) in sorted(
            no_hit_queries.items(), key=lambda x: -x[1][1]):
        f.write(f"{extendor}\t{fc:.4f}\t{tumor}\t{normal}\n")
print(f"No-hit list saved: {no_hit_path}")

# ── Per-database breakdown ────────────────────────────────────────────────────
print(f"\n{'Database':<30} {'Queries_with_hits':>18} {'Total_raw_hits':>15}")
print(f"{'-'*65}")
print(f"{'human_genome':<30} {len(hit_in_genome):>18,} {n_lines:>15,}")
print(f"{'refseq_select_rna':<30} {df_select['qseqid'].nunique() if len(df_select) else 0:>18,} {len(df_select):>15,}")
print(f"{'refseq_rna':<30} {df_rna['qseqid'].nunique() if len(df_rna) else 0:>18,} {len(df_rna):>15,}")

# ── Note on genome hit count ─────────────────────────────────────────────────
print(f"""
NOTE: The genome produced {n_lines:,} raw hits from {len(all_queries)} queries
({n_lines/len(all_queries):.0f} hits/query on average). This is because many 
extendors match repetitive elements (Alu, LINE, SINE) that occur thousands of 
times in the genome. The -max_target_seqs 10 flag limits unique SUBJECTS 
(chromosomes), not HSPs per subject, so one chromosome can contribute many rows.

For downstream analysis, use ge2_human_genome_best_hits.tsv (one row per query,
highest bitscore) instead of the raw 18GB file.
""")
print("Done.")
