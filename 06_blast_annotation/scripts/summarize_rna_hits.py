"""
Summarise BLAST results against RefSeq RNA databases only (ge3 and ge2).
Skips the genome hits (no gene names there).
Outputs:
  - {tag}_rna_best_hits.tsv   : best hit per query (highest bitscore), with gene name
  - {tag}_rna_no_hit.txt      : queries with no RNA hit (intronic / non-human candidates)
"""

import os
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

def load_fasta_queries(fasta_path):
    queries = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                q = line[1:].strip()
                queries[q] = parse_header(q)
    return queries

def load_blast(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=COLS)
    df = pd.read_csv(path, sep='\t', names=COLS)
    df['bitscore'] = pd.to_numeric(df['bitscore'], errors='coerce')
    return df

def extract_gene_name(stitle):
    """Pull a short gene label from the RefSeq stitle string."""
    s = str(stitle)
    # e.g. "Homo sapiens MMP1 (MMP1), mRNA" -> "MMP1"
    import re
    m = re.search(r'\(([A-Z][A-Z0-9\-]+)\)', s)
    if m:
        return m.group(1)
    # fallback: first token after "Homo sapiens "
    s2 = s.replace('Homo sapiens ', '')
    return s2[:40]

for tag, label in [('ge3', 'log10FC ≥ 3'), ('ge2', 'log10FC ≥ 2')]:
    print(f"\n{'='*72}")
    print(f"  {tag.upper()} — tumor-specific extendors ({label}, count_normal ≥ 10)")
    print(f"{'='*72}")

    # Load all queries from FASTA
    fasta = os.path.join(ODIR, f'tumor_extendors_log10FC_{tag}.fasta')
    all_queries = load_fasta_queries(fasta)
    print(f"Total queries: {len(all_queries)}")

    # Load RNA results
    df_sel = load_blast(os.path.join(ODIR, f'{tag}_refseq_select_rna_human.txt'))
    df_rna = load_blast(os.path.join(ODIR, f'{tag}_refseq_rna_human.txt'))

    df_sel['database'] = 'refseq_select_rna'
    df_rna['database'] = 'refseq_rna'

    print(f"refseq_select_rna  raw hits: {len(df_sel):>6,}   "
          f"unique queries: {df_sel['qseqid'].nunique():>5,}")
    print(f"refseq_rna         raw hits: {len(df_rna):>6,}   "
          f"unique queries: {df_rna['qseqid'].nunique():>5,}")

    # Combine and take best hit per query (across both RNA dbs)
    df_all = pd.concat([df_sel, df_rna], ignore_index=True)

    if len(df_all) == 0:
        print("  No RNA hits found.")
        continue

    # Parse metadata
    parsed = [parse_header(q) for q in df_all['qseqid']]
    df_all['extendor_id']  = [x[0] for x in parsed]
    df_all['log10FC']      = [x[1] for x in parsed]
    df_all['count_tumor']  = [x[2] for x in parsed]
    df_all['count_normal'] = [x[3] for x in parsed]
    df_all['gene_name']    = df_all['stitle'].apply(extract_gene_name)

    # Best hit per query (highest bitscore, prefer refseq_rna for broader coverage)
    df_best = (df_all.sort_values(['bitscore', 'database'], ascending=[False, True])
                     .drop_duplicates('qseqid')
                     .reset_index(drop=True))

    # Queries with no RNA hit at all
    hit_queries  = set(df_best['qseqid'])
    no_hit       = {q: v for q, v in all_queries.items() if q not in hit_queries}

    print(f"\nQueries with RNA hit:    {len(hit_queries):>5,}")
    print(f"Queries with NO RNA hit: {len(no_hit):>5,}  "
          f"(intronic / non-human candidates)")

    # ── Save best-hit table ───────────────────────────────────────────────────
    out_cols = ['extendor_id', 'log10FC', 'count_tumor', 'count_normal',
                'gene_name', 'pident', 'length', 'evalue', 'bitscore',
                'sseqid', 'stitle', 'sscinames', 'scomnames', 'database']
    best_path = os.path.join(ODIR, f'{tag}_rna_best_hits.tsv')
    df_best[out_cols].sort_values('log10FC', ascending=False).to_csv(
        best_path, sep='\t', index=False)
    print(f"Best RNA hits saved:     {best_path}")

    # ── Save no-hit list ──────────────────────────────────────────────────────
    no_hit_path = os.path.join(ODIR, f'{tag}_rna_no_hit.txt')
    with open(no_hit_path, 'w') as f:
        f.write('extendor\tlog10FC\tcount_tumor\tcount_normal\n')
        for q, (extendor, fc, tumor, normal) in sorted(
                no_hit.items(), key=lambda x: -x[1][1]):
            f.write(f'{extendor}\t{fc:.4f}\t{tumor}\t{normal}\n')
    print(f"No-RNA-hit list saved:   {no_hit_path}")

    # ── Gene hit summary ─────────────────────────────────────────────────────
    gene_counts = (df_best['gene_name']
                   .value_counts()
                   .reset_index()
                   .rename(columns={'index': 'gene_name', 'gene_name': 'n_extendors'}))

    print(f"\n  Top 20 matched genes (by number of extendors):")
    print(f"  {'Gene':<35} {'N_extendors':>12} {'Median_log10FC':>15}")
    print(f"  {'-'*65}")
    for _, row in gene_counts.head(20).iterrows():
        g = row.iloc[0]
        n = row.iloc[1]
        med_fc = df_best[df_best['gene_name'] == g]['log10FC'].median()
        print(f"  {g:<35} {n:>12} {med_fc:>15.2f}")

    # ── Top extendors with RNA hit (by log10FC) ───────────────────────────────
    print(f"\n  Top 20 extendors with RNA hit (by log10FC):")
    print(f"  {'Extendor':<52} {'log10FC':>8} {'tumor':>7} {'normal':>7} "
          f"{'pident':>7}  Gene")
    print(f"  {'-'*110}")
    for _, row in df_best.sort_values('log10FC', ascending=False).head(20).iterrows():
        print(f"  {row['extendor_id']:<52} {row['log10FC']:>8.2f} "
              f"{int(row['count_tumor']):>7} {int(row['count_normal']):>7} "
              f"{row['pident']:>7.1f}  {row['gene_name']}")

    # ── Top no-hit extendors (by log10FC) ─────────────────────────────────────
    if no_hit:
        print(f"\n  Top 20 extendors with NO RNA hit (by log10FC):")
        print(f"  {'Extendor':<52} {'log10FC':>8} {'tumor':>7} {'normal':>7}")
        print(f"  {'-'*80}")
        sorted_no_hit = sorted(no_hit.items(), key=lambda x: -x[1][1])
        for q, (extendor, fc, tumor, normal) in sorted_no_hit[:20]:
            print(f"  {extendor:<52} {fc:>8.2f} {tumor:>7} {normal:>7}")

print(f"\n{'='*72}")
print("Done.")
