#!/usr/bin/env python3
"""
Summarise BLAST results for tumor-specific extendors (RNA databases only).

Directory structure expected:
  blast_tumor_specific/
  ├── tumor_extendors_log10FC_ge{2,3}.fasta
  └── BLAST_results/
      └── RNA/
          ├── refseq_rna/         ge{2,3}_refseq_rna_human.txt
          ├── refseq_rna_select/  ge{2,3}_refseq_select_rna_human.txt
          ├── ge{2,3}_rna_best_hits.tsv       (combined, written here)
          ├── ge{2,3}_rna_no_hit.txt          (written here)
          ├── ge{2,3}_gene_summary.tsv        (NEW, written here)
          ├── refseq_rna/ge{2,3}_refseq_rna_best_hits.tsv        (NEW)
          └── refseq_rna_select/ge{2,3}_refseq_select_best_hits.tsv  (NEW)

The human_genome database is intentionally skipped: it provides no gene names
and generates massive hit counts from pericentromeric repetitive elements.

Usage:
  python3 summarize_results.py [--outdir DIR]
  DIR defaults to the blast_tumor_specific directory (parent of BLAST_results).
"""

import argparse
import os
import re
import pandas as pd

# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument(
    '--outdir',
    default=os.path.dirname(os.path.abspath(__file__)),
    help='blast_tumor_specific directory (default: script dir)'
)
args = parser.parse_args()

BASE     = args.outdir
RNA_BASE = os.path.join(BASE, "BLAST_results", "RNA")
RNA_DIR  = os.path.join(RNA_BASE, "refseq_rna")
SEL_DIR  = os.path.join(RNA_BASE, "refseq_rna_select")

# ── Constants ─────────────────────────────────────────────────────────────────
BLAST_COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
              'stitle', 'sscinames', 'scomnames']

OUT_COLS = ['extendor_id', 'log10FC', 'count_tumor', 'count_normal',
            'gene_name', 'pident', 'length', 'evalue', 'bitscore',
            'sseqid', 'stitle', 'database']

# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_header(qseqid: str):
    parts = qseqid.split('|')
    extendor = parts[0]
    meta = {}
    for p in parts[1:]:
        if '=' in p:
            k, v = p.split('=', 1)
            meta[k] = v
    return (extendor,
            float(meta.get('log10FC', 0)),
            int(meta.get('tumor', 0)),
            int(meta.get('normal', 0)))


def load_fasta_queries(fasta_path: str) -> dict:
    queries = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                q = line[1:].strip()
                queries[q] = parse_header(q)
    return queries


def load_blast(path: str) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=BLAST_COLS)
    df = pd.read_csv(path, sep='\t', names=BLAST_COLS)
    for col in ('bitscore', 'pident', 'evalue', 'length'):
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def extract_gene_name(stitle: str) -> str:
    s = str(stitle)
    m = re.search(r'\(([A-Z][A-Z0-9\-\.]+)\)', s)
    if m:
        return m.group(1)
    return s.replace('Homo sapiens ', '')[:40].strip()


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    parsed = [parse_header(q) for q in df['qseqid']]
    df = df.copy()
    df['extendor_id']  = [x[0] for x in parsed]
    df['log10FC']      = [x[1] for x in parsed]
    df['count_tumor']  = [x[2] for x in parsed]
    df['count_normal'] = [x[3] for x in parsed]
    df['gene_name']    = df['stitle'].apply(extract_gene_name)
    return df


def best_per_query(df: pd.DataFrame) -> pd.DataFrame:
    """Best hit per query by bitscore; refseq_rna preferred on ties."""
    return (df.sort_values(['bitscore', 'database'], ascending=[False, True])
              .drop_duplicates('qseqid')
              .reset_index(drop=True))


def gene_summary(df_best: pd.DataFrame) -> pd.DataFrame:
    return (df_best
            .groupby('gene_name')
            .agg(n_extendors=('qseqid', 'count'),
                 median_log10FC=('log10FC', 'median'),
                 max_log10FC=('log10FC', 'max'),
                 median_pident=('pident', 'median'))
            .reset_index()
            .sort_values('n_extendors', ascending=False))


# ── Per-tag processing ────────────────────────────────────────────────────────
for tag, label in [('ge3', 'log10FC ≥ 3'), ('ge2', 'log10FC ≥ 2')]:
    print(f"\n{'='*72}")
    print(f"  {tag.upper()} — tumor-specific extendors ({label}, count_normal ≥ 10)")
    print(f"{'='*72}")

    fasta_path = os.path.join(BASE, f'tumor_extendors_log10FC_{tag}.fasta')
    if not os.path.exists(fasta_path):
        print(f"  FASTA not found: {fasta_path} — skipping.")
        continue
    all_queries = load_fasta_queries(fasta_path)
    print(f"  Total queries: {len(all_queries):,}")

    # ── Load raw BLAST results ────────────────────────────────────────────────
    df_rna = load_blast(os.path.join(RNA_DIR, f'{tag}_refseq_rna_human.txt'))
    df_sel = load_blast(os.path.join(SEL_DIR, f'{tag}_refseq_select_rna_human.txt'))

    df_rna['database'] = 'refseq_rna'
    df_sel['database'] = 'refseq_select_rna'

    print(f"  refseq_rna         raw hits: {len(df_rna):>6,}  "
          f"unique queries: {df_rna['qseqid'].nunique():>5,}")
    print(f"  refseq_select_rna  raw hits: {len(df_sel):>6,}  "
          f"unique queries: {df_sel['qseqid'].nunique():>5,}")

    # ── Per-database annotated best-hit tables ────────────────────────────────
    for db_df, db_dir, db_label in [
        (df_rna, RNA_DIR, 'refseq_rna'),
        (df_sel, SEL_DIR, 'refseq_select_rna'),
    ]:
        if len(db_df) == 0:
            continue
        db_ann  = annotate(db_df)
        db_best = best_per_query(db_ann)
        db_path = os.path.join(db_dir, f'{tag}_{db_label}_best_hits.tsv')
        db_best[OUT_COLS].sort_values('log10FC', ascending=False).to_csv(
            db_path, sep='\t', index=False)
        print(f"  Per-DB best hits ({db_label}): {db_path}")

    # ── Combined best-hit table (one row per query, best across both DBs) ─────
    df_all = pd.concat([df_rna, df_sel], ignore_index=True)

    if len(df_all) == 0:
        print("  No RNA hits found.")
        no_hit_path = os.path.join(RNA_BASE, f'{tag}_rna_no_hit.txt')
        with open(no_hit_path, 'w') as f:
            f.write('extendor\tlog10FC\tcount_tumor\tcount_normal\n')
            for q, (ext, fc, t, n) in sorted(all_queries.items(),
                                              key=lambda x: -x[1][1]):
                f.write(f'{ext}\t{fc:.4f}\t{t}\t{n}\n')
        print(f"  No-hit list: {no_hit_path}")
        continue

    df_ann  = annotate(df_all)
    df_best = best_per_query(df_ann)

    # Queries with no RNA hit
    hit_queries = set(df_best['qseqid'])
    no_hit = {q: v for q, v in all_queries.items() if q not in hit_queries}

    print(f"\n  Queries with RNA hit:    {len(hit_queries):>5,}")
    print(f"  Queries with NO RNA hit: {len(no_hit):>5,}  "
          "(intronic / pericentromeric / non-human candidates)")

    # Save combined best-hit table
    best_path = os.path.join(RNA_BASE, f'{tag}_rna_best_hits.tsv')
    df_best[OUT_COLS].sort_values('log10FC', ascending=False).to_csv(
        best_path, sep='\t', index=False)
    print(f"  Combined best hits:      {best_path}")

    # Save no-hit list
    no_hit_path = os.path.join(RNA_BASE, f'{tag}_rna_no_hit.txt')
    with open(no_hit_path, 'w') as f:
        f.write('extendor\tlog10FC\tcount_tumor\tcount_normal\n')
        for q, (ext, fc, t, n) in sorted(no_hit.items(),
                                          key=lambda x: -x[1][1]):
            f.write(f'{ext}\t{fc:.4f}\t{t}\t{n}\n')
    print(f"  No-RNA-hit list:         {no_hit_path}")

    # Save gene summary
    gs = gene_summary(df_best)
    gene_path = os.path.join(RNA_BASE, f'{tag}_gene_summary.tsv')
    gs.to_csv(gene_path, sep='\t', index=False)
    print(f"  Gene summary:            {gene_path}")

    # ── Print top 20 genes ────────────────────────────────────────────────────
    print(f"\n  Top 20 matched genes (by n_extendors):")
    print(f"  {'Gene':<35} {'N':>5}  {'median_FC':>10}  {'max_FC':>7}  "
          f"{'median_pid':>10}")
    print(f"  {'-'*72}")
    for _, row in gs.head(20).iterrows():
        print(f"  {row['gene_name']:<35} {int(row['n_extendors']):>5}  "
              f"{row['median_log10FC']:>10.2f}  {row['max_log10FC']:>7.2f}  "
              f"{row['median_pident']:>10.1f}")

    # ── Print top 20 extendors with RNA hit (by log10FC) ─────────────────────
    print(f"\n  Top 20 extendors with RNA hit (by log10FC):")
    print(f"  {'Extendor':<52} {'FC':>6} {'tumor':>7} {'norm':>5} "
          f"{'pid':>5}  Gene")
    print(f"  {'-'*100}")
    for _, row in df_best.sort_values('log10FC', ascending=False).head(20).iterrows():
        print(f"  {row['extendor_id']:<52} {row['log10FC']:>6.2f} "
              f"{int(row['count_tumor']):>7} {int(row['count_normal']):>5} "
              f"{row['pident']:>5.1f}  {row['gene_name']}")

    # ── Print top 20 no-hit extendors (by log10FC) ───────────────────────────
    if no_hit:
        print(f"\n  Top 20 extendors with NO RNA hit (by log10FC):")
        print(f"  {'Extendor':<52} {'FC':>6} {'tumor':>7} {'norm':>5}")
        print(f"  {'-'*75}")
        for q, (ext, fc, t, n) in sorted(
                no_hit.items(), key=lambda x: -x[1][1])[:20]:
            print(f"  {ext:<52} {fc:>6.2f} {t:>7} {n:>5}")

print(f"\n{'='*72}")
print(f"All output files written under: {RNA_BASE}")
print("Done.")
