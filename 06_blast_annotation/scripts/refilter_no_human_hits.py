#!/usr/bin/env python3
"""
Re-filter no_human_hits for ge3 and ge2 using the new directory structure.

Logic (same as summarize_ge2.py):
  A query has NO human hit if it is absent from ALL three databases:
    1. human_genome  (GRCh38.p13) -- streamed line-by-line for large files
    2. refseq_rna    (human)       -- pandas
    3. refseq_select_rna (human)   -- pandas

Outputs (to BLAST_results/human_database_GRCh38.p13/):
  ge3_no_human_hits.tsv   (header + TSV, sorted by log10FC desc)
  ge2_no_human_hits.tsv   (header + TSV, sorted by log10FC desc)

Note: ge3_human_genome.txt is ~95MB and ge2_human_genome.txt is ~18GB;
both are streamed to keep memory usage low.
"""

import os
import sys
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
GENOME_DIR = os.path.join(BASE, "BLAST_results", "human_database_GRCh38.p13")
RNA_DIR    = os.path.join(BASE, "BLAST_results", "RNA", "refseq_rna")
SEL_DIR    = os.path.join(BASE, "BLAST_results", "RNA", "refseq_rna_select")

BLAST_COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
              'stitle', 'sscinames', 'scomnames']


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


def stream_genome_hits(genome_path: str) -> set:
    """Stream genome BLAST file, return set of qseqid strings with any hit."""
    print(f"  Streaming {genome_path} "
          f"({os.path.getsize(genome_path)/1e9:.2f} GB)...", flush=True)
    hit_set = set()
    n = 0
    with open(genome_path) as f:
        for line in f:
            n += 1
            if n % 10_000_000 == 0:
                print(f"    ...{n/1e6:.0f}M lines, {len(hit_set)} unique queries",
                      flush=True)
            parts = line.split('\t', 1)
            if parts:
                hit_set.add(parts[0])
    print(f"  Done: {n:,} lines, {len(hit_set):,} unique queries with genome hits",
          flush=True)
    return hit_set


def load_rna_hits(path: str) -> set:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return set()
    df = pd.read_csv(path, sep='\t', names=BLAST_COLS, usecols=['qseqid'])
    return set(df['qseqid'])


def process_tag(tag: str):
    print(f"\n{'='*60}")
    print(f"Processing {tag.upper()} ...")
    print(f"{'='*60}", flush=True)

    fasta_path   = os.path.join(BASE, f"tumor_extendors_log10FC_{tag}.fasta")
    genome_path  = os.path.join(GENOME_DIR, f"{tag}_human_genome.txt")
    rna_path     = os.path.join(RNA_DIR,    f"{tag}_refseq_rna_human.txt")
    sel_path     = os.path.join(SEL_DIR,    f"{tag}_refseq_select_rna_human.txt")
    out_path     = os.path.join(GENOME_DIR, f"{tag}_no_human_hits.tsv")

    # Load all query IDs
    all_queries = load_fasta_queries(fasta_path)
    print(f"Total queries: {len(all_queries):,}", flush=True)

    # Collect hits from each database
    genome_hits = stream_genome_hits(genome_path)

    print(f"  Loading RNA hits...", flush=True)
    rna_hits  = load_rna_hits(rna_path)
    sel_hits  = load_rna_hits(sel_path)
    print(f"  refseq_rna hits:        {len(rna_hits):,}", flush=True)
    print(f"  refseq_select_rna hits: {len(sel_hits):,}", flush=True)

    # Union of all hits across all 3 databases
    all_hits = genome_hits | rna_hits | sel_hits
    print(f"  Queries with ANY human hit: {len(all_hits):,}", flush=True)

    # Subtract from full query set
    no_hit_queries = {q: v for q, v in all_queries.items()
                      if q not in all_hits}
    print(f"  Queries with NO human hit:  {len(no_hit_queries):,}", flush=True)

    # Write output: clean TSV with header, sorted by log10FC descending
    with open(out_path, 'w') as f:
        f.write('extendor\tlog10FC\tcount_tumor\tcount_normal\n')
        for q, (extendor, fc, tumor, normal) in sorted(
                no_hit_queries.items(), key=lambda x: -x[1][1]):
            f.write(f'{extendor}\t{fc:.4f}\t{tumor}\t{normal}\n')

    print(f"  Saved: {out_path}", flush=True)
    print(f"  Entries: {len(no_hit_queries):,}", flush=True)

    # Quick check: show top 10 by log10FC
    print(f"\n  Top 10 no-human-hit extendors (log10FC desc):")
    print(f"  {'Extendor':<52} {'log10FC':>8} {'tumor':>7} {'normal':>7}")
    print(f"  {'-'*78}")
    for q, (ext, fc, t, n) in sorted(
            no_hit_queries.items(), key=lambda x: -x[1][1])[:10]:
        print(f"  {ext:<52} {fc:>8.4f} {t:>7} {n:>7}")


if __name__ == '__main__':
    for tag in ['ge3', 'ge2']:
        process_tag(tag)

    print(f"\n{'='*60}")
    print("Done.")
