#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J blast_no_human_hits
#SBATCH -o blast_no_human_hits_%j.out
#SBATCH -t 3:00:00
#SBATCH -p shared
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL

set -euo pipefail
ml PDC
ml blast+

export BLASTDB="/sw/data/blast_databases"

QUERY="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/2026-02-20_min_50_samples/extendor_analysis/blast_tumor_specific/BLAST_results/no_human_hits_blast/ge2_no_human_hits.fasta"
OUT_DIR="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/2026-02-20_min_50_samples/extendor_analysis/blast_tumor_specific/BLAST_results/no_human_hits_blast"
THREADS=16

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sscinames scomnames"

echo "=== BLAST no_human_hits (71 ge2 extendors) ==="
echo "Query:   ${QUERY}"
echo "Start:   $(date)"
echo "Threads: ${THREADS}"
echo ""

# ── Step 1: core_nt (broad nucleotide database) ───────────────────────────────
echo "--- Step 1: blastn vs core_nt ---"
echo "Start: $(date)"
blastn \
  -query        ${QUERY} \
  -db           core_nt \
  -out          ${OUT_DIR}/ge2_no_human_hits_core_nt.txt \
  -outfmt       "${OUTFMT}" \
  -num_threads  ${THREADS} \
  -evalue       1e-5 \
  -word_size    11 \
  -dust         no \
  -perc_identity 80 \
  -qcov_hsp_perc 50 \
  -max_target_seqs 20
echo "Done core_nt: $(date)"
echo "Hits: $(wc -l < ${OUT_DIR}/ge2_no_human_hits_core_nt.txt)"

echo ""

# ── Step 2: refseq_rna ────────────────────────────────────────────────────────
echo "--- Step 2: blastn vs refseq_rna ---"
echo "Start: $(date)"
blastn \
  -query        ${QUERY} \
  -db           refseq_rna \
  -out          ${OUT_DIR}/ge2_no_human_hits_refseq_rna.txt \
  -outfmt       "${OUTFMT}" \
  -num_threads  ${THREADS} \
  -evalue       1e-3 \
  -word_size    11 \
  -dust         no \
  -perc_identity 80 \
  -qcov_hsp_perc 50 \
  -max_target_seqs 20
echo "Done refseq_rna: $(date)"
echo "Hits: $(wc -l < ${OUT_DIR}/ge2_no_human_hits_refseq_rna.txt)"

echo ""

# ── Step 3: Quick summary ─────────────────────────────────────────────────────
echo "--- Step 3: Summary ---"
python3 - <<'PYEOF'
import os, collections

OUT_DIR = "/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched/2026-02-20_min_50_samples/extendor_analysis/blast_tumor_specific/BLAST_results/no_human_hits_blast"

for db, fname in [
    ("core_nt",    "ge2_no_human_hits_core_nt.txt"),
    ("refseq_rna", "ge2_no_human_hits_refseq_rna.txt"),
]:
    path = os.path.join(OUT_DIR, fname)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"\n[{db}] No hits found.")
        continue

    hits = {}       # qseqid → best (bitscore, stitle, sscinames)
    with open(path) as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 15:
                continue
            qseqid, sseqid, pident, length, mm, gap, qs, qe, ss, se, ev, bs, stitle, sscinames, scomnames = cols[:15]
            bs = float(bs)
            if qseqid not in hits or bs > hits[qseqid][0]:
                hits[qseqid] = (bs, pident, stitle[:80], sscinames)

    print(f"\n{'='*72}")
    print(f"[{db}]  Queries with hit: {len(hits)}/71")
    print(f"{'='*72}")
    print(f"{'Extendor (short)':<52} {'pident':>6}  {'Species / Title'}")
    print(f"{'-'*100}")
    for q, (bs, pid, title, sp) in sorted(hits.items(), key=lambda x: -x[1][0]):
        short = q.split('|')[0][:50]
        print(f"{short:<52} {float(pid):>6.1f}  [{sp}] {title[:50]}")
PYEOF

echo ""
echo "=== All done: $(date) ==="
