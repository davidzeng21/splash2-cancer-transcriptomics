#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J blast_tumor_extendors
#SBATCH -o blast_tumor_extendors_%j.out
#SBATCH -t 6:00:00
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH --begin=now
#SBATCH --mail-type=ALL

# Load BLAST+
ml PDC blast+
source /cfs/klemming/home/j/jlzeng/myenv/bin/activate
# Set BLASTDB so taxonomy lookups (sscinames, scomnames) work
export BLASTDB=/sw/data/blast_databases

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched"
EXTENDOR_DIR="${BASE_DIR}/2026-02-20_min_50_samples/extendor_analysis"
OUT_DIR="${EXTENDOR_DIR}/blast_tumor_specific"

TSV_GE3="${EXTENDOR_DIR}/extendor_proportions_log10FC_ge3_count_normal_ge10.tsv"
TSV_GE2="${EXTENDOR_DIR}/extendor_proportions_log10FC_ge2_count_normal_ge10.tsv"

FASTA_GE3="${OUT_DIR}/tumor_extendors_log10FC_ge3.fasta"
FASTA_GE2="${OUT_DIR}/tumor_extendors_log10FC_ge2.fasta"

# ── Databases ────────────────────────────────────────────────────────────────
# human_genome: dedicated GRCh38 database, no taxid filter needed
HUMAN_GENOME="/sw/data/blast_databases/human_genome"
# refseq_select_rna: one representative transcript per gene, filtered to human
REFSEQ_SELECT_RNA="/sw/data/blast_databases/refseq_select_rna"
# refseq_rna: full RefSeq RNA, filtered to human (broader coverage)
REFSEQ_RNA="/sw/data/blast_databases/refseq_rna"
HUMAN_TAXID=9606

THREADS=16

# ── BLAST parameter rationale ─────────────────────────────────────────────────
# Queries are 50bp k-mers (extendors = anchor25 + target25, underscore removed).
#
# -dust no       : DUST masks low-complexity (AT-rich) regions; for short k-mers
#                  this can suppress real hits. Disabled to maximise sensitivity.
# -word_size 11  : blastn default; safe for 50bp (smaller = slower but no gain here).
# -qcov_hsp_perc : Require ≥80% of the 50bp query (≥40bp) to align; prevents
#                  spurious partial-end hits.
# -perc_identity : 90% allows ~5 mismatches over 50bp (SNP-level variation).
# -evalue        : Database-size dependent:
#                    genome  (~3 Gbp) → 1e-5  (large search space, strict)
#                    RNA dbs (~few hundred Mbp) → 1e-3  (smaller, more lenient)
EVALUE_GENOME=1e-5
EVALUE_RNA=1e-3

# ── Output format ─────────────────────────────────────────────────────────────
# Standard tabular fields + subject title + scientific name + common name
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle sscinames scomnames"

# ── Setup ────────────────────────────────────────────────────────────────────
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "=== BLAST Analysis: Tumor-Specific Extendors vs Human Genome & Transcriptome ==="
echo "Date: $(date)"
echo "Output dir: ${OUT_DIR}"
echo ""

# Verify taxonomy database files (required for sscinames/scomnames)
echo "Checking taxonomy database files..."
for f in "${BLASTDB}/taxonomy4blast.sqlite3" "${BLASTDB}/taxdb.btd" "${BLASTDB}/taxdb.bti"; do
    if [ -f "$f" ]; then
        echo "  OK: $f"
    else
        echo "  MISSING: $f  <-- sscinames/scomnames may not work!"
    fi
done
echo ""

# ── Step 1: Generate FASTA files from TSV ────────────────────────────────────
echo "=== Step 1: Generating FASTA query files ==="
echo "Start: $(date)"

python3 << 'PYEOF'
import os, sys

def tsv_to_fasta(tsv_path, fasta_path):
    """
    Convert extendor TSV to FASTA.
    Extendor column format: ANCHOR_TARGET (25bp_25bp).
    The underscore is removed to form the 50bp query sequence.
    FASTA header: >extendor|log10FC=X.XXXX|tumor=XXXXX|normal=XX
    """
    written = 0
    with open(tsv_path) as fin, open(fasta_path, 'w') as fout:
        header = fin.readline().strip().split('\t')
        # Expected columns: extendor anchor target count_normal count_tumor
        #                   prop_tumor prop_normal log10FC
        idx = {col: i for i, col in enumerate(header)}
        for line in fin:
            parts = line.strip().split('\t')
            if not parts or len(parts) < len(header):
                continue
            extendor_id = parts[idx['extendor']]       # e.g. ACGT..._TGCA...
            log10fc     = float(parts[idx['log10FC']])
            count_tumor = int(parts[idx['count_tumor']])
            count_normal= int(parts[idx['count_normal']])
            # Remove underscore to get clean 50bp sequence
            sequence    = extendor_id.replace('_', '')
            fasta_name  = (f"{extendor_id}|"
                           f"log10FC={log10fc:.4f}|"
                           f"tumor={count_tumor}|"
                           f"normal={count_normal}")
            fout.write(f">{fasta_name}\n{sequence}\n")
            written += 1
    return written

base = "/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/CPTAC_Lung_Adeno_matched"
edir = f"{base}/2026-02-20_min_50_samples/extendor_analysis"
odir = f"{edir}/blast_tumor_specific"

tsv_ge3   = f"{edir}/extendor_proportions_log10FC_ge3_count_normal_ge10.tsv"
tsv_ge2   = f"{edir}/extendor_proportions_log10FC_ge2_count_normal_ge10.tsv"
fasta_ge3 = f"{odir}/tumor_extendors_log10FC_ge3.fasta"
fasta_ge2 = f"{odir}/tumor_extendors_log10FC_ge2.fasta"

n3 = tsv_to_fasta(tsv_ge3, fasta_ge3)
print(f"ge3 FASTA: {n3} sequences -> {fasta_ge3}")

n2 = tsv_to_fasta(tsv_ge2, fasta_ge2)
print(f"ge2 FASTA: {n2} sequences -> {fasta_ge2}")
PYEOF

echo "FASTA generation done: $(date)"
echo ""

# Quick sanity check
echo "ge3 sequences: $(grep -c '^>' ${FASTA_GE3})"
echo "ge2 sequences: $(grep -c '^>' ${FASTA_GE2})"
echo "ge3 first 3 headers:"
grep '^>' ${FASTA_GE3} | head -3
echo ""

# ── Step 2: BLAST ge3 (log10FC ≥ 3) extendors ───────────────────────────────
echo "=== Step 2: BLAST log10FC ≥ 3 extendors (${FASTA_GE3}) ==="

# 2a. vs human genome
echo "--- [ge3] blastn vs human_genome ---"
echo "Start: $(date)"
blastn -query ${FASTA_GE3} \
    -db ${HUMAN_GENOME} \
    -task blastn \
    -word_size 11 \
    -evalue ${EVALUE_GENOME} \
    -perc_identity 90 \
    -qcov_hsp_perc 80 \
    -dust no \
    -strand both \
    -max_target_seqs 10 \
    -num_threads ${THREADS} \
    -outfmt "${OUTFMT}" \
    -out ge3_human_genome.txt
echo "Hits: $(wc -l < ge3_human_genome.txt)"
echo "End: $(date)"
echo ""

# 2b. vs RefSeq select RNA (human)
echo "--- [ge3] blastn vs refseq_select_rna (human) ---"
echo "Start: $(date)"
blastn -query ${FASTA_GE3} \
    -db ${REFSEQ_SELECT_RNA} \
    -taxids ${HUMAN_TAXID} \
    -task blastn \
    -word_size 11 \
    -evalue ${EVALUE_RNA} \
    -perc_identity 90 \
    -qcov_hsp_perc 80 \
    -dust no \
    -strand both \
    -max_target_seqs 10 \
    -num_threads ${THREADS} \
    -outfmt "${OUTFMT}" \
    -out ge3_refseq_select_rna_human.txt
echo "Hits: $(wc -l < ge3_refseq_select_rna_human.txt)"
echo "End: $(date)"
echo ""

# 2c. vs RefSeq RNA full (human) — broader coverage for noncoding RNA
echo "--- [ge3] blastn vs refseq_rna (human) ---"
echo "Start: $(date)"
blastn -query ${FASTA_GE3} \
    -db ${REFSEQ_RNA} \
    -taxids ${HUMAN_TAXID} \
    -task blastn \
    -word_size 11 \
    -evalue ${EVALUE_RNA} \
    -perc_identity 90 \
    -qcov_hsp_perc 80 \
    -dust no \
    -strand both \
    -max_target_seqs 10 \
    -num_threads ${THREADS} \
    -outfmt "${OUTFMT}" \
    -out ge3_refseq_rna_human.txt
echo "Hits: $(wc -l < ge3_refseq_rna_human.txt)"
echo "End: $(date)"
echo ""

# ── Step 3: BLAST ge2 (log10FC ≥ 2) extendors ───────────────────────────────
echo "=== Step 3: BLAST log10FC ≥ 2 extendors (${FASTA_GE2}) ==="

# 3a. vs human genome
echo "--- [ge2] blastn vs human_genome ---"
echo "Start: $(date)"
blastn -query ${FASTA_GE2} \
    -db ${HUMAN_GENOME} \
    -task blastn \
    -word_size 11 \
    -evalue ${EVALUE_GENOME} \
    -perc_identity 90 \
    -qcov_hsp_perc 80 \
    -dust no \
    -strand both \
    -max_target_seqs 10 \
    -num_threads ${THREADS} \
    -outfmt "${OUTFMT}" \
    -out ge2_human_genome.txt
echo "Hits: $(wc -l < ge2_human_genome.txt)"
echo "End: $(date)"
echo ""

# 3b. vs RefSeq select RNA (human)
echo "--- [ge2] blastn vs refseq_select_rna (human) ---"
echo "Start: $(date)"
blastn -query ${FASTA_GE2} \
    -db ${REFSEQ_SELECT_RNA} \
    -taxids ${HUMAN_TAXID} \
    -task blastn \
    -word_size 11 \
    -evalue ${EVALUE_RNA} \
    -perc_identity 90 \
    -qcov_hsp_perc 80 \
    -dust no \
    -strand both \
    -max_target_seqs 10 \
    -num_threads ${THREADS} \
    -outfmt "${OUTFMT}" \
    -out ge2_refseq_select_rna_human.txt
echo "Hits: $(wc -l < ge2_refseq_select_rna_human.txt)"
echo "End: $(date)"
echo ""

# 3c. vs RefSeq RNA full (human)
echo "--- [ge2] blastn vs refseq_rna (human) ---"
echo "Start: $(date)"
blastn -query ${FASTA_GE2} \
    -db ${REFSEQ_RNA} \
    -taxids ${HUMAN_TAXID} \
    -task blastn \
    -word_size 11 \
    -evalue ${EVALUE_RNA} \
    -perc_identity 90 \
    -qcov_hsp_perc 80 \
    -dust no \
    -strand both \
    -max_target_seqs 10 \
    -num_threads ${THREADS} \
    -outfmt "${OUTFMT}" \
    -out ge2_refseq_rna_human.txt
echo "Hits: $(wc -l < ge2_refseq_rna_human.txt)"
echo "End: $(date)"
echo ""

# ── Step 4: Summarize results ─────────────────────────────────────────────────
echo "=== Step 4: Summarizing results (RNA databases only) ==="
echo "Start: $(date)"

SUMMARIZE_SCRIPT="${BASE_DIR}/2026-02-20_min_50_samples/extendor_analysis/blast_tumor_specific/summarize_results.py"
python3 ${SUMMARIZE_SCRIPT} --outdir ${OUT_DIR}

echo ""
echo "=== All done ==="
echo "Date: $(date)"
