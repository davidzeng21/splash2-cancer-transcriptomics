#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J bam2fq_reprocess
#SBATCH -o %x_%A_%a.out
#SBATCH --array=1-24
#SBATCH -t 30:00
#SBATCH -p shared
#SBATCH -c 20
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL

set -euo pipefail
module load samtools

# ---------- USER EDITABLE ----------
BAM_DIR="/cfs/klemming/projects/supr/ref_free_cancer/jialin/tcga_endometrial/bams"
OUT_DIR="/cfs/klemming/projects/supr/ref_free_cancer/jialin/tcga_endometrial/fastqs"
TMP_BASE="/cfs/klemming/projects/supr/ref_free_cancer/jialin/tcga_endometrial/tmp"
USE_COLLATE=true
COMPRESSION_LEVEL=6
THREADS=${SLURM_CPUS_PER_TASK:-20}
SAMTOOLS_THREADS=$((THREADS > 1 ? THREADS - 1 : 1))
# -----------------------------------

mkdir -p "${OUT_DIR}"
mkdir -p "${TMP_BASE}"
cd "${TMP_BASE}"

# List of failed samples (24 files with blank output)
FAILED_SAMPLES=(
"TCGA-AJ-A2QL-01A"
"TCGA-AJ-A2QL-11A"
"TCGA-AX-A05Y-01A"
"TCGA-AX-A0IZ-01A"
"TCGA-AX-A0J0-01A"
"TCGA-AX-A1CF-01A"
"TCGA-AX-A1CF-11A"
"TCGA-AX-A1CI-01A"
"TCGA-AX-A1CI-11A"
"TCGA-AX-A1CK-01A"
"TCGA-AX-A1CK-11A"
"TCGA-AX-A2H8-01A"
"TCGA-AX-A2H8-11A"
"TCGA-AX-A2HA-01A"
"TCGA-AX-A2HA-11A"
"TCGA-AX-A2HC-01A"
"TCGA-AX-A2HD-01A"
"TCGA-BG-A2AD-01A"
"TCGA-BK-A0CB-01A"
"TCGA-BK-A0CB-11A"
"TCGA-BK-A13C-01A"
"TCGA-BK-A13C-11A"
"TCGA-E6-A1M0-01A"
"TCGA-E6-A1M0-11A"
)

# Get the sample for this array task
sample="${FAILED_SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}"
BAM="${BAM_DIR}/${sample}.bam"

echo "[$(date)] Reprocessing sample: ${sample}"
echo "Input BAM: ${BAM}"

if [ ! -f "${BAM}" ]; then
  echo "Error: BAM file not found: ${BAM}"
  exit 1
fi

# Temporary filenames
COLL_BAM="${TMP_BASE}/${sample}.collated.bam"
FINAL="${OUT_DIR}/${sample}.fastq.gz"
BACKUP="${OUT_DIR}/${sample}.fastq.gz.blank_backup"

# Backup the blank file
if [ -f "${FINAL}" ]; then
  echo "[$(date)] Backing up blank file to ${BACKUP}"
  mv "${FINAL}" "${BACKUP}"
fi

# 1) Optionally collate to group read pairs
if [ "${USE_COLLATE}" = true ]; then
  echo "[$(date)] Running samtools collate (threads=${SAMTOOLS_THREADS})..."
  samtools collate -@ "${SAMTOOLS_THREADS}" -o "${COLL_BAM}" "${BAM}"
  INPUT_FOR_FASTQ="${COLL_BAM}"
else
  echo "[$(date)] Skipping collate; using BAM directly"
  INPUT_FOR_FASTQ="${BAM}"
fi

# 2) Run samtools fastq with MODIFIED command to handle unpaired reads
#    The key change: output unpaired/singleton reads to stdout as well (no -0 /dev/null)
#    Use -s flag to capture singletons/unpaired to stdout along with any paired reads
echo "[$(date)] Running samtools fastq -O (capturing all reads including unpaired)..."
samtools fastq -@ "${SAMTOOLS_THREADS}" -O \
    "${INPUT_FOR_FASTQ}" | pigz -p "${THREADS}" -${COMPRESSION_LEVEL} -c > "${FINAL}"

# Validate output
if [ ! -s "${FINAL}" ]; then
  echo "Error: Output file ${FINAL} is empty or doesn't exist!"
  exit 1
fi

FILE_SIZE=$(stat -c%s "${FINAL}")
if [ "${FILE_SIZE}" -le 100 ]; then
  echo "Error: Output file ${FINAL} is suspiciously small (${FILE_SIZE} bytes)!"
  exit 1
fi

echo "[$(date)] SUCCESS! Output size: $(ls -lh ${FINAL} | awk '{print $5}')"

# 3) Cleanup
echo "[$(date)] Cleaning up temporary files..."
if [ "${USE_COLLATE}" = true ]; then
  rm -f "${COLL_BAM}"
fi

echo "[$(date)] Done for ${sample}. Output: ${FINAL}"

