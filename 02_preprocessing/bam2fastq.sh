#!/bin/bash -l
#SBATCH -A naiss2025-22-738
#SBATCH -J bam2fq
#SBATCH -o bam2fq_%j.out
#SBATCH --array=1-6
#SBATCH -t 1:00:00
#SBATCH -p shared
#SBATCH -c 20
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL
#SBATCH --begin=now

set -euo pipefail
module load samtools    # adjust module names/versions to your cluster

# ---------- USER EDITABLE ----------
BAM_DIR="/cfs/klemming/projects/supr/ref_free_cancer/jialin/tcga_test_run/bams"           # dir with input .bam files
OUT_DIR="/cfs/klemming/projects/supr/ref_free_cancer/jialin/tcga_test_run/fastqs"         # where final .fastq.gz go
TMP_BASE="/cfs/klemming/projects/supr/ref_free_cancer/jialin/tcga_test_run/tmp"  # per-job tmp area (change if needed)
USE_COLLATE=true                  # set to false to skip collate (see notes)
COMPRESSION_LEVEL=6               # 1=fastest, 9=best compression, 6=default balance
THREADS=${SLURM_CPUS_PER_TASK:-20}
SAMTOOLS_THREADS=$((THREADS > 1 ? THREADS - 1 : 1))  # Leave 1 core for I/O
# -----------------------------------

mkdir -p "${OUT_DIR}"
mkdir -p "${TMP_BASE}"
cd "${TMP_BASE}"

# Build list of BAMs (stable ordering) and pick the file for this array index
mapfile -t BAMFILES < <(find "${BAM_DIR}" -maxdepth 1 -type f -name "*.bam" | sort)

# Check if any BAM files were found
if [ "${#BAMFILES[@]}" -eq 0 ]; then
  echo "Error: No BAM files found in ${BAM_DIR}"
  exit 1
fi

# Check if this array task ID is within bounds
if [ ${SLURM_ARRAY_TASK_ID:-1} -gt "${#BAMFILES[@]}" ]; then
  echo "Warning: Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of BAMs (${#BAMFILES[@]}). Exiting."
  exit 0  # Exit cleanly so SLURM doesn't mark as failed
fi

BAM="${BAMFILES[$((SLURM_ARRAY_TASK_ID-1))]}"
sample="$(basename "${BAM}" .bam)"
echo "[$(date)] Processing sample: ${sample}"
echo "Input BAM: ${BAM}"

# Temporary filenames
COLL_BAM="${TMP_BASE}/${sample}.collated.bam"
FINAL="${OUT_DIR}/${sample}.fastq.gz"

# 1) Optionally collate to group read pairs (lighter than full name-sort)
if [ "${USE_COLLATE}" = true ]; then
  echo "[$(date)] Running samtools collate (threads=${SAMTOOLS_THREADS})..."
  samtools collate -@ "${SAMTOOLS_THREADS}" -o "${COLL_BAM}" "${BAM}"
  INPUT_FOR_FASTQ="${COLL_BAM}"
else
  echo "[$(date)] Skipping collate; using BAM directly"
  INPUT_FOR_FASTQ="${BAM}"
fi

# 2) Run samtools fastq with -O to preserve original qualities
#    Output interleaved paired+singleton reads, discard only supplementary/secondary
#    This keeps all primary alignments including unmapped reads (good for reference-free analysis)
echo "[$(date)] Running samtools fastq -O (outputting interleaved format with unmapped reads)..."
samtools fastq -@ "${SAMTOOLS_THREADS}" -O \
    -0 /dev/null \
    "${INPUT_FOR_FASTQ}" | pigz -p "${THREADS}" -"${COMPRESSION_LEVEL}" -c > "${FINAL}"

# Validate output
if [ ! -s "${FINAL}" ]; then
  echo "Error: Output file ${FINAL} is empty or doesn't exist!"
  exit 1
fi

# 3) cleanup
echo "[$(date)] Cleaning up temporary files..."
if [ "${USE_COLLATE}" = true ]; then
  rm -f "${COLL_BAM}"
fi

echo "[$(date)] Done for ${sample}. Output: ${FINAL}"
