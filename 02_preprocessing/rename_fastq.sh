#!/bin/bash
# Script to rename FASTQ files based on sample IDs from GDC sample sheet

set -uo pipefail

# ---------- USER EDITABLE ----------
FASTQ_DIR="/cfs/klemming/projects/supr/ref_free_cancer/jialin/CPTAC_Lung_Adeno_matched/fastqs"
SAMPLE_SHEET="gdc_sample_sheet.filtered.tsv"
DRY_RUN=${1:-true}  # Set to false to actually rename, true for dry-run
# -----------------------------------

if [ ! -d "${FASTQ_DIR}" ]; then
    echo "ERROR: FASTQ directory ${FASTQ_DIR} does not exist!"
    exit 1
fi

if [ ! -f "${SAMPLE_SHEET}" ]; then
    echo "ERROR: Sample sheet ${SAMPLE_SHEET} does not exist!"
    exit 1
fi

echo "========================================="
echo "FASTQ File Renaming Script"
echo "========================================="
echo "FASTQ Directory: ${FASTQ_DIR}"
echo "Sample Sheet:    ${SAMPLE_SHEET}"
echo "Mode:            $([ "${DRY_RUN}" = "true" ] && echo "DRY RUN" || echo "ACTUAL RENAME")"
echo ""

# Check if running in dry-run mode
if [ "${DRY_RUN}" = "true" ]; then
    echo "⚠️  DRY RUN MODE - No files will be renamed"
    echo "   Run with 'false' as argument to actually rename files"
    echo ""
fi

# Counter variables
TOTAL=0
SUCCESS=0
MISSING=0
ALREADY_RENAMED=0

# Create a temporary file with the mappings
TEMP_FILE=$(mktemp)
tail -n +2 "${SAMPLE_SHEET}" | awk -F'\t' '{if (NF >= 8 && $2 != "") print $2"\t"$6"\t"$8}' > "${TEMP_FILE}"

# Read the mappings
while IFS=$'\t' read -r file_name case_id tissue_type; do
    # Trim whitespace
    file_name=$(echo "${file_name}" | xargs)
    # Extract first value from comma-separated list and trim whitespace
    case_id=$(echo "${case_id}" | cut -d',' -f1 | xargs)
    tissue_type=$(echo "${tissue_type}" | cut -d',' -f1 | xargs | tr '[:upper:]' '[:lower:]')
    
    # Skip empty lines
    if [ -z "${file_name}" ] || [ -z "${case_id}" ] || [ -z "${tissue_type}" ]; then
        continue
    fi
    
    TOTAL=$((TOTAL + 1))
    
    # Convert BAM filename to FASTQ filename (remove .bam, add .fastq.gz)
    # Sample sheet has BAM filenames, but we need to find corresponding FASTQ files
    fastq_base=$(echo "${file_name}" | sed 's/\.bam$//')
    fastq_name="${fastq_base}.fastq.gz"
    
    # Construct full paths
    OLD_PATH="${FASTQ_DIR}/${fastq_name}"
    NEW_NAME="${case_id}_${tissue_type}.fastq.gz"
    NEW_PATH="${FASTQ_DIR}/${NEW_NAME}"
    
    # Check if source file exists
    if [ ! -f "${OLD_PATH}" ]; then
        echo "❌ MISSING: ${fastq_name}"
        MISSING=$((MISSING + 1))
        continue
    fi
    
    # Check if file is already renamed
    if [ "${fastq_name}" = "${NEW_NAME}" ]; then
        echo "✓ Already renamed: ${fastq_name}"
        ALREADY_RENAMED=$((ALREADY_RENAMED + 1))
        continue
    fi
    
    # Check if target name already exists (collision)
    if [ -f "${NEW_PATH}" ] && [ "${OLD_PATH}" != "${NEW_PATH}" ]; then
        echo "⚠️  COLLISION: ${NEW_NAME} already exists! Skipping ${fastq_name}"
        MISSING=$((MISSING + 1))
        continue
    fi
    
    # Perform rename (or dry-run)
    if [ "${DRY_RUN}" = "true" ]; then
        echo "WOULD RENAME: ${fastq_name} → ${NEW_NAME}"
        SUCCESS=$((SUCCESS + 1))
    else
        if mv "${OLD_PATH}" "${NEW_PATH}"; then
            echo "✓ RENAMED: ${fastq_name} → ${NEW_NAME}"
            SUCCESS=$((SUCCESS + 1))
        else
            echo "❌ FAILED: Could not rename ${fastq_name}"
            MISSING=$((MISSING + 1))
        fi
    fi
    
done < "${TEMP_FILE}"

# Cleanup temp file
rm -f "${TEMP_FILE}"

echo ""
echo "========================================="
echo "Summary:"
echo "========================================="
echo "Total entries in TSV:    ${TOTAL}"
echo "Successfully processed:  ${SUCCESS}"
echo "Already renamed:         ${ALREADY_RENAMED}"
echo "Missing or failed:       ${MISSING}"
echo ""

if [ "${DRY_RUN}" = "true" ]; then
    echo "========================================="
    echo "⚠️  This was a DRY RUN - no changes made"
    echo "To actually rename files, run:"
    echo "  $0 false"
    echo "========================================="
else
    echo "========================================="
    echo "✓ Renaming complete!"
    echo "========================================="
fi
