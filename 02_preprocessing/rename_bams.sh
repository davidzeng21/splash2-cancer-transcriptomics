#!/bin/bash
# Script to rename BAM files based on sample IDs from GDC sample sheet

set -uo pipefail

# ---------- USER EDITABLE ----------
BAM_DIR="/cfs/klemming/projects/supr/ref_free_cancer/jialin/tcga_endometrial/bams"
SAMPLE_SHEET="/cfs/klemming/home/j/jlzeng/cancer_proj/jialin/tcga_endometrial/gdc_sample_sheet.2025-11-18.tsv"
DRY_RUN=${1:-true}  # Set to false to actually rename, true for dry-run
# -----------------------------------

if [ ! -d "${BAM_DIR}" ]; then
    echo "ERROR: BAM directory ${BAM_DIR} does not exist!"
    exit 1
fi

if [ ! -f "${SAMPLE_SHEET}" ]; then
    echo "ERROR: Sample sheet ${SAMPLE_SHEET} does not exist!"
    exit 1
fi

echo "========================================="
echo "BAM File Renaming Script"
echo "========================================="
echo "BAM Directory: ${BAM_DIR}"
echo "Sample Sheet:  ${SAMPLE_SHEET}"
echo "Mode:          $([ "${DRY_RUN}" = "true" ] && echo "DRY RUN" || echo "ACTUAL RENAME")"
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
tail -n +2 "${SAMPLE_SHEET}" | awk -F'\t' '{if (NF >= 7 && $2 != "") print $2"\t"$7}' > "${TEMP_FILE}"

# Read the mappings
while IFS=$'\t' read -r file_name sample_id; do
    # Trim whitespace
    file_name=$(echo "${file_name}" | xargs)
    sample_id=$(echo "${sample_id}" | xargs)
    
    # Skip empty lines
    if [ -z "${file_name}" ] || [ -z "${sample_id}" ]; then
        continue
    fi
    
    TOTAL=$((TOTAL + 1))
    
    # Construct full paths
    OLD_PATH="${BAM_DIR}/${file_name}"
    NEW_NAME="${sample_id}.bam"
    NEW_PATH="${BAM_DIR}/${NEW_NAME}"
    
    # Check if source file exists
    if [ ! -f "${OLD_PATH}" ]; then
        echo "❌ MISSING: ${file_name}"
        MISSING=$((MISSING + 1))
        continue
    fi
    
    # Check if file is already renamed
    if [ "${file_name}" = "${NEW_NAME}" ]; then
        echo "✓ Already renamed: ${file_name}"
        ALREADY_RENAMED=$((ALREADY_RENAMED + 1))
        continue
    fi
    
    # Check if target name already exists (collision)
    if [ -f "${NEW_PATH}" ] && [ "${OLD_PATH}" != "${NEW_PATH}" ]; then
        echo "⚠️  COLLISION: ${NEW_NAME} already exists! Skipping ${file_name}"
        MISSING=$((MISSING + 1))
        continue
    fi
    
    # Perform rename (or dry-run)
    if [ "${DRY_RUN}" = "true" ]; then
        echo "WOULD RENAME: ${file_name} → ${NEW_NAME}"
        SUCCESS=$((SUCCESS + 1))
    else
        if mv "${OLD_PATH}" "${NEW_PATH}"; then
            echo "✓ RENAMED: ${file_name} → ${NEW_NAME}"
            SUCCESS=$((SUCCESS + 1))
        else
            echo "❌ FAILED: Could not rename ${file_name}"
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
