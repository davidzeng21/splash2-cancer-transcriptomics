SOURCE_DIR="./"
DEST_DIR="download/"
MANIFEST_FILE="gdc_manifest.2025-11-28.184605.txt"

SUCCESS_COUNT=0
ERROR_LOG="move_uuid_folder_errors.log"
> "$ERROR_LOG" # Clear previous log

echo "--- Starting Move Process ---"

# Read the manifest, skipping the header (tail -n +2)
tail -n +2 "$MANIFEST_FILE" | while IFS=$'\t' read -r UUID FILENAME REST; do
    UUID_FOLDER="${SOURCE_DIR}${UUID}"
    
    # Check if the folder exists
    if [ -d "$UUID_FOLDER" ]; then
        # Use mv to move the entire non-empty UUID folder to the DEST_DIR
        mv "$UUID_FOLDER" "$DEST_DIR"
        
        if [ $? -eq 0 ]; then
            SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
            echo "Moved folder: $UUID"
        else
            echo "❌ ERROR: Failed to move $UUID_FOLDER (Check permissions)." >> "$ERROR_LOG"
        fi
    else
        echo "Skipping: $UUID_FOLDER (Folder does not exist in SOURCE_DIR)." >> "$ERROR_LOG"
    fi
done

# --- Summary ---
echo "--- Process Summary ---"
echo "Total UUID folders moved: **${SUCCESS_COUNT}**"

if [ -s "$ERROR_LOG" ]; then
    echo "❌ **PROBLEM DETECTED:** Errors occurred. Please check the **${ERROR_LOG}** file for details."
else
    echo "✅ **SUCCESS:** All existing UUID folders were moved to ${DEST_DIR}"
fi
