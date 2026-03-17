#!/bin/bash

# ============================================================
# Snippy Pipeline: Run → Rename → Collect CSV
# ============================================================

# ---------- Configuration ----------
IN_DIR="split_output"
REF="ref.gbk"
THREADS=8
RESULTS_DIR="snippy_results"
MAP_FILE="name_mapping.csv"
OUTPUT_DIR="all_csv_files"
# -----------------------------------

set -euo pipefail

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ============================================================
# STEP 1: Run Snippy on all samples
# ============================================================
log "========== STEP 1: Running Snippy =========="
mkdir -p "$RESULTS_DIR"

for R1 in "$IN_DIR"/*_1.fastq.gz; do
    filename=$(basename "$R1")
    sample_name="${filename%_1.fastq.gz}"
    R2="$IN_DIR/${sample_name}_2.fastq.gz"

    if [ ! -f "$R2" ]; then
        log "WARNING: R2 not found for $sample_name, skipping."
        continue
    fi

    log "Processing: $sample_name"

    snippy --cpus "$THREADS" \
           --outdir "$RESULTS_DIR/$sample_name" \
           --ref "$REF" \
           --R1 "$R1" \
           --R2 "$R2" \
           --mincov 10 \
           --minfrac 0.9 \
           --basequal 20 \
           --mapqual 30
done

log "STEP 1 complete."

# ============================================================
# STEP 2: Rename folders and files using mapping CSV
# ============================================================
log "========== STEP 2: Renaming samples =========="

if [ ! -f "$MAP_FILE" ]; then
    log "ERROR: Cannot find $MAP_FILE — skipping rename step."
else
    tail -n +2 "$MAP_FILE" | while IFS=, read -r old_id new_name || [ -n "$old_id" ]; do
        new_name=$(echo "$new_name" | tr -d '\r')
        old_id=$(echo "$old_id"   | tr -d '\r')

        [[ -z "$old_id" || -z "$new_name" ]] && continue

        OLD_PATH="$RESULTS_DIR/$old_id"
        NEW_PATH="$RESULTS_DIR/$new_name"

        if [ -d "$OLD_PATH" ]; then
            if [ -d "$NEW_PATH" ]; then
                log "WARNING: $NEW_PATH already exists, skipping $old_id."
                continue
            fi
            mv "$OLD_PATH" "$NEW_PATH"
            find "$NEW_PATH" -maxdepth 1 -type f -name "snps.*" | while read -r file; do
                dir=$(dirname "$file")
                base=$(basename "$file")
                mv "$file" "$dir/${base/snps/$new_name}"
            done
            log "Renamed: $old_id -> $new_name"
        elif [ -d "$NEW_PATH" ]; then
            log "Skip: $new_name already exists (old ID $old_id not found)."
        else
            log "Skip: Folder not found for ID: $old_id"
        fi
    done
    log "STEP 2 complete."
fi

# ============================================================
# STEP 3: Collect all CSV files into one directory
# ============================================================
log "========== STEP 3: Collecting CSV files =========="
mkdir -p "$OUTPUT_DIR"
count=0

for folder in "$RESULTS_DIR"/*/; do
    sample_name=$(basename "$folder")
    target_file="$folder/${sample_name}.csv"

    if [ -f "$target_file" ]; then
        cp "$target_file" "$OUTPUT_DIR/"
        log "Copied: ${sample_name}.csv"
        ((count++))
    else
        log "WARNING: No ${sample_name}.csv found in $sample_name"
    fi
done

log "STEP 3 complete. Total $count CSV files collected."
log "========== Pipeline Finished =========="
