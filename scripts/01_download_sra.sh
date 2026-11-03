#!/usr/bin/env bash
# Step 1 - Download SRA runs and convert them to gzipped paired-end FASTQ.
# Env: atacama-core (sra-tools, pigz)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd prefetch fasterq-dump pigz
mkdir -p "$SRA_DIR" "$RAW_DIR" "$TMP_DIR"

# prefetch dislikes comments/blank lines, so give it a cleaned list
CLEAN_LIST="$TMP_DIR/accessions.clean.txt"
printf '%s\n' "${SAMPLES[@]}" > "$CLEAN_LIST"

log "Prefetching ${#SAMPLES[@]} run(s)"
prefetch --option-file "$CLEAN_LIST" --output-directory "$SRA_DIR"

for SAMPLE in "${SAMPLES[@]}"; do
  if [[ -s "$RAW_DIR/${SAMPLE}_1.fastq.gz" && -s "$RAW_DIR/${SAMPLE}_2.fastq.gz" ]]; then
    log "$SAMPLE: FASTQ already present, skipping"
    continue
  fi

  SRA_FILE="$SRA_DIR/$SAMPLE/$SAMPLE.sra"
  [[ -f "$SRA_FILE" ]] || SRA_FILE="$SRA_DIR/$SAMPLE/$SAMPLE.sralite"

  log "$SAMPLE: fasterq-dump"
  fasterq-dump "$SRA_FILE" --split-files -m "$SORT_MEM" -e "$THREADS" \
    --temp "$TMP_DIR" --outdir "$RAW_DIR"

  # fasterq-dump writes plain FASTQ; compress so downstream steps can use *.fastq.gz
  log "$SAMPLE: compressing"
  pigz -p "$THREADS" "$RAW_DIR/${SAMPLE}"_*.fastq
done
log "Done. FASTQ files are in $RAW_DIR"
