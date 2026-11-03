#!/usr/bin/env bash
# Step 5 - Per-sample assembly of trimmed reads with MEGAHIT.
# Env: atacama-core (megahit)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd megahit
mkdir -p "$MEGAHIT_DIR" "$ASSEMBLY_DIR"

for SAMPLE in "${SAMPLES[@]}"; do
  R1="$TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz"
  R2="$TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz"
  OUTDIR="$MEGAHIT_DIR/$SAMPLE"      # MEGAHIT refuses to write into an existing folder...

  if [[ -s "$ASSEMBLY_DIR/${SAMPLE}.fa" ]]; then
    log "$SAMPLE: assembly exists, skipping"
    continue
  fi

  have_files "$R1" "$R2" || { log "$SAMPLE: trimmed reads not found in $TRIM_DIR (run step 2), skipping"; continue; }

  EXTRA=()
  [[ -d "$OUTDIR" ]] && EXTRA=(--continue)   # ...unless we resume an interrupted run

  log "$SAMPLE: MEGAHIT"
  megahit -1 "$R1" -2 "$R2" -t "$THREADS" -o "$OUTDIR" ${EXTRA[@]+"${EXTRA[@]}"}

  cp "$OUTDIR/final.contigs.fa" "$ASSEMBLY_DIR/${SAMPLE}.fa"
done
