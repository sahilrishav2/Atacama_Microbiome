#!/usr/bin/env bash
# Step 2 - Read QC (FastQC + MultiQC) and trimming (fastp), then QC again on trimmed reads.
# Env: atacama-core (fastqc, multiqc, fastp)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd fastqc multiqc fastp
compgen -G "$RAW_DIR/*.fastq.gz" > /dev/null || die "no FASTQ files in $RAW_DIR (run step 1 first)"
mkdir -p "$QC_RAW_DIR" "$QC_TRIM_DIR" "$FASTP_REPORT_DIR" "$TRIM_DIR"

log "FastQC on raw reads"
fastqc "$RAW_DIR"/*.fastq.gz --threads "$THREADS" --outdir "$QC_RAW_DIR"
multiqc "$QC_RAW_DIR" --outdir "$QC_RAW_DIR" --filename multiqc_raw --force

for SAMPLE in "${SAMPLES[@]}"; do
  R1="$RAW_DIR/${SAMPLE}_1.fastq.gz";           R2="$RAW_DIR/${SAMPLE}_2.fastq.gz"
  O1="$TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz";  O2="$TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz"

  if have_files "$O1" "$O2"; then
    log "$SAMPLE: already trimmed, skipping"
    continue
  fi
  have_files "$R1" "$R2" || { log "$SAMPLE: raw reads not found in $RAW_DIR, skipping"; continue; }

  log "$SAMPLE: fastp"
  fastp -i "$R1" -I "$R2" -o "$O1" -O "$O2" \
    --detect_adapter_for_pe --thread "$FASTP_THREADS" \
    --cut_front --cut_tail --cut_mean_quality 20 \
    --html "$FASTP_REPORT_DIR/${SAMPLE}.fastp.html" \
    --json "$FASTP_REPORT_DIR/${SAMPLE}.fastp.json"
done

log "FastQC on trimmed reads"
fastqc "$TRIM_DIR"/*.fastq.gz --threads "$THREADS" --outdir "$QC_TRIM_DIR"
multiqc "$QC_TRIM_DIR" --outdir "$QC_TRIM_DIR" --filename multiqc_trimmed --force
