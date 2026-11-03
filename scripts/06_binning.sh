#!/usr/bin/env bash
# Step 6 - Map reads back to contigs (BWA + samtools), compute depth, bin with MetaBAT2.
# Env: atacama-core (bwa, samtools, metabat2)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd bwa samtools jgi_summarize_bam_contig_depths metabat2

for SAMPLE in "${SAMPLES[@]}"; do
  CONTIGS="$ASSEMBLY_DIR/${SAMPLE}.fa"
  R1="$TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz"
  R2="$TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz"
  WORK="$BINNING_DIR/$SAMPLE"
  BAM="$WORK/${SAMPLE}.sorted.bam"
  DEPTH="$WORK/${SAMPLE}.depth.txt"
  BINS="$WORK/bins"
  mkdir -p "$BINS"

  if compgen -G "$BINS/bin.*.fa" > /dev/null; then
    log "$SAMPLE: bins exist, skipping"
    continue
  fi

  have_files "$CONTIGS" "$R1" "$R2" || { log "$SAMPLE: contigs or trimmed reads missing (run steps 2 and 5), skipping"; continue; }

  [[ -s "${CONTIGS}.bwt" ]] || { log "$SAMPLE: bwa index"; bwa index "$CONTIGS"; }

  if [[ ! -s "${BAM}.bai" ]]; then
    log "$SAMPLE: mapping + sorting"
    bwa mem -t "$THREADS" "$CONTIGS" "$R1" "$R2" \
      | samtools view -@ "$THREADS" -b - \
      | samtools sort -@ "$THREADS" -o "$BAM"
    samtools index "$BAM"
  fi

  log "$SAMPLE: contig depths"
  jgi_summarize_bam_contig_depths --outputDepth "$DEPTH" --referenceFasta "$CONTIGS" "$BAM"

  log "$SAMPLE: MetaBAT2"
  metabat2 -i "$CONTIGS" -a "$DEPTH" -o "$BINS/bin" -t "$THREADS"
done
