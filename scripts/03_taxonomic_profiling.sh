#!/usr/bin/env bash
# Step 3 - Taxonomic profiling of trimmed reads with MetaPhlAn.
# The marker database is pinned with METAPHLAN_INDEX (config/config.sh) because step 4 (HUMAnN)
# only accepts profiles made with the vJun23 (or MetaPhlAn 3) database.
# Env: atacama-metaphlan   (run 'metaphlan --install --index <METAPHLAN_INDEX>' once, see README)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd metaphlan merge_metaphlan_tables.py
mkdir -p "$METAPHLAN_DIR"

DB_ARGS=()
[[ -n "$METAPHLAN_DB_DIR" ]] && DB_ARGS=(--bowtie2db "$METAPHLAN_DB_DIR")

for SAMPLE in "${SAMPLES[@]}"; do
  R1="$TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz"
  R2="$TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz"
  BT2="$METAPHLAN_DIR/${SAMPLE}.bowtie2.bz2"
  OUT="$METAPHLAN_DIR/${SAMPLE}_metagenome.txt"

  if [[ -s "$OUT" ]]; then
    log "$SAMPLE: profile exists, skipping"
    continue
  fi
  have_files "$R1" "$R2" || { log "$SAMPLE: trimmed reads not found in $TRIM_DIR (run step 2), skipping"; continue; }

  log "$SAMPLE: MetaPhlAn"
  metaphlan "$R1,$R2" --input_type fastq \
    --index "$METAPHLAN_INDEX" ${DB_ARGS[@]+"${DB_ARGS[@]}"} \
    --bowtie2out "$BT2" --nproc "$THREADS" \
    -o "$OUT" --stat_q 0.15 \
    --ignore_eukaryotes --ignore_archaea
done

compgen -G "$METAPHLAN_DIR/*_metagenome.txt" > /dev/null || die "no MetaPhlAn profiles found in $METAPHLAN_DIR"
log "Merging per-sample profiles"
merge_metaphlan_tables.py "$METAPHLAN_DIR"/*_metagenome.txt > "$METAPHLAN_DIR/merged_abundance_table.txt"
