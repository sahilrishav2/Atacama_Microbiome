#!/usr/bin/env bash
# Step 14 - Bacteriophage abundance (TPM) per sample from DIAMOND hits against the viral proteins of step 8.
#   1. every trimmed read file is searched with diamond blastx (best hit per read)
#   2. analysis/compute_phage_tpm.py assigns each read pair to the phage of its best hit, divides the count by the
#      phage's total protein length and scales every sample to 1e6 (TPM); see the header of that script
# Result: results/14_phage_abundance/bacteriophage_abundance_tpm.csv (rows = phages, columns = samples),
#         used by step 15.
# Env: atacama-core (diamond, python)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd diamond python3
[[ -s "$VIRAL_DMND.dmnd" || -s "$VIRAL_DMND" ]] || die "viral DIAMOND database '$VIRAL_DMND' not found (run step 8 first)"
have_files "$VIRAL_FAA" || die "viral protein FASTA '$VIRAL_FAA' not found; it is needed for the phage lengths (step 8 creates it)"

HITS="$PHAGE_DIR/diamond"
mkdir -p "$HITS"

for SAMPLE in "${SAMPLES[@]}"; do
  for MATE in 1 2; do
    READS="$TRIM_DIR/${SAMPLE}_${MATE}.trimmed.fastq.gz"
    OUT="$HITS/${SAMPLE}_R${MATE}.tsv"
    if [[ -e "$OUT" ]]; then
      log "$SAMPLE R$MATE: DIAMOND output exists, skipping"
      continue
    fi
    have_files "$READS" || { log "$SAMPLE R$MATE: trimmed reads not found in $TRIM_DIR (run step 2), skipping"; continue; }
    log "$SAMPLE R$MATE: DIAMOND blastx"
    diamond blastx -d "$VIRAL_DMND" -q "$READS" -o "$OUT.partial" \
      --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
      --max-target-seqs 1 --evalue "$PHAGE_EVALUE" --threads "$DIAMOND_THREADS"
    mv "$OUT.partial" "$OUT"
  done
done

log "Counting read pairs per phage and calculating TPM"
python3 analysis/compute_phage_tpm.py \
  --hits-dir "$HITS" --reference-faa "$VIRAL_FAA" --outdir "$PHAGE_DIR" \
  --name-filter "$PHAGE_NAME_FILTER" --min-identity "$PHAGE_MIN_IDENTITY" --min-fragments "$PHAGE_MIN_FRAGMENTS" \
  --samples "${SAMPLES[@]}"
log "Done. Phage abundance table: $PHAGE_DIR/bacteriophage_abundance_tpm.csv"
