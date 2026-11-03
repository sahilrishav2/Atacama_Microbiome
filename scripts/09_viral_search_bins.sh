#!/usr/bin/env bash
# Step 9 - Screen high-quality bins (>= MIN_COMPLETENESS, <= MAX_CONTAMINATION) against viral proteins.
# Protein predictions come from CheckM2 (<checkm2_out>/protein_files/<bin>.faa).
# Env: atacama-core (diamond)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd diamond
[[ -s "$VIRAL_DMND.dmnd" || -s "$VIRAL_DMND" ]] || die "viral DIAMOND database '$VIRAL_DMND' not found (run step 8 first)"
mkdir -p "$VIRAL_HITS_DIR"

for SAMPLE in "${SAMPLES[@]}"; do
  REPORT="$CHECKM2_DIR/$SAMPLE/quality_report.tsv"
  PROT_DIR="$CHECKM2_DIR/$SAMPLE/protein_files"
  [[ -s "$REPORT" ]] || { log "$SAMPLE: no CheckM2 report, skipping"; continue; }

  # column 1 = bin name, 2 = completeness, 3 = contamination
  mapfile -t PASS_BINS < <(awk -F'\t' -v c="$MIN_COMPLETENESS" -v x="$MAX_CONTAMINATION" \
                             'NR>1 && $2>=c && $3<=x {print $1}' "$REPORT")
  log "$SAMPLE: ${#PASS_BINS[@]} bin(s) pass (>=${MIN_COMPLETENESS}% complete, <=${MAX_CONTAMINATION}% contaminated)"

  mkdir -p "$VIRAL_HITS_DIR/$SAMPLE"
  for BIN in ${PASS_BINS[@]+"${PASS_BINS[@]}"}; do
    FAA="$PROT_DIR/${BIN}.faa"
    OUT="$VIRAL_HITS_DIR/$SAMPLE/${BIN}.diamond.tsv"
    if [[ ! -s "$FAA" ]]; then
      log "$SAMPLE/$BIN: protein file $FAA not found, skipping"
      continue
    fi
    log "$SAMPLE/$BIN: DIAMOND blastp"
    diamond blastp -d "$VIRAL_DMND" -q "$FAA" -o "$OUT" --sensitive \
      --threads "$DIAMOND_THREADS" \
      --outfmt 6 qseqid sseqid stitle pident length evalue bitscore
  done
done

# Combined table with a header
COMBINED="$VIRAL_HITS_DIR/all_viral_hits.tsv"
printf 'sample\tbin\tqseqid\tsseqid\tstitle\tpident\tlength\tevalue\tbitscore\n' > "$COMBINED"
shopt -s nullglob
for F in "$VIRAL_HITS_DIR"/*/*.diamond.tsv; do
  S="$(basename "$(dirname "$F")")"
  B="$(basename "$F" .diamond.tsv)"
  awk -v s="$S" -v b="$B" 'BEGIN{OFS="\t"} {print s, b, $0}' "$F" >> "$COMBINED"
done
log "Combined hits: $COMBINED"
