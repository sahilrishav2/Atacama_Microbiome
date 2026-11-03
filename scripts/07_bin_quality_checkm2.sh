#!/usr/bin/env bash
# Step 7 - Completeness / contamination of every bin with CheckM2.
# Env: atacama-checkm2   (run 'checkm2 database --download' once)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
load_samples
require_cmd checkm2
mkdir -p "$CHECKM2_DIR"

for SAMPLE in "${SAMPLES[@]}"; do
  BINS="$BINNING_DIR/$SAMPLE/bins"
  OUT="$CHECKM2_DIR/$SAMPLE"

  if ! compgen -G "$BINS/bin.*.fa" > /dev/null; then
    log "$SAMPLE: no bins found in $BINS, skipping"
    continue
  fi
  if [[ -s "$OUT/quality_report.tsv" ]]; then
    log "$SAMPLE: CheckM2 report exists, skipping"
    continue
  fi

  log "$SAMPLE: CheckM2"
  checkm2 predict --threads "$CHECKM2_THREADS" --input "$BINS" -x fa --output-directory "$OUT"
done

# One combined table across samples
COMBINED="$CHECKM2_DIR/all_samples_quality_report.tsv"
header_written=false
for SAMPLE in "${SAMPLES[@]}"; do
  REPORT="$CHECKM2_DIR/$SAMPLE/quality_report.tsv"
  [[ -s "$REPORT" ]] || continue
  if ! $header_written; then
    { printf 'Sample\t'; head -n 1 "$REPORT"; } > "$COMBINED"
    header_written=true
  fi
  tail -n +2 "$REPORT" | awk -v s="$SAMPLE" 'BEGIN{OFS="\t"} {print s, $0}' >> "$COMBINED"
done
if $header_written; then log "Combined report: $COMBINED"; fi
