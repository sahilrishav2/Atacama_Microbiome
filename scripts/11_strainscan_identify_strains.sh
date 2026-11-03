#!/usr/bin/env bash
# Step 11 - Strain identification with StrainScan. Every read pair in STRAIN_FASTQ_DIR (default: the trimmed
# reads of step 2) is searched against every database listed in STRAINSCAN_DB_LIST; if the list is empty, the
# database built in step 10 is used.
# Results:  results/11_strainscan_identify/<sample>/<database>/
# Run details go to processing.log and per-run time/memory to resource_usage.tsv (same folder).
# Needs a local clone of StrainScan (STRAINSCAN_DIR) with its dependencies, as in step 10.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd python realpath
[[ -f "$STRAINSCAN_DIR/StrainScan.py" ]] \
  || die "StrainScan.py not found in STRAINSCAN_DIR='$STRAINSCAN_DIR' (see README: StrainScan)"

load_read_pairs "$STRAIN_FASTQ_DIR"
(( ${#PAIR_NAMES[@]} > 0 )) || die "no read pairs found in $STRAIN_FASTQ_DIR"
load_db_list "$STRAINSCAN_DB_LIST" "$STRAINSCAN_OUT_DIR/$STRAIN_DB_NAME"
for i in "${!DBS[@]}"; do
  [[ -d "${DBS[i]}" ]] || die "StrainScan database folder not found: ${DBS[i]} (build it with step 10 or fix $STRAINSCAN_DB_LIST)"
  DBS[i]="$(realpath "${DBS[i]}")"      # absolute, because StrainScan is run from its own folder
done

OUT_ROOT="$(realpath -m "$STRAINSCAN_ID_DIR")"
LOG="$(realpath -m "${STRAINSCAN_LOG:-$OUT_ROOT/processing.log}")"
USAGE="$OUT_ROOT/resource_usage.tsv"
SUFFIX="${STRAIN_OUTPUT_SUFFIX:+_$STRAIN_OUTPUT_SUFFIX}"
mkdir -p "$OUT_ROOT" "$(dirname "$LOG")"
usage_table_init "$USAGE"
{
  echo
  echo "===== $(date '+%F %T') StrainScan: ${#PAIR_NAMES[@]} read pair(s), ${#DBS[@]} database(s)"
  printf '  database: %s\n' "${DBS[@]}"
} >> "$LOG"

FAILED=0
for i in "${!PAIR_NAMES[@]}"; do
  SAMPLE="${PAIR_NAMES[i]}"
  R1="$(realpath "${PAIR_R1[i]}")"
  R2="$(realpath "${PAIR_R2[i]}")"
  for DB in "${DBS[@]}"; do
    DB_NAME="$(basename "$DB")"
    OUT="$OUT_ROOT/$SAMPLE/${DB_NAME}${SUFFIX}"
    if [[ -e "$OUT/.done" ]]; then
      log "$SAMPLE x $DB_NAME: already done, skipping"
      continue
    fi
    rm -rf "$OUT"                          # leftovers of an interrupted run
    mkdir -p "$(dirname "$OUT")"
    log "$SAMPLE x $DB_NAME: StrainScan"
    if (cd "$STRAINSCAN_DIR" && run_logged "$LOG" "$USAGE" StrainScan "$SAMPLE" "$DB_NAME" scan -- \
          python StrainScan.py -i "$R1" -j "$R2" -d "$DB" -o "$OUT"); then
      mkdir -p "$OUT"
      touch "$OUT/.done"
    else
      log "$SAMPLE x $DB_NAME: FAILED (see $LOG)"
      FAILED=$((FAILED + 1))
    fi
  done
done

summarize_usage "$USAGE" "$LOG"
(( FAILED == 0 )) || die "$FAILED run(s) failed; details in $LOG"
log "Done. Results are in $OUT_ROOT"
