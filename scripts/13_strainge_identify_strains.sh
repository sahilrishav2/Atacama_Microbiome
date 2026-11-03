#!/usr/bin/env bash
# Step 13 - Strain identification with StrainGE (StrainGST). Every read pair in STRAIN_FASTQ_DIR (default: the
# trimmed reads of step 2) is k-merized once and then searched against every database listed in
# STRAINGE_DB_LIST; if the list is empty, the database built in step 12 is used.
# Results:  results/13_strainge_identify/<sample>/<sample>_<database>.tsv  (+ the k-merized reads <sample>.hdf5)
# Run details go to processing.log and per-run time/memory to resource_usage.tsv (same folder).
# Env: atacama-strainge
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd straingst realpath

load_read_pairs "$STRAIN_FASTQ_DIR"
(( ${#PAIR_NAMES[@]} > 0 )) || die "no read pairs found in $STRAIN_FASTQ_DIR"
load_db_list "$STRAINGE_DB_LIST" "$STRAINGE_OUT_DIR/$STRAINGE_DB_NAME/$STRAINGE_DB_NAME"

# An entry is either a database prefix (<prefix>_pan-genome-db.hdf5, as produced by step 12) or an *.hdf5 file.
DB_FILES=(); DB_NAMES=()
for ENTRY in "${DBS[@]}"; do
  case "$ENTRY" in
    *.hdf5) FILE="$ENTRY" ;;
    *)      FILE="${ENTRY}_pan-genome-db.hdf5" ;;
  esac
  [[ -s "$FILE" ]] || die "StrainGE database not found: $FILE (build it with step 12 or fix $STRAINGE_DB_LIST)"
  DB_FILES+=("$(realpath "$FILE")")
  NAME="$(basename "$FILE")"
  NAME="${NAME%_pan-genome-db.hdf5}"
  NAME="${NAME%.hdf5}"
  DB_NAMES+=("${NAME//[ .]/_}")
done

OUT_ROOT="$(realpath -m "$STRAINGE_ID_DIR")"
LOG="$(realpath -m "${STRAINGE_LOG:-$OUT_ROOT/processing.log}")"
USAGE="$OUT_ROOT/resource_usage.tsv"
SUFFIX="${STRAIN_OUTPUT_SUFFIX:+_$STRAIN_OUTPUT_SUFFIX}"
mkdir -p "$OUT_ROOT" "$(dirname "$LOG")"
usage_table_init "$USAGE"
{
  echo
  echo "===== $(date '+%F %T') StrainGE: ${#PAIR_NAMES[@]} read pair(s), ${#DB_FILES[@]} database(s), k=$STRAINGE_KMER"
  printf '  database: %s\n' "${DB_FILES[@]}"
} >> "$LOG"

FAILED=0
for i in "${!PAIR_NAMES[@]}"; do
  SAMPLE="${PAIR_NAMES[i]}"
  R1="$(realpath "${PAIR_R1[i]}")"
  R2="$(realpath "${PAIR_R2[i]}")"
  OUT_DIR="$OUT_ROOT/$SAMPLE"
  KMERS="$OUT_DIR/${SAMPLE}${SUFFIX}.hdf5"
  mkdir -p "$OUT_DIR"

  # k-merize the reads once; the result is reused for every database
  if [[ ! -s "$KMERS" ]]; then
    log "$SAMPLE: k-merize reads (k=$STRAINGE_KMER)"
    TMP="$OUT_DIR/.tmp_${SAMPLE}.hdf5"
    rm -f "$TMP"
    if run_logged "$LOG" "$USAGE" StrainGE "$SAMPLE" - kmerize -- \
         straingst kmerize -k "$STRAINGE_KMER" -o "$TMP" "$R1" "$R2"; then
      mv "$TMP" "$KMERS"
    else
      log "$SAMPLE: k-merizing FAILED (see $LOG)"
      rm -f "$TMP"
      FAILED=$((FAILED + 1))
      continue
    fi
  fi

  for j in "${!DB_FILES[@]}"; do
    RESULT="$OUT_DIR/${SAMPLE}${SUFFIX}_${DB_NAMES[j]}.tsv"
    if [[ -s "$RESULT" ]]; then
      log "$SAMPLE x ${DB_NAMES[j]}: already done, skipping"
      continue
    fi
    log "$SAMPLE x ${DB_NAMES[j]}: straingst run"
    rm -f "$RESULT.partial"
    if run_logged "$LOG" "$USAGE" StrainGE "$SAMPLE" "${DB_NAMES[j]}" run -- \
         straingst run -o "$RESULT.partial" "${DB_FILES[j]}" "$KMERS"; then
      mv "$RESULT.partial" "$RESULT"
    else
      log "$SAMPLE x ${DB_NAMES[j]}: FAILED (see $LOG)"
      rm -f "$RESULT.partial"
      FAILED=$((FAILED + 1))
    fi
  done
done

summarize_usage "$USAGE" "$LOG"
(( FAILED == 0 )) || die "$FAILED run(s) failed; details in $LOG"
log "Done. Results are in $OUT_ROOT"
