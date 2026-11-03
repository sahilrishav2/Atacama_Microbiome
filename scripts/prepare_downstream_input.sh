#!/usr/bin/env bash
# Builds the input tables of the downstream documents. For every metadata file you put in DOWNSTREAM_DIR the
# matching table is written, with its columns in the order of that metadata file:
#   metadata.csv        -> countdata.csv and relative_abundance.csv
#   metadata_plant.csv  -> countdata_plant.csv
#   metadata_soil.csv   -> countdata_soil.csv
# DOWNSTREAM_SOURCE=shotgun (default): the merged MetaPhlAn table of step 3 is converted (rank DOWNSTREAM_RANK,
#   default species). The values are MetaPhlAn relative abundances (percent), not counts, so countdata*.csv and
#   relative_abundance.csv hold the same kind of numbers.
# DOWNSTREAM_SOURCE=amplicon: the prevalence-filtered count table and its relative abundance table, made by
#   scripts/amplicon_analysis.sh, are subset and ordered.
# The metadata files are yours (see downstream_input/README.md); their sample IDs must match those of the tables.
# Env: atacama-core (python)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd python3
mkdir -p "$DOWNSTREAM_DIR"

case "$DOWNSTREAM_SOURCE" in
  shotgun)
    MERGED="$METAPHLAN_DIR/merged_abundance_table.txt"
    have_files "$MERGED" || die "$MERGED not found (run step 3 first)"
    # build_table <output table> <metadata file> <source: counts|relative> [sample column]
    build_table() {
      local out="$DOWNSTREAM_DIR/$1" meta="$DOWNSTREAM_DIR/$2" column="${4:-}" extra=()
      [[ -n "$column" ]] && extra=(--sample-column "$column")
      log "$1 <- $MERGED ($DOWNSTREAM_RANK level, samples of $2)"
      python3 analysis/metaphlan_to_table.py --metaphlan "$MERGED" --rank "$DOWNSTREAM_RANK" \
        --out "$out" --keep-samples-from "$meta" ${extra[@]+"${extra[@]}"}
    }
    ;;
  amplicon)
    PCT="$(awk -v t="$AMPLICON_PREVALENCE" 'BEGIN { printf "%d", 100 * t + 0.5 }')"
    COUNTS="$AMPLICON_DIR/merged_bacteria_${PCT}%_prevalence.csv"
    RELATIVE="$AMPLICON_DIR/relative_abundance_merged_bacteria_${PCT}%_prevalence.csv"
    have_files "$COUNTS" "$RELATIVE" || die "$COUNTS or $RELATIVE not found (run scripts/amplicon_analysis.sh first)"
    build_table() {
      local out="$DOWNSTREAM_DIR/$1" meta="$DOWNSTREAM_DIR/$2" source_table="$COUNTS" column="${4:-}" extra=()
      [[ "$3" == relative ]] && source_table="$RELATIVE"
      [[ -n "$column" ]] && extra=(--sample-column "$column")
      log "$1 <- $source_table (samples of $2)"
      python3 analysis/subset_table.py --table "$source_table" --out "$out" --keep-samples-from "$meta" \
        ${extra[@]+"${extra[@]}"}
    }
    ;;
  *) die "DOWNSTREAM_SOURCE must be 'shotgun' or 'amplicon' (got '$DOWNSTREAM_SOURCE')" ;;
esac

MADE=0
if have_files "$DOWNSTREAM_DIR/metadata.csv"; then
  build_table countdata.csv metadata.csv counts
  build_table relative_abundance.csv metadata.csv relative
  MADE=$((MADE + 1))
fi
if have_files "$DOWNSTREAM_DIR/metadata_plant.csv"; then
  build_table countdata_plant.csv metadata_plant.csv counts Sample
  MADE=$((MADE + 1))
fi
if have_files "$DOWNSTREAM_DIR/metadata_soil.csv"; then
  build_table countdata_soil.csv metadata_soil.csv counts Sample
  MADE=$((MADE + 1))
fi
(( MADE > 0 )) || die "no metadata file (metadata.csv, metadata_plant.csv, metadata_soil.csv) found in $DOWNSTREAM_DIR; see downstream_input/README.md"
log "Done. Tables are in $DOWNSTREAM_DIR"
