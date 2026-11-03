#!/usr/bin/env bash
# Step 10 - Build a StrainScan database from many genomes of one species (e.g. C. acnes).
# Needs a local clone of StrainScan (see README); set STRAINSCAN_DIR in config/config.sh.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd python realpath
[[ -d "$STRAINSCAN_DIR" ]] || { echo "ERROR: STRAINSCAN_DIR '$STRAINSCAN_DIR' not found" >&2; exit 1; }
[[ -d "$STRAIN_GENOMES_DIR" ]] || { echo "ERROR: genome folder '$STRAIN_GENOMES_DIR' not found" >&2; exit 1; }

# Resolve to absolute paths before changing directory
GENOMES="$(realpath "$STRAIN_GENOMES_DIR")"
SUBSAMPLE_OUT="$(realpath -m "$STRAINSCAN_OUT_DIR/${STRAIN_DB_NAME}_Subsample_Output_Dir")"
DB_OUT="$(realpath -m "$STRAINSCAN_OUT_DIR/$STRAIN_DB_NAME")"
mkdir -p "$STRAINSCAN_OUT_DIR"

cd "$STRAINSCAN_DIR"
log "StrainScan: subsample genomes"
python StrainScan_subsample.py -i "$GENOMES" -o "$SUBSAMPLE_OUT"

log "StrainScan: build database"
python StrainScan_build.py -i "$SUBSAMPLE_OUT/Rep_ref" -o "$DB_OUT" -t "$THREADS"
