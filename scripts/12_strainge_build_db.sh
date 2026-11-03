#!/usr/bin/env bash
# Step 12 - Build a custom StrainGE database (StrainGST) from many reference genomes of one species.
# Input: STRAINGE_GENOMES_DIR with the genomes (*.fa) AND their assembly report (*.txt).
# Steps: prepare genomes -> k-merize each -> all-vs-all similarity -> cluster near-identical references ->
#        build the pan-genome database
# Result: results/12_strainge_db/<name>/<name>_pan-genome-db.hdf5
# Env: atacama-strainge; prepare_strainge_db.py (separate download, see README) must be in STRAINGE_DIR.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd python3 straingst
PREPARE="$STRAINGE_DIR/prepare_strainge_db.py"
[[ -f "$PREPARE" ]] || die "prepare_strainge_db.py not found in STRAINGE_DIR='$STRAINGE_DIR' (download it, see README: StrainGE)"
[[ -d "$STRAINGE_GENOMES_DIR" ]] || die "genome folder '$STRAINGE_GENOMES_DIR' not found"
compgen -G "$STRAINGE_GENOMES_DIR/*.fa" > /dev/null || compgen -G "$STRAINGE_GENOMES_DIR/*.fa.gz" > /dev/null \
  || die "no genomes (*.fa) in $STRAINGE_GENOMES_DIR"
compgen -G "$STRAINGE_GENOMES_DIR/*.txt" > /dev/null \
  || die "no assembly report (*.txt) in $STRAINGE_GENOMES_DIR; prepare_strainge_db.py needs it next to the genomes"

NAME="$STRAINGE_DB_NAME"
OUT="$STRAINGE_OUT_DIR/$NAME"
PREP="$OUT/prepared_genomes"
SIMILARITIES="$OUT/${NAME}_similarities.tsv"
CLUSTERS="$OUT/${NAME}_clusters.tsv"
KEEP="$OUT/${NAME}_references_to_keep.txt"
DB="$OUT/${NAME}_pan-genome-db.hdf5"
read -ra CLUSTER_ARGS <<< "$STRAINGE_CLUSTER_ARGS"

if [[ -s "$DB" ]]; then
  log "Database already exists: $DB"
  exit 0
fi
mkdir -p "$OUT"

if [[ ! -e "$OUT/.prepared" ]]; then
  log "StrainGE: prepare genomes"
  rm -rf "$PREP"
  mkdir -p "$PREP"
  python3 "$PREPARE" "$STRAINGE_GENOMES_DIR" -o "$PREP"
  touch "$OUT/.prepared"
fi
compgen -G "$PREP/*.fa.gz" > /dev/null || die "prepare_strainge_db.py produced no *.fa.gz files in $PREP"

log "StrainGE: k-merize genomes (k=$STRAINGE_KMER)"
for GENOME in "$PREP"/*.fa.gz; do
  KMERS="${GENOME%.fa.gz}.hdf5"
  [[ -s "$KMERS" ]] && continue
  rm -f "$OUT/.tmp_kmerize.hdf5"
  straingst kmerize -k "$STRAINGE_KMER" -o "$OUT/.tmp_kmerize.hdf5" "$GENOME"
  mv "$OUT/.tmp_kmerize.hdf5" "$KMERS"
done

if [[ ! -s "$SIMILARITIES" ]]; then
  log "StrainGE: all-vs-all k-mer similarity"
  straingst kmersim --all-vs-all -t "$THREADS" -S jaccard -S subset "$PREP"/*.hdf5 > "$SIMILARITIES.partial"
  mv "$SIMILARITIES.partial" "$SIMILARITIES"
fi

if [[ ! -s "$KEEP" ]]; then
  log "StrainGE: cluster similar references (${CLUSTER_ARGS[*]})"
  straingst cluster -i "$SIMILARITIES" "${CLUSTER_ARGS[@]}" --clusters-out "$CLUSTERS" "$PREP"/*.hdf5 > "$KEEP.partial"
  mv "$KEEP.partial" "$KEEP"
fi

log "StrainGE: create pan-genome database"
rm -f "$OUT/.tmp_db.hdf5"
straingst createdb -f "$KEEP" -o "$OUT/.tmp_db.hdf5"
mv "$OUT/.tmp_db.hdf5" "$DB"
log "Done. Database: $DB"
