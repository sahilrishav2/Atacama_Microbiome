#!/usr/bin/env bash
# Step 8 - Download NCBI RefSeq viral proteins and build a DIAMOND database (run once).
# Env: atacama-core (wget, diamond)
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd wget diamond
mkdir -p "$VIRAL_DB_DIR"

if [[ ! -s "$VIRAL_FAA" ]]; then
  log "Downloading RefSeq viral protein FASTA files"
  # RefSeq viral proteins may be split over several files (viral.1.protein.faa.gz, viral.2...).
  # All are merged here. To reproduce a run that used only viral.1, point VIRAL_FAA in
  # config/config.sh at that single file instead.
  wget --no-verbose --recursive --no-parent --no-directories \
       --accept 'viral.*.protein.faa.gz' --directory-prefix "$VIRAL_DB_DIR" \
       https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
  zcat "$VIRAL_DB_DIR"/viral.*.protein.faa.gz > "$VIRAL_FAA"
fi

log "diamond makedb"
diamond makedb --in "$VIRAL_FAA" -d "$VIRAL_DMND" --threads "$DIAMOND_THREADS"
