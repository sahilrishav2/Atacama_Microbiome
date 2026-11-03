#!/usr/bin/env bash
# One-time download of the SILVA 138.1 training set formatted for DADA2 (silva_nr99_v138.1_wSpecies_train_set.fa.gz,
# about 140 MB) used by Amplicon_data_processing.Rmd for the taxonomy assignment.
# SILVA is licensed CC-BY 4.0; please cite SILVA (Quast et al. 2013) and DADA2 (Callahan et al. 2016).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd wget gzip
if have_files "$SILVA_DB"; then
  log "SILVA training set already present: $SILVA_DB"
  exit 0
fi
mkdir -p "$(dirname "$SILVA_DB")"
log "Downloading $SILVA_URL"
wget -O "$SILVA_DB.partial" "$SILVA_URL"
gzip -t "$SILVA_DB.partial" || die "the downloaded file is not a valid gzip file; delete it and try again, or download it from https://benjjneb.github.io/dada2/training.html"
mv "$SILVA_DB.partial" "$SILVA_DB"
log "Done: $SILVA_DB"
