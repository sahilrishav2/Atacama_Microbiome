#!/usr/bin/env bash
# One-time download of the HUMAnN databases used by step 4:
#   ChocoPhlAn pangenomes, UniRef90 (DIAMOND format) and the utility mapping files (for HUMANN_REGROUP).
# These are large downloads (tens of GB in total); put them somewhere with space via HUMANN_DB_DIR.
# Env: atacama-humann
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd humann_databases
mkdir -p "$HUMANN_DB_DIR"

# fetch <database> <build> <expected folder>
fetch() {
  local db="$1" build="$2" target="$3"
  if [[ -d "$target" && -n "$(ls -A "$target" 2> /dev/null)" ]]; then
    log "$db: $target already populated, skipping"
    return 0
  fi
  log "Downloading $db ($build)"
  humann_databases --download "$db" "$build" "$HUMANN_DB_DIR" --update-config yes
}

fetch chocophlan       full             "$HUMANN_NUC_DB"
fetch uniref           uniref90_diamond "$HUMANN_PROT_DB"
fetch utility_mapping  full             "$HUMANN_UTIL_DB"

log "HUMAnN databases ready in $HUMANN_DB_DIR"
