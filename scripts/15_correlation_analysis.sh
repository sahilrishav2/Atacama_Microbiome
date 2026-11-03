#!/usr/bin/env bash
# Step 15 - Correlation analyses in R:
#   (a) abundance of one bacterium (default Cutibacterium acnes, from MetaPhlAn) vs bacteriophages
#   (b) bacteriophages vs the most abundant bacterial pathways (from HUMAnN)
# Needs a bacteriophage abundance (TPM) table: made by step 14, or your own (PHAGE_TPM, see README "Step 15 inputs").
# Env: atacama-r
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd Rscript

MERGED_TAXA="$METAPHLAN_DIR/merged_abundance_table.txt"
MERGED_PATHWAYS="$HUMANN_DIR/merged/merged_pathabundance_cpm.tsv"
have_files "$MERGED_TAXA"     || die "$MERGED_TAXA not found (run step 3 first)"
have_files "$MERGED_PATHWAYS" || die "$MERGED_PATHWAYS not found (run step 4 first)"
have_files "$PHAGE_TPM"       || die "bacteriophage abundance table '$PHAGE_TPM' not found.
       Run step 14 first, or set PHAGE_TPM to your own CSV (one phage per row, one sample per column,
       see README, 'Step 15 inputs')."

INPUT_DIR="$CORR_DIR/input"
mkdir -p "$INPUT_DIR"

log "Preparing correlation inputs"
Rscript analysis/prepare_correlation_inputs.R \
  --metaphlan "$MERGED_TAXA" --pathways "$MERGED_PATHWAYS" \
  --species "$CORR_SPECIES" --top-n "$CORR_TOP_PATHWAYS" --outdir "$INPUT_DIR"

log "${CORR_LABEL} vs bacteriophages"
Rscript analysis/correlate_cacnes_phages.R \
  --bacteria "$INPUT_DIR/cutibacterium_abundance.csv" --phage "$PHAGE_TPM" \
  --label "$CORR_LABEL" --outdir "$CORR_DIR"

log "Bacteriophages vs pathways"
Rscript analysis/correlate_phages_pathways.R \
  --phage "$PHAGE_TPM" --pathways "$INPUT_DIR/abundant_pathways.csv" --outdir "$CORR_DIR"

log "Done. Results are in $CORR_DIR"
