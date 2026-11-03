#!/usr/bin/env bash
# Processes ONE 16S amplicon dataset with DADA2 (Amplicon_data/Amplicon_data_processing.Rmd): quality filtering,
# ASVs, chimera removal, SILVA taxonomy, and a genus table for the merging step. Run it once per dataset.
#
#   bash scripts/amplicon_process_dataset.sh <dataset> <fastq_dir>                       # only the quality plots
#   bash scripts/amplicon_process_dataset.sh <dataset> <fastq_dir> <trunc_f> <trunc_r>   # full processing
#
# <dataset>   a short name (letters, digits, . _ -), e.g. the BioProject; it names the results and the genus table
# <fastq_dir> folder with <sample>_1.fastq[.gz] and <sample>_2.fastq[.gz]; step 1 produces such files from SRA runs
# <trunc_f>/<trunc_r> truncation lengths of the forward/reverse reads, chosen from the quality plots of the first call
# Results: results/amplicon/<dataset>/ and results/amplicon/genus_tables/<dataset>.tsv, report in results/amplicon/reports/
# Env: atacama-amplicon (and scripts/setup_silva_database.sh once)
set -euo pipefail
USAGE="usage: bash scripts/amplicon_process_dataset.sh <dataset> <fastq_dir> [<trunc_f> <trunc_r>]"
if [[ $# -ne 2 && $# -ne 4 ]]; then echo "$USAGE" >&2; exit 1; fi
DATASET="$1"
FASTQ_DIR="$(cd "$2" 2> /dev/null && pwd)" || { echo "ERROR: read folder '$2' not found" >&2; exit 1; }
[[ "$DATASET" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "ERROR: dataset name '$DATASET' may only contain letters, digits, . _ -" >&2; exit 1; }
QUALITY_ONLY=TRUE; TRUNC_F=-1; TRUNC_R=-1
if [[ $# -eq 4 ]]; then
  [[ "$3" =~ ^[0-9]+$ && "$4" =~ ^[0-9]+$ ]] || { echo "ERROR: truncation lengths must be whole numbers" >&2; exit 1; }
  QUALITY_ONLY=FALSE; TRUNC_F="$3"; TRUNC_R="$4"
fi

cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd Rscript pandoc
Rscript -e 'pk <- c("dada2", "phyloseq", "rmarkdown", "knitr")
            miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
            if (length(miss) > 0) { cat("Missing R package(s):", miss, "\n"); quit(status = 1) }' \
  || die "install the missing R packages or activate the atacama-amplicon environment (envs/amplicon.yml)"
if [[ "$QUALITY_ONLY" == FALSE ]]; then
  have_files "$SILVA_DB" || die "SILVA training set '$SILVA_DB' not found: run scripts/setup_silva_database.sh"
  REPORT="${DATASET}_processing.html"
else
  REPORT="${DATASET}_quality.html"
  log "No truncation lengths given: only the read-quality plots are made (see the report), then run again with them"
fi

log "DADA2 processing of '$DATASET' ($FASTQ_DIR)"
RMDPARAM_dataset="$DATASET" RMDPARAM_fastq_dir="$FASTQ_DIR" \
RMDPARAM_trunc_len_f="$TRUNC_F" RMDPARAM_trunc_len_r="$TRUNC_R" RMDPARAM_quality_only="$QUALITY_ONLY" \
RMDPARAM_out_dir="$AMPLICON_DIR" RMDPARAM_genus_dir="$AMPLICON_DIR/genus_tables" \
RMDPARAM_silva_db="$SILVA_DB" RMDPARAM_multithread="$THREADS" \
  Rscript analysis/render_rmd.R Amplicon_data/Amplicon_data_processing.Rmd "$AMPLICON_DIR/reports" "$REPORT"
log "Done. Report: $AMPLICON_DIR/reports/$REPORT"
