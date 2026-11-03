#!/usr/bin/env bash
# 16S amplicon analysis across datasets, from the genus tables made by amplicon_process_dataset.sh:
#   1. merging_multiple_count_abundance_tables.Rmd   -> merged_bacteria.csv
#   2. filtering_prevalent_bacteria.Rmd              -> merged_bacteria_<N>%_prevalence.csv   (N = AMPLICON_PREVALENCE x 100)
#   3. Differential_abundance_analysis.Rmd           -> relative_abundance_*.csv and diff_<A>_<B>.csv (limma-voom)
# Input: results/amplicon/genus_tables/*.tsv and the metadata file AMPLICON_METADATA (columns Sample and Origin;
# the samples of all datasets, see downstream_input/README.md). Output: AMPLICON_DIR, reports in AMPLICON_DIR/reports.
# Env: atacama-amplicon
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd Rscript pandoc
compgen -G "$AMPLICON_DIR/genus_tables/*.tsv" > /dev/null \
  || die "no genus tables in $AMPLICON_DIR/genus_tables (run scripts/amplicon_process_dataset.sh for each dataset first)"
have_files "$AMPLICON_METADATA" || die "metadata file '$AMPLICON_METADATA' not found (columns Sample and $AMPLICON_GROUP_COLUMN; set AMPLICON_METADATA)"
Rscript -e 'pk <- c("dplyr", "microbiome", "limma", "rmarkdown", "knitr")
            miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
            if (length(miss) > 0) { cat("Missing R package(s):", miss, "\n"); quit(status = 1) }' \
  || die "install the missing R packages or activate the atacama-amplicon environment (envs/amplicon.yml)"

PCT="$(awk -v t="$AMPLICON_PREVALENCE" 'BEGIN { printf "%d", 100 * t + 0.5 }')"
MERGED="$AMPLICON_DIR/merged_bacteria.csv"
FILTERED="$AMPLICON_DIR/merged_bacteria_${PCT}%_prevalence.csv"
REPORTS="$AMPLICON_DIR/reports"

log "Merging the genus tables (collapse: $AMPLICON_COLLAPSE)"
RMDPARAM_input_dir="$AMPLICON_DIR/genus_tables" RMDPARAM_out_dir="$AMPLICON_DIR" RMDPARAM_collapse="$AMPLICON_COLLAPSE" \
  Rscript analysis/render_rmd.R Amplicon_data/merging_multiple_count_abundance_tables.Rmd "$REPORTS"

log "Prevalence filtering (>= ${PCT}% of the samples)"
RMDPARAM_input="$MERGED" RMDPARAM_metadata="$AMPLICON_METADATA" RMDPARAM_threshold="$AMPLICON_PREVALENCE" \
RMDPARAM_out_dir="$AMPLICON_DIR" \
  Rscript analysis/render_rmd.R Amplicon_data/filtering_prevalent_bacteria.Rmd "$REPORTS"

log "Differential abundance: $AMPLICON_GROUP_A vs $AMPLICON_GROUP_B"
RMDPARAM_input="$FILTERED" RMDPARAM_metadata="$AMPLICON_METADATA" RMDPARAM_group_column="$AMPLICON_GROUP_COLUMN" \
RMDPARAM_group_a="$AMPLICON_GROUP_A" RMDPARAM_group_b="$AMPLICON_GROUP_B" RMDPARAM_fdr="$AMPLICON_FDR" \
RMDPARAM_out_dir="$AMPLICON_DIR" \
  Rscript analysis/render_rmd.R Amplicon_data/Differential_abundance_analysis.Rmd "$REPORTS"
log "Done. Tables are in $AMPLICON_DIR, reports in $REPORTS"
