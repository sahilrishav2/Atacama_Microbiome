#!/usr/bin/env bash
# Downstream community analysis for amplicon or shotgun data: renders Microbial_community_analysis.Rmd
# (alpha and beta diversity) and Co-occurrence_network.Rmd (Jaccard similarity network) from the tables in
# DOWNSTREAM_DIR (see downstream_input/README.md).
# Output (DOWNSTREAM_OUT_DIR): one HTML report per document, plus network.graphml.
# Independent of the numbered shotgun steps. Env: atacama-downstream
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
# shellcheck source=config/config.sh
source config/config.sh
require_cmd Rscript pandoc
[[ -d "$DOWNSTREAM_DIR" ]] || die "input folder '$DOWNSTREAM_DIR' not found (set DOWNSTREAM_DIR)"

# stop early, with the list of what is missing, instead of failing halfway through a report
Rscript -e 'pk <- c("phyloseq", "microbiome", "ggplot2", "RColorBrewer", "ade4", "factoextra", "dplyr", "vegan", "igraph", "rmarkdown", "knitr")
            miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
            if (length(miss) > 0) { cat("Missing R package(s):", miss, "\n"); quit(status = 1) }' \
  || die "install the missing R packages or activate the atacama-downstream environment (envs/downstream.yml)"

mkdir -p "$DOWNSTREAM_OUT_DIR"
for RMD in Microbial_community_analysis Co-occurrence_network; do
  [[ -s "$RMD.Rmd" ]] || die "$RMD.Rmd not found in the repository root"
  log "Rendering $RMD.Rmd"
  RMDPARAM_data_dir="$DOWNSTREAM_DIR" RMDPARAM_out_dir="$DOWNSTREAM_OUT_DIR" \
    Rscript analysis/render_rmd.R "$RMD.Rmd" "$DOWNSTREAM_OUT_DIR"
done
log "Done. Reports are in $DOWNSTREAM_OUT_DIR"
