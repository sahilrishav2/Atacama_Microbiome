#!/usr/bin/env Rscript
# Spearman correlation of one bacterium (default: Cutibacterium acnes) with bacteriophage abundance,
# drawn as a one-row corrplot with red (negative) - white - blue (positive) circles.
#
# Inputs (CSV, first column = feature names, one column per sample, same sample IDs in both files):
#   --bacteria  table with exactly ONE row: the abundance of the bacterium across samples
#   --phage     table with one row per bacteriophage (e.g. TPM)
# Outputs in --outdir:  <prefix>_correlations.csv  (rho, p, BH-adjusted q per phage)  and  <prefix>_corrplot.png

.f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(if (length(.f) > 0) dirname(normalizePath(.f[1])) else "analysis", "utils.R"))

usage <- "Usage: Rscript analysis/correlate_cacnes_phages.R --bacteria FILE.csv --phage FILE.csv [options]

  --bacteria FILE   one-row abundance table of the bacterium (required)
  --phage FILE      bacteriophage abundance table, one row per phage (required)
  --label TEXT      name of the bacterium in the plot    [Cutibacterium acnes]
  --outdir DIR      output folder                        [results/15_correlation]
  --prefix TEXT     output file prefix                   [cacnes_phage]
  --title TEXT      plot title                           [auto]
  --width INCHES    figure width                         [auto, grows with the number of phages]
  --height INCHES   figure height                        [5]
  -h, --help        show this message"

opts <- parse_args(list(bacteria = NA_character_, phage = NA_character_,
                        label = "Cutibacterium acnes", outdir = "results/15_correlation",
                        prefix = "cacnes_phage", title = "auto",
                        width = NA_real_, height = NA_real_), usage)
need_pkgs("corrplot")
suppressPackageStartupMessages(library(corrplot))
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

bac <- read_matrix(opts$bacteria, "bacteria table")
phage <- read_matrix(opts$phage, "phage table")
if (nrow(bac) != 1) {
  stop("The bacteria table must contain exactly one row but has ", nrow(bac), ": ", opts$bacteria, call. = FALSE)
}
rownames(bac) <- opts$label
if (opts$label %in% rownames(phage)) {
  stop("'", opts$label, "' also appears in the phage table; choose another --label.", call. = FALSE)
}

samples <- common_samples(bac, phage, "bacteria table", "phage table")
bac <- drop_constant(bac[, samples, drop = FALSE], "bacterium")
if (nrow(bac) == 0) stop(opts$label, " has the same value in all shared samples; no correlation possible.", call. = FALSE)
phage <- drop_constant(phage[, samples, drop = FALSE], "phage(s)")
if (nrow(phage) == 0) stop("No phage with variable abundance across the shared samples.", call. = FALSE)

res <- spearman_one_vs_many(as.numeric(bac[1, ]), phage)
names(res) <- c("phage", "spearman_rho", "p_value", "n_samples")
res$q_value_BH <- p.adjust(res$p_value, method = "BH")
res <- res[, c("phage", "spearman_rho", "p_value", "q_value_BH", "n_samples")]
table_path <- file.path(opts$outdir, paste0(opts$prefix, "_correlations.csv"))
write.csv(res, table_path, row.names = FALSE)

# one-row matrix for corrplot (the bacterium's correlation with itself is not shown)
bac_phage_cor <- matrix(res$spearman_rho, nrow = 1, dimnames = list(opts$label, res$phage))
title <- if (opts$title == "auto") paste("Correlation of", opts$label, "with phages") else opts$title
w <- if (is.na(opts$width)) max(7, 0.45 * ncol(bac_phage_cor) + 3) else opts$width
h <- if (is.na(opts$height)) 5 else opts$height

png_path <- file.path(opts$outdir, paste0(opts$prefix, "_corrplot.png"))
png(png_path, width = w, height = h, units = "in", res = 300)
corrplot(bac_phage_cor,
         method = "circle",
         col = colorRampPalette(c("red", "white", "blue"))(200),
         tl.col = "black",
         tl.cex = 0.9,
         cl.cex = 0.8,
         mar = c(2, 2, 2, 2),
         title = title)
invisible(dev.off())

message("Wrote ", table_path)
message("Wrote ", png_path)
