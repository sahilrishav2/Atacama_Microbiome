#!/usr/bin/env Rscript
# Spearman correlation of every bacteriophage with every (abundant) bacterial pathway,
# drawn with ggcorrplot (squares, correlation printed in each cell).
#
# Inputs (CSV, first column = feature names, one column per sample, same sample IDs in both files):
#   --phage     bacteriophage abundance table, one row per phage (e.g. TPM)
#   --pathways  pathway abundance table, one row per pathway (e.g. from prepare_correlation_inputs.R)
# Outputs in --outdir:
#   <prefix>_rho.csv, <prefix>_pvalue.csv   phage x pathway matrices
#   <prefix>_long.csv                       one row per pair with rho, p and Benjamini-Hochberg q
#   <prefix>_corrplot.png

.f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(if (length(.f) > 0) dirname(normalizePath(.f[1])) else "analysis", "utils.R"))

usage <- "Usage: Rscript analysis/correlate_phages_pathways.R --phage FILE.csv --pathways FILE.csv [options]

  --phage FILE       bacteriophage abundance table, one row per phage (required)
  --pathways FILE    pathway abundance table, one row per pathway (required)
  --outdir DIR       output folder                         [results/15_correlation]
  --prefix TEXT      output file prefix                    [phage_pathway]
  --width PIXELS     figure width                          [8800]
  --height PIXELS    figure height                         [5400]
  --res DPI          figure resolution                     [600]
  --lab-size N       size of the numbers inside the cells  [3]
  --axis-size N      size of the axis labels               [14]
  -h, --help         show this message"

opts <- parse_args(list(phage = NA_character_, pathways = NA_character_,
                        outdir = "results/15_correlation", prefix = "phage_pathway",
                        width = 8800, height = 5400, res = 600,
                        lab_size = 3, axis_size = 14), usage)
need_pkgs(c("psych", "ggcorrplot", "ggplot2"))
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

t1 <- read_matrix(opts$phage, "phage table")
t2 <- read_matrix(opts$pathways, "pathway table")
overlap <- intersect(rownames(t1), rownames(t2))
if (length(overlap) > 0) {
  stop("Feature names occur in both tables (e.g. '", overlap[1], "'); they must be distinct.", call. = FALSE)
}

# the same samples, in the same order, in both tables
samples <- common_samples(t1, t2, "phage table", "pathway table")
t1 <- drop_constant(t1[, samples, drop = FALSE], "phage(s)")
t2 <- drop_constant(t2[, samples, drop = FALSE], "pathway(s)")
if (nrow(t1) == 0 || nrow(t2) == 0) stop("Nothing left to correlate after removing constant features.", call. = FALSE)

# samples x (phages + pathways), then all pairwise Spearman correlations with p-values.
# adjust = "none": multiple-testing correction is done below, on the phage x pathway block only.
combined <- as.data.frame(t(rbind(t1, t2)))
ct <- psych::corr.test(combined, method = "spearman", adjust = "none", ci = FALSE)

# cross-correlation block: rows = phages, columns = pathways
phage_names <- rownames(t1)
pathway_names <- rownames(t2)
cross_cor <- ct$r[phage_names, pathway_names, drop = FALSE]
cross_p <- ct$p[phage_names, pathway_names, drop = FALSE]
cross_q <- matrix(p.adjust(as.vector(cross_p), method = "BH"), nrow = nrow(cross_p), dimnames = dimnames(cross_p))

base <- file.path(opts$outdir, opts$prefix)
write.csv(cross_cor, paste0(base, "_rho.csv"))
write.csv(cross_p, paste0(base, "_pvalue.csv"))
long <- data.frame(phage = rep(rownames(cross_cor), times = ncol(cross_cor)),
                   pathway = rep(colnames(cross_cor), each = nrow(cross_cor)),
                   spearman_rho = as.vector(cross_cor),
                   p_value = as.vector(cross_p),
                   q_value_BH = as.vector(cross_q),
                   stringsAsFactors = FALSE)
long <- long[order(long$p_value), ]
write.csv(long, paste0(base, "_long.csv"), row.names = FALSE)

p <- ggcorrplot::ggcorrplot(cross_cor, method = "square", type = "full",
                            lab = TRUE, lab_size = opts$lab_size) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold", size = opts$axis_size),
                 axis.text.x = ggplot2::element_text(face = "bold", size = opts$axis_size, angle = 45, hjust = 1))

png_path <- paste0(base, "_corrplot.png")
png(png_path, height = opts$height, width = opts$width, res = opts$res)
print(p)
invisible(dev.off())

message("Wrote ", base, "_rho.csv, _pvalue.csv, _long.csv")
message("Wrote ", png_path)
