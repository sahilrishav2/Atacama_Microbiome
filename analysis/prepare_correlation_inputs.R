#!/usr/bin/env Rscript
# Builds the two tables the correlation scripts need from the outputs of steps 3 and 4:
#   cutibacterium_abundance.csv   one row: relative abundance of one species (default C. acnes) from MetaPhlAn
#   abundant_pathways.csv         the most abundant HUMAnN pathways (unstratified), one row per pathway
# Both have the samples as columns and the sample names cleaned to the plain run accession.
#
# The bacteriophage abundance (TPM) table is not made here: it comes from step 14 or is your own (see README, "Step 15 inputs").

.f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(if (length(.f) > 0) dirname(normalizePath(.f[1])) else "analysis", "utils.R"))

usage <- "Usage: Rscript analysis/prepare_correlation_inputs.R [options]

  --metaphlan FILE         merged MetaPhlAn table (step 3)   [results/03_metaphlan/merged_abundance_table.txt]
  --pathways FILE          merged HUMAnN pathway table in CPM (step 4)
                                                             [results/04_humann/merged/merged_pathabundance_cpm.tsv]
  --species NAME           species as written in MetaPhlAn, without 's__'   [Cutibacterium_acnes]
  --top-n N                number of most abundant pathways to keep         [20]
  --min-prevalence X       keep only pathways detected (> 0) in at least this fraction of samples, 0-1   [0]
  --metaphlan-suffix REGEX text removed from the end of MetaPhlAn sample names   [_metagenome$]
  --outdir DIR             output folder                     [results/15_correlation/input]
  --bacteria-file NAME     output file name for the species  [cutibacterium_abundance.csv]
  --pathways-file NAME     output file name for the pathways [abundant_pathways.csv]
  -h, --help               show this message"

opts <- parse_args(list(metaphlan = "results/03_metaphlan/merged_abundance_table.txt",
                        pathways = "results/04_humann/merged/merged_pathabundance_cpm.tsv",
                        species = "Cutibacterium_acnes", top_n = 20, min_prevalence = 0,
                        metaphlan_suffix = "_metagenome$",
                        outdir = "results/15_correlation/input",
                        bacteria_file = "cutibacterium_abundance.csv",
                        pathways_file = "abundant_pathways.csv"), usage)
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- MetaPhlAn: abundance of one species ----------------------------------------
if (!file.exists(opts$metaphlan)) stop("MetaPhlAn table not found: ", opts$metaphlan, call. = FALSE)
lines <- readLines(opts$metaphlan, warn = FALSE)
header_line <- grep("^#?clade_name", lines)[1]
if (is.na(header_line)) {
  stop("No 'clade_name' header found in ", opts$metaphlan, ". Is this a merged MetaPhlAn table?", call. = FALSE)
}
lines <- lines[header_line:length(lines)]
lines[1] <- sub("^#", "", lines[1])
lines <- lines[!grepl("^#", lines)]
mp <- read.delim(text = lines, check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
clades <- mp[[1]]
mp <- mp[, !(names(mp) %in% c("clade_name", "NCBI_tax_id", "clade_taxid", "additional_species")), drop = FALSE]
mp[] <- lapply(mp, as.numeric)
names(mp) <- sub(opts$metaphlan_suffix, "", names(mp))

# species-level rows end in s__<name>; strain/SGB rows below them end in t__<id> and are not used
last_label <- sub(".*\\|", "", clades)
target <- paste0("s__", opts$species)
hit <- which(last_label == target)
if (length(hit) == 0) {
  present <- unique(last_label[startsWith(last_label, "s__")])
  stop("Species '", target, "' not found in ", opts$metaphlan, ". Species present include: ",
       paste(head(present, 8), collapse = ", "), if (length(present) > 8) ", ..." else "", call. = FALSE)
}
bac <- mp[hit[1], , drop = FALSE]
rownames(bac) <- gsub("_", " ", opts$species)

# ---- HUMAnN: most abundant unstratified pathways ------------------------------------
if (!file.exists(opts$pathways)) stop("HUMAnN pathway table not found: ", opts$pathways, call. = FALSE)
hum <- read.delim(opts$pathways, check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
ids <- hum[[1]]
pw <- hum[-1]
names(pw) <- sub("_(Abundance|Coverage)(-[A-Za-z]+)?$", "", names(pw))
keep <- !grepl("|", ids, fixed = TRUE) & !(ids %in% c("UNMAPPED", "UNINTEGRATED"))   # community-level pathways only
pw <- pw[keep, , drop = FALSE]
pw[] <- lapply(pw, as.numeric)
rownames(pw) <- ids[keep]

prevalence <- rowMeans(pw > 0, na.rm = TRUE)
pw <- pw[prevalence >= opts$min_prevalence, , drop = FALSE]
pw <- pw[order(rowMeans(pw, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
pw <- pw[seq_len(min(opts$top_n, nrow(pw))), , drop = FALSE]

# ---- report and write ----------------------------------------------------------------
shared <- intersect(colnames(bac), colnames(pw))
message("Species: ", rownames(bac), " | pathways kept: ", nrow(pw), " | samples: ",
        ncol(bac), " (MetaPhlAn), ", ncol(pw), " (HUMAnN), ", length(shared), " shared")
if (length(shared) == 0) {
  warning("No sample names are shared between the MetaPhlAn and HUMAnN tables. ",
          "Check --metaphlan-suffix (MetaPhlAn names: ", paste(head(colnames(bac), 2), collapse = ", "),
          "; HUMAnN names: ", paste(head(colnames(pw), 2), collapse = ", "), ").", call. = FALSE)
}

bac_path <- file.path(opts$outdir, opts$bacteria_file)
pw_path <- file.path(opts$outdir, opts$pathways_file)
write.csv(bac, bac_path)
write.csv(pw, pw_path)
message("Wrote ", bac_path)
message("Wrote ", pw_path)
