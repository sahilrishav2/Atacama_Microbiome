#!/usr/bin/env Rscript
# Renders an R Markdown document with parameters taken from environment variables named RMDPARAM_<parameter>.
# Usage: Rscript analysis/render_rmd.R DOCUMENT.Rmd OUTPUT_DIR [OUTPUT_FILE.html]
# Values "TRUE"/"FALSE" become logical and plain numbers become numeric; everything else stays text.
# The working directory (where the documents look for relative paths) is the current directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: Rscript analysis/render_rmd.R DOCUMENT.Rmd OUTPUT_DIR [OUTPUT_FILE.html]", call. = FALSE)
}
env <- Sys.getenv()
keys <- grep("^RMDPARAM_", names(env), value = TRUE)
params <- as.list(env[keys])
names(params) <- sub("^RMDPARAM_", "", keys)
params <- lapply(params, function(v) {
  if (v %in% c("TRUE", "FALSE")) {
    as.logical(v)
  } else if (grepl("^-?[0-9]+(\\.[0-9]+)?$", v)) {
    as.numeric(v)
  } else {
    v
  }
})
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
rmarkdown::render(args[1], params = params, output_dir = args[2],
                  output_file = if (length(args) == 3) args[3] else NULL,
                  knit_root_dir = getwd(), quiet = TRUE)
