# Helpers shared by the scripts in analysis/ (base R only, plus the packages checked by need_pkgs()).

# ---- command line ------------------------------------------------------------
# spec : named list of defaults.
#          - character default NA_character_  -> option is required
#          - logical default                  -> option is a flag (present = TRUE)
#          - numeric default                  -> value is converted to numeric
# "--some-option value" ends up in opts$some_option.
parse_args <- function(spec, usage, args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("-h", "--help"))) {
    cat(usage, "\n", sep = "")
    quit(save = "no", status = 0)
  }
  opts <- spec
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (!grepl("^--", a)) stop("Unexpected argument '", a, "'. Use --help.", call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", a))
    if (!(key %in% names(spec))) stop("Unknown option '", a, "'. Use --help.", call. = FALSE)
    if (is.logical(spec[[key]])) {
      opts[[key]] <- TRUE
      i <- i + 1
    } else {
      if (i == length(args)) stop("Option '", a, "' needs a value.", call. = FALSE)
      val <- args[i + 1]
      opts[[key]] <- if (is.numeric(spec[[key]])) as.numeric(val) else val
      i <- i + 2
    }
  }
  is_required <- vapply(spec, function(x) is.character(x) && length(x) == 1 && is.na(x), logical(1))
  unset <- names(spec)[is_required][vapply(opts[is_required], is.na, logical(1))]
  if (length(unset) > 0) {
    stop("Missing required option(s): ", paste0("--", gsub("_", "-", unset), collapse = ", "),
         ". Use --help.", call. = FALSE)
  }
  opts
}

need_pkgs <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing R package(s): ", paste(missing_pkgs, collapse = ", "),
         ". Activate the 'atacama-r' conda environment (envs/r-analysis.yml) or install them.",
         call. = FALSE)
  }
}

# ---- input tables --------------------------------------------------------------
# CSV with feature names in the first column and one numeric column per sample.
read_matrix <- function(path, what) {
  if (!file.exists(path)) stop(what, " not found: ", path, call. = FALSE)
  df <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  bad <- names(df)[!vapply(df, is.numeric, logical(1))]
  if (length(bad) > 0) {
    stop(what, " (", path, ") has non-numeric column(s): ", paste(head(bad, 5), collapse = ", "),
         ". Expected: first column = feature names, then one numeric column per sample.", call. = FALSE)
  }
  if (anyDuplicated(colnames(df)) > 0) stop(what, " (", path, ") has duplicated sample names.", call. = FALSE)
  df
}

# Sample names present in both tables; stops if there are too few to correlate.
common_samples <- function(a, b, name_a, name_b, min_n = 3) {
  common <- intersect(colnames(a), colnames(b))
  message(length(common), " samples shared between the ", name_a, " (", ncol(a), ") and the ",
          name_b, " (", ncol(b), ")")
  if (length(common) < min_n) {
    stop("Need at least ", min_n, " shared samples but found ", length(common),
         ". Check that sample IDs match between the tables (", name_a, ": ",
         paste(head(colnames(a), 3), collapse = ", "), " ... ", name_b, ": ",
         paste(head(colnames(b), 3), collapse = ", "), " ...).", call. = FALSE)
  }
  if (length(common) < 10) {
    warning("Only ", length(common), " shared samples: Spearman correlations are unreliable at this size.",
            call. = FALSE)
  }
  common
}

# Rows with zero variance have no defined correlation, so they are removed (and reported).
drop_constant <- function(df, what) {
  if (nrow(df) == 0) return(df)
  sds <- apply(df, 1, sd, na.rm = TRUE)
  constant <- is.na(sds) | sds == 0
  if (any(constant)) {
    message("Dropping ", sum(constant), " constant ", what, ": ",
            paste(head(rownames(df)[constant], 10), collapse = ", "),
            if (sum(constant) > 10) ", ..." else "")
  }
  df[!constant, , drop = FALSE]
}

# ---- statistics ----------------------------------------------------------------
# Spearman correlation of the vector x with every row of the matrix/data frame mat
# (pairwise-complete observations; p-values from the t approximation, as in psych::corr.test).
spearman_one_vs_many <- function(x, mat) {
  rows <- lapply(seq_len(nrow(mat)), function(i) {
    y <- as.numeric(mat[i, ])
    ok <- is.finite(x) & is.finite(y)
    out <- list(feature = rownames(mat)[i], rho = NA_real_, p_value = NA_real_, n = sum(ok))
    if (out$n >= 3) {
      ct <- tryCatch(suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)),
                     error = function(e) NULL)
      if (!is.null(ct)) {
        out$rho <- unname(ct$estimate)
        out$p_value <- ct$p.value
      }
    }
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
