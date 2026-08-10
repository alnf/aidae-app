#!/usr/bin/env Rscript
# Convert study count matrices referenced in config.yaml to sibling .rds files.
#
# Usage (from repository root):
#   Rscript scripts/convert_study_counts.R --study <study_id>
#   Rscript scripts/convert_study_counts.R --study <study_id> --overwrite
#
# Runs check_study_config.R first; aborts on config errors.

args_all <- commandArgs(trailingOnly = TRUE)

parse_cli <- function(args) {
  out <- list(study = NA_character_, overwrite = FALSE, show_help = FALSE, skip_check = FALSE)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      out$show_help <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--overwrite")) {
      out$overwrite <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--skip-check")) {
      out$skip_check <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--study")) {
      if (i >= length(args)) stop("--study requires a value", call. = FALSE)
      out$study <- trimws(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (startsWith(a, "-")) stop("Unknown option: ", a, call. = FALSE)
    stop("Unexpected argument: ", a, call. = FALSE)
  }
  out
}

print_help <- function() {
  message(
    "Usage:\n",
    "  Rscript scripts/convert_study_counts.R --study <study_id> [--overwrite]\n",
    "\n",
    "Wraps convert_study_counts_to_rds() after check_study_config.R.\n"
  )
}

parsed <- parse_cli(args_all)
if (isTRUE(parsed$show_help)) {
  print_help()
  quit(status = 0L)
}
if (is.na(parsed$study) || !nzchar(parsed$study)) {
  print_help()
  message("ERROR: --study <study_id> is required")
  quit(status = 1L)
}

if (!dir.exists("scripts") || !dir.exists("data")) {
  message("ERROR: Run from repository root. getwd()=", getwd())
  quit(status = 1L)
}

sid <- parsed$study

if (!isTRUE(parsed$skip_check)) {
  message("Running check_study_config.R --study ", sid, " ...")
  st <- system2(
    "Rscript",
    c("scripts/check_study_config.R", "--study", sid),
    stdout = "",
    stderr = ""
  )
  if (!identical(as.integer(st), 0L)) {
    message("ERROR: check_study_config.R failed; aborting convert")
    quit(status = 1L)
  }
}

source("scripts/rds_utils.R")
paths <- tryCatch(
  convert_study_counts_to_rds(sid, overwrite = isTRUE(parsed$overwrite)),
  error = function(e) e
)
if (inherits(paths, "error")) {
  message("ERROR: ", conditionMessage(paths))
  quit(status = 1L)
}
message("OK: wrote ", length(paths), " RDS file(s):\n  ", paste(paths, collapse = "\n  "))
quit(status = 0L)
