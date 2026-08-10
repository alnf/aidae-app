#!/usr/bin/env Rscript
# Build data/<study_id>/<gdegs_file> from deg_lists + comps.tsv via gdf_utils.R.
#
# Usage (from repository root):
#   Rscript scripts/build_gdegs_long.R --study <study_id>
#
# Runs check_study_config.R first; aborts on config errors.
# Exits non-zero if the long table has 0 rows or any deg_lists entry was skipped.

args_all <- commandArgs(trailingOnly = TRUE)

parse_cli <- function(args) {
  out <- list(study = NA_character_, show_help = FALSE, skip_check = FALSE)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      out$show_help <- TRUE
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
    "  Rscript scripts/build_gdegs_long.R --study <study_id>\n",
    "\n",
    "Wraps build_gene_deg_long() + write_gene_deg_long() after check_study_config.R.\n",
    "Fails if any deg_lists entry is skipped or the result has 0 rows.\n"
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
    message("ERROR: check_study_config.R failed; aborting build_gdegs_long")
    quit(status = 1L)
  }
}

source("scripts/gdf_utils.R")

cfg <- read_study_config(sid)
deg_lists <- normalize_deg_lists(sid, cfg)
if (length(deg_lists) < 1L) {
  message("ERROR: no deg_lists / deg_file for study ", sid)
  quit(status = 1L)
}

# Preflight: every entry must be buildable (mirror build_gene_deg_long skip conditions)
comp <- read_comparison_table(sid, cfg)
req_comp <- c("joint", "group1", "group2", "name")
if (!all(req_comp %in% colnames(comp))) {
  message("ERROR: comparison table missing columns: ", paste(setdiff(req_comp, colnames(comp)), collapse = ", "))
  quit(status = 1L)
}

skipped <- character()
for (entry in deg_lists) {
  label <- if (!is.null(entry$label)) as.character(entry$label) else ""
  deg_rel <- if (!is.null(entry$deg_file)) as.character(entry$deg_file) else ""
  if (!nzchar(label) || !nzchar(deg_rel)) {
    skipped <- c(skipped, paste0("(empty label/deg_file)"))
    next
  }
  deg_path <- file.path("data", sid, deg_rel)
  if (!file.exists(deg_path)) {
    skipped <- c(skipped, paste0(label, ": missing file ", deg_path))
    next
  }
  if (!any(comp$joint == label)) {
    skipped <- c(skipped, paste0(label, ": no comps.tsv joint match"))
    next
  }
  res <- tryCatch(
    read.table(deg_path, sep = "\t", header = TRUE, check.names = FALSE, nrows = 1L),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    skipped <- c(skipped, paste0(label, ": unreadable DEG file"))
    next
  }
  req_deg <- c("ens_gene", "symbol", "padj", "log2FoldChange")
  if (!all(req_deg %in% colnames(res))) {
    skipped <- c(skipped, paste0(label, ": missing DEG columns"))
    next
  }
}

if (length(skipped) > 0L) {
  message("ERROR: would skip ", length(skipped), " deg_lists entr(y/ies):")
  for (s in skipped) message("  - ", s)
  quit(status = 1L)
}

warns <- character()
withCallingHandlers(
  {
    df <- build_gene_deg_long(sid)
  },
  warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
if (length(warns) > 0L) {
  message("ERROR: build_gene_deg_long emitted warning(s) (treated as failure):")
  for (w in warns) message("  - ", w)
  quit(status = 1L)
}
if (is.null(df) || nrow(df) < 1L) {
  message("ERROR: build_gene_deg_long returned 0 rows for study ", sid)
  quit(status = 1L)
}

out_path <- tryCatch(write_gene_deg_long(df, sid), error = function(e) e)
if (inherits(out_path, "error")) {
  message("ERROR: ", conditionMessage(out_path))
  quit(status = 1L)
}
message("OK: wrote ", nrow(df), " rows to ", out_path)
quit(status = 0L)
