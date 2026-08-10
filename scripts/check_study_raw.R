#!/usr/bin/env Rscript
# Validate raw study files (DEG / metadata / counts) before or without a complete config.
#
# Usage (from repository root):
#   Rscript scripts/check_study_raw.R --study <study_id>
#   Rscript scripts/check_study_raw.R --study <id> --apply   # apply SUGGEST_RENAME after user approval
#
# Exit 0 if no errors (warnings/infos allowed). Exit 1 on schema/alignment failures.

args_all <- commandArgs(trailingOnly = TRUE)

parse_cli <- function(args) {
  out <- list(study = NA_character_, apply = FALSE, show_help = FALSE)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      out$show_help <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--apply")) {
      out$apply <- TRUE
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
    "  Rscript scripts/check_study_raw.R --study <study_id>\n",
    "  Rscript scripts/check_study_raw.R --study <study_id> --apply\n",
    "\n",
    "Discovers degs/*.tsv (except gdegs_long), exprs/*, metadata/* under data/<study_id>/.\n",
    "Suggests column renames as SUGGEST_RENAME lines; --apply writes them (only after user OK).\n"
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

source("scripts/study_check_utils.R")
require_repo_root()

sid <- parsed$study
root <- study_root(sid)
if (!dir.exists(root)) {
  message("ERROR: study directory not found: ", root)
  quit(status = 1L)
}

disc <- discover_study_files(sid)
message("[", sid, "] checking raw files under ", root)
message(
  "  deg_files=", length(disc$deg_files),
  " counts_files=", length(disc$counts_files),
  " meta_files=", length(disc$meta_files)
)

errors <- character()
warnings <- character()
infos <- character()
all_suggestions <- list()  # path -> list of renames

if (length(disc$deg_files) < 1L) {
  errors <- c(errors, "no DEG .tsv/.txt/.csv files found under degs/ (excluding gdegs_long)")
}

for (dp in disc$deg_files) {
  res <- check_deg_file_schema(dp, prefix = "")
  errors <- c(errors, res$errors)
  warnings <- c(warnings, res$warnings)
  infos <- c(infos, res$infos)
  if (length(res$suggestions) > 0L) {
    # Keep only suggestions that fix missing required cols
    missing <- setdiff(DEG_REQUIRED_COLS, res$colnames)
    keep <- Filter(function(s) s$to %in% missing, res$suggestions)
    if (length(keep) > 0L) all_suggestions[[dp]] <- keep
  }
}

if (length(disc$meta_files) < 1L) {
  errors <- c(errors, "no metadata tables found under metadata/ (excluding comps.*)")
}

meta_sample_lists <- list()
for (mp in disc$meta_files) {
  res <- check_metadata_schema(mp)
  errors <- c(errors, res$errors)
  warnings <- c(warnings, res$warnings)
  infos <- c(infos, res$infos)
  meta_sample_lists[[mp]] <- res$sample_ids
}

if (length(disc$counts_files) < 1L) {
  errors <- c(errors, "no counts files found under exprs/")
} else if (length(meta_sample_lists) > 0L) {
  # Prefer pairing: for each counts file, check against union of all metadata SampleNumbers,
  # and also each metadata file individually if only one meta file.
  all_meta_ids <- unique(unlist(meta_sample_lists, use.names = FALSE))
  # Deduplicate counts: prefer .tsv over .rds for same basename stem
  counts_to_check <- disc$counts_files
  stems <- sub("\\.(tsv|txt|csv|rds|eds)$", "", basename(counts_to_check), ignore.case = TRUE)
  keep <- logical(length(counts_to_check))
  for (st in unique(stems)) {
    idx <- which(stems == st)
    tabs <- idx[grepl("\\.(tsv|txt|csv)$", counts_to_check[idx], ignore.case = TRUE)]
    if (length(tabs) > 0L) {
      keep[tabs[[1L]]] <- TRUE
    } else {
      keep[idx[[1L]]] <- TRUE
    }
  }
  counts_to_check <- counts_to_check[keep]

  for (cp in counts_to_check) {
    res <- check_counts_vs_metadata(cp, all_meta_ids)
    errors <- c(errors, res$errors)
    warnings <- c(warnings, res$warnings)
    infos <- c(infos, res$infos)
  }
}

print_check_messages(errors, warnings, infos)

if (isTRUE(parsed$apply)) {
  if (length(all_suggestions) < 1L) {
    message("No rename suggestions to apply.")
  } else {
    message("--apply: writing column renames (ensure the user already approved these SUGGEST_RENAME lines)")
    for (path in names(all_suggestions)) {
      ren <- all_suggestions[[path]]
      message("  renaming in ", path, ": ",
              paste(vapply(ren, function(s) paste0(s$from, "->", s$to), character(1L)), collapse = ", "))
      apply_column_renames(path, ren)
    }
    message("Re-run without --apply to verify.")
  }
}

if (length(errors) > 0L) {
  message("\n[", sid, "] FAILED with ", length(errors), " error(s)")
  quit(status = 1L)
}
message("\n[", sid, "] OK (raw checks passed",
        if (length(warnings) > 0L) paste0("; ", length(warnings), " warning(s)") else "",
        ")")
quit(status = 0L)
