#!/usr/bin/env Rscript
# Sanity checks for Gene tab inputs: load_gene_tab_study() / load_gene_tab_study_panels()
# and load_gene_tab_gdegs().
#
# Usage (from repository root):
#   Rscript scripts/check_gene_tab_data.R
#   Rscript scripts/check_gene_tab_data.R --study <study_id>
#
# Without --study: checks every study listed in config.yaml.
# With --study: checks that study only (need not be registered in config.yaml yet).
#
# Exits with status 1 if any study fails required checks.

args_all <- commandArgs(trailingOnly = TRUE)

parse_cli <- function(args) {
  out <- list(study = NA_character_, show_help = FALSE)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      out$show_help <- TRUE
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
    "  Rscript scripts/check_gene_tab_data.R\n",
    "  Rscript scripts/check_gene_tab_data.R --study <study_id>\n"
  )
}

parsed <- parse_cli(args_all)
if (isTRUE(parsed$show_help)) {
  print_help()
  quit(status = 0L)
}

if (!file.exists("config.yaml") && (is.na(parsed$study) || !nzchar(parsed$study))) {
  message("ERROR: Run from repository root (config.yaml not found). getwd()=", getwd())
  quit(status = 1L)
}
if (!dir.exists("scripts") || !dir.exists("data")) {
  message("ERROR: Run from repository root (scripts/ and data/ not found). getwd()=", getwd())
  quit(status = 1L)
}

source("scripts/perf_utils.R")
source("scripts/study_data.R")

check_one_study <- function(sid) {
  prefix <- paste0("[", sid, "] ")
  failed <- FALSE

  panel_set <- load_gene_tab_study_panels(sid)
  if (is.null(panel_set) || length(panel_set$panels) < 1L) {
    message(prefix, "FAIL: load_gene_tab_study_panels returned NULL / no panels")
    return(TRUE)
  }

  mode <- panel_set$mode
  total_samples <- 0L
  for (panel in panel_set$panels) {
    mm <- panel$mm
    metadata <- panel$metadata
    if (is.null(mm) || is.null(metadata)) {
      message(prefix, "FAIL: panel missing mm/metadata (key=", panel$key, ")")
      return(TRUE)
    }
    if (ncol(mm) != nrow(metadata)) {
      message(
        prefix, "FAIL: panel key=", panel$key,
        " ncol(mm) ", ncol(mm), " != nrow(metadata) ", nrow(metadata)
      )
      return(TRUE)
    }
    cn <- colnames(mm)
    sn <- as.character(metadata$SampleNumber)
    if (!identical(cn, sn)) {
      message(
        prefix, "FAIL: panel key=", panel$key,
        " colnames(mm) not identical to metadata$SampleNumber (aligned order)"
      )
      return(TRUE)
    }
    total_samples <- total_samples + ncol(mm)
  }

  scfg <- yaml::read_yaml(file.path("data", sid, "config.yaml"))
  gtf <- scfg$gene_tab_facet
  if (!is.null(gtf) && nzchar(as.character(gtf))) {
    gtf <- as.character(gtf)
    # facet column must exist in each panel's metadata
    for (panel in panel_set$panels) {
      if (!gtf %in% colnames(panel$metadata)) {
        message(
          prefix, "FAIL: gene_tab_facet \"", gtf,
          "\" not in metadata columns (panel key=", panel$key, ")"
        )
        return(TRUE)
      }
    }
  }

  gdegs_n <- NA_integer_
  rel <- scfg$gdegs_file
  if (!is.null(rel) && nzchar(as.character(rel))) {
    gpath <- file.path("data", sid, rel)
    gd <- load_gene_tab_gdegs(sid)
    if (is.null(gd)) {
      message(
        prefix, "FAIL: gdegs_file set but load_gene_tab_gdegs returned NULL",
        " (path exists: ", file.exists(gpath), ")"
      )
      return(TRUE)
    }
    gdegs_n <- nrow(gd)
    if ("study_id" %in% names(gd) && !all(gd$study_id == sid)) {
      message(prefix, "FAIL: gdegs column study_id is not all \"", sid, "\"")
      return(TRUE)
    }
    if (!is.null(gtf) && nzchar(as.character(gtf))) {
      gtf <- as.character(gtf)
      if (!gtf %in% colnames(gd)) {
        message(
          prefix, "FAIL: gdegs table missing column \"", gtf,
          "\" (rebuild with build_gene_deg_long after setting gene_tab_facet)"
        )
        return(TRUE)
      }
    }
    if (nrow(gd) > 0L && "ens_gene" %in% names(gd)) {
      u <- unique(gd$ens_gene)
      n_try <- min(50L, length(u))
      eg_sample <- u[seq_len(n_try)]
      # ens_gene should appear in at least one panel matrix
      found_any <- FALSE
      for (panel in panel_set$panels) {
        if (any(eg_sample %in% rownames(panel$mm))) {
          found_any <- TRUE
          break
        }
      }
      if (!found_any) {
        # Also try all rownames union
        all_ids <- unique(unlist(lapply(panel_set$panels, function(p) rownames(p$mm)), use.names = FALSE))
        miss <- eg_sample[!eg_sample %in% all_ids]
        if (length(miss) == length(eg_sample)) {
          message(
            prefix, "FAIL: ens_gene from gdegs not found in any panel counts rownames (e.g. ",
            eg_sample[[1L]], ")"
          )
          return(TRUE)
        }
      }
    }
  }

  msg <- paste0(
    prefix, "OK: mode=", mode,
    ", panels=", length(panel_set$panels),
    ", samples=", total_samples
  )
  if (!is.na(gdegs_n)) msg <- paste0(msg, ", gdegs_rows=", gdegs_n)
  message(msg)
  FALSE
}

main <- function() {
  studies <- character()
  if (!is.na(parsed$study) && nzchar(parsed$study)) {
    studies <- parsed$study
  } else {
    if (!file.exists("config.yaml")) {
      message("ERROR: config.yaml not found and --study not given. getwd()=", getwd())
      quit(status = 1L)
    }
    cfg <- yaml::read_yaml("config.yaml")
    studies <- cfg$studies
    if (is.null(studies) || length(studies) == 0L) {
      message("No studies in config.yaml — nothing to check.")
      quit(status = 0L)
    }
    studies <- as.character(studies)
  }

  failed <- character(0L)
  for (sid in studies) {
    if (isTRUE(check_one_study(sid))) {
      failed <- c(failed, sid)
    }
  }

  if (length(failed) > 0L) {
    message("\nFAILED: ", paste(unique(failed), collapse = ", "))
    quit(status = 1L)
  }
  message("\nAll checks passed.")
  quit(status = 0L)
}

main()
