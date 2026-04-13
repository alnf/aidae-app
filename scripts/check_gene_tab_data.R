#!/usr/bin/env Rscript
# Sanity checks for Gene tab inputs: load_gene_tab_study() and load_gene_tab_gdegs()
# for each study listed in config.yaml.
#
# Usage (from repository root):
#   Rscript scripts/check_gene_tab_data.R
#
# Exits with status 1 if any study fails required checks.

if (!file.exists("config.yaml")) {
  message("ERROR: Run from repository root (config.yaml not found). getwd()=", getwd())
  quit(status = 1L)
}

source("scripts/perf_utils.R")
source("scripts/study_data.R")

main <- function() {

  cfg <- yaml::read_yaml("config.yaml")
  studies <- cfg$studies
  if (is.null(studies) || length(studies) == 0L) {
    message("No studies in config.yaml — nothing to check.")
    quit(status = 0L)
  }

  failed <- character(0L)

  for (sid in studies) {
    prefix <- paste0("[", sid, "] ")
    gt <- load_gene_tab_study(sid)
    if (is.null(gt)) {
      message(prefix, "FAIL: load_gene_tab_study returned NULL")
      failed <- c(failed, sid)
      next
    }

    if (ncol(gt$mm) != nrow(gt$metadata)) {
      message(
        prefix, "FAIL: ncol(mm) ", ncol(gt$mm),
        " != nrow(metadata) ", nrow(gt$metadata)
      )
      failed <- c(failed, sid)
      next
    }

    cn <- colnames(gt$mm)
    sn <- as.character(gt$metadata$SampleNumber)
    if (!identical(cn, sn)) {
      message(prefix, "FAIL: colnames(mm) not identical to metadata$SampleNumber (aligned order)")
      failed <- c(failed, sid)
      next
    }

    scfg <- yaml::read_yaml(file.path("data", sid, "config.yaml"))
    gtf <- scfg$gene_tab_facet
    if (!is.null(gtf) && nzchar(as.character(gtf))) {
      gtf <- as.character(gtf)
      if (!gtf %in% colnames(gt$metadata)) {
        message(prefix, "FAIL: gene_tab_facet \"", gtf, "\" not in metadata columns")
        failed <- c(failed, sid)
        next
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
        failed <- c(failed, sid)
        next
      }
      gdegs_n <- nrow(gd)
      if ("study_id" %in% names(gd) && !all(gd$study_id == sid)) {
        message(prefix, "FAIL: gdegs column study_id is not all \"", sid, "\"")
        failed <- c(failed, sid)
        next
      }
      if (!is.null(gtf) && nzchar(as.character(gtf))) {
        gtf <- as.character(gtf)
        if (!gtf %in% colnames(gd)) {
          message(
            prefix, "FAIL: gdegs table missing column \"", gtf,
            "\" (rebuild with build_gene_deg_long after setting gene_tab_facet)"
          )
          failed <- c(failed, sid)
          next
        }
      }
      if (nrow(gd) > 0L && "ens_gene" %in% names(gd)) {
        u <- unique(gd$ens_gene)
        n_try <- min(50L, length(u))
        eg_sample <- u[seq_len(n_try)]
        miss <- eg_sample[!eg_sample %in% rownames(gt$mm)]
        if (length(miss) > 0L) {
          message(
            prefix, "FAIL: ens_gene from gdegs not found in counts rownames (e.g. ",
            miss[[1L]], ")"
          )
          failed <- c(failed, sid)
          next
        }
      }
    }

    msg <- paste0(prefix, "OK: samples=", ncol(gt$mm))
    if (!is.na(gdegs_n)) msg <- paste0(msg, ", gdegs_rows=", gdegs_n)
    message(msg)
  }

  if (length(failed) > 0L) {
    message("\nFAILED: ", paste(unique(failed), collapse = ", "))
    quit(status = 1L)
  }
  message("\nAll checks passed.")
  quit(status = 0L)
}

main()
