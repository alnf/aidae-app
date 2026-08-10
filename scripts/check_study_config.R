#!/usr/bin/env Rscript
# Formal consistency check for data/<study_id>/config.yaml and referenced files.
#
# Usage (from repository root):
#   Rscript scripts/check_study_config.R --study <study_id>
#
# Exit 0 if no errors (warnings/infos allowed). Exit 1 on consistency failures.
# Missing derived artifacts (gdegs_file, ora_file) are WARN only.

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
    "  Rscript scripts/check_study_config.R --study <study_id>\n",
    "\n",
    "Checks config.yaml paths, deg_lists, comps.tsv, metadata_filter, thresholds,\n",
    "and schema of referenced DEG/metadata/counts files.\n"
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
cfg_path <- file.path(root, "config.yaml")
if (!file.exists(cfg_path)) {
  message("ERROR: config not found: ", cfg_path)
  quit(status = 1L)
}

cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) e)
if (inherits(cfg, "error")) {
  message("ERROR: cannot parse YAML: ", cfg_path, " — ", conditionMessage(cfg))
  quit(status = 1L)
}

errors <- character()
warnings <- character()
infos <- character()

message("[", sid, "] checking config ", cfg_path)

# --- required keys ---
if (is.null(cfg$name) || !nzchar(as.character(cfg$name)[[1L]])) {
  errors <- c(errors, "missing required key: name")
}

deg_lists <- normalize_deg_lists_cfg(cfg, sid)
if (length(deg_lists) < 1L) {
  errors <- c(errors, "missing deg_lists (or legacy deg_file)")
}

for (key in c("metadata_file", "comparison_file", "gdegs_file", "ora_file")) {
  if (is.null(cfg[[key]]) || !nzchar(as.character(cfg[[key]])[[1L]])) {
    errors <- c(errors, paste0("missing required key: ", key))
  }
}

has_top_counts <- !is.null(cfg$counts_file) && nzchar(as.character(cfg$counts_file)[[1L]])
has_per_list_counts <- FALSE
if (length(deg_lists) > 0L) {
  for (entry in deg_lists) {
    if (!is.null(entry$counts_file) && nzchar(as.character(entry$counts_file)[[1L]])) {
      has_per_list_counts <- TRUE
      break
    }
  }
}
if (!has_top_counts && !has_per_list_counts) {
  errors <- c(errors, "no counts_file at study level or on any deg_lists[] entry")
}
if (!has_top_counts && has_per_list_counts) {
  warnings <- c(
    warnings,
    paste0(
      "no top-level counts_file: Gene tab uses per-list matrices via ",
      "load_gene_tab_study_panels() (grouped mode). ",
      "load_gene_tab_study() returns NULL without a study-level counts_file."
    )
  )
}

# --- resolve paths ---
resolve_study_path <- function(rel) {
  if (is.null(rel) || !nzchar(as.character(rel)[[1L]])) return(NA_character_)
  file.path(root, as.character(rel)[[1L]])
}

meta_path <- resolve_study_path(cfg$metadata_file)
comp_path <- resolve_study_path(cfg$comparison_file)
gdegs_path <- resolve_study_path(cfg$gdegs_file)
ora_path <- resolve_study_path(cfg$ora_file)

if (!is.na(meta_path) && !file.exists(meta_path)) {
  errors <- c(errors, paste0("metadata_file not found: ", meta_path))
}
if (!is.na(comp_path) && !file.exists(comp_path)) {
  errors <- c(errors, paste0("comparison_file not found: ", comp_path))
}
if (!is.na(gdegs_path) && !file.exists(gdegs_path)) {
  warnings <- c(warnings, paste0("derived missing (build later): gdegs_file=", gdegs_path))
}
if (!is.na(ora_path) && !file.exists(ora_path)) {
  warnings <- c(warnings, paste0("derived missing (precompute later): ora_file=", ora_path))
}

# --- thresholds ---
check_thresholds <- function(thr, where) {
  if (is.null(thr)) return(invisible(NULL))
  if (!is.list(thr)) {
    errors <<- c(errors, paste0(where, ": thresholds must be a mapping"))
    return(invisible(NULL))
  }
  bad <- setdiff(names(thr), THRESHOLD_KEYS)
  if (length(bad) > 0L) {
    errors <<- c(
      errors,
      paste0(where, ": unknown thresholds keys [", paste(bad, collapse = ", "),
             "]; allowed: ", paste(THRESHOLD_KEYS, collapse = ", "))
    )
  }
  invisible(NULL)
}
invisible(check_thresholds(cfg$thresholds, "study"))

# --- metadata ---
meta <- NULL
meta_sample_ids <- character()
meta_pheno <- character()
if (!is.na(meta_path) && file.exists(meta_path)) {
  mres <- check_metadata_schema(meta_path)
  errors <- c(errors, mres$errors)
  warnings <- c(warnings, mres$warnings)
  meta <- mres$meta
  meta_sample_ids <- mres$sample_ids
  meta_pheno <- unique(mres$pheno)
}

# --- gene_tab_facet ---
gtf <- cfg$gene_tab_facet
if (!is.null(gtf) && nzchar(as.character(gtf)[[1L]])) {
  gtf <- as.character(gtf)[[1L]]
  if (!is.null(meta) && !(gtf %in% colnames(meta))) {
    errors <- c(errors, paste0("gene_tab_facet \"", gtf, "\" not in metadata columns"))
  }
}

# --- comps ---
comp <- NULL
if (!is.na(comp_path) && file.exists(comp_path)) {
  peek <- read_tsv_header_and_nrow(comp_path)
  if (!isTRUE(peek$ok)) {
    errors <- c(errors, paste0("comparison_file unreadable (", peek$error, "): ", comp_path))
  } else {
    missing_c <- setdiff(COMP_REQUIRED_COLS, peek$colnames)
    if (length(missing_c) > 0L) {
      errors <- c(
        errors,
        paste0("comparison_file missing columns [", paste(missing_c, collapse = ", "), "]")
      )
    }
    comp <- tryCatch(
      read.table(comp_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) e
    )
    if (inherits(comp, "error")) {
      errors <- c(errors, paste0("comparison_file read failed: ", conditionMessage(comp)))
      comp <- NULL
    }
  }
}

# --- deg_lists ---
labels <- character()
if (length(deg_lists) > 0L) {
  for (i in seq_along(deg_lists)) {
    entry <- deg_lists[[i]]
    where <- paste0("deg_lists[", i, "]")
    lab <- if (!is.null(entry$label)) as.character(entry$label)[[1L]] else ""
    drel <- if (!is.null(entry$deg_file)) as.character(entry$deg_file)[[1L]] else ""
    if (!nzchar(lab)) {
      errors <- c(errors, paste0(where, ": missing label"))
      next
    }
    if (lab %in% labels) {
      errors <- c(errors, paste0(where, ": duplicate label \"", lab, "\""))
    }
    labels <- c(labels, lab)
    if (!nzchar(drel)) {
      errors <- c(errors, paste0(where, " label=", lab, ": missing deg_file"))
      next
    }
    dpath <- file.path(root, drel)
    if (!file.exists(dpath)) {
      errors <- c(errors, paste0(where, " label=", lab, ": deg_file not found: ", dpath))
    } else {
      dres <- check_deg_file_schema(dpath, prefix = paste0(where, " "))
      errors <- c(errors, dres$errors)
      warnings <- c(warnings, dres$warnings)
      infos <- c(infos, dres$infos)
    }
    check_thresholds(entry$thresholds, paste0(where, " label=", lab))

    # counts for this list
    crel <- if (!is.null(entry$counts_file) && nzchar(as.character(entry$counts_file)[[1L]])) {
      as.character(entry$counts_file)[[1L]]
    } else if (has_top_counts) {
      as.character(cfg$counts_file)[[1L]]
    } else {
      NA_character_
    }
    if (is.na(crel) || !nzchar(crel)) {
      errors <- c(errors, paste0(where, " label=", lab, ": no counts_file (list or study-level)"))
    } else {
      cpath <- file.path(root, crel)
      if (!file.exists(cpath)) {
        # also accept sibling .rds
        cpath_rds <- sub("\\.(tsv|txt|csv)$", ".rds", cpath, ignore.case = TRUE)
        if (!file.exists(cpath_rds)) {
          errors <- c(errors, paste0(where, " label=", lab, ": counts_file not found: ", cpath))
        } else {
          cpath <- cpath_rds
        }
      }
      if (file.exists(cpath) && length(meta_sample_ids) > 0L) {
        cres <- check_counts_vs_metadata(cpath, meta_sample_ids, prefix = paste0(where, " "))
        errors <- c(errors, cres$errors)
        warnings <- c(warnings, cres$warnings)
      }
    }

    # metadata_filter
    mf <- entry$metadata_filter
    if (!is.null(mf) && is.list(mf) && length(mf) > 0L && !is.null(meta)) {
      for (col in names(mf)) {
        if (!col %in% colnames(meta)) {
          errors <- c(
            errors,
            paste0(where, " label=", lab, ": metadata_filter key \"", col, "\" not in metadata")
          )
          next
        }
        vals <- as.character(unlist(mf[[col]]))
        keep <- meta[[col]] %in% vals
        n_keep <- sum(keep, na.rm = TRUE)
        if (n_keep < 1L) {
          errors <- c(
            errors,
            paste0(where, " label=", lab, ": metadata_filter yields 0 samples (", col, ")")
          )
        } else if (length(meta_sample_ids) > 0L && !is.na(crel) && nzchar(crel)) {
          cpath2 <- file.path(root, crel)
          if (!file.exists(cpath2)) {
            cpath2 <- sub("\\.(tsv|txt|csv)$", ".rds", cpath2, ignore.case = TRUE)
          }
          if (file.exists(cpath2)) {
            cc <- read_counts_sample_ids(cpath2)
            if (isTRUE(cc$ok)) {
              filt_ids <- as.character(meta$SampleNumber[keep])
              in_counts <- intersect(filt_ids, cc$sample_ids)
              miss <- setdiff(filt_ids, cc$sample_ids)
              if (length(in_counts) < 1L) {
                errors <- c(
                  errors,
                  paste0(
                    where, " label=", lab,
                    ": metadata_filter samples have zero overlap with counts colnames"
                  )
                )
              } else if (length(miss) > 0L) {
                warnings <- c(
                  warnings,
                  paste0(
                    where, " label=", lab,
                    ": ", length(miss), "/", length(filt_ids),
                    " filtered metadata samples absent from counts (e.g. ", miss[[1L]],
                    "); app uses intersection"
                  )
                )
              }
            }
          }
        }
      }
    }

    # comps joint match
    if (!is.null(comp) && all(COMP_REQUIRED_COLS %in% colnames(comp))) {
      n_joint <- sum(comp$joint == lab, na.rm = TRUE)
      if (n_joint == 0L) {
        errors <- c(errors, paste0(where, " label=", lab, ": no comps.tsv row with joint == label"))
      } else if (n_joint > 1L) {
        errors <- c(
          errors,
          paste0(where, " label=", lab, ": ", n_joint, " comps.tsv rows with joint == label (need exactly 1)")
        )
      } else {
        row <- comp[comp$joint == lab, , drop = FALSE][1L, ]
        g1 <- as.character(row$group1)
        g2 <- as.character(row$group2)
        if (length(meta_pheno) > 0L) {
          if (!(g1 %in% meta_pheno)) {
            warnings <- c(
              warnings,
              paste0(where, " label=", lab, ": group1 \"", g1, "\" not in metadata$PhenoNames")
            )
          }
          if (!(g2 %in% meta_pheno)) {
            warnings <- c(
              warnings,
              paste0(where, " label=", lab, ": group2 \"", g2, "\" not in metadata$PhenoNames")
            )
          }
        }
      }
    }
  }
}

# orphan joints
if (!is.null(comp) && "joint" %in% colnames(comp) && length(labels) > 0L) {
  orphans <- setdiff(unique(as.character(comp$joint)), labels)
  if (length(orphans) > 0L) {
    warnings <- c(
      warnings,
      paste0("comps.tsv joint values with no deg_lists label: ", paste(orphans, collapse = ", "))
    )
  }
}

# study-level counts if present
if (has_top_counts) {
  cpath <- file.path(root, as.character(cfg$counts_file)[[1L]])
  if (!file.exists(cpath)) {
    cpath_rds <- sub("\\.(tsv|txt|csv)$", ".rds", cpath, ignore.case = TRUE)
    if (!file.exists(cpath_rds)) {
      errors <- c(errors, paste0("study counts_file not found: ", cpath))
    } else {
      cpath <- cpath_rds
    }
  }
  if (file.exists(cpath) && length(meta_sample_ids) > 0L) {
    cres <- check_counts_vs_metadata(cpath, meta_sample_ids, prefix = "study ")
    errors <- c(errors, cres$errors)
    warnings <- c(warnings, cres$warnings)
  }
}

print_check_messages(errors, warnings, infos)

if (length(errors) > 0L) {
  message("\n[", sid, "] FAILED with ", length(errors), " error(s)")
  quit(status = 1L)
}
message(
  "\n[", sid, "] OK (config checks passed",
  if (length(warnings) > 0L) paste0("; ", length(warnings), " warning(s)") else "",
  ")"
)
quit(status = 0L)
