#!/usr/bin/env Rscript
# Precompute ORA into data/<study_id>/<ora_file> (per-study RDS v2, bundled pathway_file column).
#
# Usage (from project root):
#   Rscript scripts/precompute_ora.R
#     All studies, all *.txt under databases/pathways/.
#   Rscript scripts/precompute_ora.R [pathway.txt ...]
#     All studies, only listed pathway files (names under databases/pathways/).
#   Rscript scripts/precompute_ora.R --study <study_id> [pathway.txt ...]
#     One study; optional pathway file names limit work like above.
#   Rscript scripts/precompute_ora.R --study <id> --deg-file <deg_file> [pathway.txt ...]
#   Rscript scripts/precompute_ora.R --study <id> --deg-label <label> [pathway.txt ...]
#     One DEG list from one study (path as in config, or comparison label). Merges into
#     existing enrichment.rds when present (replaces rows for that comparison × pathway files run).

args_all <- commandArgs(trailingOnly = TRUE)

normalize_pathway_arg <- function(s) {
  pathway_file <- s
  if (grepl("^databases/pathways/", s)) {
    pathway_file <- sub("^databases/pathways/", "", s)
  }
  pathway_file
}

parse_precompute_cli <- function(args) {
  out <- list(
    study = NA_character_,
    deg_file = "",
    deg_label = "",
    pathway_positional = character(),
    show_help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      out$show_help <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--study")) {
      if (i >= length(args)) {
        stop("--study requires a value", call. = FALSE)
      }
      out$study <- trimws(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (identical(a, "--deg-file")) {
      if (i >= length(args)) {
        stop("--deg-file requires a value", call. = FALSE)
      }
      out$deg_file <- trimws(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (identical(a, "--deg-label")) {
      if (i >= length(args)) {
        stop("--deg-label requires a value", call. = FALSE)
      }
      out$deg_label <- trimws(args[[i + 1L]])
      i <- i + 2L
      next
    }
    if (startsWith(a, "-")) {
      stop("Unknown option: ", a, call. = FALSE)
    }
    out$pathway_positional <- c(out$pathway_positional, a)
    i <- i + 1L
  }

  has_df <- nzchar(out$deg_file)
  has_lb <- nzchar(out$deg_label)
  if (has_df && has_lb) {
    stop("Use only one of --deg-file or --deg-label", call. = FALSE)
  }
  if ((has_df || has_lb) && (is.na(out$study) || !nzchar(out$study))) {
    stop("--deg-file / --deg-label require --study", call. = FALSE)
  }

  deg_filter <- NULL
  if (has_df) {
    deg_filter <- list(deg_file = out$deg_file)
  } else if (has_lb) {
    deg_filter <- list(label = out$deg_label)
  }

  list(
    study = out$study,
    deg_filter = deg_filter,
    pathway_positional = vapply(out$pathway_positional, normalize_pathway_arg, character(1L)),
    show_help = out$show_help
  )
}

print_usage <- function() {
  message(paste0(
    "Usage:\n",
    "  Rscript scripts/precompute_ora.R\n",
    "    Precompute ORA for all studies (all pathway *.txt under databases/pathways/).\n",
    "  Rscript scripts/precompute_ora.R [pathway_filename ...]\n",
    "    All studies; only the listed pathway files (under databases/pathways/). Multiple allowed.\n",
    "  Rscript scripts/precompute_ora.R --study <study_id> [pathway_filename ...]\n",
    "    One study; optional pathway file list limits work.\n",
    "  Rscript scripts/precompute_ora.R --study <id> --deg-file <deg_file_in_config> [pathway ...]\n",
    "  Rscript scripts/precompute_ora.R --study <id> --deg-label <comparison_label> [pathway ...]\n",
    "    One DEG list only; merges into existing per-study RDS (v2) when present.\n",
    "Environment:\n",
    "  EXPRS_APP_ROOT   If set, working directory is set to this path before sourcing.\n"
  ))
}

parsed <- parse_precompute_cli(args_all)

if (isTRUE(parsed$show_help)) {
  print_usage()
  quit(status = 0L)
}

root <- Sys.getenv("EXPRS_APP_ROOT", unset = "")
if (nzchar(root)) {
  setwd(root)
}

source("scripts/heatmap_utils.R")
source("scripts/result_table_indices.R")
source("scripts/study_data.R")
source("scripts/perf_utils.R")
source("scripts/pathway_signatures.R")
source("scripts/ora_cache.R")

resolve_pathway_files <- function(pathway_positional) {
  if (length(pathway_positional) < 1L || !any(nzchar(pathway_positional))) {
    return(list_pathway_txt_files())
  }
  pathway_positional[nzchar(pathway_positional)]
}

discover_study_ids <- function() {
  main_config <- if (file.exists("config.yaml")) yaml::read_yaml("config.yaml") else list(studies = character(0))
  if (!is.null(main_config$studies) && length(main_config$studies) > 0L) {
    return(as.character(main_config$studies))
  }
  data_dirs <- list.dirs("data", full.names = TRUE, recursive = FALSE)
  keep <- vapply(data_dirs, function(d) file.exists(file.path(d, "config.yaml")), logical(1L))
  basename(data_dirs[keep])
}

study_label_from_config <- function(sid) {
  cfg_path <- file.path("data", sid, "config.yaml")
  if (!file.exists(cfg_path)) {
    return(sid)
  }
  sc <- yaml::read_yaml(cfg_path)
  if (!is.null(sc$name)) as.character(sc$name) else sid
}

#' Run ORA for one study; returns list(long_df, pathway_files_ok) or NULL if skipped.
precompute_one_study_chunks <- function(sid, slbl, pathway_files_to_run, deg_filter) {
  if (length(study_deg_lists(sid)) == 0L) {
    message("[skip] ", sid, " — no deg_lists")
    return(NULL)
  }

  n_pf <- length(pathway_files_to_run)
  n_lists <- if (is.null(deg_filter)) length(study_deg_lists(sid)) else 1L
  den <- max(1L, n_lists * n_pf)

  parts <- list()
  pathway_files_ok <- character(0)

  for (pi in seq_len(n_pf)) {
    pf <- pathway_files_to_run[[pi]]
    t2g <- parse_pathway_file_to_term2gene(pf)
    if (is.null(t2g) || nrow(t2g) == 0L) {
      message("[skip pathway] ", sid, " | ", pf, " — missing or empty")
      next
    }

    t2g_u <- data.frame(
      term = t2g$term,
      gene = toupper(trimws(as.character(t2g$gene))),
      stringsAsFactors = FALSE
    )

    chunk <- ora_enrichment_long_df_single_study(
      sid,
      slbl,
      t2g_u,
      min_gs_size = 1L,
      max_gs_size = 50000L,
      progress = NULL,
      den = den,
      pathway_rel_file = pf,
      deg_filter = deg_filter
    )

    if (is.null(chunk)) {
      chunk <- ora_long_df_empty_rows(pf)
    }

    pathway_files_ok <- c(pathway_files_ok, pf)
    parts[[length(parts) + 1L]] <- chunk
  }

  if (length(parts) < 1L) {
    message("[skip] ", sid, " — no usable pathway files")
    return(NULL)
  }

  long_df <- do.call(rbind, parts)
  list(long_df = long_df, pathway_files_ok = pathway_files_ok)
}

save_study_ora_rds <- function(
    sid,
    long_df,
    pathway_files_ok,
    deg_filter,
    out_path) {
  merge_partial <- !is.null(deg_filter)

  if (merge_partial) {
    lists_f <- ora_filter_deg_lists(study_deg_lists(sid), deg_filter)
    cmp_lbl <- ora_deg_entry_comparison_label(lists_f[[1L]])

    existing <- if (file.exists(out_path)) {
      tryCatch(readRDS(out_path), error = function(e) NULL)
    } else {
      NULL
    }

    if (!is.null(existing)) {
      ver <- suppressWarnings(as.integer(existing$version))
      if (length(ver) != 1L || is.na(ver) || ver != ORA_CACHE_VERSION) {
        stop(
          "Existing ", out_path, " is not ORA cache version ", ORA_CACHE_VERSION,
          "; merge requires v2. Remove the file or run a full --study precompute.",
          call. = FALSE
        )
      }
      long_df <- ora_long_df_merge_replace_comparison_pathways(
        existing$long_df,
        long_df,
        cmp_lbl,
        pathway_files_ok
      )
      pathway_files_union <- unique(c(
        as.character(existing$pathway_files),
        pathway_files_ok
      ))
    } else {
      message(
        "[note] No existing RDS at ", out_path,
        "; saving rows for this DEG list only (run full --study ", sid, " to fill others)."
      )
      pathway_files_union <- pathway_files_ok
    }
  } else {
    pathway_files_union <- pathway_files_ok
  }

  obj <- list(
    version = ORA_CACHE_VERSION,
    study_id = sid,
    pathway_files = pathway_files_union,
    created = Sys.time(),
    long_df = long_df
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, out_path)
  message(
    "[saved] ", sid, " -> ", out_path,
    " (", nrow(long_df), " rows; ", length(pathway_files_union), " pathway files in index)"
  )
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop("Install clusterProfiler: BiocManager::install(\"clusterProfiler\")", call. = FALSE)
}

pathway_files_to_run <- resolve_pathway_files(parsed$pathway_positional)
if (length(pathway_files_to_run) < 1L) {
  stop("No pathway .txt files to process (check databases/pathways/).", call. = FALSE)
}

message(
  "Pathway file(s): ",
  length(pathway_files_to_run),
  if (length(pathway_files_to_run) <= 5L) {
    paste0(" — ", paste(pathway_files_to_run, collapse = ", "))
  } else {
    paste0(" (first 5: ", paste(head(pathway_files_to_run, 5L), collapse = ", "), ", …)")
  }
)

deg_filter <- parsed$deg_filter

if (!is.null(deg_filter)) {
  message(
    "Single DEG list: ",
    if (!is.null(deg_filter$deg_file)) paste0("deg_file=", deg_filter$deg_file) else paste0("label=", deg_filter$label)
  )
}

run_batch <- function(study_ids) {
  study_ids <- as.character(study_ids)
  message("Studies: ", paste(study_ids, collapse = ", "))

  for (sid in study_ids) {
    out_path <- study_ora_rds_abs_path(sid)
    if (is.na(out_path) || !nzchar(out_path)) {
      message("[skip] ", sid, " — no config")
      next
    }

    slbl <- study_label_from_config(sid)
    if (length(study_deg_lists(sid)) == 0L) {
      message("[skip] ", sid, " — no deg_lists")
      next
    }

    res <- precompute_one_study_chunks(sid, slbl, pathway_files_to_run, deg_filter)
    if (is.null(res)) {
      next
    }

    save_study_ora_rds(sid, res$long_df, res$pathway_files_ok, deg_filter, out_path)
  }
}

if (is.na(parsed$study) || !nzchar(parsed$study)) {
  if (!is.null(deg_filter)) {
    stop("Internal: deg_filter without --study", call. = FALSE)
  }
  run_batch(discover_study_ids())
} else {
  sid <- parsed$study
  cfg_path <- file.path("data", sid, "config.yaml")
  if (!file.exists(cfg_path)) {
    stop("Study config not found: ", cfg_path, call. = FALSE)
  }
  message("Study: ", sid, " (", study_label_from_config(sid), ")")
  run_batch(c(sid))
}

message("Done.")
