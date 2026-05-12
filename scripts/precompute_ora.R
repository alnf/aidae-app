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
#   Rscript scripts/precompute_ora.R [--cores N] ...
#     Parallel ORA: --cores overrides EXPRS_ORA_WORKERS for this process (use 1 for sequential).
#     At most one parallel layer runs at a time:
#     - ≥2 pathway/ontology files: parallel across those files; DEG lists run sequentially per file.
#     - Exactly 1 pathway file: parallel across DEG lists when there are ≥2 (like the Shiny ORA tab).
#     - Full batch, one pathway file, every study has ≤1 DEG list: parallel across studies.
#   --incremental / EXPRS_ORA_INCREMENTAL: merge-save RDS after each DEG list (resume on crash).
#     With ontology-parallel: one RDS write per DEG list after all ontologies finish that list (no locking).

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
    cores = NA_integer_,
    incremental = NA,
    pathway_mode = "all",
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
    if (identical(a, "--incremental")) {
      out$incremental <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--no-incremental")) {
      out$incremental <- FALSE
      i <- i + 1L
      next
    }
    if (identical(a, "--cores")) {
      if (i >= length(args)) {
        stop("--cores requires a positive integer", call. = FALSE)
      }
      nc <- suppressWarnings(as.integer(trimws(args[[i + 1L]])))
      if (length(nc) != 1L || is.na(nc) || nc < 1L) {
        stop("--cores must be a positive integer", call. = FALSE)
      }
      out$cores <- nc
      i <- i + 2L
      next
    }
    if (identical(a, "--light")) {
      if (!identical(out$pathway_mode, "all")) {
        stop("Use only one of --light or --heavy", call. = FALSE)
      }
      out$pathway_mode <- "light"
      i <- i + 1L
      next
    }
    if (identical(a, "--heavy")) {
      if (!identical(out$pathway_mode, "all")) {
        stop("Use only one of --light or --heavy", call. = FALSE)
      }
      out$pathway_mode <- "heavy"
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
    cores = out$cores,
    incremental = out$incremental,
    pathway_mode = out$pathway_mode,
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
    "  Rscript scripts/precompute_ora.R [--cores N] ...\n",
    "    Optional: parallel ORA (one layer: ontologies vs DEG lists vs studies — see file header).\n",
    "    --cores 1 is sequential. Overrides EXPRS_ORA_WORKERS for this run when N >= 1.\n",
    "  Rscript scripts/precompute_ora.R [--incremental | --no-incremental] ...\n",
    "    --incremental: merge-save RDS for crash-resume (per-DEG after all ontologies when parallel).\n",
    "  Rscript scripts/precompute_ora.R [--light | --heavy] ...\n",
    "    Optional ontology bucket from databases/pathways_list.yaml.\n",
    "    --light: run pathways_light or setdiff(pathways_list, pathways_heavy).\n",
    "    --heavy: run pathways_heavy only.\n",
    "    Default: use EXPRS_ORA_INCREMENTAL=1/true if neither flag is given.\n",
    "Environment:\n",
    "  EXPRS_APP_ROOT         If set, working directory is set to this path before sourcing.\n",
    "  EXPRS_ORA_WORKERS      Default parallel worker count when --cores is omitted (see README).\n",
    "  EXPRS_ORA_INCREMENTAL  If 1/true/yes and no --incremental/--no-incremental flag, enable incremental RDS.\n"
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

resolve_pathway_files <- function(pathway_positional, pathway_mode = "all") {
  mode <- as.character(pathway_mode)[[1L]]
  if (!mode %in% c("all", "light", "heavy")) {
    stop("Internal: unknown pathway_mode=", mode, call. = FALSE)
  }
  pathway_positional <- as.character(pathway_positional)
  pathway_positional <- pathway_positional[nzchar(pathway_positional)]
  if (length(pathway_positional) > 0L && !identical(mode, "all")) {
    stop("Do not combine explicit pathway files with --light/--heavy.", call. = FALSE)
  }
  list_cfg_path <- "databases/pathways_list.yaml"
  list_cfg <- if (file.exists(list_cfg_path)) yaml::read_yaml(list_cfg_path) else list()
  heavy_cfg <- read_pathways_list_value(list_cfg$pathways_heavy)
  if (length(heavy_cfg) > 0L) {
    heavy_cfg <- unique(heavy_cfg)
  }
  all_cfg <- read_pathways_list_value(list_cfg$pathways_list)
  light_cfg <- read_pathways_list_value(list_cfg$pathways_light)
  if (length(light_cfg) < 1L && length(all_cfg) > 0L && length(heavy_cfg) > 0L) {
    light_cfg <- setdiff(all_cfg, heavy_cfg)
  }

  if (identical(mode, "heavy")) {
    if (length(heavy_cfg) < 1L) {
      stop("No pathways_heavy configured in databases/pathways_list.yaml.", call. = FALSE)
    }
    return(heavy_cfg)
  }
  if (identical(mode, "light")) {
    if (length(light_cfg) < 1L) {
      stop("No light ontology bucket available (define pathways_light or pathways_list/pathways_heavy).", call. = FALSE)
    }
    return(light_cfg)
  }
  if (length(pathway_positional) < 1L || !any(nzchar(pathway_positional))) {
    return(list_pathway_txt_files())
  }
  pathway_positional[nzchar(pathway_positional)]
}

#' Order pathway files by decreasing on-disk size.
#'
#' This reduces stragglers in ontology-parallel runs (largest files start first)
#' without paying a long single-threaded pre-parse at startup.
order_pathway_files_largest_first <- function(pathway_files) {
  pfs <- as.character(pathway_files)
  if (length(pfs) < 2L) {
    return(pfs)
  }
  paths <- file.path("databases/pathways", pfs)
  finfo <- file.info(paths)
  bytes <- as.numeric(finfo$size)
  bytes[is.na(bytes)] <- 0
  ord <- order(bytes, decreasing = TRUE, na.last = TRUE)
  pfs_ord <- pfs[ord]
  message(
    "[ORA precompute] ontology scheduling: largest-first by file size",
    " (top 5: ",
    paste(head(pfs_ord, 5L), collapse = ", "),
    ")."
  )
  pfs_ord
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

#' One pathway index for precompute: parse file, run ORA for this study, return chunk or skip.
#' @param ora_workers Passed to `ora_enrichment_long_df_single_study` (use `1L` when multi-ontology `mclapply` is active).
#' @param quiet_skip If TRUE, do not `message()` on missing pathway (parent may log in order).
#' @param out_path Per-study ORA RDS path; required when `incremental` is TRUE.
#' @param incremental If TRUE, run DEG lists sequentially and merge-save after each list (crash-resume).
precompute_pathway_one_index <- function(
    sid,
    slbl,
    pathway_files_to_run,
    pi,
    den,
    deg_filter,
    ora_workers,
    quiet_skip = FALSE,
    out_path = "",
    incremental = FALSE) {
  pf <- pathway_files_to_run[[pi]]
  if (length(pathway_files_to_run) >= 2L) {
    message(
      "[ontology-start] ", sid,
      " | ", pi, "/", length(pathway_files_to_run),
      " | ", pf,
      " | pid=", Sys.getpid()
    )
  }
  t2g <- parse_pathway_file_to_term2gene(pf)
  if (is.null(t2g) || nrow(t2g) == 0L) {
    if (!isTRUE(quiet_skip)) {
      message("[skip pathway] ", sid, " | ", pf, " — missing or empty")
    }
    return(list(skip = TRUE, pf = pf))
  }

  t2g_u <- data.frame(
    term = t2g$term,
    gene = toupper(trimws(as.character(t2g$gene))),
    stringsAsFactors = FALSE
  )

  use_incremental <- isTRUE(incremental) && nzchar(trimws(as.character(out_path)))

  if (use_incremental) {
    tasks <- ora_build_deg_tasks(
      c(sid),
      stats::setNames(list(as.character(slbl)), sid),
      deg_filter
    )
    if (length(tasks) < 1L) {
      return(list(skip = TRUE, pf = pf))
    }
    cmp_done_for_pf <- character(0)
    if (file.exists(out_path)) {
      obj_pf <- tryCatch(readRDS(out_path), error = function(e) NULL)
      if (!is.null(obj_pf) && is.data.frame(obj_pf$long_df) && nrow(obj_pf$long_df) > 0L &&
          all(c("comparison", "pathway_file") %in% colnames(obj_pf$long_df))) {
        ld_pf <- obj_pf$long_df
        mm <- as.character(ld_pf$pathway_file) == as.character(pf)
        if (any(mm)) {
          cmp_done_for_pf <- unique(as.character(ld_pf$comparison[mm]))
        }
      }
    }
    for (ti in seq_along(tasks)) {
      tk <- tasks[[ti]]
      cmp <- ora_deg_entry_comparison_label(tk$entry)
      if (cmp %in% cmp_done_for_pf) {
        message("[resume-skip] ", sid, " — ", cmp, " × ", pf)
        next
      }
      fr <- ora_enrichment_one_deg_entry(
        tk$sid,
        tk$study_label,
        tk$entry,
        t2g_u,
        min_gs_size = 1L,
        max_gs_size = 50000L
      )
      ora_precompute_incremental_merge_save(sid, out_path, fr, cmp, pf, slbl)
    }
    chk <- if (file.exists(out_path)) {
      obj <- tryCatch(readRDS(out_path), error = function(e) NULL)
      if (!is.null(obj) && is.data.frame(obj$long_df) && nrow(obj$long_df) > 0L &&
          "pathway_file" %in% colnames(obj$long_df)) {
        obj$long_df[as.character(obj$long_df$pathway_file) == as.character(pf), , drop = FALSE]
      } else {
        NULL
      }
    } else {
      NULL
    }
    if (is.null(chk) || nrow(chk) < 1L) {
      chk <- ora_long_df_empty_rows(pf)
    }
    return(list(skip = FALSE, pf = pf, chunk = chk))
  }

  chunk <- ora_enrichment_long_df_single_study(
    sid,
    slbl,
    t2g_u,
    min_gs_size = 1L,
    max_gs_size = 50000L,
    progress = NULL,
    den = den,
    pathway_rel_file = pf,
    deg_filter = deg_filter,
    workers = ora_workers
  )

  if (is.null(chunk)) {
    chunk <- ora_long_df_empty_rows(pf)
  }

  list(skip = FALSE, pf = pf, chunk = chunk)
}

#' Run ORA for one study; returns list(long_df, pathway_files_ok) or NULL if skipped.
#'
#' Parallelism (at most one layer uses `mclapply` with `mc.cores > 1`):
#' - **Ontology / pathway-file parallel** when there are ≥2 pathway files: those run in parallel; DEG lists are always sequential within each file (`workers = 1` per job).
#' - **DEG-list parallel** only when there is **exactly one** pathway file and ≥2 DEG lists (same idea as Shiny ORA for one ontology).
#' @param workers Passed to `ora_enrichment_long_df_single_study` when `n_pf == 1L`; for multiple pathway files inner jobs use `workers = 1L`.
#' @param out_path Target per-study ORA RDS (for `--incremental` flush after each DEG list).
#' @param incremental If TRUE, flush RDS after each DEG list (single ontology: each list; multi-ontology: one write per DEG after parallel round).
precompute_one_study_chunks <- function(
    sid,
    slbl,
    pathway_files_to_run,
    deg_filter,
    workers = NULL,
    out_path = "",
    incremental = FALSE) {
  if (length(study_deg_lists(sid)) == 0L) {
    message("[skip] ", sid, " — no deg_lists")
    return(NULL)
  }

  n_pf <- length(pathway_files_to_run)
  n_lists <- if (is.null(deg_filter)) length(study_deg_lists(sid)) else 1L
  den <- max(1L, n_lists * n_pf)

  deg_w <- if (n_pf == 1L) ora_resolve_ora_workers(workers, n_lists) else 1L
  path_w <- if (n_pf >= 2L) ora_resolve_ora_workers(workers, n_pf) else 1L

  incremental_effective <- isTRUE(incremental) && nzchar(trimws(as.character(out_path)))
  if (incremental_effective && deg_w >= 2L) {
    message(
      "[ORA precompute] ", sid,
      ": --incremental: DEG-list parallel disabled (sequential lists + RDS flush)."
    )
  }

  parts <- list()
  pathway_files_ok <- character(0)

  if (path_w >= 2L) {
    if (incremental_effective) {
      message(
        "[ORA precompute] ", sid, ": incremental + ontology-parallel (", path_w, " workers) — ",
        "each DEG list: parallel `enricher` across ", n_pf,
        " ontologies, then **one** RDS write (no file locking)."
      )
      tasks_deg <- ora_build_deg_tasks(
        c(sid),
        stats::setNames(list(as.character(slbl)), sid),
        deg_filter
      )
      pathway_files_ok <- character(0)
      for (ki in seq_along(tasks_deg)) {
        deg_t0 <- unname(as.numeric(Sys.time()))
        tk <- tasks_deg[[ki]]
        cmp <- ora_deg_entry_comparison_label(tk$entry)
        deg_progress <- paste0(ki, "/", length(tasks_deg))
        # One read per DEG round — not one readRDS per ontology (that pegged I/O and hid CPU use).
        pfs_done_cmp <- ora_precompute_resume_pathways_done_for_comparison(out_path, cmp)
        pfs_todo <- as.character(pathway_files_to_run)
        all_resume <- length(pfs_todo) > 0L && all(pfs_todo %in% pfs_done_cmp)
        if (all_resume) {
          message("[resume-skip] ", sid, " — full round ", cmp, " (all ontologies already in RDS)")
          next
        }
        prep <- ora_prepare_deg_gene_sets(tk$sid, tk$entry)
        if (is.null(prep)) {
          message("[skip] ", sid, " — ", cmp, " has no usable DEG genes after filters")
          next
        }
        rnd <- parallel::mclapply(
          seq_len(n_pf),
          function(pi) {
            t0 <- unname(as.numeric(Sys.time()))
            pf <- pathway_files_to_run[[pi]]
            message(
              "[ontology-start] ", sid,
              " | deg ", deg_progress,
              " | ", cmp,
              " | ", pi, "/", n_pf,
              " | ", pf,
              " | pid=", Sys.getpid()
            )
            t2g <- parse_pathway_file_to_term2gene(pf)
            if (is.null(t2g) || nrow(t2g) == 0L) {
              dt <- unname(as.numeric(Sys.time())) - t0
              message(
                "[ontology-done] ", sid,
                " | deg ", deg_progress,
                " | ", cmp,
                " | ", pi, "/", n_pf,
                " | ", pf,
                " | file-empty",
                " | elapsed=", sprintf("%.2f", dt), "s",
                " | pid=", Sys.getpid()
              )
              return(list(kind = "file", pf = pf, fr = NULL, elapsed_sec = dt))
            }
            if (pf %in% pfs_done_cmp) {
              dt <- unname(as.numeric(Sys.time())) - t0
              message(
                "[ontology-done] ", sid,
                " | deg ", deg_progress,
                " | ", cmp,
                " | ", pi, "/", n_pf,
                " | ", pf,
                " | resume-skip",
                " | elapsed=", sprintf("%.2f", dt), "s",
                " | pid=", Sys.getpid()
              )
              return(list(kind = "resume", pf = pf, fr = NULL, elapsed_sec = dt))
            }
            t2g_u <- data.frame(
              term = t2g$term,
              gene = toupper(trimws(as.character(t2g$gene))),
              stringsAsFactors = FALSE
            )
            fr <- ora_enrichment_one_deg_entry(
              tk$sid,
              tk$study_label,
              tk$entry,
              t2g_u,
              min_gs_size = 1L,
              max_gs_size = 50000L,
              prepared = prep
            )
            dt <- unname(as.numeric(Sys.time())) - t0
            nfr <- if (is.null(fr) || !is.data.frame(fr)) 0L else nrow(fr)
            message(
              "[ontology-done] ", sid,
              " | deg ", deg_progress,
              " | ", cmp,
              " | ", pi, "/", n_pf,
              " | ", pf,
              " | rows=", nfr,
              " | elapsed=", sprintf("%.2f", dt), "s",
              " | pid=", Sys.getpid()
            )
            list(kind = "merge", pf = pf, fr = fr, elapsed_sec = dt)
          },
          mc.cores = path_w,
          mc.preschedule = FALSE
        )
        rd <- vapply(
          rnd,
          function(x) {
            if (is.null(x$elapsed_sec) || length(x$elapsed_sec) != 1L || is.na(x$elapsed_sec)) {
              return(NA_real_)
            }
            as.numeric(x$elapsed_sec)
          },
          numeric(1L)
        )
        rd <- rd[!is.na(rd)]
        mr <- ora_precompute_incremental_merge_round(sid, out_path, rnd, cmp, slbl)
        if (is.list(mr) && isTRUE(mr$wrote) && length(mr$pathway_files_merged) > 0L) {
          pathway_files_ok <- unique(c(pathway_files_ok, mr$pathway_files_merged))
        }
        deg_dt <- unname(as.numeric(Sys.time())) - deg_t0
        message(
          "[deg-round] ", sid,
          " | ", cmp, " (", deg_progress, ")",
          " | elapsed=", sprintf("%.2f", deg_dt), "s",
          " | ontology_elapsed_sum=", sprintf("%.2f", sum(rd)), "s",
          " | ontology_elapsed_max=", sprintf("%.2f", if (length(rd) > 0L) max(rd) else 0), "s"
        )
      }
      if (!file.exists(out_path)) {
        message("[skip] ", sid, " — incremental rounds produced no RDS (all pathways skipped?)")
        return(NULL)
      }
      obj_final <- tryCatch(readRDS(out_path), error = function(e) NULL)
      if (is.null(obj_final) || !is.data.frame(obj_final$long_df)) {
        message("[skip] ", sid, " — could not read incremental RDS")
        return(NULL)
      }
      if (length(pathway_files_ok) < 1L) {
        pathway_files_ok <- unique(as.character(obj_final$long_df$pathway_file))
      }
      return(list(
        long_df = obj_final$long_df,
        pathway_files_ok = pathway_files_ok,
        incremental_flushed = TRUE
      ))
    }

    message(
      "[ORA precompute] ", sid, ": ontology-parallel (", path_w, " workers) over ",
      n_pf, " pathway file(s); DEG lists sequential within each ontology."
    )
    pl_raw <- parallel::mclapply(
      seq_len(n_pf),
      function(pi) {
        precompute_pathway_one_index(
          sid, slbl, pathway_files_to_run, pi, den, deg_filter,
          ora_workers = 1L,
          quiet_skip = TRUE,
          out_path = "",
          incremental = FALSE
        )
      },
      mc.cores = path_w,
      mc.preschedule = FALSE
    )
    for (pi in seq_len(n_pf)) {
      r <- pl_raw[[pi]]
      if (isTRUE(r$skip)) {
        message("[skip pathway] ", sid, " | ", r$pf, " — missing or empty")
        next
      }
      pathway_files_ok <- c(pathway_files_ok, r$pf)
      parts[[length(parts) + 1L]] <- r$chunk
    }
  } else {
    if (deg_w >= 2L && !incremental_effective) {
      message(
        "[ORA precompute] ", sid, ": DEG-list-parallel (", deg_w, " workers); ",
        "single ontology (", n_pf, " pathway file)."
      )
    }
    if (incremental_effective) {
      message("[ORA precompute] ", sid, ": incremental RDS → ", out_path)
    }
    for (pi in seq_len(n_pf)) {
      r <- precompute_pathway_one_index(
        sid, slbl, pathway_files_to_run, pi, den, deg_filter,
        ora_workers = if (n_pf == 1L) workers else 1L,
        quiet_skip = FALSE,
        out_path = out_path,
        incremental = incremental_effective
      )
      if (isTRUE(r$skip)) {
        next
      }
      pathway_files_ok <- c(pathway_files_ok, r$pf)
      parts[[length(parts) + 1L]] <- r$chunk
    }
  }

  if (length(parts) < 1L) {
    message("[skip] ", sid, " — no usable pathway files")
    return(NULL)
  }

  long_df <- do.call(rbind, parts)
  list(
    long_df = long_df,
    pathway_files_ok = pathway_files_ok,
    incremental_flushed = isTRUE(incremental_effective)
  )
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

pathway_files_to_run <- resolve_pathway_files(parsed$pathway_positional, parsed$pathway_mode)
if (length(pathway_files_to_run) < 1L) {
  stop("No pathway .txt files to process (check databases/pathways/).", call. = FALSE)
}
pathway_files_to_run <- order_pathway_files_largest_first(pathway_files_to_run)
if (!identical(parsed$pathway_mode, "all")) {
  message("Ontology bucket mode: ", parsed$pathway_mode, " (--", parsed$pathway_mode, ")")
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

precompute_workers <- if (is.na(parsed$cores)) NULL else as.integer(parsed$cores)
if (!is.na(parsed$cores)) {
  message("ORA workers: ", parsed$cores, " (--cores)")
} else {
  ev <- Sys.getenv("EXPRS_ORA_WORKERS", unset = "")
  if (nzchar(trimws(ev))) {
    message("ORA workers: ", trimws(ev), " (EXPRS_ORA_WORKERS)")
  }
}

precompute_incremental <- if (is.na(parsed$incremental)) {
  v <- tolower(trimws(Sys.getenv("EXPRS_ORA_INCREMENTAL", unset = "")))
  v %in% c("1", "true", "yes")
} else {
  isTRUE(parsed$incremental)
}
if (isTRUE(precompute_incremental)) {
  message("ORA incremental RDS: ON (--incremental or EXPRS_ORA_INCREMENTAL)")
}

#' One study: precompute chunks + save RDS (used sequential or inside `mclapply`).
precompute_run_one_study <- function(sid, pathway_files_to_run, deg_filter, workers, incremental = FALSE) {
  out_path <- study_ora_rds_abs_path(sid)
  if (is.na(out_path) || !nzchar(out_path)) {
    message("[skip] ", sid, " — no config")
    return(invisible(NULL))
  }

  slbl <- study_label_from_config(sid)
  if (length(study_deg_lists(sid)) == 0L) {
    message("[skip] ", sid, " — no deg_lists")
    return(invisible(NULL))
  }

  res <- precompute_one_study_chunks(
    sid,
    slbl,
    pathway_files_to_run,
    deg_filter,
    workers = workers,
    out_path = out_path,
    incremental = incremental
  )
  if (is.null(res)) {
    return(invisible(NULL))
  }

  if (isTRUE(res$incremental_flushed)) {
    message(
      "[final] ", sid, ": RDS updated incrementally during run (",
      nrow(res$long_df), " rows in memory view) →\n  ", out_path
    )
    return(invisible(TRUE))
  }

  save_study_ora_rds(sid, res$long_df, res$pathway_files_ok, deg_filter, out_path)
  invisible(TRUE)
}

#' Study-batch parallel only when inner precompute would not use DEG- or pathway-level `mclapply`
#' (avoids nested parallel). Requires exactly one pathway file in this run and ≤1 DEG list per study.
precompute_study_batch_parallel_plan <- function(study_ids, pathway_files_to_run, deg_filter, workers) {
  study_ids <- as.character(study_ids)
  if (length(study_ids) < 2L) {
    return(list(ok = FALSE, w = 1L))
  }
  if (!is.null(deg_filter)) {
    return(list(ok = FALSE, w = 1L))
  }
  if (length(pathway_files_to_run) != 1L) {
    return(list(ok = FALSE, w = 1L))
  }
  for (sid in study_ids) {
    if (length(study_deg_lists(sid)) > 1L) {
      return(list(ok = FALSE, w = 1L))
    }
  }
  w <- ora_resolve_ora_workers(workers, length(study_ids))
  if (w < 2L) {
    return(list(ok = FALSE, w = 1L))
  }
  list(ok = TRUE, w = w)
}

run_batch <- function(study_ids, workers = NULL) {
  study_ids <- as.character(study_ids)
  message("Studies: ", paste(study_ids, collapse = ", "))

  sp <- precompute_study_batch_parallel_plan(study_ids, pathway_files_to_run, deg_filter, workers)
  if (isTRUE(sp$ok)) {
    message(
      "[ORA precompute] study-parallel: ", sp$w, " workers × ", length(study_ids),
      " studies (one pathway file each, ≤1 DEG list per study)."
    )
    parallel::mclapply(
      study_ids,
      function(sid) {
        tryCatch(
          precompute_run_one_study(sid, pathway_files_to_run, deg_filter, workers, incremental = precompute_incremental),
          error = function(e) {
            message("[error] ", sid, ": ", conditionMessage(e))
            NULL
          }
        )
      },
      mc.cores = sp$w,
      mc.preschedule = FALSE
    )
    return(invisible(NULL))
  }

  for (sid in study_ids) {
    precompute_run_one_study(sid, pathway_files_to_run, deg_filter, workers, incremental = precompute_incremental)
  }
  invisible(NULL)
}

if (is.na(parsed$study) || !nzchar(parsed$study)) {
  if (!is.null(deg_filter)) {
    stop("Internal: deg_filter without --study", call. = FALSE)
  }
  run_batch(discover_study_ids(), workers = precompute_workers)
} else {
  sid <- parsed$study
  cfg_path <- file.path("data", sid, "config.yaml")
  if (!file.exists(cfg_path)) {
    stop("Study config not found: ", cfg_path, call. = FALSE)
  }
  message("Study: ", sid, " (", study_label_from_config(sid), ")")
  run_batch(c(sid), workers = precompute_workers)
}

message("Done.")
