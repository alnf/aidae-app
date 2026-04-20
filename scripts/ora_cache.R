# ORA long-form cache: precompute enrichment per study × DEG list × pathway DB, save RDS, filter when plotting.
# Per-study RDS v2: one `long_df` with columns `pathway_file` and `comparison`; v1 legacy (single pathway) still loads.
# Depends: scripts/study_data.R, scripts/pathway_signatures.R, scripts/heatmap_utils.R (filter_heatmap_row_index, deg_list_threshold_defaults).
# Optional: Bioconductor clusterProfiler.

ORA_CACHE_VERSION <- 2L
ORA_CACHE_VERSION_LEGACY <- 1L

ora_gene_ratio_numeric <- function(gr) {
  vapply(seq_along(gr), function(i) {
    parts <- strsplit(trimws(as.character(gr[i])), "/", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) return(NA_real_)
    a <- suppressWarnings(as.numeric(parts[[1L]]))
    b <- suppressWarnings(as.numeric(parts[[2L]]))
    if (is.na(a) || is.na(b) || b <= 0) NA_real_ else a / b
  }, numeric(1L))
}

ora_pick_top_pathway_ids <- function(long_df, n_show) {
  if (is.null(long_df) || nrow(long_df) < 1L || !"ID" %in% colnames(long_df)) {
    return(character(0))
  }
  if (!"p_adj" %in% colnames(long_df)) {
    return(head(unique(as.character(long_df$ID)), n_show))
  }
  uids <- unique(as.character(long_df$ID))
  best <- vapply(uids, function(uid) {
    sub <- long_df[as.character(long_df$ID) == uid, "p_adj", drop = TRUE]
    suppressWarnings(min(as.numeric(sub), na.rm = TRUE))
  }, numeric(1L))
  ord <- order(best, na.last = TRUE)
  head(uids[ord], n_show)
}

ora_pathway_labels_for_ids <- function(long_df, ids) {
  if (length(ids) < 1L) {
    return(list(levels = character(0), id_to_label = character(0)))
  }
  desc_per_id <- vapply(ids, function(id) {
    sub <- long_df[as.character(long_df$ID) == id, , drop = FALSE]
    if (nrow(sub) < 1L) {
      return(NA_character_)
    }
    as.character(sub$Description[[1L]])
  }, character(1L))
  ok <- !is.na(desc_per_id)
  ids <- ids[ok]
  desc_per_id <- desc_per_id[ok]
  if (length(ids) < 1L) {
    return(list(levels = character(0), id_to_label = character(0)))
  }
  dup <- duplicated(desc_per_id)
  lab <- desc_per_id
  lab[dup] <- paste0(lab[dup], " (", ids[dup], ")")
  names(lab) <- ids
  list(levels = rev(lab), id_to_label = lab)
}

#' Default path for cached ORA RDS for a pathway file under `databases/pathways/` (legacy all-in-one cache).
ora_cache_rds_path <- function(pathway_rel_file) {
  if (is.null(pathway_rel_file) || !nzchar(as.character(pathway_rel_file))) {
    return(NA_character_)
  }
  base <- tools::file_path_sans_ext(basename(as.character(pathway_rel_file)))
  file.path("data", "ora_cache", paste0(base, ".rds"))
}

#' Absolute path to per-study ORA RDS from `data/<study_id>/config.yaml` key `ora_file` (default `ora/enrichment.rds`).
#'
#' All ORA-tab cache reads go through this path (`ora_try_load_per_study_caches`); precompute writes the same path.
study_ora_rds_abs_path <- function(study_id) {
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    return(NA_character_)
  }
  cfg <- yaml::read_yaml(cfg_path)
  rel <- if (!is.null(cfg$ora_file) && nzchar(as.character(cfg$ora_file))) {
    as.character(cfg$ora_file)
  } else {
    "ora/enrichment.rds"
  }
  file.path("data", study_id, rel)
}

study_ora_rds_compatible <- function(obj, pathway_rel_file, study_id) {
  if (is.null(obj) || !is.list(obj)) return(FALSE)
  if (!identical(as.character(obj$study_id), as.character(study_id))) return(FALSE)
  if (!is.data.frame(obj$long_df)) return(FALSE)
  pf <- as.character(pathway_rel_file)
  ver <- suppressWarnings(as.integer(obj$version))
  if (length(ver) != 1L || is.na(ver)) return(FALSE)

  if (ver == ORA_CACHE_VERSION_LEGACY) {
    if (is.null(obj$pathway_file)) return(FALSE)
    if (!identical(as.character(obj$pathway_file), pf)) return(FALSE)
    return(TRUE)
  }

  if (ver == ORA_CACHE_VERSION) {
    pfs <- obj$pathway_files
    if (is.null(pfs) || length(pfs) < 1L) {
      if ("pathway_file" %in% colnames(obj$long_df)) {
        pfs <- unique(as.character(obj$long_df$pathway_file))
      } else {
        return(FALSE)
      }
    } else {
      pfs <- as.character(pfs)
    }
    return(pf %in% pfs)
  }

  FALSE
}

#' Try loading per-study RDS files (one `long_df` per study, all DEG lists in `comparison` column).
#' Paths come from each study’s `config.yaml` via `study_ora_rds_abs_path()`. Skips live `enricher` when every
#' study with DEG lists has a compatible RDS for `pathway_rel_file` (version, pathway membership).
#' @return `list(ok_all = logical, long_by_sid = named list)`
ora_try_load_per_study_caches <- function(study_ids, pathway_rel_file) {
  long_by_sid <- list()
  any_need <- FALSE
  for (sid in study_ids) {
    if (length(study_deg_lists(sid)) > 0L) {
      any_need <- TRUE
      break
    }
  }
  if (!any_need) {
    return(list(ok_all = TRUE, long_by_sid = long_by_sid))
  }

  ok_all <- TRUE
  for (sid in study_ids) {
    if (length(study_deg_lists(sid)) == 0L) {
      next
    }
    abs_path <- study_ora_rds_abs_path(sid)
    if (is.na(abs_path) || !nzchar(abs_path) || !file.exists(abs_path)) {
      ok_all <- FALSE
      next
    }
    obj <- tryCatch(readRDS(abs_path), error = function(e) NULL)
    if (is.null(obj) || !study_ora_rds_compatible(obj, pathway_rel_file, sid)) {
      ok_all <- FALSE
      next
    }
    ld <- obj$long_df
    ver <- suppressWarnings(as.integer(obj$version))
    if (length(ver) == 1L && !is.na(ver) && ver == ORA_CACHE_VERSION) {
      if (!"pathway_file" %in% colnames(ld)) {
        ok_all <- FALSE
        next
      }
      ld <- ld[as.character(ld$pathway_file) == as.character(pathway_rel_file), , drop = FALSE]
      ld$pathway_file <- NULL
    }
    long_by_sid[[sid]] <- ld
  }
  list(ok_all = ok_all, long_by_sid = long_by_sid)
}

#' One table per study: rows are pathway × DEG list; `comparison` is the DEG list label.
#'
#' Uses `minGSSize` as given (typically 1) so filters apply later.
#'
#' @param pathway_rel_file If non-NULL, adds column `pathway_file` (for bundled per-study RDS).
#' @param progress Optional `function(amount, detail)`.
#' @param den Denominator for progress increments (number of study×DEG steps).
ora_long_df_empty_rows <- function(pathway_rel_file = NULL) {
  df <- data.frame(
    comparison = character(0),
    ID = character(0),
    Description = character(0),
    gene_ratio = numeric(0),
    Count = integer(0),
    p_adj = numeric(0),
    setSize = integer(0),
    study_label = character(0),
    stringsAsFactors = FALSE
  )
  if (!is.null(pathway_rel_file) && nzchar(as.character(pathway_rel_file))) {
    df$pathway_file <- character(0)
  }
  df
}

ora_deg_entry_comparison_label <- function(entry) {
  if (!is.null(entry$label) && nzchar(as.character(entry$label))) {
    as.character(entry$label)
  } else {
    basename(as.character(entry$deg_file))
  }
}

#' Keep one DEG list entry; `deg_filter` is `list(deg_file = "...")` or `list(label = "...")` (paths as in config).
ora_filter_deg_lists <- function(lists, deg_filter) {
  if (is.null(deg_filter)) {
    return(lists)
  }
  if (!is.null(deg_filter$deg_file) && nzchar(as.character(deg_filter$deg_file))) {
    target <- as.character(deg_filter$deg_file)
    hit <- vapply(lists, function(e) {
      identical(as.character(e$deg_file), target)
    }, logical(1L))
  } else if (!is.null(deg_filter$label) && nzchar(as.character(deg_filter$label))) {
    target <- as.character(deg_filter$label)
    hit <- vapply(lists, function(e) {
      lbl <- if (!is.null(e$label) && nzchar(as.character(e$label))) {
        as.character(e$label)
      } else {
        basename(as.character(e$deg_file))
      }
      identical(lbl, target)
    }, logical(1L))
  } else {
    stop("deg_filter must contain deg_file or label", call. = FALSE)
  }
  n <- sum(hit)
  if (n < 1L) {
    stop("No deg_lists entry matches --deg-file / --deg-label", call. = FALSE)
  }
  if (n > 1L) {
    stop("Multiple deg_lists entries match --deg-file / --deg-label", call. = FALSE)
  }
  lists[hit]
}

#' Replace rows for `comparison` × `pathway_file` in `pathway_files_ok` with `new_df` (partial ORA refresh).
ora_long_df_merge_replace_comparison_pathways <- function(
    existing_df,
    new_df,
    comparison_label,
    pathway_files_ok) {
  cmp <- as.character(comparison_label)
  pf_ok <- unique(as.character(pathway_files_ok))
  base <- existing_df
  if (!is.null(base) && nrow(base) > 0L && "pathway_file" %in% colnames(base) && "comparison" %in% colnames(base)) {
    ex_cmp <- as.character(base$comparison)
    ex_pf <- as.character(base$pathway_file)
    keep <- !(ex_cmp == cmp & ex_pf %in% pf_ok)
    base <- base[keep, , drop = FALSE]
  } else if (is.null(base) && !is.null(new_df) && nrow(new_df) > 0L) {
    base <- new_df[integer(0), , drop = FALSE]
  } else if (is.null(base)) {
    base <- ora_long_df_empty_rows(pathway_rel_file = if (length(pf_ok) > 0L) pf_ok[[1L]] else NULL)
  }
  if (is.null(new_df) || nrow(new_df) < 1L) {
    return(base)
  }
  rbind(base, new_df)
}

ora_enrichment_long_df_single_study <- function(
    sid,
    study_label,
    t2g_u,
    min_gs_size,
    max_gs_size,
    progress,
    den,
    pathway_rel_file = NULL,
    deg_filter = NULL) {
  lists <- study_deg_lists(sid)
  if (length(lists) == 0L) {
    return(NULL)
  }
  lists <- ora_filter_deg_lists(lists, deg_filter)

  long_rows <- list()
  study_lbl <- study_label
  if (is.null(study_lbl) || !nzchar(as.character(study_lbl))) study_lbl <- sid

  for (k in seq_along(lists)) {
    entry <- lists[[k]]
    deg_rel <- entry$deg_file
    if (is.null(deg_rel) || !nzchar(as.character(deg_rel))) next
    cmp_lbl <- if (!is.null(entry$label) && nzchar(as.character(entry$label))) {
      as.character(entry$label)
    } else {
      basename(as.character(deg_rel))
    }

    if (is.function(progress)) {
      try(progress(1 / den, paste(as.character(study_lbl), "—", cmp_lbl)), silent = TRUE)
    }

    thr <- deg_list_threshold_defaults(sid, deg_rel)
    loaded <- load_study_data(sid, deg_rel)
    res <- loaded$res
    mm <- loaded$mm
    if (is.null(res) || is.null(mm) || nrow(res) == 0L) next

    idx <- filter_heatmap_row_index(
      res, mm,
      thr$fdr, thr$base_mean, thr$log2fc, thr$svalue
    )
    if (is.null(idx) || length(idx) == 0L) next

    query <- toupper(trimws(as.character(res$symbol[idx])))
    query <- query[!is.na(query) & nzchar(query)]
    query <- unique(query)

    uni <- toupper(trimws(as.character(res$symbol)))
    uni <- uni[!is.na(uni) & nzchar(uni)]
    uni <- unique(uni)

    if (length(query) < 1L) next

    e_ora <- tryCatch(
      list(ok = TRUE, val = clusterProfiler::enricher(
        gene = query,
        universe = uni,
        TERM2GENE = t2g_u,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size
      )),
      error = function(err) list(ok = FALSE, msg = conditionMessage(err))
    )
    if (!isTRUE(e_ora$ok)) next

    er <- e_ora$val
    df <- as.data.frame(er)
    if (is.null(df) || nrow(df) == 0L) next

    gr <- ora_gene_ratio_numeric(df$GeneRatio)
    p_adj <- if ("p.adjust" %in% colnames(df)) {
      as.numeric(df$p.adjust)
    } else if ("pvalue" %in% colnames(df)) {
      as.numeric(df$pvalue)
    } else {
      rep(NA_real_, nrow(df))
    }

    set_sz <- if ("setSize" %in% colnames(df)) {
      as.integer(df$setSize)
    } else {
      rep(NA_integer_, nrow(df))
    }

    for (ir in seq_len(nrow(df))) {
      long_rows[[length(long_rows) + 1L]] <- data.frame(
        comparison = cmp_lbl,
        ID = as.character(df$ID[ir]),
        Description = as.character(df$Description[ir]),
        gene_ratio = gr[ir],
        Count = as.integer(df$Count[ir]),
        p_adj = p_adj[ir],
        setSize = set_sz[ir],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(long_rows) == 0L) {
    return(NULL)
  }

  long_df <- do.call(rbind, long_rows)
  long_df$study_label <- as.character(study_lbl)
  if (!is.null(pathway_rel_file) && nzchar(as.character(pathway_rel_file))) {
    long_df$pathway_file <- as.character(pathway_rel_file)
  }
  long_df
}

#' Run enricher for every study × DEG list; no Count / top-N filtering (store full enricher rows).
#'
#' @param progress Optional `function(amount, detail)` for Shiny (e.g. `shiny::incProgress`).
#' @return `list(error = NULL, long_by_sid = named list)` or `list(error = character(1), long_by_sid = NULL)`.
ora_run_enrichment_long_by_study <- function(
    study_ids,
    study_labels,
    pathway_rel_file,
    min_gs_size = 1L,
    max_gs_size = 50000L,
    progress = NULL) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    return(list(
      error = paste0("Install Bioconductor package: BiocManager::install(\"clusterProfiler\")"),
      long_by_sid = NULL
    ))
  }

  t2g <- parse_pathway_file_to_term2gene(pathway_rel_file)
  if (is.null(t2g) || nrow(t2g) == 0L) {
    return(list(error = "No pathways in selected file or file missing.", long_by_sid = NULL))
  }

  t2g_u <- data.frame(
    term = t2g$term,
    gene = toupper(trimws(as.character(t2g$gene))),
    stringsAsFactors = FALSE
  )

  n_steps <- 0L
  for (sid2 in study_ids) {
    n_steps <- n_steps + length(study_deg_lists(sid2))
  }
  den <- max(1L, n_steps)

  long_by_sid <- list()

  for (sid in study_ids) {
    slbl <- study_labels[[sid]]
    df <- ora_enrichment_long_df_single_study(
      sid,
      slbl,
      t2g_u,
      min_gs_size,
      max_gs_size,
      progress,
      den,
      pathway_rel_file = NULL
    )
    if (!is.null(df)) {
      long_by_sid[[sid]] <- df
    }
  }

  list(error = NULL, long_by_sid = long_by_sid)
}

#' Filter cached / fresh long tables and build the same payload as `oraTabServer` needs for plotting.
ora_build_plot_payload <- function(
    long_by_sid,
    study_ids,
    study_labels,
    min_ol,
    min_ct,
    min_gene_ratio,
    n_show) {
  out <- stats::setNames(vector("list", length(study_ids)), study_ids)
  lbs <- list()

  for (sid in study_ids) {
    if (length(study_deg_lists(sid)) == 0L) {
      out[[sid]] <- list(plot_df = NULL, msg = "No DEG lists in study config.")
      next
    }

    ld <- long_by_sid[[sid]]
    if (is.null(ld) || !is.data.frame(ld) || nrow(ld) < 1L) {
      out[[sid]] <- list(
        plot_df = NULL,
        msg = "No enrichment results (check DEG lists, thresholds in config, and overlap filters)."
      )
      next
    }

    df <- ld
    if ("Count" %in% colnames(df)) {
      df <- df[df$Count >= min_ct, , drop = FALSE]
    }
    if ("setSize" %in% colnames(df)) {
      ok <- is.na(df$setSize) | df$setSize >= min_ol
      df <- df[ok, , drop = FALSE]
    }
    if ("gene_ratio" %in% colnames(df) && min_gene_ratio > 0) {
      gr <- as.numeric(df$gene_ratio)
      ok <- !is.na(gr) & gr >= min_gene_ratio
      df <- df[ok, , drop = FALSE]
    }

    if (nrow(df) < 1L) {
      out[[sid]] <- list(
        plot_df = NULL,
        msg = "No pathways left after overlap / count / gene ratio filters (adjust sidebar thresholds)."
      )
      next
    }
    lbs[[sid]] <- df
  }

  if (length(lbs) == 0L) {
    return(list(error = NULL, by_study = out, pathway_level_order = NULL))
  }

  combined_long <- do.call(rbind, lbs)
  global_ids <- ora_pick_top_pathway_ids(combined_long, n_show)
  if (length(global_ids) == 0L) {
    global_ids <- unique(as.character(combined_long$ID))
    global_ids <- head(global_ids, n_show)
  }
  pl <- ora_pathway_labels_for_ids(combined_long, global_ids)
  pathway_levels <- pl$levels
  id_to_label <- pl$id_to_label
  pathway_level_order <- pathway_levels

  for (sid in study_ids) {
    if (is.null(lbs[[sid]])) {
      next
    }
    long_df <- lbs[[sid]]
    plot_df <- long_df[as.character(long_df$ID) %in% global_ids, , drop = FALSE]
    n_comp <- length(unique(long_df$comparison))
    slbl <- unique(as.character(long_df$study_label))[[1L]]

    if (nrow(plot_df) == 0L) {
      cmp_u <- unique(as.character(long_df$comparison))
      cmp0 <- if (length(cmp_u)) cmp_u[[1L]] else "—"
      plot_df <- data.frame(
        comparison = cmp0,
        ID = global_ids[[1L]],
        gene_ratio = NA_real_,
        Count = NA_integer_,
        p_adj = NA_real_,
        pathway_rank = NA_real_,
        study_label = slbl,
        stringsAsFactors = FALSE
      )
      if ("setSize" %in% colnames(long_df)) {
        plot_df$setSize <- NA_integer_
      }
      plot_df$Description <- factor(
        unname(id_to_label[as.character(plot_df$ID)]),
        levels = pathway_levels
      )
    } else {
      best_p <- vapply(split(plot_df$p_adj, plot_df$ID), function(z) {
        suppressWarnings(min(as.numeric(z), na.rm = TRUE))
      }, numeric(1L))
      plot_df$pathway_rank <- unname(best_p[as.character(plot_df$ID)])
      plot_df$study_label <- slbl
      lab <- unname(id_to_label[as.character(plot_df$ID)])
      plot_df$Description <- factor(lab, levels = pathway_levels)
    }

    out[[sid]] <- list(plot_df = plot_df, msg = NULL, n_comp = n_comp)
  }

  list(error = NULL, by_study = out, pathway_level_order = pathway_level_order)
}

#' @param obj List with `version`, `pathway_file`, `study_ids`, `long_by_sid` (named list per study_id).
ora_save_cache <- function(obj, rds_path) {
  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, rds_path)
}

ora_read_cache <- function(rds_path) {
  readRDS(rds_path)
}

#' Align cached `long_by_sid` with current `study_ids` (drop extras, NULL for missing).
ora_cache_subset_studies <- function(cached_long_by_sid, study_ids) {
  out <- stats::setNames(vector("list", length(study_ids)), study_ids)
  for (sid in study_ids) {
    out[[sid]] <- cached_long_by_sid[[sid]]
  }
  out
}

cache_compatible_with_app <- function(cached, pathway_rel_file, study_ids) {
  if (is.null(cached) || !is.list(cached)) return(FALSE)
  if (is.null(cached$version)) return(FALSE)
  cv <- suppressWarnings(as.integer(cached$version))
  if (!cv %in% c(ORA_CACHE_VERSION_LEGACY, ORA_CACHE_VERSION)) return(FALSE)
  if (!identical(as.character(cached$pathway_file), as.character(pathway_rel_file))) return(FALSE)
  if (is.null(cached$study_ids) || is.null(cached$long_by_sid)) return(FALSE)
  if (!setequal(as.character(cached$study_ids), as.character(study_ids))) return(FALSE)
  TRUE
}
