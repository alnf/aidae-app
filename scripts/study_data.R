 # Caches to avoid repeated disk I/O within a session
.study_data_cache <- new.env(parent = emptyenv())
.metadata_cache <- new.env(parent = emptyenv())
.comparison_cache <- new.env(parent = emptyenv())
# Gene tab: full-sample counts + metadata; global DEGs long table (separate from heatmap cache keys)
.gene_tab_cache <- new.env(parent = emptyenv())

study_deg_lists_from_cfg <- function(cfg, study_id) {
  if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    return(cfg$deg_lists)
  }
  if (!is.null(cfg$deg_file)) {
    return(list(list(label = if (!is.null(cfg$name)) cfg$name else study_id, deg_file = cfg$deg_file)))
  }
  list()
}

resolve_counts_file_for_deg <- function(cfg, study_id, deg_file = NULL, label = NULL) {
  lists <- study_deg_lists_from_cfg(cfg, study_id)
  idx <- NA_integer_
  if (!is.null(deg_file) && nzchar(as.character(deg_file))) {
    idx <- match(as.character(deg_file), vapply(lists, function(x) x$deg_file, character(1L)))
  } else if (!is.null(label) && nzchar(as.character(label))) {
    idx <- match(as.character(label), vapply(lists, function(x) x$label, character(1L)))
  }
  if (!is.na(idx)) {
    cf <- lists[[idx]]$counts_file
    if (!is.null(cf) && nzchar(as.character(cf))) return(as.character(cf))
  }
  cf <- cfg$counts_file
  if (!is.null(cf) && nzchar(as.character(cf))) return(as.character(cf))
  NULL
}

study_deg_lists <- function(study_id) {
  if (is.null(study_id) || study_id == "") return(list())
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(list())
  cfg <- yaml::read_yaml(cfg_path)
  study_deg_lists_from_cfg(cfg, study_id)
}

#' Catalog of studies and DEG lists for the Config tab (available scope).
build_study_deg_catalog <- function(study_ids, study_labels) {
  out <- list()
  for (sid in as.character(study_ids)) {
    lists <- study_deg_lists(sid)
    if (length(lists) < 1L) next
    slbl <- study_labels[[sid]]
    if (is.null(slbl) || !nzchar(as.character(slbl))) slbl <- sid
    out[[length(out) + 1L]] <- list(
      study_id = sid,
      study_label = as.character(slbl),
      deg_entries = lists
    )
  }
  out
}

#' All studies and deg_file paths enabled (session default).
default_visibility_state <- function(catalog) {
  study_ids <- vapply(catalog, function(x) x$study_id, character(1L))
  deg_by_study <- stats::setNames(vector("list", length(study_ids)), study_ids)
  for (item in catalog) {
    sid <- item$study_id
    deg_by_study[[sid]] <- vapply(item$deg_entries, function(e) as.character(e$deg_file), character(1L))
  }
  list(study_ids = study_ids, deg_by_study = deg_by_study)
}

#' Subset deg list entries to visible deg_file paths.
#' When `visible_deg_files` is NULL, returns `lists` unchanged (no filter).
filter_deg_lists_visible <- function(lists, visible_deg_files) {
  if (is.null(lists) || length(lists) < 1L) return(list())
  if (is.null(visible_deg_files)) return(lists)
  if (length(visible_deg_files) < 1L) return(list())
  vis <- as.character(visible_deg_files)
  keep <- vapply(lists, function(e) {
    df <- as.character(if (is.null(e$deg_file)) "" else e$deg_file)
    nzchar(df) && df %in% vis
  }, logical(1L))
  lists[keep]
}

#' Parse Config tab checkbox inputs into visible study_ids and deg_by_study.
parse_visibility_from_inputs <- function(catalog, input, ns_prefix = "") {
  study_ids <- character(0)
  deg_by_study <- list()
  for (item in catalog) {
    sid <- item$study_id
    study_in <- input[[paste0(ns_prefix, "study_", sid)]]
    if (!isTRUE(study_in)) next
    all_deg <- vapply(item$deg_entries, function(e) as.character(e$deg_file), character(1L))
    deg_in <- input[[paste0(ns_prefix, "deg_", sid)]]
    if (is.null(deg_in)) deg_in <- character(0)
    deg_in <- intersect(as.character(deg_in), all_deg)
    if (length(deg_in) < 1L) next
    study_ids <- c(study_ids, sid)
    deg_by_study[[sid]] <- deg_in
  }
  list(study_ids = study_ids, deg_by_study = deg_by_study)
}

# ColorBrewer qualitative "Set3" (n <= 12); recycle in order if more studies.
.study_palette_set3 <- function(n) {
  n <- max(1L, as.integer(n))
  base <- c(
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
    "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
  )
  if (n <= length(base)) return(base[seq_len(n)])
  rep(base, length.out = n)
}

# Named vector: display study label -> hex for UpSet (Set3 by main-config study order).
study_color_map_for_labels <- function(study_ids, study_labels) {
  if (is.null(study_ids) || length(study_ids) < 1L) return(character(0))
  pal <- .study_palette_set3(length(study_ids))
  names(pal) <- study_ids
  out <- character(0)
  for (sid in study_ids) {
    lab <- study_labels[[sid]]
    if (is.null(lab) || !nzchar(trimws(as.character(lab)))) lab <- as.character(sid)
    lab <- trimws(as.character(lab))
    if (lab %in% names(out)) next
    out[[lab]] <- pal[[sid]]
  }
  out
}

samples_for_comparison <- function(study_id, label) {
  message("samples_for_comparison entered: study_id=", study_id, " label=", label)
  if (is.null(study_id) || study_id == "" || is.null(label) || label == "") {
    message("samples_for_comparison: return NULL (null/empty study_id or label)")
    return(NULL)
  }
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    message("samples_for_comparison: return NULL (config file not found)")
    return(NULL)
  }
  cfg <- yaml::read_yaml(cfg_path)
  if (is.null(cfg$comparison_file) || is.null(cfg$metadata_file)) {
    message("samples_for_comparison: return NULL (comparison_file or metadata_file missing in config)")
    return(NULL)
  }
  comp_path <- file.path("data", study_id, cfg$comparison_file)
  meta_path <- file.path("data", study_id, cfg$metadata_file)
  if (!file.exists(comp_path) || !file.exists(meta_path)) {
    message("samples_for_comparison: return NULL (file not found: comp=", file.exists(comp_path), " meta=", file.exists(meta_path), " meta_path=", meta_path, ")")
    return(NULL)
  }
  comp_key <- paste(study_id, cfg$comparison_file, sep = "|")
  if (exists(comp_key, envir = .comparison_cache, inherits = FALSE)) {
    comp <- get(comp_key, envir = .comparison_cache, inherits = FALSE)
  } else {
    comp <- read.table(comp_path, sep = "\t", header = TRUE, check.names = FALSE)
    assign(comp_key, comp, envir = .comparison_cache)
  }
  if (!all(c("joint", "group1", "group2") %in% colnames(comp))) {
    message("samples_for_comparison: return NULL (comparison file missing joint/group1/group2)")
    return(NULL)
  }
  idx <- which(comp$joint == label)
  if (length(idx) == 0L) {
    message("samples_for_comparison: return NULL (no row with joint == label)")
    return(NULL)
  }
  group1 <- comp$group1[idx[1L]]
  group2 <- comp$group2[idx[1L]]
  meta_key <- paste(study_id, cfg$metadata_file, sep = "|")
  if (exists(meta_key, envir = .metadata_cache, inherits = FALSE)) {
    meta <- get(meta_key, envir = .metadata_cache, inherits = FALSE)
  } else {
    meta <- read.table(meta_path, sep = "\t", header = TRUE, check.names = FALSE)
    assign(meta_key, meta, envir = .metadata_cache)
  }

  if (!all(c("PhenoNames", "SampleNumber") %in% colnames(meta))) {
    message("samples_for_comparison: return NULL (metadata missing PhenoNames or SampleNumber)")
    return(NULL)
  }
  keep <- meta$PhenoNames %in% c(group1, group2)
  deg_entry <- if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    idx <- match(label, vapply(cfg$deg_lists, function(x) x$label, character(1L)))
    if (!is.na(idx)) cfg$deg_lists[[idx]] else NULL
  } else if (!is.null(cfg$deg_file) && (if (!is.null(cfg$name)) cfg$name else study_id) == label) {
    list(metadata_filter = NULL)
  } else NULL
  mf <- if (!is.null(deg_entry) && !is.null(deg_entry$metadata_filter)) deg_entry$metadata_filter else NULL
  if (!is.null(mf) && is.list(mf) && length(mf) > 0L) {
    for (col in names(mf)) {
      if (col %in% colnames(meta)) {
        vals <- as.character(unlist(mf[[col]]))
        keep <- keep & (meta[[col]] %in% vals)
      }
    }
  }
  out <- meta$SampleNumber[keep]
  message("Samples for heatmap (SampleNumber): ", paste(out, collapse = ", "))
  out
}

load_study_data <- function(study_id, deg_file) {
  if (is.null(study_id) || study_id == "" || is.null(deg_file) || deg_file == "") return(list(res = NULL, mm = NULL, col_annot = NULL))

  cache_key <- paste(study_id, deg_file, sep = "|")
  if (exists(cache_key, envir = .study_data_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .study_data_cache, inherits = FALSE))
  }

  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(list(res = NULL, mm = NULL, col_annot = NULL))
  cfg <- yaml::read_yaml(cfg_path)

  deg_path <- file.path("data", study_id, deg_file)
  counts_rel <- resolve_counts_file_for_deg(cfg, study_id, deg_file = deg_file)
  if (is.null(counts_rel) || !nzchar(counts_rel)) {
    return(list(res = NULL, mm = NULL, col_annot = NULL))
  }
  counts_path <- file.path("data", study_id, counts_rel)
  if (!file.exists(deg_path)) return(list(res = NULL, mm = NULL, col_annot = NULL))
  counts_eds_path <- sub("\\.(tsv|txt|csv|rds)$", ".eds", counts_path, ignore.case = TRUE)
  if (identical(counts_eds_path, counts_path)) {
    counts_eds_path <- paste0(counts_path, ".eds")
  }
  counts_rds_path <- sub("\\.(tsv|txt|csv)$", ".rds", counts_path, ignore.case = TRUE)
  if (identical(counts_rds_path, counts_path)) {
    counts_rds_path <- paste0(counts_path, ".rds")
  }
  if (!file.exists(counts_eds_path) && !file.exists(counts_rds_path) && !file.exists(counts_path)) {
    return(list(res = NULL, mm = NULL, col_annot = NULL))
  }

  res <- perf_time(
    sprintf("load_study_data[%s|%s]: read_deg", study_id, deg_file),
    read.table(deg_path, sep = "\t", header = TRUE, check.names = FALSE)
  )
  message(
    "load_study_data: raw nrow(res) = ", nrow(res),
    " (study=", study_id, ", deg_file=", deg_file, ")"
  )
  mm <- perf_time(
    sprintf("load_study_data[%s|%s]: read_counts", study_id, deg_file),
    {
      if (file.exists(counts_eds_path)) {
        message("Using counts EDS file at ", counts_eds_path)
        readRDS(counts_eds_path)
      } else if (file.exists(counts_rds_path)) {
        message("Using counts RDS file at ", counts_rds_path)
        readRDS(counts_rds_path)
      } else {
        read.table(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
      }
    }
  )

  lists <- study_deg_lists_from_cfg(cfg, study_id)
  idx <- match(deg_file, vapply(lists, function(x) x$deg_file, character(1L)))
  label <- if (!is.na(idx) && !is.null(lists[[idx]]$label)) lists[[idx]]$label else NULL

  sample_names <- NULL
  if (!is.null(label)) {
    sample_names <- perf_time(
      sprintf("load_study_data[%s|%s]: samples_for_comparison", study_id, deg_file),
      samples_for_comparison(study_id, label)
    )
  }

  if (length(sample_names) > 0L) {
    keep <- intersect(sample_names, colnames(mm))
    if (length(keep) > 0L) mm <- mm[, keep, drop = FALSE]
  }

  required_base <- c("padj", "log2FoldChange", "symbol")
  if (!all(required_base %in% colnames(res))) return(list(res = NULL, mm = NULL, col_annot = NULL))

  gene_id_col <- if ("ens_gene" %in% colnames(res)) "ens_gene" else "symbol"
  gene_ids <- as.character(res[[gene_id_col]])
  keep <- !is.na(gene_ids) & nzchar(gene_ids) & (gene_ids %in% rownames(mm))
  res <- res[keep, , drop = FALSE]
  gene_ids_kept <- as.character(res[[gene_id_col]])
  mm <- mm[gene_ids_kept, , drop = FALSE]
  message(
    "load_study_data: filtered nrow(res) = ", nrow(res),
    " (kept ", sum(keep), " / ", length(keep), " rows; matched by ", gene_id_col, ")"
  )

  col_annot <- NULL
  if (!is.null(cfg$metadata_file)) {
    meta_path <- file.path("data", study_id, cfg$metadata_file)
    if (file.exists(meta_path)) {
      meta_key <- paste(study_id, cfg$metadata_file, sep = "|")
      if (exists(meta_key, envir = .metadata_cache, inherits = FALSE)) {
        meta <- get(meta_key, envir = .metadata_cache, inherits = FALSE)
      } else {
        meta <- perf_time(
          sprintf("load_study_data[%s|%s]: read_metadata", study_id, deg_file),
          read.table(meta_path, sep = "\t", header = TRUE, check.names = FALSE)
        )
        assign(meta_key, meta, envir = .metadata_cache)
      }
      if (all(c("SampleNumber", "PhenoNames") %in% colnames(meta))) {
        rownames(meta) <- meta$SampleNumber
        sample_ids <- colnames(mm)
        if (all(sample_ids %in% rownames(meta))) {
          col_annot <- meta[sample_ids, "PhenoNames", drop = FALSE]
          colnames(col_annot) <- "PhenoNames"
        }
      }
    }
  }

  out <- list(res = res, mm = mm, col_annot = col_annot)
  assign(cache_key, out, envir = .study_data_cache)
  out
}

# Merge thresholds from data/<study_id>/config.yaml into defaults from make_heatmap().
# 1) Study-level cfg$thresholds (fdr, base_mean, log2fc, svalue — all optional) applies to every DEG list.
# 2) If deg_lists[].thresholds exists for the selected list, it overrides those keys only.
# Legacy: single top-level deg_file without deg_lists uses cfg$thresholds only.
deg_list_threshold_defaults <- function(study_id, deg_file) {
  d <- default_heatmap_thresholds()
  if (is.null(study_id) || study_id == "" || is.null(deg_file) || deg_file == "") return(d)
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(d)
  cfg <- yaml::read_yaml(cfg_path)

  merge_thr <- function(d, thr) {
    if (is.null(thr) || !is.list(thr)) return(d)
    if (!is.null(thr$fdr)) d$fdr <- as.numeric(thr$fdr)
    if (!is.null(thr$base_mean)) d$base_mean <- as.numeric(thr$base_mean)
    if (!is.null(thr$log2fc)) d$log2fc <- as.numeric(thr$log2fc)
    if (!is.null(thr$svalue)) d$svalue <- as.numeric(thr$svalue)
    d
  }

  if (!is.null(cfg$thresholds) && is.list(cfg$thresholds)) {
    d <- merge_thr(d, cfg$thresholds)
  }

  if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    idx <- match(deg_file, vapply(cfg$deg_lists, function(x) x$deg_file, character(1L)))
    if (!is.na(idx) && is.list(cfg$deg_lists[[idx]]) && !is.null(cfg$deg_lists[[idx]]$thresholds)) {
      d <- merge_thr(d, cfg$deg_lists[[idx]]$thresholds)
    }
  }

  d
}

# Full expression matrix + metadata for the Gene tab (all samples in both tables; no DEG metadata_filter).
# Cached under gene_tab_study|<study_id>. Returns NULL on failure, else list(mm, metadata, study_id, cfg).
load_gene_tab_study <- function(study_id) {
  if (is.null(study_id) || study_id == "") return(NULL)
  cache_key <- paste0("gene_tab_study|", study_id)
  if (exists(cache_key, envir = .gene_tab_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gene_tab_cache, inherits = FALSE))
  }
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(NULL)
  cfg <- yaml::read_yaml(cfg_path)
  if (is.null(cfg$counts_file) || is.null(cfg$metadata_file)) return(NULL)
  counts_path <- file.path("data", study_id, cfg$counts_file)
  meta_path <- file.path("data", study_id, cfg$metadata_file)
  counts_eds_path <- sub("\\.(tsv|txt|csv|rds)$", ".eds", counts_path, ignore.case = TRUE)
  if (identical(counts_eds_path, counts_path)) {
    counts_eds_path <- paste0(counts_path, ".eds")
  }
  counts_rds_path <- sub("\\.(tsv|txt|csv)$", ".rds", counts_path, ignore.case = TRUE)
  if (identical(counts_rds_path, counts_path)) {
    counts_rds_path <- paste0(counts_path, ".rds")
  }
  if (!file.exists(meta_path)) return(NULL)
  if (!file.exists(counts_eds_path) && !file.exists(counts_rds_path) && !file.exists(counts_path)) return(NULL)

  mm <- perf_time(
    sprintf("load_gene_tab_study[%s]: read_counts", study_id),
    {
      if (file.exists(counts_eds_path)) {
        message("Using counts EDS file at ", counts_eds_path)
        readRDS(counts_eds_path)
      } else if (file.exists(counts_rds_path)) {
        message("Using counts RDS file at ", counts_rds_path)
        readRDS(counts_rds_path)
      } else {
        read.table(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
      }
    }
  )

  meta_key <- paste(study_id, cfg$metadata_file, sep = "|")
  if (exists(meta_key, envir = .metadata_cache, inherits = FALSE)) {
    meta <- get(meta_key, envir = .metadata_cache, inherits = FALSE)
  } else {
    meta <- perf_time(
      sprintf("load_gene_tab_study[%s]: read_metadata", study_id),
      read.table(meta_path, sep = "\t", header = TRUE, check.names = FALSE)
    )
    assign(meta_key, meta, envir = .metadata_cache)
  }

  if (!all(c("SampleNumber", "PhenoNames") %in% colnames(meta))) return(NULL)

  sample_cols <- colnames(mm)
  meta_ids <- as.character(meta$SampleNumber)
  common <- intersect(sample_cols, meta_ids)
  if (length(common) == 0L) return(NULL)

  mm <- mm[, common, drop = FALSE]
  idx <- match(colnames(mm), meta_ids)
  ok <- !is.na(idx)
  if (!all(ok)) {
    mm <- mm[, ok, drop = FALSE]
    idx <- idx[ok]
  }
  metadata <- meta[idx, , drop = FALSE]
  rownames(metadata) <- NULL

  mm_num <- suppressWarnings(apply(mm, 2, as.numeric))
  if (is.null(dim(mm_num))) {
    mm_num <- matrix(mm_num, ncol = 1L)
    rownames(mm_num) <- rownames(mm)
    colnames(mm_num) <- colnames(mm)
  }
  finite_vals <- as.vector(mm_num)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  is_count_like <- length(finite_vals) > 0L && all(abs(finite_vals - round(finite_vals)) < 1e-8)

  out <- list(
    mm = mm,
    metadata = metadata,
    study_id = study_id,
    cfg = cfg,
    is_count_like = is_count_like
  )
  assign(cache_key, out, envir = .gene_tab_cache)
  out
}

load_gene_tab_matrix_for_counts <- function(study_id, cfg, counts_rel) {
  if (is.null(counts_rel) || !nzchar(as.character(counts_rel))) return(NULL)
  if (is.null(cfg$metadata_file) || !nzchar(as.character(cfg$metadata_file))) return(NULL)

  cache_key <- paste("gene_tab_matrix", study_id, counts_rel, cfg$metadata_file, sep = "|")
  if (exists(cache_key, envir = .gene_tab_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gene_tab_cache, inherits = FALSE))
  }

  counts_path <- file.path("data", study_id, counts_rel)
  meta_path <- file.path("data", study_id, cfg$metadata_file)
  counts_eds_path <- sub("\\.(tsv|txt|csv|rds)$", ".eds", counts_path, ignore.case = TRUE)
  if (identical(counts_eds_path, counts_path)) {
    counts_eds_path <- paste0(counts_path, ".eds")
  }
  counts_rds_path <- sub("\\.(tsv|txt|csv)$", ".rds", counts_path, ignore.case = TRUE)
  if (identical(counts_rds_path, counts_path)) {
    counts_rds_path <- paste0(counts_path, ".rds")
  }
  if (!file.exists(meta_path)) return(NULL)
  if (!file.exists(counts_eds_path) && !file.exists(counts_rds_path) && !file.exists(counts_path)) return(NULL)

  mm <- perf_time(
    sprintf("load_gene_tab_matrix_for_counts[%s|%s]: read_counts", study_id, counts_rel),
    {
      if (file.exists(counts_eds_path)) {
        message("Using counts EDS file at ", counts_eds_path)
        readRDS(counts_eds_path)
      } else if (file.exists(counts_rds_path)) {
        message("Using counts RDS file at ", counts_rds_path)
        readRDS(counts_rds_path)
      } else {
        read.table(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
      }
    }
  )

  meta_key <- paste(study_id, cfg$metadata_file, sep = "|")
  if (exists(meta_key, envir = .metadata_cache, inherits = FALSE)) {
    meta <- get(meta_key, envir = .metadata_cache, inherits = FALSE)
  } else {
    meta <- perf_time(
      sprintf("load_gene_tab_matrix_for_counts[%s|%s]: read_metadata", study_id, counts_rel),
      read.table(meta_path, sep = "\t", header = TRUE, check.names = FALSE)
    )
    assign(meta_key, meta, envir = .metadata_cache)
  }

  if (!all(c("SampleNumber", "PhenoNames") %in% colnames(meta))) return(NULL)

  sample_cols <- colnames(mm)
  meta_ids <- as.character(meta$SampleNumber)
  common <- intersect(sample_cols, meta_ids)
  if (length(common) == 0L) return(NULL)
  mm <- mm[, common, drop = FALSE]
  idx <- match(colnames(mm), meta_ids)
  ok <- !is.na(idx)
  if (!all(ok)) {
    mm <- mm[, ok, drop = FALSE]
    idx <- idx[ok]
  }
  metadata <- meta[idx, , drop = FALSE]
  rownames(metadata) <- NULL

  mm_num <- suppressWarnings(apply(mm, 2, as.numeric))
  if (is.null(dim(mm_num))) {
    mm_num <- matrix(mm_num, ncol = 1L)
    rownames(mm_num) <- rownames(mm)
    colnames(mm_num) <- colnames(mm)
  }
  finite_vals <- as.vector(mm_num)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  is_count_like <- length(finite_vals) > 0L && all(abs(finite_vals - round(finite_vals)) < 1e-8)

  out <- list(mm = mm, metadata = metadata, is_count_like = is_count_like, counts_file = counts_rel)
  assign(cache_key, out, envir = .gene_tab_cache)
  out
}

load_gene_tab_study_panels <- function(study_id) {
  if (is.null(study_id) || study_id == "") return(NULL)
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(NULL)
  cfg <- yaml::read_yaml(cfg_path)

  global <- load_gene_tab_study(study_id)
  if (!is.null(global)) {
    return(list(
      mode = "global",
      cfg = cfg,
      panels = list(list(
        key = "global",
        title = NULL,
        labels = character(0),
        counts_file = cfg$counts_file,
        mm = global$mm,
        metadata = global$metadata,
        is_count_like = isTRUE(global$is_count_like)
      ))
    ))
  }

  lists <- study_deg_lists_from_cfg(cfg, study_id)
  if (length(lists) < 1L) return(NULL)
  labels_all <- vapply(lists, function(x) as.character(x$label), character(1L))
  dup_labels <- unique(labels_all[duplicated(labels_all)])
  if (length(dup_labels) > 0L) {
    message("Duplicate DEG labels in study ", study_id, ": ", paste(dup_labels, collapse = ", "))
  }

  entries <- list()
  for (entry in lists) {
    lbl <- if (!is.null(entry$label)) as.character(entry$label) else ""
    dfile <- if (!is.null(entry$deg_file)) as.character(entry$deg_file) else ""
    if (!nzchar(lbl) || !nzchar(dfile)) next
    rel <- resolve_counts_file_for_deg(cfg, study_id, deg_file = dfile, label = lbl)
    if (is.null(rel) || !nzchar(rel)) next
    entries[[length(entries) + 1L]] <- list(label = lbl, deg_file = dfile, counts_file = rel)
  }
  if (length(entries) < 1L) return(NULL)

  keys <- vapply(entries, function(x) x$counts_file, character(1L))
  uniq_keys <- unique(keys)
  panels <- list()
  for (k in uniq_keys) {
    idx <- which(keys == k)
    ent <- entries[idx]
    mat <- load_gene_tab_matrix_for_counts(study_id, cfg, k)
    if (is.null(mat)) next
    panel_labels <- vapply(ent, function(x) x$label, character(1L))
    panel_deg_files <- vapply(ent, function(x) x$deg_file, character(1L))
    panels[[length(panels) + 1L]] <- list(
      key = paste0("matrix_", length(panels) + 1L),
      title = paste(panel_labels, collapse = ", "),
      labels = panel_labels,
      deg_files = panel_deg_files,
      counts_file = k,
      mm = mat$mm,
      metadata = mat$metadata,
      is_count_like = isTRUE(mat$is_count_like)
    )
  }
  if (length(panels) < 1L) return(NULL)
  list(mode = "grouped", cfg = cfg, panels = panels)
}

# Precomputed global DEGs long table (cfg$gdegs_file). NULL if unset, missing, or unreadable. Cached as gene_tab_gdegs|<study_id>.
load_gene_tab_gdegs <- function(study_id) {
  if (is.null(study_id) || study_id == "") return(NULL)
  cache_key <- paste0("gene_tab_gdegs|", study_id)
  if (exists(cache_key, envir = .gene_tab_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .gene_tab_cache, inherits = FALSE))
  }
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(NULL)
  cfg <- yaml::read_yaml(cfg_path)
  rel <- cfg$gdegs_file
  if (is.null(rel) || !nzchar(as.character(rel))) {
    assign(cache_key, NULL, envir = .gene_tab_cache)
    return(NULL)
  }
  path <- file.path("data", study_id, rel)
  if (!file.exists(path)) {
    assign(cache_key, NULL, envir = .gene_tab_cache)
    return(NULL)
  }
  ext <- tolower(tools::file_ext(path))
  tab <- if (ext %in% c("rds")) {
    readRDS(path)
  } else {
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
  assign(cache_key, tab, envir = .gene_tab_cache)
  tab
}

.read_deg_table_for_symbol_map <- function(study_id, deg_rel) {
  cache_key <- paste0("deg_symbol_map|", study_id, "|", deg_rel)
  if (exists(cache_key, envir = .study_data_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .study_data_cache, inherits = FALSE))
  }
  path <- file.path("data", study_id, deg_rel)
  if (!file.exists(path)) return(NULL)
  tab <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  assign(cache_key, tab, envir = .study_data_cache)
  tab
}

.add_symbol_id_pairs <- function(sym_to_id, sym_vec, id_vec) {
  sym_vec <- toupper(trimws(as.character(sym_vec)))
  id_vec <- as.character(id_vec)
  ok <- !is.na(sym_vec) & nzchar(sym_vec) & !is.na(id_vec) & nzchar(id_vec)
  if (!any(ok)) return(sym_to_id)
  sym_vec <- sym_vec[ok]
  id_vec <- id_vec[ok]
  new_mask <- !sym_vec %in% names(sym_to_id)
  if (!any(new_mask)) return(sym_to_id)
  sym_new <- sym_vec[new_mask]
  id_new <- id_vec[new_mask]
  keep <- !duplicated(sym_new)
  sym_to_id[sym_new[keep]] <- id_new[keep]
  sym_to_id
}

#' Cached symbol -> matrix row ID map from DEG tables and optional `gdegs_file`.
study_symbol_to_matrix_id_map <- function(study_id) {
  if (is.null(study_id) || !nzchar(study_id)) return(character(0))
  cache_key <- paste0("symbol_to_id|", study_id)
  if (exists(cache_key, envir = .study_data_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .study_data_cache, inherits = FALSE))
  }
  sym_to_id <- character(0)
  lists <- study_deg_lists(study_id)
  for (entry in lists) {
    deg_rel <- entry$deg_file
    if (is.null(deg_rel) || !nzchar(deg_rel)) next
    res <- .read_deg_table_for_symbol_map(study_id, deg_rel)
    if (is.null(res)) next
    if (all(c("symbol", "ens_gene") %in% colnames(res))) {
      sym_to_id <- .add_symbol_id_pairs(sym_to_id, res$symbol, res$ens_gene)
    } else if ("symbol" %in% colnames(res) && !is.null(rownames(res)) && nzchar(rownames(res)[[1L]])) {
      sym_to_id <- .add_symbol_id_pairs(sym_to_id, res$symbol, rownames(res))
    }
  }
  gdegs <- load_gene_tab_gdegs(study_id)
  if (!is.null(gdegs) && all(c("symbol", "ens_gene") %in% colnames(gdegs))) {
    sym_to_id <- .add_symbol_id_pairs(sym_to_id, gdegs$symbol, gdegs$ens_gene)
  }
  assign(cache_key, sym_to_id, envir = .study_data_cache)
  sym_to_id
}

#' Map gene symbols to expression-matrix row IDs for one study.
#'
#' Uses DEG tables, optional `gdegs_file`, then case-insensitive rowname match on `mm`.
#' @param matrix_ids Optional character vector of matrix row IDs (avoids passing full `mm`).
#' @return data.frame with columns `symbol`, `matrix_id`, `found` (logical).
resolve_symbols_in_study_matrix <- function(study_id, symbols, mm = NULL, matrix_ids = NULL) {
  symbols <- unique(toupper(trimws(as.character(symbols))))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) < 1L) {
    return(data.frame(symbol = character(0), matrix_id = character(0), found = logical(0), stringsAsFactors = FALSE))
  }

  sym_to_id <- if (!is.null(study_id) && nzchar(study_id)) {
    study_symbol_to_matrix_id_map(study_id)
  } else {
    character(0)
  }

  matrix_id <- sym_to_id[symbols]
  names(matrix_id) <- symbols
  matrix_id[is.na(matrix_id)] <- NA_character_

  if (!is.null(matrix_ids)) {
    rn <- as.character(matrix_ids)
  } else if (!is.null(mm) && nrow(mm) > 0L) {
    rn <- rownames(mm)
  } else {
    rn <- NULL
  }
  if (!is.null(rn)) {
    rn_map <- stats::setNames(rn, tolower(rn))
    miss <- is.na(matrix_id) | !nzchar(matrix_id)
    if (any(miss)) {
      for (sym in symbols[miss]) {
        key <- tolower(sym)
        if (!key %in% names(rn_map)) next
        hit <- rn_map[[key]]
        if (!is.null(hit) && nzchar(hit)) matrix_id[[sym]] <- hit
      }
    }
  }

  found <- !is.na(matrix_id) & nzchar(matrix_id)
  if (!is.null(rn)) {
    found <- found & matrix_id %in% rn
  }
  data.frame(symbol = symbols, matrix_id = unname(matrix_id), found = found, stringsAsFactors = FALSE)
}

