 # Caches to avoid repeated disk I/O within a session
.study_data_cache <- new.env(parent = emptyenv())
.metadata_cache <- new.env(parent = emptyenv())
.comparison_cache <- new.env(parent = emptyenv())
# Gene tab: full-sample counts + metadata; global DEGs long table (separate from heatmap cache keys)
.gene_tab_cache <- new.env(parent = emptyenv())

study_deg_lists <- function(study_id) {
  if (is.null(study_id) || study_id == "") return(list())
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(list())
  cfg <- yaml::read_yaml(cfg_path)
  if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    return(cfg$deg_lists)
  }
  if (!is.null(cfg$deg_file)) {
    return(list(list(label = if (!is.null(cfg$name)) cfg$name else study_id, deg_file = cfg$deg_file)))
  }
  list()
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
  counts_path <- file.path("data", study_id, cfg$counts_file)
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

  lists <- if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) cfg$deg_lists else
    if (!is.null(cfg$deg_file)) list(list(label = if (!is.null(cfg$name)) cfg$name else study_id, deg_file = cfg$deg_file)) else list()
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

  required <- c("padj", "log2FoldChange", "symbol", "ens_gene")
  if (!all(required %in% colnames(res))) return(list(res = NULL, mm = NULL, col_annot = NULL))

  keep <- res$ens_gene %in% rownames(mm)
  res <- res[keep, , drop = FALSE]
  mm <- mm[res$ens_gene, , drop = FALSE]
  message(
    "load_study_data: filtered nrow(res) = ", nrow(res),
    " (kept ", sum(keep), " / ", length(keep), " rows)"
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

