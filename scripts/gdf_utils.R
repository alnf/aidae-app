# Offline builder for long DE tables used by the Gene tab (p-value brackets).
# Callable from Rscript with setwd() to the project root.

# Infer token from deg file path (e.g. degs_LV_*.tsv -> "LV") — only used as fallback when gene_tab_facet is "Region".
infer_region_from_deg_path <- function(deg_rel_path) {
  bn <- basename(deg_rel_path)
  m <- regexec("^degs_([A-Za-z0-9]+)_", bn)
  reg <- regmatches(bn, m)[[1L]]
  if (length(reg) >= 2L) return(reg[[2L]])
  NA_character_
}

# Values for cfg$gene_tab_facet: per-row from DE table, else single level from metadata_filter, else optional path heuristic for Region, else NA.
resolve_facet_vector <- function(res, deg_entry, deg_rel_path, facet_col) {
  if (is.null(facet_col) || !nzchar(facet_col)) {
    return(rep(NA_character_, nrow(res)))
  }
  if (facet_col %in% colnames(res)) {
    return(as.character(res[[facet_col]]))
  }
  mf <- deg_entry$metadata_filter
  if (is.list(mf) && !is.null(mf[[facet_col]])) {
    rv <- unlist(mf[[facet_col]])
    if (length(rv) == 1L) {
      return(rep(as.character(rv[[1L]]), nrow(res)))
    }
  }
  if (identical(facet_col, "Region")) {
    ir <- infer_region_from_deg_path(deg_rel_path)
    if (!is.na(ir) && nzchar(ir)) {
      return(rep(ir, nrow(res)))
    }
  }
  rep(NA_character_, nrow(res))
}

read_study_config <- function(study_id, data_root = "data") {
  cfg_path <- file.path(data_root, study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    stop("config not found: ", cfg_path, call. = FALSE)
  }
  yaml::read_yaml(cfg_path)
}

read_comparison_table <- function(study_id, cfg, data_root = "data") {
  if (is.null(cfg$comparison_file) || cfg$comparison_file == "") {
    stop("comparison_file missing in config for study ", study_id, call. = FALSE)
  }
  comp_path <- file.path(data_root, study_id, cfg$comparison_file)
  if (!file.exists(comp_path)) {
    stop("comparison file not found: ", comp_path, call. = FALSE)
  }
  read.table(comp_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

normalize_deg_lists <- function(study_id, cfg) {
  if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    return(cfg$deg_lists)
  }
  if (!is.null(cfg$deg_file)) {
    return(list(list(
      label = if (!is.null(cfg$name)) cfg$name else study_id,
      deg_file = cfg$deg_file
    )))
  }
  list()
}

# Optional facet column name from config (`gene_tab_facet`); NULL/empty => no facet column in output.
parse_gene_tab_facet <- function(cfg) {
  fc <- cfg$gene_tab_facet
  if (is.null(fc)) return(NULL)
  fc <- as.character(fc)
  if (length(fc) != 1L || !nzchar(fc)) return(NULL)
  fc
}

#' Build a long DE table for one study (all comparisons in deg_lists).
#'
#' Base columns: study_id, ens_gene, symbol, pvalue, padj, log2FC, group1, group2, comp, joint.
#' If `gene_tab_facet` is set in config, one extra column with that name (values aligned with metadata / DE tables for bracket placement).
build_gene_deg_long <- function(study_id, data_root = "data") {
  cfg <- read_study_config(study_id, data_root)
  facet_col <- parse_gene_tab_facet(cfg)
  deg_lists <- normalize_deg_lists(study_id, cfg)
  if (length(deg_lists) == 0L) {
    warning("No deg_lists or deg_file for study ", study_id, "; returning empty table", call. = FALSE)
    return(empty_gene_deg_long(facet_col))
  }
  comp <- read_comparison_table(study_id, cfg, data_root)
  req_comp <- c("joint", "group1", "group2", "name")
  if (!all(req_comp %in% colnames(comp))) {
    stop("comparison table must contain columns: ", paste(req_comp, collapse = ", "), call. = FALSE)
  }

  pieces <- list()
  for (entry in deg_lists) {
    label <- entry$label
    deg_rel <- entry$deg_file
    if (is.null(deg_rel) || deg_rel == "" || is.null(label) || label == "") next
    deg_path <- file.path(data_root, study_id, deg_rel)
    if (!file.exists(deg_path)) {
      warning("Skipping missing DEG file: ", deg_path, call. = FALSE)
      next
    }
    rows <- which(comp$joint == label)
    if (length(rows) == 0L) {
      warning("No comps.tsv row with joint == ", label, " for study ", study_id, call. = FALSE)
      next
    }
    r <- rows[[1L]]
    group1 <- as.character(comp$group1[[r]])
    group2 <- as.character(comp$group2[[r]])
    comp_name <- as.character(comp$name[[r]])

    res <- read.table(deg_path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    req_deg <- c("ens_gene", "symbol", "padj", "log2FoldChange")
    if (!all(req_deg %in% colnames(res))) {
      warning("DEG file missing required columns (ens_gene, symbol, padj, log2FoldChange): ", deg_path, call. = FALSE)
      next
    }
    facet_vec <- resolve_facet_vector(res, entry, deg_rel, facet_col)
    n <- nrow(res)
    row <- data.frame(
      study_id = rep(study_id, n),
      ens_gene = as.character(res$ens_gene),
      symbol = as.character(res$symbol),
      pvalue = if ("pvalue" %in% colnames(res)) as.numeric(res$pvalue) else rep(NA_real_, n),
      padj = as.numeric(res$padj),
      log2FC = as.numeric(res$log2FoldChange),
      group1 = rep(group1, n),
      group2 = rep(group2, n),
      comp = rep(comp_name, n),
      joint = rep(label, n),
      stringsAsFactors = FALSE
    )
    if (!is.null(facet_col)) {
      row[[facet_col]] <- facet_vec
    }
    pieces[[length(pieces) + 1L]] <- row
  }
  if (length(pieces) == 0L) {
    return(empty_gene_deg_long(facet_col))
  }
  do.call(rbind, pieces)
}

empty_gene_deg_long <- function(facet_col = NULL) {
  d <- data.frame(
    study_id = character(),
    ens_gene = character(),
    symbol = character(),
    pvalue = numeric(),
    padj = numeric(),
    log2FC = numeric(),
    group1 = character(),
    group2 = character(),
    comp = character(),
    joint = character(),
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_col) && nzchar(as.character(facet_col))) {
    d[[facet_col]] <- character()
  }
  d
}

#' Write `build_gene_deg_long` output using `gdegs_file` (global DEGs) from the study config.
#'
#' Creates parent directories as needed. Format is chosen by file extension (`.tsv` / `.txt` or `.rds`).
write_gene_deg_long <- function(df, study_id, data_root = "data") {
  cfg <- read_study_config(study_id, data_root)
  if (is.null(cfg$gdegs_file) || !nzchar(as.character(cfg$gdegs_file))) {
    stop("gdegs_file is not set in config for study ", study_id, call. = FALSE)
  }
  out_path <- file.path(data_root, study_id, cfg$gdegs_file)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(out_path))
  if (ext %in% c("rds")) {
    saveRDS(df, out_path)
  } else {
    write.table(df, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  invisible(out_path)
}
