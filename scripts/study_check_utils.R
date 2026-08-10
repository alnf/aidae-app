# Shared helpers for study raw/config consistency checks.
# Sourced by check_study_raw.R and check_study_config.R (repo root as cwd).

DEG_REQUIRED_COLS <- c("ens_gene", "symbol", "padj", "log2FoldChange")
DEG_OPTIONAL_COLS <- c("pvalue", "baseMean")
META_REQUIRED_COLS <- c("SampleNumber", "PhenoNames")
COMP_REQUIRED_COLS <- c("group1", "group2", "name", "joint")
THRESHOLD_KEYS <- c("fdr", "base_mean", "log2fc", "svalue")

# Near-miss column aliases -> canonical name (case-insensitive match on keys).
DEG_COL_ALIASES <- list(
  ens_gene = c(
    "ens_gene", "ensembl", "ensembl_id", "ensembl_gene_id", "geneid", "gene_id"
  ),
  symbol = c(
    "symbol", "gene_symbol", "gene_name", "genename", "hgnc_symbol", "external_gene_name",
    "gene_name_clean", "gene"
  ),
  padj = c("padj", "adj.p.val", "adj_p_val", "fdr", "qvalue", "q.value", "p.adjust", "padjust"),
  log2FoldChange = c(
    "log2foldchange", "log2fc", "logfc", "log2_fold_change", "log2.fold.change", "lfc"
  ),
  pvalue = c("pvalue", "p.value", "p_value", "pval"),
  baseMean = c("basemean", "base_mean", "avg_expr", "aveexpr")
)

study_root <- function(study_id, data_root = "data") {
  file.path(data_root, study_id)
}

read_tsv_header_and_nrow <- function(path) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, error = "file_missing", colnames = character(), nrow = 0L))
  }
  info <- file.info(path)
  if (!is.na(info$size) && info$size == 0) {
    return(list(ok = FALSE, error = "empty_file", colnames = character(), nrow = 0L))
  }
  tab <- tryCatch(
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
               nrows = 5L),
    error = function(e) e
  )
  if (inherits(tab, "error")) {
    return(list(ok = FALSE, error = conditionMessage(tab), colnames = character(), nrow = 0L))
  }
  if (ncol(tab) < 1L) {
    return(list(ok = FALSE, error = "no_columns", colnames = character(), nrow = 0L))
  }
  # Count data rows without loading whole file when possible
  nlines <- tryCatch(
    length(readLines(path, warn = FALSE)) - 1L,
    error = function(e) NA_integer_
  )
  if (is.na(nlines) || nlines < 0L) nlines <- nrow(tab)
  list(ok = TRUE, error = NULL, colnames = colnames(tab), nrow = as.integer(nlines), peek = tab)
}

read_counts_sample_ids <- function(path) {
  if (!file.exists(path)) {
    return(list(ok = FALSE, error = "file_missing", sample_ids = character()))
  }
  info <- file.info(path)
  if (!is.na(info$size) && info$size == 0) {
    return(list(ok = FALSE, error = "empty_file", sample_ids = character()))
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("rds", "eds")) {
    mm <- tryCatch(readRDS(path), error = function(e) e)
    if (inherits(mm, "error")) {
      return(list(ok = FALSE, error = conditionMessage(mm), sample_ids = character()))
    }
    ids <- colnames(mm)
    if (is.null(ids) || length(ids) < 1L) {
      return(list(ok = FALSE, error = "no_colnames", sample_ids = character()))
    }
    return(list(ok = TRUE, error = NULL, sample_ids = as.character(ids), nrow = nrow(mm)))
  }
  hdr <- tryCatch(
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE, nrows = 1L),
    error = function(e) e
  )
  if (inherits(hdr, "error")) {
    return(list(ok = FALSE, error = conditionMessage(hdr), sample_ids = character()))
  }
  ids <- colnames(hdr)
  if (length(ids) < 1L) {
    return(list(ok = FALSE, error = "no_columns", sample_ids = character()))
  }
  list(ok = TRUE, error = NULL, sample_ids = as.character(ids), nrow = NA_integer_)
}

normalize_col_key <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

suggest_deg_renames <- function(colnames_vec) {
  cn <- as.character(colnames_vec)
  keys <- vapply(cn, normalize_col_key, character(1L))
  suggestions <- list()
  for (canon in names(DEG_COL_ALIASES)) {
    if (canon %in% cn) next
    aliases <- DEG_COL_ALIASES[[canon]]
    alias_keys <- vapply(aliases, normalize_col_key, character(1L))
    hit <- which(keys %in% alias_keys)
    if (length(hit) == 1L) {
      from <- cn[[hit]]
      if (!identical(from, canon)) {
        suggestions[[length(suggestions) + 1L]] <- list(from = from, to = canon)
      }
    } else if (length(hit) > 1L) {
      # Prefer exact alias order match: first alias that hits
      for (a in aliases) {
        ak <- normalize_col_key(a)
        idx <- which(keys == ak)
        if (length(idx) == 1L && !identical(cn[[idx]], canon)) {
          suggestions[[length(suggestions) + 1L]] <- list(from = cn[[idx]], to = canon)
          break
        }
      }
    }
  }
  suggestions
}

check_deg_file_schema <- function(path, prefix = "") {
  # Returns list(ok, errors, warnings, infos, suggestions)
  errors <- character()
  warnings <- character()
  infos <- character()
  suggestions <- list()

  peek <- read_tsv_header_and_nrow(path)
  if (!isTRUE(peek$ok)) {
    errors <- c(errors, paste0(prefix, "DEG file unreadable (", peek$error, "): ", path))
    return(list(ok = FALSE, errors = errors, warnings = warnings, infos = infos, suggestions = suggestions))
  }
  if (peek$nrow < 1L) {
    errors <- c(errors, paste0(prefix, "DEG file has no data rows: ", path))
  }
  missing <- setdiff(DEG_REQUIRED_COLS, peek$colnames)
  if (length(missing) > 0L) {
    suggestions <- suggest_deg_renames(peek$colnames)
    sug_to <- vapply(suggestions, function(s) s$to, character(1L))
    still <- setdiff(missing, sug_to)
    if (length(still) > 0L) {
      errors <- c(
        errors,
        paste0(
          prefix, "DEG missing required columns [", paste(still, collapse = ", "),
          "]: ", path
        )
      )
    } else {
      # All missing columns have rename suggestions — still a fail until applied
      errors <- c(
        errors,
        paste0(
          prefix, "DEG missing required columns [", paste(missing, collapse = ", "),
          "] (rename suggestions available): ", path
        )
      )
    }
    for (s in suggestions) {
      if (s$to %in% missing) {
        infos <- c(
          infos,
          paste0("SUGGEST_RENAME\t", path, "\t", s$from, "\t", s$to)
        )
      }
    }
  } else {
    suggestions <- suggest_deg_renames(peek$colnames)
  }

  for (opt in DEG_OPTIONAL_COLS) {
    if (opt %in% peek$colnames) {
      infos <- c(infos, paste0(prefix, "optional column present: ", opt, " (", basename(path), ")"))
    } else {
      infos <- c(infos, paste0(prefix, "optional column absent: ", opt, " (", basename(path), ")"))
    }
  }

  list(
    ok = length(errors) == 0L,
    errors = errors,
    warnings = warnings,
    infos = infos,
    suggestions = suggestions,
    colnames = peek$colnames,
    nrow = peek$nrow
  )
}

check_metadata_schema <- function(path, prefix = "") {
  errors <- character()
  warnings <- character()
  infos <- character()
  peek <- read_tsv_header_and_nrow(path)
  if (!isTRUE(peek$ok)) {
    errors <- c(errors, paste0(prefix, "metadata unreadable (", peek$error, "): ", path))
    return(list(ok = FALSE, errors = errors, warnings = warnings, infos = infos, sample_ids = character(), pheno = character()))
  }
  if (peek$nrow < 1L) {
    errors <- c(errors, paste0(prefix, "metadata has no data rows: ", path))
  }
  missing <- setdiff(META_REQUIRED_COLS, peek$colnames)
  if (length(missing) > 0L) {
    errors <- c(
      errors,
      paste0(prefix, "metadata missing columns [", paste(missing, collapse = ", "), "]: ", path)
    )
  }
  # Load SampleNumber / PhenoNames for alignment checks
  sample_ids <- character()
  pheno <- character()
  meta <- tryCatch(
    read.table(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) e
  )
  if (!inherits(meta, "error")) {
    if ("SampleNumber" %in% colnames(meta)) {
      sample_ids <- as.character(meta$SampleNumber)
    }
    if ("PhenoNames" %in% colnames(meta)) {
      pheno <- as.character(meta$PhenoNames)
    }
  }
  list(
    ok = length(errors) == 0L,
    errors = errors,
    warnings = warnings,
    infos = infos,
    sample_ids = sample_ids,
    pheno = pheno,
    colnames = peek$colnames,
    meta = if (inherits(meta, "error")) NULL else meta
  )
}

check_counts_vs_metadata <- function(counts_path, meta_sample_ids, prefix = "") {
  errors <- character()
  warnings <- character()
  infos <- character()
  cc <- read_counts_sample_ids(counts_path)
  if (!isTRUE(cc$ok)) {
    errors <- c(errors, paste0(prefix, "counts unreadable (", cc$error, "): ", counts_path))
    return(list(ok = FALSE, errors = errors, warnings = warnings, infos = infos, sample_ids = character()))
  }
  if (length(cc$sample_ids) < 1L) {
    errors <- c(errors, paste0(prefix, "counts have no sample columns: ", counts_path))
    return(list(ok = FALSE, errors = errors, warnings = warnings, infos = infos, sample_ids = character()))
  }
  common <- intersect(cc$sample_ids, meta_sample_ids)
  if (length(common) < 1L) {
    errors <- c(
      errors,
      paste0(
        prefix, "zero overlap between counts colnames and metadata$SampleNumber: ",
        counts_path
      )
    )
  } else {
    only_counts <- setdiff(cc$sample_ids, meta_sample_ids)
    only_meta <- setdiff(meta_sample_ids, cc$sample_ids)
    if (length(only_counts) > 0L || length(only_meta) > 0L) {
      warnings <- c(
        warnings,
        paste0(
          prefix, "counts vs SampleNumber set mismatch (",
          "only_in_counts=", length(only_counts),
          ", only_in_metadata=", length(only_meta),
          ", common=", length(common), "): ", basename(counts_path)
        )
      )
    }
    # Order: compare shared samples in counts column order vs metadata order of those samples
    meta_common_order <- meta_sample_ids[meta_sample_ids %in% common]
    counts_common_order <- cc$sample_ids[cc$sample_ids %in% common]
    if (!identical(counts_common_order, meta_common_order)) {
      warnings <- c(
        warnings,
        paste0(
          prefix, "counts colnames order differs from metadata$SampleNumber order ",
          "(app will reorder): ", basename(counts_path)
        )
      )
    }
  }
  list(
    ok = length(errors) == 0L,
    errors = errors,
    warnings = warnings,
    infos = infos,
    sample_ids = cc$sample_ids
  )
}

discover_study_files <- function(study_id, data_root = "data") {
  root <- study_root(study_id, data_root)
  deg_dir <- file.path(root, "degs")
  exprs_dir <- file.path(root, "exprs")
  meta_dir <- file.path(root, "metadata")

  deg_files <- character()
  if (dir.exists(deg_dir)) {
    deg_files <- list.files(deg_dir, pattern = "\\.(tsv|txt|csv)$", full.names = TRUE, ignore.case = TRUE)
    deg_files <- deg_files[!grepl("gdegs_long", basename(deg_files), ignore.case = TRUE)]
  }

  counts_files <- character()
  if (dir.exists(exprs_dir)) {
    # Prefer tabular; include rds if no tabular sibling needed for discovery list
    tab <- list.files(exprs_dir, pattern = "\\.(tsv|txt|csv)$", full.names = TRUE, ignore.case = TRUE)
    rds <- list.files(exprs_dir, pattern = "\\.(rds|eds)$", full.names = TRUE, ignore.case = TRUE)
    counts_files <- c(tab, rds)
  }

  meta_files <- character()
  comps_files <- character()
  if (dir.exists(meta_dir)) {
    all_meta <- list.files(meta_dir, pattern = "\\.(tsv|txt|csv)$", full.names = TRUE, ignore.case = TRUE)
    is_comps <- grepl("^comps\\.", basename(all_meta), ignore.case = TRUE)
    comps_files <- all_meta[is_comps]
    meta_files <- all_meta[!is_comps]
  }

  list(
    root = root,
    deg_files = deg_files,
    counts_files = counts_files,
    meta_files = meta_files,
    comps_files = comps_files
  )
}

apply_column_renames <- function(path, renames) {
  # renames: list of list(from=, to=)
  if (length(renames) < 1L) return(invisible(FALSE))
  tab <- read.table(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  cn <- colnames(tab)
  for (r in renames) {
    if (r$from %in% cn && !(r$to %in% cn)) {
      cn[cn == r$from] <- r$to
    }
  }
  colnames(tab) <- cn
  write.table(tab, path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(TRUE)
}

normalize_deg_lists_cfg <- function(cfg, study_id) {
  if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    return(cfg$deg_lists)
  }
  if (!is.null(cfg$deg_file) && nzchar(as.character(cfg$deg_file))) {
    return(list(list(
      label = if (!is.null(cfg$name) && nzchar(as.character(cfg$name))) as.character(cfg$name) else study_id,
      deg_file = as.character(cfg$deg_file)
    )))
  }
  list()
}

print_check_messages <- function(errors, warnings, infos) {
  for (m in infos) message("INFO: ", m)
  for (m in warnings) message("WARN: ", m)
  for (m in errors) message("ERROR: ", m)
}

require_repo_root <- function() {
  if (!dir.exists("scripts") || !dir.exists("data")) {
    message("ERROR: Run from repository root (scripts/ and data/ not found). getwd()=", getwd())
    quit(status = 1L)
  }
}
