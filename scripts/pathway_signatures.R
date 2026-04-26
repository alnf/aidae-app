# Parse Enrichr-style tab-separated pathway files into TERM2GENE for clusterProfiler::enricher().
# Depends: none (base R only).

if (!exists(".pathway_t2g_cache", inherits = FALSE)) {
  .pathway_t2g_cache <- new.env(parent = emptyenv())
}

#' List *.txt pathway files under databases/pathways (no dotfiles).
#' Excludes editor/checkpoint artifacts (e.g. *checkpoint*.txt).
list_pathway_txt_files <- function(dir = "databases/pathways") {
  if (!dir.exists(dir)) return(character(0))
  files <- list.files(dir, pattern = "\\.txt$", full.names = FALSE)
  files <- files[!grepl("^\\.", files)]
  files <- files[!grepl("checkpoint", files, ignore.case = TRUE)]
  sort(files)
}

#' ORA pathway file choices: prefer configured `pathways_list`; otherwise scan *.txt.
#'
#' `pathways_list` can be defined in root `config.yaml` and/or `data/<study_id>/config.yaml`.
#' Values should contain pathway filenames relative to `databases/pathways/`.
read_pathways_list_value <- function(raw_val) {
  if (is.null(raw_val)) {
    return(character(0))
  }
  vals <- trimws(as.character(unlist(raw_val)))
  vals <- vals[nzchar(vals)]
  if (length(vals) < 1L) {
    return(character(0))
  }
  vals
}

read_pathways_list_from_file <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("yaml", "yml")) {
    obj <- yaml::read_yaml(path)
    if (is.list(obj) && !is.null(obj$pathways_list)) {
      return(read_pathways_list_value(obj$pathways_list))
    }
    return(read_pathways_list_value(obj))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) < 1L) {
    return(character(0))
  }
  # Allow either one filename per line or comma-separated lines.
  vals <- trimws(unlist(strsplit(lines, ",", fixed = TRUE)))
  vals <- vals[nzchar(vals)]
  vals
}

resolve_main_pathways_list <- function(main_cfg) {
  raw <- main_cfg$pathways_list
  vals <- read_pathways_list_value(raw)
  if (length(vals) == 1L && file.exists(vals[[1L]])) {
    return(read_pathways_list_from_file(vals[[1L]]))
  }
  vals
}

list_ora_pathway_files <- function(study_ids = character(0), dir = "databases/pathways", main_config_path = "config.yaml") {
  configured <- character(0)
  if (file.exists(main_config_path)) {
    main_cfg <- yaml::read_yaml(main_config_path)
    vals <- resolve_main_pathways_list(main_cfg)
    if (length(vals) > 0L) {
      configured <- c(configured, vals)
    }
  }
  for (sid in as.character(study_ids)) {
    cfg_path <- file.path("data", sid, "config.yaml")
    if (!file.exists(cfg_path)) next
    cfg <- yaml::read_yaml(cfg_path)
    vals <- read_pathways_list_value(cfg$pathways_list)
    if (length(vals) > 0L) {
      configured <- c(configured, vals)
    }
  }
  configured <- sort(unique(configured))
  if (length(configured) > 0L) {
    return(configured)
  }
  list_pathway_txt_files(dir = dir)
}

#' Custom ontology from Excel (long format), unlike Enrichr `.txt` files (one pathway per line).
#'
#' **Input:** column 1 = gene symbol, column 2 = category (pathway). Headers `symbol`/`category`
#' preferred (`gene`/`pathway` allowed); otherwise the first two columns are used.
#'
#' **Logic:** trim both columns. Categories are taken **as the user typed them** (after trimming);
#' level order follows **first occurrence** in the sheet. For each category, the gene set is the
#' union of all symbols in rows with that category. **Gene symbols are uppercased on read** so
#' overlap with DEG tables is case-insensitive regardless of how the user typed symbols.
#'
#' The result is a `TERM2GENE` table (`term`, `gene` rows) suitable for `clusterProfiler::enricher`,
#' same as built from `.txt`, but built explicitly from per-row symbol–category assignments.
#'
#' @return `list(ok = logical, error = character(1) or NULL, t2g = data.frame(term, gene) or NULL)`
parse_ontology_xlsx_for_ora <- function(path) {
  if (is.null(path) || !nzchar(as.character(path)) || !file.exists(as.character(path))) {
    return(list(ok = FALSE, error = "File not found.", t2g = NULL))
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    return(list(
      ok = FALSE,
      error = "Install readxl: install.packages(\"readxl\")",
      t2g = NULL
    ))
  }
  raw <- tryCatch(readxl::read_xlsx(path), error = function(e) NULL)
  if (is.null(raw) || nrow(raw) < 1L) {
    return(list(ok = FALSE, error = "Empty or unreadable .xlsx.", t2g = NULL))
  }
  if (ncol(raw) < 2L) {
    return(list(ok = FALSE, error = "Need at least two columns (symbol, category).", t2g = NULL))
  }
  nm <- tolower(trimws(as.character(names(raw))))
  has_sym_name <- any(nm %in% c("symbol", "gene"))
  has_cat_name <- any(nm %in% c("category", "pathway"))
  sym_idx <- if (any(nm == "symbol")) which(nm == "symbol")[[1L]] else if (any(nm == "gene")) which(nm == "gene")[[1L]] else 1L
  cat_idx <- if (any(nm == "category")) which(nm == "category")[[1L]] else if (any(nm == "pathway")) which(nm == "pathway")[[1L]] else 2L
  if (sym_idx > ncol(raw) || cat_idx > ncol(raw) || sym_idx == cat_idx) {
    sym_idx <- 1L
    cat_idx <- 2L
  }
  syms <- as.character(raw[[sym_idx]])
  cats <- as.character(raw[[cat_idx]])
  # Normalise hidden whitespace copied from Excel/web tables.
  syms <- gsub("\u00A0", " ", syms, fixed = TRUE)
  syms <- gsub("\u200B", "", syms, fixed = TRUE)
  cats <- gsub("\u00A0", " ", cats, fixed = TRUE)
  cats <- gsub("\u200B", "", cats, fixed = TRUE)
  syms <- trimws(syms)
  cats <- trimws(cats)

  swapped_by_heuristic <- FALSE
  if (!has_sym_name && !has_cat_name) {
    # Heuristic for unnamed columns: gene symbols are usually more symbol-like and more unique.
    score_symbol_like <- function(v) {
      v <- trimws(as.character(v))
      v <- v[!is.na(v) & nzchar(v)]
      if (length(v) < 1L) return(0)
      mean(grepl("^[A-Za-z0-9._-]+$", v) & nchar(v) <= 30)
    }
    s1 <- score_symbol_like(syms)
    s2 <- score_symbol_like(cats)
    u1 <- length(unique(syms[!is.na(syms) & nzchar(syms)]))
    u2 <- length(unique(cats[!is.na(cats) & nzchar(cats)]))
    if ((s2 > s1 + 0.18 && u2 > u1) || (u1 <= 3L && u2 > u1)) {
      tmp <- syms
      syms <- cats
      cats <- tmp
      swapped_by_heuristic <- TRUE
    }
  }

  ok_row <- !is.na(syms) & nzchar(syms) & !is.na(cats) & nzchar(cats)
  syms <- syms[ok_row]
  cats <- cats[ok_row]
  if (length(syms) < 1L) {
    return(list(ok = FALSE, error = "No valid symbol/category rows.", t2g = NULL, swapped = swapped_by_heuristic))
  }
  gene <- toupper(syms)
  # Remove any remaining whitespace in symbols (gene symbols should not contain spaces).
  gene <- gsub("\\s+", "", gene, perl = TRUE)
  ok_g <- nzchar(gene)
  gene <- gene[ok_g]
  cats <- cats[ok_g]
  if (length(gene) < 1L) {
    return(list(ok = FALSE, error = "No valid gene symbols after normalisation.", t2g = NULL, swapped = swapped_by_heuristic))
  }
  # Categories as provided (trimmed only); first-occurrence order via factor then character.
  term <- as.character(factor(cats, levels = unique(cats)))
  t2g <- data.frame(term = term, gene = gene, stringsAsFactors = FALSE)
  t2g <- t2g[!duplicated(paste(t2g$term, t2g$gene, sep = "\t")), , drop = FALSE]
  rownames(t2g) <- NULL
  if (nrow(t2g) < 1L) {
    return(list(ok = FALSE, error = "No gene–category pairs after cleaning.", t2g = NULL, swapped = swapped_by_heuristic))
  }
  list(ok = TRUE, error = NULL, t2g = t2g, swapped = swapped_by_heuristic)
}

#' Read a pathway file into a long data.frame with columns term, gene.
#'
#' Lines: tab-separated; first field = pathway name; remaining non-empty fields = gene symbols
#' (Enrichr often leaves an empty second column after the name).
parse_pathway_file_to_term2gene <- function(rel_filename) {
  full <- file.path("databases/pathways", rel_filename)
  key <- paste0("t2g|", rel_filename)
  if (exists(key, envir = .pathway_t2g_cache, inherits = FALSE)) {
    return(get(key, envir = .pathway_t2g_cache, inherits = FALSE))
  }
  if (!file.exists(full)) return(NULL)
  lines <- readLines(full, warn = FALSE, encoding = "UTF-8")
  term <- character(0)
  gene <- character(0)
  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) next
    tn <- trimws(parts[[1L]])
    if (!nzchar(tn)) next
    genes <- trimws(parts[-1L])
    genes <- genes[nzchar(genes)]
    genes <- unique(genes)
    if (length(genes) < 1L) next
    term <- c(term, rep(tn, length(genes)))
    gene <- c(gene, genes)
  }
  if (length(term) == 0L) {
    out <- data.frame(term = character(0), gene = character(0), stringsAsFactors = FALSE)
  } else {
    out <- data.frame(term = term, gene = gene, stringsAsFactors = FALSE)
  }
  assign(key, out, envir = .pathway_t2g_cache)
  out
}


pathway_terms_from_file <- function(rel_filename) {
  t2g <- parse_pathway_file_to_term2gene(rel_filename)
  if (is.null(t2g) || nrow(t2g) == 0L) return(character(0))
  sort(unique(as.character(t2g$term)))
}

#' Subset TERM2GENE to one pathway term.
term2gene_single_term <- function(rel_filename, term) {
  t2g <- parse_pathway_file_to_term2gene(rel_filename)
  if (is.null(t2g) || nrow(t2g) == 0L) return(NULL)
  sub <- t2g[t2g$term == term, , drop = FALSE]
  if (nrow(sub) == 0L) return(NULL)
  sub
}
