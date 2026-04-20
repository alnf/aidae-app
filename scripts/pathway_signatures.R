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
