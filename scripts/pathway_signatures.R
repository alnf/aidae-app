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
