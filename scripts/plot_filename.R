# Contextual plot download basenames (no extension). Used with ggiraph opts_toolbar(pngname=...).

sanitize_plot_part <- function(x, max_len = 40L) {
  x <- as.character(x %||% "")
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!nzchar(x)) x <- "na"
  if (nchar(x) > max_len) x <- substr(x, 1L, max_len)
  x
}

format_param_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) < 1L || is.na(x[[1L]])) return("na")
  x <- x[[1L]]
  if (abs(x - round(x)) < 1e-9) return(as.character(as.integer(round(x))))
  gsub("\\.?0+$", "", format(x, scientific = FALSE, trim = TRUE))
}

build_plot_filename <- function(..., max_total = 200L) {
  parts <- list(...)
  parts <- parts[!vapply(parts, function(p) is.null(p) || !nzchar(as.character(p)), logical(1L))]
  if (length(parts) < 1L) return("plot")
  out <- paste(vapply(parts, sanitize_plot_part, character(1L)), collapse = "-")
  if (nchar(out) > max_total) out <- gsub("-+$", "", substr(out, 1L, max_total))
  out
}

ontology_label <- function(pathway_file, custom_ontology = NULL) {
  if (!is.null(custom_ontology) && is.list(custom_ontology) &&
      isTRUE(custom_ontology$active) && !is.null(custom_ontology$source)) {
    src <- as.character(custom_ontology$source)
    if (nzchar(src)) return(sanitize_plot_part(tools::file_path_sans_ext(basename(src))))
  }
  pf <- as.character(pathway_file %||% "")
  if (nzchar(pf)) return(sanitize_plot_part(tools::file_path_sans_ext(basename(pf))))
  "no-ontology"
}

deg_list_label <- function(deg_file) {
  df <- as.character(deg_file %||% "")
  if (!nzchar(df)) return("no-deg")
  sanitize_plot_part(tools::file_path_sans_ext(basename(df)))
}

degs_threshold_slug <- function(
    fdr, log2fc, base_mean = NULL, svalue = NULL, custom_genes = NULL,
    has_svalue_col = FALSE, has_base_mean_col = FALSE) {
  parts <- c(
    paste0("padj", format_param_num(fdr)),
    paste0("lfc", format_param_num(log2fc))
  )
  bm <- suppressWarnings(as.numeric(base_mean))
  if (isTRUE(has_base_mean_col) && length(bm) > 0L && !is.na(bm) && bm > 0) {
    parts <- c(parts, paste0("bm", format_param_num(bm)))
  }
  sv <- suppressWarnings(as.numeric(svalue))
  if (isTRUE(has_svalue_col) && length(sv) > 0L && !is.na(sv) && sv > 0) {
    parts <- c(parts, paste0("sv", format_param_num(sv)))
  }
  if (!is.null(custom_genes) && length(custom_genes) > 0L) {
    parts <- c(parts, paste0("custom", length(custom_genes), "genes"))
  }
  paste(parts, collapse = "-")
}

degs_plot_basename <- function(study_id, deg_file, plot_kind, ...) {
  build_plot_filename(
    "degs", sanitize_plot_part(study_id), deg_list_label(deg_file),
    degs_threshold_slug(...), plot_kind
  )
}

info_pca_basename <- function(study_id, panel_key, color_by) {
  build_plot_filename(
    "info", sanitize_plot_part(study_id), sanitize_plot_part(panel_key),
    "pca", "colorby", sanitize_plot_part(color_by)
  )
}

ora_filter_slug <- function(oi) {
  if (is.null(oi) || !is.list(oi)) return("filters-na")
  paste(
    c(
      paste0("ol", format_param_num(oi$min_overlap)),
      paste0("ct", format_param_num(oi$min_count)),
      paste0("gr", format_param_num(oi$min_gene_ratio)),
      paste0("padj", format_param_num(oi$max_p_adj)),
      paste0("pv", format_param_num(oi$max_p_value)),
      paste0("top", format_param_num(oi$show_category))
    ),
    collapse = "-"
  )
}

ora_dotplot_basename <- function(pathway_file, custom_ontology, oi) {
  build_plot_filename("ora", ontology_label(pathway_file, custom_ontology), ora_filter_slug(oi), "dotplot")
}

ora_heatmap_pathway_basename <- function(pathway_file, custom_ontology, oi, pathway_desc) {
  build_plot_filename(
    "ora", ontology_label(pathway_file, custom_ontology), ora_filter_slug(oi),
    "pathway", sanitize_plot_part(pathway_desc), "heatmap"
  )
}

ora_heatmap_comparison_basename <- function(pathway_file, custom_ontology, oi, study_id, comparison) {
  build_plot_filename(
    "ora", ontology_label(pathway_file, custom_ontology), ora_filter_slug(oi),
    sanitize_plot_part(study_id), sanitize_plot_part(comparison), "heatmap"
  )
}

upset_plot_basename <- function(kind, max_degree, min_intersection, fdr, log2fc, base_mean = NULL, svalue = NULL) {
  slug <- paste(
    c(
      paste0("deg", format_param_num(max_degree)),
      paste0("min", format_param_num(min_intersection)),
      paste0("padj", format_param_num(fdr)),
      paste0("lfc", format_param_num(log2fc))
    ),
    collapse = "-"
  )
  bm <- suppressWarnings(as.numeric(base_mean))
  if (length(bm) > 0L && !is.na(bm) && bm > 0) slug <- paste(slug, paste0("bm", format_param_num(bm)), sep = "-")
  sv <- suppressWarnings(as.numeric(svalue))
  if (length(sv) > 0L && !is.na(sv) && sv > 0) slug <- paste(slug, paste0("sv", format_param_num(sv)), sep = "-")
  build_plot_filename("upset", slug, kind)
}

ggiraph_toolbar_pngname <- function(pngname) {
  if (!requireNamespace("ggiraph", quietly = TRUE)) return(NULL)
  nm <- as.character(pngname %||% "plot")
  if (!nzchar(nm)) nm <- "plot"
  if (nchar(nm) > 180L) nm <- substr(nm, 1L, 180L)
  ggiraph::opts_toolbar(pngname = nm)
}
