# Category-annotated per-study expression heatmaps (Panels tab).
# Depends: ComplexHeatmap, circlize; pheno_colors() from heatmap_utils.R.

`%||%` <- function(x, y) if (is.null(x)) y else x

# mm per sample column: wider for small panels, tapers as sample count grows.
.panel_mm_per_column <- function(n_samples) {
  ns <- suppressWarnings(as.integer(n_samples))
  if (is.na(ns) || ns < 1L) ns <- 1L
  nf <- as.numeric(ns)
  mm <- 2.15 + 3.9 * (nf^(-0.31))
  max(2.4, min(5.5, mm))
}

# Layout for sample columns (no column names shown).
.panel_column_layout <- function(n_samples, cluster_columns = TRUE) {
  ns <- suppressWarnings(as.integer(n_samples))
  if (is.na(ns) || ns < 1L) ns <- 1L
  mm_per_col <- .panel_mm_per_column(ns)
  list(
    body_width = grid::unit(ns * mm_per_col, "mm"),
    column_gap = grid::unit(0.35, "mm"),
    column_dend_height = if (isTRUE(cluster_columns)) {
      grid::unit(7, "mm")
    } else {
      grid::unit(0, "mm")
    },
    mm_per_col = mm_per_col
  )
}

#' Plot width (px) for Shiny `plotOutput`, from a built heatmap object.
#' Uses [ComplexHeatmap::ht_size()] for the full layout (left annotations,
#' matrix body, row names). Extra chrome must not be added on top — an oversized
#' canvas centers the heatmap and creates a large empty band on the left.
panel_heatmap_plot_width_px <- function(
    ht,
    show_row_names = TRUE,
    n_categories = 1L,
    n_samples = NULL,
    res = 96L) {
  if (is.null(ht)) return(620L)
  px_per_mm <- res / 25.4
  layout_mm <- as.numeric(grid::convertWidth(
    ComplexHeatmap::ht_size(ht)$width,
    "mm",
    valueOnly = TRUE
  ))
  layout_px <- layout_mm * px_per_mm
  w <- ceiling(layout_px * 1.04)
  as.integer(min(12000L, max(480L, w)))
}

panel_category_colors <- function(terms) {
  terms <- as.character(terms)
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) < 1L) return(character(0))
  base <- c(
    "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69",
    "#FCCDE5", "#BC80BD", "#CCEBC5", "#FFED6F", "#8DD3C7", "#D9D9D9",
    "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
  )
  cols <- base[seq_along(terms) %% length(base) + 1L]
  stats::setNames(cols, terms)
}

#' Build row index + labels from ontology TERM2GENE and symbol resolution table.
prepare_panel_heatmap_rows <- function(t2g, resolved, deduplicate_genes = FALSE) {
  if (is.null(t2g) || nrow(t2g) < 1L || is.null(resolved) || nrow(resolved) < 1L) {
    return(NULL)
  }
  t2g <- t2g[, c("term", "gene"), drop = FALSE]
  t2g$gene <- toupper(trimws(as.character(t2g$gene)))
  t2g$term <- trimws(as.character(t2g$term))
  ok <- !is.na(t2g$gene) & nzchar(t2g$gene) & !is.na(t2g$term) & nzchar(t2g$term)
  t2g <- t2g[ok, , drop = FALSE]
  if (nrow(t2g) < 1L) return(NULL)

  resolved <- resolved[resolved$found %in% TRUE, , drop = FALSE]
  if (nrow(resolved) < 1L) return(NULL)
  resolved$symbol <- toupper(trimws(as.character(resolved$symbol)))
  resolved$matrix_id <- as.character(resolved$matrix_id)

  rows <- merge(t2g, resolved[, c("symbol", "matrix_id"), drop = FALSE],
    by.x = "gene", by.y = "symbol", all.x = FALSE, sort = FALSE)
  if (nrow(rows) < 1L) return(NULL)

  term_levels <- unique(as.character(t2g$term))
  rows$term <- factor(as.character(rows$term), levels = term_levels)

  if (isTRUE(deduplicate_genes)) {
    ord <- order(rows$term, rows$gene)
    rows <- rows[ord, , drop = FALSE]
    keep <- !duplicated(rows$gene)
    rows <- rows[keep, , drop = FALSE]
  } else {
    ord <- order(rows$term, rows$gene)
    rows <- rows[ord, , drop = FALSE]
  }

  term_levels <- unique(as.character(rows$term))
  rows$term <- factor(as.character(rows$term), levels = term_levels)
  rows
}

#' Keep ontology rows with at least one significant comparison in study gdegs.
#' @return list(rows, n_before, n_after, sig_unavailable_msg)
filter_panel_row_df_significant <- function(
    row_df,
    study_id,
    panel = NULL,
    sample_filter_col = NULL,
    sample_filter_value = NULL,
    fdr = 0.05,
    log2fc = 1,
    base_mean = 0,
    svalue = 0.005) {
  if (is.null(row_df) || nrow(row_df) < 1L) {
    return(list(rows = row_df, n_before = 0L, n_after = 0L, sig_msg = NULL))
  }
  gdegs <- load_gene_tab_gdegs(study_id)
  if (is.null(gdegs) || nrow(gdegs) < 1L) {
    return(list(
      rows = row_df,
      n_before = nrow(row_df),
      n_after = nrow(row_df),
      sig_msg = "Significance filter unavailable (no gdegs for this study)."
    ))
  }
  gd <- gdegs
  facet_col <- as.character(sample_filter_col %||% "")[1L]
  facet_val <- as.character(sample_filter_value %||% "")[1L]
  if (nzchar(facet_col) && facet_col %in% colnames(gd) && nzchar(facet_val)) {
    gd <- gd[as.character(gd[[facet_col]]) == facet_val, , drop = FALSE]
  }
  if (!is.null(panel) && is.list(panel$labels) && length(panel$labels) > 0L && "joint" %in% colnames(gd)) {
    gd <- gd[as.character(gd$joint) %in% as.character(panel$labels), , drop = FALSE]
  }
  if (nrow(gd) < 1L) {
    return(list(
      rows = row_df[0, , drop = FALSE],
      n_before = nrow(row_df),
      n_after = 0L,
      sig_msg = "No DEG comparisons matched the current region/matrix filter."
    ))
  }
  lfc_col <- if ("log2FC" %in% colnames(gd)) {
    "log2FC"
  } else if ("log2FoldChange" %in% colnames(gd)) {
    "log2FoldChange"
  } else {
    NULL
  }
  if (is.null(lfc_col) || !"padj" %in% colnames(gd)) {
    return(list(
      rows = row_df,
      n_before = nrow(row_df),
      n_after = nrow(row_df),
      sig_msg = "Significance filter unavailable (gdegs missing padj/log2FC)."
    ))
  }
  padj <- as.numeric(gd$padj)
  lfc <- as.numeric(gd[[lfc_col]])
  sig <- padj <= fdr & abs(lfc) >= log2fc
  if ("baseMean" %in% colnames(gd)) {
    sig <- sig & as.numeric(gd$baseMean) >= base_mean
  }
  if ("svalue" %in% colnames(gd)) {
    sig <- sig & as.numeric(gd$svalue) <= svalue
  }
  sig[is.na(sig)] <- FALSE
  if (!any(sig)) {
    return(list(
      rows = row_df[0, , drop = FALSE],
      n_before = nrow(row_df),
      n_after = 0L,
      sig_msg = NULL
    ))
  }
  sig_sym <- unique(toupper(trimws(as.character(gd$symbol[sig]))))
  sig_sym <- sig_sym[!is.na(sig_sym) & nzchar(sig_sym)]
  sig_id <- unique(as.character(gd$ens_gene[sig]))
  sig_id <- sig_id[!is.na(sig_id) & nzchar(sig_id)]
  keep <- toupper(as.character(row_df$gene)) %in% sig_sym |
    as.character(row_df$matrix_id) %in% sig_id
  filtered <- row_df[keep, , drop = FALSE]
  list(
    rows = filtered,
    n_before = nrow(row_df),
    n_after = nrow(filtered),
    sig_msg = NULL
  )
}

#' Prepare z-scored matrix and annotations for a panel heatmap (fast; no Heatmap object).
#' @return list ready for [panel_heatmap_draw()] or `list(error = ...)`.
prepare_panel_heatmap_data <- function(
    mm,
    metadata,
    row_df,
    pheno_col = "PhenoNames",
    pheno_order = NULL,
    sample_filter_col = NULL,
    sample_filter_value = NULL,
    is_count_like = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    matrix_title = NULL) {
  if (is.null(mm) || is.null(metadata) || is.null(row_df) || nrow(row_df) < 1L) {
    return(list(error = "No genes available for panel heatmap."))
  }
  if (!"SampleNumber" %in% colnames(metadata)) {
    return(list(error = "Metadata must include SampleNumber."))
  }
  pheno_col <- as.character(pheno_col %||% "PhenoNames")[1L]
  if (!nzchar(pheno_col) || !pheno_col %in% colnames(metadata)) {
    return(list(error = paste0("Metadata column not found: ", pheno_col)))
  }

  ids <- as.character(row_df$matrix_id)
  keep <- ids %in% rownames(mm)
  row_df <- row_df[keep, , drop = FALSE]
  ids <- as.character(row_df$matrix_id)
  if (length(ids) < 1L) {
    return(list(error = "No ontology genes matched the expression matrix."))
  }

  sample_ids <- colnames(mm)
  meta <- metadata
  meta$SampleNumber <- as.character(meta$SampleNumber)
  meta <- meta[meta$SampleNumber %in% sample_ids, , drop = FALSE]
  if (nrow(meta) < 1L) {
    return(list(error = "No metadata samples overlap the expression matrix."))
  }
  idx <- match(sample_ids, meta$SampleNumber)
  ok <- !is.na(idx)
  if (!all(ok)) {
    mm <- mm[, ok, drop = FALSE]
    idx <- idx[ok]
  }
  meta <- meta[idx, , drop = FALSE]
  rownames(meta) <- meta$SampleNumber

  sample_filter_col <- as.character(sample_filter_col %||% "")[1L]
  sample_filter_value <- as.character(sample_filter_value %||% "")[1L]
  if (nzchar(sample_filter_col) && sample_filter_col %in% colnames(meta) && nzchar(sample_filter_value)) {
    keep_samples <- as.character(meta[[sample_filter_col]]) == sample_filter_value
    if (!any(keep_samples)) {
      return(list(error = paste0("No samples with ", sample_filter_col, " = ", sample_filter_value)))
    }
    meta <- meta[keep_samples, , drop = FALSE]
    mm <- mm[, rownames(meta), drop = FALSE]
    if (ncol(mm) < 1L) {
      return(list(error = "No expression samples left after sample filter."))
    }
  }

  if (!pheno_col %in% colnames(meta)) {
    pheno_col <- if ("PhenoNames" %in% colnames(meta)) "PhenoNames" else colnames(meta)[1L]
  }

  pheno_chr <- as.character(meta[[pheno_col]])
  if (is.null(pheno_order) || length(pheno_order) < 1L) {
    pheno_order <- unique(pheno_chr)
  } else {
    pheno_order <- as.character(pheno_order)
    extra <- setdiff(unique(pheno_chr), pheno_order)
    pheno_order <- c(pheno_order, extra)
  }
  meta[[pheno_col]] <- factor(pheno_chr, levels = pheno_order)
  meta <- meta[order(meta[[pheno_col]], meta$SampleNumber), , drop = FALSE]
  sample_order <- as.character(meta$SampleNumber)

  # Subset genes before column reorder (critical on large full-study matrices).
  m <- mm[ids, sample_order, drop = FALSE]
  if (inherits(m, "data.frame")) {
    m <- as.matrix(m)
  }
  rn <- rownames(m)
  cn <- colnames(m)
  m <- suppressWarnings(apply(m, 2, as.numeric))
  if (is.null(dim(m))) {
    m <- matrix(m, nrow = length(ids))
  }
  rownames(m) <- rn
  colnames(m) <- cn
  if (isTRUE(is_count_like)) {
    m <- log2(m + 0.5)
  }
  m_z <- t(scale(t(m)))
  finite_rows <- rowSums(!is.finite(m_z)) == 0L
  if (!any(finite_rows)) {
    return(list(error = "No genes with finite z-scores after normalization."))
  }
  row_df <- row_df[finite_rows, , drop = FALSE]
  m_z <- m_z[finite_rows, , drop = FALSE]

  term_present <- unique(as.character(row_df$term))
  list(
    error = NULL,
    m_z = m_z,
    row_df = row_df,
    meta = meta,
    pheno_col = pheno_col,
    term_present = term_present,
    cluster_rows = isTRUE(cluster_rows),
    cluster_columns = isTRUE(cluster_columns),
    show_row_names = isTRUE(show_row_names),
    matrix_title = matrix_title,
    n_genes = nrow(m_z),
    n_samples = ncol(m_z),
    n_pheno_groups = length(levels(meta[[pheno_col]])),
    categories = term_present
  )
}

#' Build a ComplexHeatmap object from [prepare_panel_heatmap_data()] (no draw).
panel_heatmap_build <- function(plot_data) {
  if (is.null(plot_data) || !is.null(plot_data$error)) {
    return(NULL)
  }
  m_z <- plot_data$m_z
  row_df <- plot_data$row_df
  meta <- plot_data$meta
  pheno_col <- plot_data$pheno_col %||% "PhenoNames"
  pheno_vals <- meta[[pheno_col]]
  term_present <- plot_data$term_present
  cat_cols <- panel_category_colors(term_present)
  category_vals <- factor(as.character(row_df$term), levels = term_present)
  row_split <- category_vals
  pheno_cols <- pheno_colors(levels(pheno_vals))

  top_anno <- ComplexHeatmap::HeatmapAnnotation(
    df = stats::setNames(data.frame(pheno_vals, stringsAsFactors = FALSE), pheno_col),
    col = stats::setNames(list(pheno_cols), pheno_col),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 11),
    annotation_legend_param = stats::setNames(
      list(list(
        title_gp = grid::gpar(fontsize = 12),
        labels_gp = grid::gpar(fontsize = 11)
      )),
      pheno_col
    )
  )
  left_anno <- ComplexHeatmap::rowAnnotation(
    Category = category_vals,
    col = list(Category = cat_cols),
    show_annotation_name = FALSE,
    show_legend = TRUE,
    annotation_legend_param = list(
      Category = list(
        title = "Category",
        nrow = min(6L, max(1L, length(term_present))),
        title_gp = grid::gpar(fontsize = 12),
        labels_gp = grid::gpar(fontsize = 11)
      )
    ),
    width = grid::unit(4, "mm")
  )
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2c7bb6", "white", "#d7191c"))
  has_ragg <- requireNamespace("ragg", quietly = TRUE)
  raster_device <- if (has_ragg) "agg_png" else "png"
  col_layout <- .panel_column_layout(ncol(m_z), isTRUE(plot_data$cluster_columns))

  ComplexHeatmap::Heatmap(
    m_z,
    name = "z-score",
    col = col_fun,
    width = col_layout$body_width,
    column_gap = col_layout$column_gap,
    column_dend_height = col_layout$column_dend_height,
    row_dend_width = if (isTRUE(plot_data$cluster_rows)) {
      grid::unit(10, "mm")
    } else {
      grid::unit(0, "mm")
    },
    show_row_names = isTRUE(plot_data$show_row_names),
    row_names_side = "right",
    row_labels = as.character(row_df$gene),
    show_column_names = FALSE,
    top_annotation = top_anno,
    left_annotation = left_anno,
    row_split = row_split,
    cluster_row_slices = isTRUE(plot_data$cluster_rows),
    column_split = pheno_vals,
    cluster_column_slices = isTRUE(plot_data$cluster_columns),
    cluster_rows = isTRUE(plot_data$cluster_rows),
    cluster_columns = isTRUE(plot_data$cluster_columns),
    row_title = NULL,
    column_title = plot_data$matrix_title,
    use_raster = TRUE,
    raster_device = raster_device,
    raster_quality = 1,
    heatmap_legend_param = list(
      title = "Z-score",
      legend_direction = "horizontal",
      title_position = "topcenter",
      title_gp = grid::gpar(fontsize = 12),
      labels_gp = grid::gpar(fontsize = 11)
    )
  )
}

#' Build and draw a ComplexHeatmap from [prepare_panel_heatmap_data()].
panel_heatmap_draw <- function(plot_data) {
  ht <- panel_heatmap_build(plot_data)
  if (is.null(ht)) {
    return(invisible(NULL))
  }
  ComplexHeatmap::draw(ht, merge_legend = TRUE, padding = grid::unit(c(2, 6, 2, 2), "mm"))
  invisible(ht)
}

#' Normalized expression heatmap with category row blocks and phenotype column groups.
#'
#' @param mm Expression matrix (genes x samples).
#' @param metadata Sample metadata aligned to `colnames(mm)`; needs `SampleNumber` and the column named by `pheno_col`.
#' @param pheno_col Metadata column for column grouping and top annotation (default `PhenoNames`).
#' @param row_df Output of [prepare_panel_heatmap_rows()]: columns `matrix_id`, `gene`, `term`.
#' @param is_count_like If TRUE, apply `log2(x + 0.5)` before row z-scoring.
#' @return `list(ht = heatmap object)` or `list(error = character(1))`.
make_panel_heatmap <- function(
    mm,
    metadata,
    row_df,
    pheno_col = "PhenoNames",
    pheno_order = NULL,
    sample_filter_col = NULL,
    sample_filter_value = NULL,
    is_count_like = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    matrix_title = NULL) {
  pd <- prepare_panel_heatmap_data(
    mm, metadata, row_df, pheno_col, pheno_order,
    sample_filter_col, sample_filter_value, is_count_like,
    cluster_rows, cluster_columns, show_row_names, matrix_title
  )
  if (!is.null(pd$error)) {
    return(list(error = pd$error))
  }
  ht <- panel_heatmap_build(pd)
  list(
    ht = ht,
    n_genes = pd$n_genes,
    n_samples = pd$n_samples,
    categories = pd$categories
  )
}
