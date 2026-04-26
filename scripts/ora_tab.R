# ORA tab: clusterProfiler::enricher, one faceted ggplot (facet = study; x = DEG lists; color = GeneRatio).
# Depends: scripts/ora_cache.R, scripts/heatmap_utils.R, scripts/study_data.R, scripts/pathway_signatures.R;
#   Bioconductor clusterProfiler; ggplot2. Optional: ggh4x (shared y-axis labels on outer margin only).

# Compact dimensions for message-only ggplot.
.ora_msg_plot_px <- function() {
  list(width = 400L, height = 110L)
}

.ora_msg_plot <- function(msg) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 3.4) +
    ggplot2::theme_void() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

# Pixel size for Shiny renderPlot: total width = sum of per-study panel widths (matches facet_grid space = "free_x").
.ora_faceted_plot_dims <- function(comp_per_panel, max_path) {
  np <- max(1L, as.integer(max_path))
  if (length(comp_per_panel) < 1L) {
    return(list(width = 520L, height = 340L))
  }
  w_body <- 0L
  for (nc in comp_per_panel) {
    nc <- max(1L, as.integer(nc))
    w_body <- w_body + max(118L, 68L + 48L * nc)
  }
  w <- min(9000L, max(480L, 168L + w_body))
  h <- min(1600L, max(340L, 128L + 24L * np))
  list(width = w, height = h)
}

# Non-linear height scaling for row-intent (pathway genes) heatmaps.
# Tuned so small gene sets remain readable while large sets do not explode linearly.
.ora_row_heatmap_height_px <- function(n_genes) {
  n <- suppressWarnings(as.integer(n_genes))
  if (is.na(n) || n < 1L) return(440L)
  nf <- as.numeric(n)
  # Tuned non-linear scaling: keep small sets compact, expand large sets further.
  h <- 205 + 42 * (nf^0.82) + 12 * sqrt(nf) + 9 * log1p(nf)
  as.integer(min(6800L, max(440L, round(h))))
}

# Non-linear height scaling for column-intent (pathway rows) heatmaps.
# Keeps 1-2 pathways compact while allowing larger pathway sets to grow.
.ora_comparison_heatmap_height_px <- function(n_pathways) {
  n <- suppressWarnings(as.integer(n_pathways))
  if (is.na(n) || n < 1L) return(420L)
  nf <- as.numeric(n)
  h <- 220 + 52 * (nf^0.92) + 8 * log1p(nf)
  as.integer(min(1800L, max(360L, round(h))))
}

# Width scaling for heatmaps: keep small plots compact, allow wide plots with horizontal scroll.
.ora_heatmap_width_px <- function(n_cols) {
  nc <- suppressWarnings(as.integer(n_cols))
  if (is.na(nc) || nc < 1L) nc <- 1L
  if (nc <= 1L) return(380L)
  if (nc == 2L) return(500L)
  if (nc == 3L) return(620L)
  w <- 220 + 62 * nc
  as.integer(min(9000L, max(760L, round(w))))
}

#' Faceted ggplot: facet = study; x = DEG list; color = gene_ratio (blue → red).
#' Shared y (pathways) across facets; dots only where enriched. Panel width scales with
#' number of DEG lists per study (`facet_grid(..., space = "free_x")`). Y labels once on the
#' margin when package ggh4x is available (`facet_grid2(..., axes = "margins")`).
.ora_faceted_comparison_plot <- function(combined_df, study_label_levels, pathway_level_order = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(combined_df) || nrow(combined_df) < 1L) {
    return(.ora_msg_plot("No pathways to display."))
  }
  combined_df$study_label <- factor(
    as.character(combined_df$study_label),
    levels = study_label_levels
  )

  plot_df <- combined_df
  # Character (not one global factor) so each facet's discrete x length = that study's DEG lists only
  # (needed for facet_grid(space = "free_x") panel widths).
  plot_df$comparison <- as.character(plot_df$comparison)

  if (!is.null(pathway_level_order) && length(pathway_level_order) > 0L) {
    plot_df$Description <- factor(as.character(plot_df$Description), levels = pathway_level_order)
  } else {
    y_levels <- if (is.factor(plot_df$Description)) {
      levels(plot_df$Description)
    } else {
      unique(as.character(plot_df$Description))
    }
    plot_df$Description <- factor(as.character(plot_df$Description), levels = y_levels)
  }

  facet_obj <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    ggh4x::facet_grid2(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      axes = "margins",
      remove_labels = "none",
      drop = FALSE
    )
  } else {
    ggplot2::facet_grid(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    )
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = comparison,
    y = Description,
    color = gene_ratio,
    size = Count
  )) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0, add = 0.35)) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_color_gradient(
      name = "Gene ratio",
      low = "blue",
      high = "red",
      limits = c(0, 1),
      na.value = NA
    ) +
    ggplot2::scale_size_area(name = "Overlap\n(Count)", max_size = 9) +
    facet_obj +
    ggplot2::labs(
      title = "Overrepresentation analysis (ORA)",
      subtitle = "Thresholds: defaults or per study/DEG list in config.yaml",
      x = "Comparison (DEG list)",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(face = "bold", size = 12.5),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 11),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10.5),
      legend.title = ggplot2::element_text(size = 11)
    )
}

.ora_faceted_comparison_plot_girafe <- function(combined_df, study_label_levels, pathway_level_order = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggiraph", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(combined_df) || nrow(combined_df) < 1L) {
    return(.ora_msg_plot("No pathways to display."))
  }
  plot_df <- combined_df
  plot_df$study_label <- factor(
    as.character(plot_df$study_label),
    levels = study_label_levels
  )
  plot_df$comparison <- as.character(plot_df$comparison)
  if (!is.null(pathway_level_order) && length(pathway_level_order) > 0L) {
    plot_df$Description <- factor(as.character(plot_df$Description), levels = pathway_level_order)
  } else {
    y_levels <- unique(as.character(plot_df$Description))
    plot_df$Description <- factor(as.character(plot_df$Description), levels = y_levels)
  }
  plot_df$selection_id <- paste(
    as.character(plot_df$study_id),
    as.character(plot_df$study_label),
    as.character(plot_df$comparison),
    as.character(plot_df$ID),
    as.character(plot_df$Description),
    sep = "|||"
  )
  plot_df$tooltip <- paste0(
    "<b>Study:</b> ", as.character(plot_df$study_label),
    "<br><b>Comparison:</b> ", as.character(plot_df$comparison),
    "<br><b>Pathway:</b> ", as.character(plot_df$Description),
    "<br><b>Count:</b> ", as.character(plot_df$Count),
    "<br><b>Gene ratio:</b> ", formatC(as.numeric(plot_df$gene_ratio), digits = 3, format = "f")
  )

  facet_obj <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    ggh4x::facet_grid2(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      axes = "margins",
      remove_labels = "none",
      drop = FALSE
    )
  } else {
    ggplot2::facet_grid(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      drop = FALSE
    )
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = comparison,
    y = Description,
    color = gene_ratio,
    size = Count
  )) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(data_id = selection_id, tooltip = tooltip),
      na.rm = TRUE
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0, add = 0.35)) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_color_gradient(
      name = "Gene ratio",
      low = "blue",
      high = "red",
      limits = c(0, 1),
      na.value = NA
    ) +
    ggplot2::scale_size_area(name = "Overlap\n(Count)", max_size = 9) +
    facet_obj +
    ggplot2::labs(
      title = "Overrepresentation analysis (ORA)",
      subtitle = "Click a dot; choose whether it selects a pathway (row intent) or comparison (column intent).",
      x = "Comparison (DEG list)",
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(face = "bold", size = 12.5),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 11),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10.5),
      legend.title = ggplot2::element_text(size = 11)
    )
}

ora_parse_gene_set <- function(x) {
  if (is.null(x) || length(x) < 1L) return(character(0))
  vals <- trimws(as.character(unlist(x)))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) < 1L) return(character(0))
  parts <- trimws(unlist(strsplit(vals, "/", fixed = TRUE)))
  parts <- toupper(parts)
  unique(parts[!is.na(parts) & nzchar(parts)])
}

ora_find_deg_file_for_comparison <- function(study_id, comparison_label) {
  lists <- study_deg_lists(study_id)
  if (length(lists) < 1L) return(NULL)
  idx <- which(vapply(lists, function(e) {
    identical(as.character(ora_deg_entry_comparison_label(e)), as.character(comparison_label))
  }, logical(1L)))
  if (length(idx) != 1L) return(NULL)
  lists[[idx]]$deg_file
}

if (!exists(".ora_sig_gene_cache", inherits = FALSE)) {
  .ora_sig_gene_cache <- new.env(parent = emptyenv())
}

ora_significant_symbols <- function(study_id, deg_file) {
  if (is.null(deg_file) || !nzchar(as.character(deg_file))) return(character(0))
  key <- paste0(as.character(study_id), "|", as.character(deg_file))
  if (exists(key, envir = .ora_sig_gene_cache, inherits = FALSE)) {
    return(get(key, envir = .ora_sig_gene_cache, inherits = FALSE))
  }
  thr <- deg_list_threshold_defaults(study_id, deg_file)

  # Fast path: use precomputed gdegs long table when available (avoids loading counts matrices).
  gdegs <- load_gene_tab_gdegs(study_id)
  if (!is.null(gdegs) && is.data.frame(gdegs) && "symbol" %in% colnames(gdegs)) {
    lists <- study_deg_lists(study_id)
    idx_cmp <- match(as.character(deg_file), vapply(lists, function(x) as.character(x$deg_file), character(1L)))
    cmp_lbl <- if (!is.na(idx_cmp) && idx_cmp >= 1L) as.character(ora_deg_entry_comparison_label(lists[[idx_cmp]])) else NULL
    cmp_col <- if ("joint" %in% colnames(gdegs)) {
      "joint"
    } else if ("comp" %in% colnames(gdegs)) {
      "comp"
    } else {
      NULL
    }
    gsub <- if (!is.null(cmp_col) && !is.null(cmp_lbl) && nzchar(cmp_lbl)) {
      gdegs[as.character(gdegs[[cmp_col]]) == cmp_lbl, , drop = FALSE]
    } else {
      gdegs
    }
    if (nrow(gsub) > 0L) {
      keep <- rep(TRUE, nrow(gsub))
      if ("padj" %in% colnames(gsub)) {
        keep <- keep & !is.na(as.numeric(gsub$padj)) & as.numeric(gsub$padj) <= as.numeric(thr$fdr)
      } else if ("pvalue" %in% colnames(gsub)) {
        keep <- keep & !is.na(as.numeric(gsub$pvalue)) & as.numeric(gsub$pvalue) <= as.numeric(thr$fdr)
      }
      lfc_col <- if ("log2FC" %in% colnames(gsub)) {
        "log2FC"
      } else if ("log2FoldChange" %in% colnames(gsub)) {
        "log2FoldChange"
      } else {
        NULL
      }
      if (!is.null(lfc_col)) {
        keep <- keep & !is.na(as.numeric(gsub[[lfc_col]])) & abs(as.numeric(gsub[[lfc_col]])) >= as.numeric(thr$log2fc)
      }
      syms <- toupper(trimws(as.character(gsub$symbol[keep])))
      syms <- unique(syms[!is.na(syms) & nzchar(syms)])
      assign(key, syms, envir = .ora_sig_gene_cache)
      return(syms)
    }
  }

  # Fallback: full DEG + counts load to apply the exact heatmap threshold logic.
  loaded <- load_study_data(study_id, deg_file)
  res <- loaded$res
  mm <- loaded$mm
  if (is.null(res) || is.null(mm) || nrow(res) < 1L) {
    assign(key, character(0), envir = .ora_sig_gene_cache)
    return(character(0))
  }
  idx <- filter_heatmap_row_index(
    res,
    mm,
    thr$fdr,
    thr$base_mean,
    thr$log2fc,
    thr$svalue
  )
  if (is.null(idx) || length(idx) < 1L || !("symbol" %in% colnames(res))) {
    assign(key, character(0), envir = .ora_sig_gene_cache)
    return(character(0))
  }
  syms <- toupper(trimws(as.character(res$symbol[idx])))
  syms <- unique(syms[!is.na(syms) & nzchar(syms)])
  assign(key, syms, envir = .ora_sig_gene_cache)
  syms
}

ora_term_genes <- function(pathway_file, pathway_id, pathway_desc, term2gene_df = NULL) {
  t2g <- if (!is.null(term2gene_df) && is.data.frame(term2gene_df) && nrow(term2gene_df) > 0L) {
    term2gene_df
  } else {
    parse_pathway_file_to_term2gene(pathway_file)
  }
  if (is.null(t2g) || nrow(t2g) < 1L) return(character(0))
  term <- as.character(t2g$term)
  keep <- term %in% c(as.character(pathway_id), as.character(pathway_desc))
  genes <- toupper(trimws(as.character(t2g$gene[keep])))
  unique(genes[!is.na(genes) & nzchar(genes)])
}

ora_read_deg_logfc <- function(study_id, comparison_label, deg_file = NULL) {
  gdegs <- load_gene_tab_gdegs(study_id)
  if (!is.null(gdegs) && is.data.frame(gdegs) &&
      all(c("symbol", "log2FC") %in% colnames(gdegs))) {
    cmp_col <- if ("joint" %in% colnames(gdegs)) "joint" else if ("comp" %in% colnames(gdegs)) "comp" else NULL
    if (!is.null(cmp_col)) {
      gsub <- gdegs[as.character(gdegs[[cmp_col]]) == as.character(comparison_label), , drop = FALSE]
      if (nrow(gsub) > 0L) {
        out <- data.frame(
          symbol = toupper(trimws(as.character(gsub$symbol))),
          log2FC = as.numeric(gsub$log2FC),
          stringsAsFactors = FALSE
        )
        out <- out[!is.na(out$symbol) & nzchar(out$symbol), , drop = FALSE]
        out <- out[!is.na(out$log2FC), , drop = FALSE]
        if (nrow(out) > 0L) return(out)
      }
    }
  }

  if (is.null(deg_file) || !nzchar(as.character(deg_file))) return(NULL)
  loaded <- load_study_data(study_id, deg_file)
  res <- loaded$res
  if (is.null(res) || nrow(res) < 1L || !all(c("symbol", "log2FoldChange") %in% colnames(res))) {
    return(NULL)
  }
  out <- data.frame(symbol = toupper(trimws(as.character(res$symbol))), log2FC = as.numeric(res$log2FoldChange), stringsAsFactors = FALSE)
  out <- out[!is.na(out$symbol) & nzchar(out$symbol), , drop = FALSE]
  out <- out[!is.na(out$log2FC), , drop = FALSE]
  if (nrow(out) > 0L) out else NULL
}

ora_pathway_genes <- function(pathway_file, pathway_id, pathway_desc, term2gene_df = NULL) {
  ora_term_genes(pathway_file, pathway_id, pathway_desc, term2gene_df = term2gene_df)
}

ora_build_pathway_mode_matrix <- function(
    pathway_file,
    pathway_id,
    pathway_desc,
    study_ids,
    study_labels,
    combined_df,
    hide_empty_comparisons = FALSE,
    term2gene_df = NULL) {
  pathway_genes <- character(0)
  sub <- NULL
  if (!is.null(combined_df) && nrow(combined_df) > 0L) {
    sub <- combined_df[as.character(combined_df$ID) == as.character(pathway_id), , drop = FALSE]
  }
  if (!is.null(sub) && nrow(sub) > 0L && "geneID" %in% colnames(sub)) {
    if (all(!is.na(sub$geneID) & nzchar(as.character(sub$geneID)))) {
      pathway_genes <- ora_parse_gene_set(sub$geneID)
    }
  }
  if (length(pathway_genes) < 1L && !is.null(sub) && nrow(sub) > 0L) {
    term_genes <- ora_term_genes(pathway_file, pathway_id, pathway_desc, term2gene_df = term2gene_df)
    if (length(term_genes) > 0L) {
      for (i in seq_len(nrow(sub))) {
        sid <- as.character(sub$study_id[i])
        cmp <- as.character(sub$comparison[i])
        deg_file <- ora_find_deg_file_for_comparison(sid, cmp)
        sig_syms <- ora_significant_symbols(sid, deg_file)
        if (length(sig_syms) < 1L) next
        pathway_genes <- c(pathway_genes, intersect(term_genes, sig_syms))
      }
      pathway_genes <- unique(pathway_genes)
    }
  }
  if (length(pathway_genes) < 1L) {
    pathway_genes <- ora_pathway_genes(pathway_file, pathway_id, pathway_desc, term2gene_df = term2gene_df)
  }
  if (length(pathway_genes) < 1L) {
    return(list(error = "No genes found for selected pathway.", mat = NULL))
  }

  vals <- list()
  sig_vals <- list()
  for (sid in study_ids) {
    lists <- study_deg_lists(sid)
    if (length(lists) < 1L) next
    s_lbl <- study_labels[[sid]]
    if (is.null(s_lbl) || !nzchar(as.character(s_lbl))) s_lbl <- sid
    for (entry in lists) {
      deg_file <- entry$deg_file
      if (is.null(deg_file) || !nzchar(as.character(deg_file))) next
      cmp_lbl <- ora_deg_entry_comparison_label(entry)
      col_name <- paste0(as.character(s_lbl), "::", cmp_lbl)
      tab <- ora_read_deg_logfc(sid, cmp_lbl, deg_file)
      if (is.null(tab) || nrow(tab) < 1L) next
      hit <- match(pathway_genes, tab$symbol)
      x <- rep(NA_real_, length(pathway_genes))
      ok <- !is.na(hit)
      x[ok] <- tab$log2FC[hit[ok]]
      vals[[col_name]] <- x
      sig_syms <- ora_significant_symbols(sid, deg_file)
      sig_vals[[col_name]] <- pathway_genes %in% sig_syms
    }
  }
  if (length(vals) < 1L) {
    return(list(error = "No DEG tables available for selected pathway.", mat = NULL, sig_mat = NULL))
  }
  mat <- do.call(cbind, vals)
  sig_mat <- do.call(cbind, sig_vals)
  rownames(mat) <- pathway_genes
  rownames(sig_mat) <- pathway_genes
  keep_cols <- rep(TRUE, ncol(mat))
  if (isTRUE(hide_empty_comparisons)) {
    keep_cols <- rep(FALSE, ncol(mat))
    if (!is.null(sub) && nrow(sub) > 0L) {
      # Preserve comparison order from the main ORA dotplot payload.
      ord <- unique(paste0(as.character(combined_df$study_label), "::", as.character(combined_df$comparison)))
      sig_cols <- unique(paste0(as.character(sub$study_label), "::", as.character(sub$comparison)))
      keep_names <- ord[ord %in% sig_cols]
      keep_cols <- colnames(mat) %in% keep_names
      if (sum(keep_cols) > 0L) {
        mat <- mat[, keep_cols, drop = FALSE]
        sig_mat <- sig_mat[, keep_cols, drop = FALSE]
        reorder_names <- keep_names[keep_names %in% colnames(mat)]
        if (length(reorder_names) > 0L) {
          mat <- mat[, reorder_names, drop = FALSE]
          sig_mat <- sig_mat[, reorder_names, drop = FALSE]
        }
      }
    }
  }
  if (ncol(mat) < 1L) {
    return(list(error = "No comparisons remained after pathway-significance filter.", mat = NULL, sig_mat = NULL))
  }
  keep_rows <- rowSums(!is.na(mat)) > 0L
  mat <- mat[keep_rows, , drop = FALSE]
  sig_mat <- sig_mat[keep_rows, , drop = FALSE]
  if (nrow(mat) < 1L) {
    return(list(error = "Selected pathway genes are absent from all DEG tables.", mat = NULL, sig_mat = NULL))
  }
  list(error = NULL, mat = mat, sig_mat = sig_mat)
}

ora_build_comparison_mode_matrix <- function(
    pathway_file,
    selected_sid,
    selected_comparison,
    combined_df,
    term2gene_df = NULL) {
  if (is.null(combined_df) || nrow(combined_df) < 1L) {
    return(list(error = "No ORA rows available for selected comparison.", mat = NULL, sig_mat = NULL))
  }
  sub <- combined_df[
    as.character(combined_df$study_id) == as.character(selected_sid) &
      as.character(combined_df$comparison) == as.character(selected_comparison),
    ,
    drop = FALSE
  ]
  if (nrow(sub) < 1L) {
    return(list(error = "Selected comparison has no pathways in current ORA view.", mat = NULL, sig_mat = NULL))
  }
  lists <- study_deg_lists(selected_sid)
  idx <- which(vapply(lists, function(e) {
    identical(as.character(ora_deg_entry_comparison_label(e)), as.character(selected_comparison))
  }, logical(1L)))
  if (length(idx) != 1L) {
    return(list(error = "Could not map selected comparison to one DEG list.", mat = NULL, sig_mat = NULL))
  }
  tab <- ora_read_deg_logfc(selected_sid, selected_comparison, lists[[idx]]$deg_file)
  if (is.null(tab) || nrow(tab) < 1L) {
    return(list(error = "Selected comparison DEG table is empty.", mat = NULL, sig_mat = NULL))
  }
  t2g <- if (!is.null(term2gene_df) && is.data.frame(term2gene_df) && nrow(term2gene_df) > 0L) {
    term2gene_df
  } else {
    parse_pathway_file_to_term2gene(pathway_file)
  }
  if (is.null(t2g) || nrow(t2g) < 1L) {
    return(list(error = "Pathway database has no term-to-gene mappings.", mat = NULL, sig_mat = NULL))
  }
  pathway_ids <- unique(as.character(sub$ID))
  pathways <- unique(data.frame(
    ID = as.character(sub$ID),
    Description = as.character(sub$Description),
    stringsAsFactors = FALSE
  ))
  gene_by_term <- stats::setNames(vector("list", length(pathway_ids)), pathway_ids)
  genes_union <- character(0)
  sig_syms <- ora_significant_symbols(selected_sid, lists[[idx]]$deg_file)
  for (pid in pathway_ids) {
    row_pid <- sub[as.character(sub$ID) == pid, , drop = FALSE]
    g <- character(0)
    if (nrow(row_pid) > 0L && "geneID" %in% colnames(row_pid)) {
      if (all(!is.na(row_pid$geneID) & nzchar(as.character(row_pid$geneID)))) {
        g <- ora_parse_gene_set(row_pid$geneID)
      }
    }
    if (length(g) < 1L) {
      g <- toupper(trimws(as.character(t2g$gene[as.character(t2g$term) == pid])))
      g <- unique(g[!is.na(g) & nzchar(g)])
      if (length(sig_syms) > 0L) {
        g <- intersect(g, sig_syms)
      }
    }
    g <- g[g %in% tab$symbol]
    gene_by_term[[pid]] <- g
    genes_union <- c(genes_union, g)
  }
  genes_union <- unique(genes_union)
  if (length(genes_union) < 1L) {
    return(list(error = "No overlapping genes between selected pathways and DEG table.", mat = NULL, sig_mat = NULL))
  }
  mat <- matrix(NA_real_, nrow = nrow(pathways), ncol = length(genes_union))
  rownames(mat) <- pathways$Description
  colnames(mat) <- genes_union
  row_map <- stats::setNames(seq_len(nrow(pathways)), pathways$ID)
  fc_map <- stats::setNames(tab$log2FC, tab$symbol)
  for (pid in names(gene_by_term)) {
    ridx <- row_map[[pid]]
    genes <- gene_by_term[[pid]]
    if (length(genes) < 1L) next
    mat[ridx, genes] <- unname(fc_map[genes])
  }
  keep_cols <- colSums(!is.na(mat)) > 0L
  mat <- mat[, keep_cols, drop = FALSE]
  if (ncol(mat) < 1L) {
    return(list(error = "No genes remained after overlap with selected comparison.", mat = NULL, sig_mat = NULL))
  }
  sig_cols <- colnames(mat) %in% sig_syms
  sig_mat <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
  rownames(sig_mat) <- rownames(mat)
  colnames(sig_mat) <- colnames(mat)
  if (any(sig_cols)) {
    for (j in which(sig_cols)) {
      sig_mat[, j] <- !is.na(mat[, j])
    }
  }
  list(error = NULL, mat = mat, sig_mat = sig_mat)
}

ora_heatmap_plot <- function(mat, title, x_lab, y_lab, sig_mat = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.ora_msg_plot("Install ggplot2 to render heatplot."))
  }
  if (is.null(mat) || nrow(mat) < 1L || ncol(mat) < 1L) {
    return(.ora_msg_plot("No matrix values to display."))
  }
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("row_key", "col_key", "value")
  row_levels <- rev(rownames(mat))
  col_levels <- colnames(mat)
  df$row_key <- factor(as.character(df$row_key), levels = row_levels)
  df$col_key <- factor(as.character(df$col_key), levels = col_levels)
  if (!is.null(sig_mat) && all(dim(sig_mat) == dim(mat))) {
    sig_df <- as.data.frame(as.table(sig_mat), stringsAsFactors = FALSE)
    colnames(sig_df) <- c("row_key", "col_key", "is_sig")
    df <- merge(df, sig_df, by = c("row_key", "col_key"), all.x = TRUE, sort = FALSE)
    df$is_sig <- as.logical(df$is_sig)
    df$is_sig[is.na(df$is_sig)] <- FALSE
  } else {
    df$is_sig <- FALSE
  }
  n_rows <- nrow(mat)
  # Slightly larger labels overall; big sets remain readable with increased height scaling.
  row_text_size <- 18.2 - 2.95 * log1p(max(0, n_rows - 1L) / 24)
  row_text_size <- min(18.2, max(8.6, row_text_size))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_key, y = row_key, fill = value)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "white",
      high = "#d7191c",
      midpoint = 0,
      na.value = "grey35",
      name = "log2FC"
    ) +
    ggplot2::labs(title = title, x = x_lab, y = y_lab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 11.5),
      axis.text.y = ggplot2::element_text(size = row_text_size),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right"
    )
  sig_df <- df[df$is_sig, , drop = FALSE]
  if (nrow(sig_df) > 0L) {
    p <- p + ggplot2::geom_text(
      data = sig_df,
      ggplot2::aes(x = col_key, y = row_key),
      label = "*",
      size = 4.2,
      color = "black",
      inherit.aes = FALSE
    )
  }
  p
}

ora_heatmap_plot_pathway_faceted <- function(mat, title, study_label_order = NULL, sig_mat = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.ora_msg_plot("Install ggplot2 to render heatplot."))
  }
  if (is.null(mat) || nrow(mat) < 1L || ncol(mat) < 1L) {
    return(.ora_msg_plot("No matrix values to display."))
  }
  cn <- colnames(mat)
  study_part <- character(length(cn))
  comp_part <- character(length(cn))
  for (i in seq_along(cn)) {
    nm <- as.character(cn[[i]])
    pos <- regexpr("::", nm, fixed = TRUE)[1L]
    if (!is.na(pos) && pos > 0L) {
      study_part[[i]] <- substr(nm, 1L, pos - 1L)
      comp_part[[i]] <- substr(nm, pos + 2L, nchar(nm))
    } else {
      study_part[[i]] <- nm
      comp_part[[i]] <- nm
    }
  }
  if (is.null(study_label_order) || length(study_label_order) < 1L) {
    study_label_order <- unique(study_part)
  }
  present_studies <- unique(study_part)
  study_label_order <- study_label_order[study_label_order %in% present_studies]
  if (length(study_label_order) < 1L) {
    study_label_order <- present_studies
  }

  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("row_key", "col_key", "value")
  row_levels <- rev(rownames(mat))
  df$row_key <- factor(as.character(df$row_key), levels = row_levels)
  if (!is.null(sig_mat) && all(dim(sig_mat) == dim(mat))) {
    sig_df <- as.data.frame(as.table(sig_mat), stringsAsFactors = FALSE)
    colnames(sig_df) <- c("row_key", "col_key", "is_sig")
    df <- merge(df, sig_df, by = c("row_key", "col_key"), all.x = TRUE, sort = FALSE)
    df$is_sig <- as.logical(df$is_sig)
    df$is_sig[is.na(df$is_sig)] <- FALSE
  } else {
    df$is_sig <- FALSE
  }
  key_df <- data.frame(
    col_key = cn,
    study = study_part,
    comparison = comp_part,
    stringsAsFactors = FALSE
  )
  df <- merge(df, key_df, by = "col_key", all.x = TRUE, sort = FALSE)
  df$study <- factor(as.character(df$study), levels = study_label_order)
  # Keep the current matrix column order inside each study facet.
  df$comparison <- factor(as.character(df$comparison), levels = unique(comp_part))

  n_rows <- nrow(mat)
  row_text_size <- 18.2 - 2.95 * log1p(max(0, n_rows - 1L) / 24)
  row_text_size <- min(18.2, max(8.6, row_text_size))

  facet_obj <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    ggh4x::facet_grid2(
      cols = ggplot2::vars(study),
      scales = "free_x",
      space = "free_x",
      axes = "margins",
      remove_labels = "none",
      drop = TRUE
    )
  } else {
    ggplot2::facet_grid(
      cols = ggplot2::vars(study),
      scales = "free_x",
      space = "free_x",
      drop = TRUE
    )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = comparison, y = row_key, fill = value)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "white",
      high = "#d7191c",
      midpoint = 0,
      na.value = "grey35",
      name = "log2FC"
    ) +
    facet_obj +
    ggplot2::labs(title = title, x = "Comparison (DEG list)", y = "Pathway genes") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      strip.text = ggplot2::element_text(face = "bold", size = 12.5),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 11.5),
      axis.text.y = ggplot2::element_text(size = row_text_size),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right"
    )
  sig_df <- df[df$is_sig, , drop = FALSE]
  if (nrow(sig_df) > 0L) {
    p <- p + ggplot2::geom_text(
      data = sig_df,
      ggplot2::aes(x = comparison, y = row_key),
      label = "*",
      size = 4.2,
      color = "black",
      inherit.aes = FALSE
    )
  }
  p
}

oraTabUI <- function(id, study_ids, study_labels, ora_file_choices, pathway_default = "") {
  ns <- shiny::NS(id)
  choices_with_empty <- ora_file_choices
  if (!("" %in% unname(choices_with_empty))) {
    choices_with_empty <- c("Select a pathway database..." = "", choices_with_empty)
  }
  sel <- ""
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::fluidRow(
          shiny::column(
            width = 7,
            shiny::conditionalPanel(
              condition = "output.ora_hide_pathway_select == '0'",
              ns = ns,
              shiny::tags$div(
                class = "ora-pathway-db-select",
                style = "max-width: 280px;",
                shiny::selectInput(
                  ns("pathway_file"),
                  "Pathway database:",
                  choices = choices_with_empty,
                  selected = sel,
                  width = "100%"
                )
              )
            ),
            shiny::conditionalPanel(
              condition = "output.ora_hide_pathway_select == '1'",
              ns = ns,
              shiny::tags$div(
                class = "text-muted",
                style = "font-size: 0.9rem; padding-top: 28px;",
                "Pathway database selection is disabled while custom ontology is loaded."
              )
            )
          ),
          shiny::column(
            width = 5,
            shiny::uiOutput(ns("ora_pathway_custom_msg"))
          )
        ),
        shiny::tags$hr(),
        shiny::tags$div(
          class = "text-muted",
          style = "font-size: 0.9rem; margin-bottom: 8px;",
          "Each study’s precomputed table is read from the path in ",
          shiny::tags$code("config.yaml"),
          " (key ",
          shiny::tags$code("ora_file"),
          ", default ",
          shiny::tags$code("ora/enrichment.rds"),
          "). The app uses that RDS when it contains the selected pathway database (same version as the app); otherwise it falls back to legacy ",
          shiny::tags$code("data/ora_cache/<pathway_basename>.rds"),
          ", then runs ",
          shiny::tags$code("enricher"),
          ". Sidebar thresholds only ",
          shiny::tags$strong("filter"),
          " cached results. A progress bar appears only while loading or computing enrichment."
        )
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::uiOutput(ns("ora_dotplot_status")),
        shiny::tags$div(
          class = "ora-plot-wrap",
          style = "overflow-x: auto; width: 100%; max-width: 100%; min-width: 0;",
          shiny::tags$div(
            style = "display: inline-block; vertical-align: top; max-width: none;",
            shiny::uiOutput(ns("ora_plot_all_ui"))
          )
        ),
        shiny::tags$hr(),
        shiny::fluidRow(
          shiny::column(
            width = 12,
            shiny::tags$div(
              class = "text-muted",
              style = "font-size: 0.95rem;",
              shiny::textOutput(ns("ora_heatmap_selection_msg"), inline = FALSE)
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 12,
            shiny::tags$div(
              class = "ora-heatmap-wrap",
              style = "overflow-x: auto; width: 100%; max-width: 100%; min-width: 0;",
              shiny::tags$div(
                style = "display: inline-block; vertical-align: top; max-width: none;",
                shiny::uiOutput(ns("ora_heatmap_plot_ui"))
              )
            )
          )
        )
      )
    )
  )
}

#' @param ora_input reactive: list(min_overlap, min_count, min_gene_ratio, show_category)
#' @param copy_genes_trigger reactive trigger for "copy visible genes" action
#' @param on_gene_select callback(symbol) for opening Gene tab and selecting symbol
#' @param custom_ontology reactive returning `list(active = logical, t2g = data.frame(term, gene) | NULL)`
oraTabServer <- function(
    id,
    study_ids,
    study_labels,
    ora_input,
    copy_genes_trigger = NULL,
    on_gene_select = NULL,
    custom_ontology = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    .ora_debug_log <- function(...) {
      message("[ORA debug] ", paste0(..., collapse = ""))
    }

    .ora_pick_discrete_index <- function(v, n) {
      if (is.null(v) || is.na(v)) return(NA_integer_)
      idx <- suppressWarnings(as.integer(round(as.numeric(v))))
      if (is.na(idx) || idx < 1L || idx > n) return(NA_integer_)
      idx
    }
    output$ora_pathway_custom_msg <- shiny::renderUI({
      co <- if (is.function(custom_ontology)) custom_ontology() else NULL
      if (!is.null(co) && is.list(co) && isTRUE(co$active) && is.data.frame(co$t2g) && nrow(co$t2g) > 0L) {
        shiny::tags$div(
          style = "color:#c0392b;font-weight:600;padding-top:26px;",
          "Custom ontology loaded"
        )
      } else {
        NULL
      }
    })
    output$ora_hide_pathway_select <- shiny::renderText({
      co <- if (is.function(custom_ontology)) custom_ontology() else NULL
      use_custom <- !is.null(co) && is.list(co) && isTRUE(co$active) && is.data.frame(co$t2g) && nrow(co$t2g) > 0L
      if (use_custom) "1" else "0"
    })
    shiny::outputOptions(output, "ora_hide_pathway_select", suspendWhenHidden = FALSE)

    # Heavy: cache load or enricher — depends on pathway file and optional custom ontology.
    ora_long_source <- shiny::reactive({
      co <- if (is.function(custom_ontology)) custom_ontology() else NULL
      use_custom <- !is.null(co) && is.list(co) && isTRUE(co$active) && is.data.frame(co$t2g) && nrow(co$t2g) > 0L
      if (!is.null(co) && is.list(co) && isTRUE(co$active) && !use_custom) {
        .ora_debug_log("custom ontology active flag is TRUE but t2g is empty/invalid")
      }
      if (use_custom) {
        n_cat <- length(unique(as.character(co$t2g$term)))
        n_pairs <- nrow(co$t2g)
        n_genes <- length(unique(as.character(co$t2g$gene)))
        .ora_debug_log(
          "custom ontology active: categories=", n_cat,
          " pairs=", n_pairs,
          " unique_genes=", n_genes
        )
        if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
          return(list(
            error = paste0(
              "Install Bioconductor package: BiocManager::install(\"clusterProfiler\")"
            ),
            long_by_sid = NULL
          ))
        }
        return(shiny::withProgress(message = "Running ORA on custom ontology…", value = 0, {
          n_steps <- 0L
          for (sid2 in study_ids) {
            n_steps <- n_steps + length(study_deg_lists(sid2))
          }
          den <- max(1L, n_steps)
          run <- ora_run_enrichment_long_by_study_t2g(
            study_ids,
            study_labels,
            co$t2g,
            pathway_label = ORA_CUSTOM_ONTOLOGY_KEY,
            min_gs_size = 1L,
            max_gs_size = 50000L,
            progress = function(amount, detail) {
              shiny::incProgress(amount, detail = detail)
            }
          )
          if (!is.null(run$error)) {
            .ora_debug_log("custom ontology run error: ", as.character(run$error))
            return(list(error = run$error, long_by_sid = NULL))
          }
          has_rows <- FALSE
          row_summ <- character(0)
          if (!is.null(run$long_by_sid) && length(run$long_by_sid) > 0L) {
            for (sid in names(run$long_by_sid)) {
              d0 <- run$long_by_sid[[sid]]
              n0 <- if (!is.null(d0) && is.data.frame(d0)) nrow(d0) else 0L
              row_summ <- c(row_summ, paste0(sid, "=", n0))
              if (!is.null(d0) && is.data.frame(d0) && nrow(d0) > 0L) {
                has_rows <- TRUE
                break
              }
            }
          }
          .ora_debug_log("custom ontology run rows by study: ", paste(row_summ, collapse = ", "))
          if (!has_rows) {
            return(list(
              error = paste0(
                "Custom ontology loaded but enrichment returned no rows. ",
                "Likely reasons: no overlap between ontology symbols and significant DEGs ",
                "(from config thresholds), or overlap exists but enrichment is not significant."
              ),
              long_by_sid = NULL
            ))
          }
          list(error = NULL, long_by_sid = run$long_by_sid)
        }))
      }

      pathway_file <- if (is.null(input$pathway_file)) "" else as.character(input$pathway_file)
      .ora_debug_log("pathway file mode: pathway_file=", pathway_file)
      if (!nzchar(pathway_file)) {
        return(list(
          error = "Select a pathway database or load a custom ontology (.xlsx) from the sidebar.",
          long_by_sid = NULL
        ))
      }

      per <- ora_try_load_per_study_caches(study_ids, pathway_file)
      if (isTRUE(per$ok_all)) {
        return(list(error = NULL, long_by_sid = per$long_by_sid))
      }

      cache_path <- ora_cache_rds_path(pathway_file)
      use_cache <- !is.na(cache_path) && nzchar(cache_path) && file.exists(cache_path)

      if (use_cache) {
        cached <- tryCatch(ora_read_cache(cache_path), error = function(e) NULL)
        if (!is.null(cached) && cache_compatible_with_app(cached, pathway_file, study_ids)) {
          return(list(
            error = NULL,
            long_by_sid = ora_cache_subset_studies(cached$long_by_sid, study_ids)
          ))
        }
      }

      pathway_abs <- file.path("databases", "pathways", pathway_file)
      if (!file.exists(pathway_abs)) {
        return(list(
          error = paste0(
            "No cached ORA RDS found for selected pathway database and pathway file is missing: ",
            pathway_file,
            ". Precompute enrichment.rds or add the pathway .txt file."
          ),
          long_by_sid = NULL
        ))
      }

      if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        return(list(
          error = paste0(
            "Install Bioconductor package: BiocManager::install(\"clusterProfiler\")"
          ),
          long_by_sid = NULL
        ))
      }

      shiny::withProgress(message = "Running ORA (overrepresentation analysis)…", value = 0, {
        n_steps <- 0L
        for (sid2 in study_ids) {
          n_steps <- n_steps + length(study_deg_lists(sid2))
        }
        den <- max(1L, n_steps)

        run <- ora_run_enrichment_long_by_study(
          study_ids,
          study_labels,
          pathway_file,
          min_gs_size = 1L,
          max_gs_size = 50000L,
          progress = function(amount, detail) {
            shiny::incProgress(amount, detail = detail)
          }
        )
        if (!is.null(run$error)) {
          return(list(error = run$error, long_by_sid = NULL))
        }

        list(error = NULL, long_by_sid = run$long_by_sid)
      })
    })

    # Light: filter + top-N — re-runs when sidebar ORA thresholds change.
    ora_by_study <- shiny::reactive({
      shiny::req(ora_input())
      oi <- ora_input()

      min_ol <- suppressWarnings(as.integer(oi$min_overlap))
      if (is.na(min_ol) || min_ol < 1L) min_ol <- 1L
      min_ct <- suppressWarnings(as.integer(oi$min_count))
      if (is.na(min_ct) || min_ct < 1L) min_ct <- 1L
      min_gr <- suppressWarnings(as.numeric(oi$min_gene_ratio))
      if (is.na(min_gr) || min_gr < 0) min_gr <- 0
      if (min_gr > 1) min_gr <- 1
      max_p_adj <- suppressWarnings(as.numeric(oi$max_p_adj))
      if (is.na(max_p_adj) || max_p_adj <= 0) max_p_adj <- 1
      if (max_p_adj > 1) max_p_adj <- 1
      n_show <- suppressWarnings(as.integer(oi$show_category))
      if (is.na(n_show) || n_show < 1L) n_show <- 20L
      n_show <- min(100L, max(1L, n_show))

      src <- ora_long_source()
      if (!is.null(src$error)) {
        .ora_debug_log("ora_long_source error: ", as.character(src$error))
        return(list(error = src$error, by_study = NULL, pathway_level_order = NULL))
      }
      ob <- ora_build_plot_payload(src$long_by_sid, study_ids, study_labels, min_ol, min_ct, min_gr, n_show, max_p_adj = max_p_adj)
      vis_summ <- character(0)
      for (sid in study_ids) {
        b <- ob$by_study[[sid]]
        nvis <- if (!is.null(b) && !is.null(b$plot_df) && is.data.frame(b$plot_df)) nrow(b$plot_df) else 0L
        vis_summ <- c(vis_summ, paste0(sid, "=", nvis))
      }
      .ora_debug_log(
        "plot payload after filters: min_count=", min_ct,
        " min_overlap=", min_ol,
        " min_gene_ratio=", format(min_gr, digits = 3),
        " max_p_adj=", format(max_p_adj, digits = 3),
        " n_show=", n_show,
        " rows_by_study=", paste(vis_summ, collapse = ", ")
      )
      ob
    })

    study_label_for <- function(sid) {
      lbl <- study_labels[[sid]]
      if (is.null(lbl) || !nzchar(as.character(lbl))) as.character(sid) else as.character(lbl)
    }

    output$ora_dotplot_status <- shiny::renderUI({
      d <- tryCatch(ora_combined_plot_df(), error = function(e) NULL)
      if (is.null(d)) return(NULL)
      co <- if (is.function(custom_ontology)) custom_ontology() else NULL
      use_custom <- !is.null(co) && is.list(co) && isTRUE(co$active) && is.data.frame(co$t2g) && nrow(co$t2g) > 0L
      if (!is.null(d$error) && nzchar(as.character(d$error))) {
        hint <- if (use_custom) {
          shiny::tags$div(
            style = "margin-top:10px; font-size:0.9rem; line-height:1.35; color:#333;",
            shiny::tags$p(
              style = "margin:0 0 6px 0;",
              shiny::tags$strong("If you expected results: "),
              "the run may have completed with ",
              shiny::tags$strong("no overlapping enrichment"),
              " (no categories pass statistical overlap), or ",
              shiny::tags$strong("sidebar filters removed all rows"),
              " (try lowering minimum overlap count, pathway size, and gene ratio)."
            ),
            shiny::tags$p(
              style = "margin:0;",
              shiny::tags$strong("Note: "),
              "ORA uses significant genes from each study’s ",
              shiny::tags$code("config.yaml"),
              " thresholds — ",
              shiny::tags$strong("not"),
              " the DEGs-tab sliders. Symbols in the .xlsx must match DEG gene symbols (matching is case-insensitive)."
            )
          )
        } else {
          NULL
        }
        shiny::tags$div(
          class = "alert alert-warning",
          style = "margin-bottom:12px;",
          shiny::tags$div(
            style = "font-weight:700; margin-bottom:6px;",
            "No ORA dot plot to show"
          ),
          shiny::tags$div(
            style = "white-space:pre-wrap;",
            htmltools::htmlEscape(as.character(d$error), attribute = FALSE)
          ),
          hint
        )
      } else {
        NULL
      }
    })

    ora_combined_plot_df <- shiny::reactive({
      ob <- ora_by_study()
      if (!is.null(ob$error)) {
        return(list(error = ob$error, combined = NULL, ob = ob))
      }
      dfs <- list()
      first_msg <- NULL
      for (sid in study_ids) {
        b <- ob$by_study[[sid]]
        if (is.null(b)) next
        if (!is.null(b$msg)) {
          if (is.null(first_msg)) first_msg <- b$msg
          next
        }
        if (!is.null(b$plot_df) && nrow(b$plot_df) > 0L) {
          df_sid <- b$plot_df
          df_sid$study_id <- sid
          dfs[[length(dfs) + 1L]] <- df_sid
        }
      }
      if (length(dfs) == 0L) {
        return(list(
          error = if (!is.null(first_msg)) first_msg else "No pathways to display.",
          combined = NULL,
          ob = ob
        ))
      }
      combined <- do.call(rbind, dfs)
      if (!is.null(ob$pathway_level_order) && length(ob$pathway_level_order) > 0L) {
        combined$Description <- factor(
          as.character(combined$Description),
          levels = ob$pathway_level_order
        )
      }
      list(error = NULL, combined = combined, ob = ob)
    })

    output$ora_plot_all_ui <- shiny::renderUI({
      d <- ora_combined_plot_df()
      comp_per_panel <- integer(0)
      max_path <- 1L
      if (!is.null(d$combined) && nrow(d$combined) > 0L) {
        for (sid in study_ids) {
          sub <- d$combined[as.character(d$combined$study_id) == sid, , drop = FALSE]
          if (nrow(sub) < 1L) next
          comp_per_panel <- c(comp_per_panel, length(unique(as.character(sub$comparison))))
          max_path <- max(max_path, length(unique(as.character(sub$ID))))
        }
      }
      dims <- if (length(comp_per_panel) > 0L) {
        .ora_faceted_plot_dims(comp_per_panel, max_path)
      } else {
        .ora_msg_plot_px()
      }
      if (requireNamespace("ggiraph", quietly = TRUE)) {
        ggiraph::girafeOutput(session$ns("ora_plot_all"), width = paste0(dims$width, "px"), height = paste0(dims$height, "px"))
      } else {
        shiny::plotOutput(session$ns("ora_plot_all_fallback"), width = paste0(dims$width, "px"), height = paste0(dims$height, "px"))
      }
    })

    if (requireNamespace("ggiraph", quietly = TRUE)) {
      output$ora_plot_all <- ggiraph::renderGirafe({
        d <- ora_combined_plot_df()
        if (!is.null(d$error)) {
          gp <- .ora_msg_plot(d$error)
        } else {
          all_lbls <- vapply(study_ids, study_label_for, character(1L))
          present <- unique(as.character(d$combined$study_label))
          study_label_levels <- all_lbls[all_lbls %in% present]
          gp <- tryCatch(
            .ora_faceted_comparison_plot_girafe(d$combined, study_label_levels, d$ob$pathway_level_order),
            error = function(err) {
              .ora_msg_plot(paste("Plot:", conditionMessage(err)))
            }
          )
        }
        if (is.null(gp)) {
          gp <- .ora_msg_plot("ggiraph is unavailable.")
        }
        ggiraph::girafe(
          ggobj = gp,
          width_svg = 14,
          height_svg = 8,
          options = list(
            ggiraph::opts_selection(
              type = "single",
              only_shiny = TRUE,
              css = "stroke:#39ff14;fill:#39ff14;stroke-width:2.2px;"
            ),
            ggiraph::opts_hover(css = "stroke:#000;stroke-width:1.2px;"),
            ggiraph::opts_sizing(rescale = FALSE)
          )
        )
      })
    }

    output$ora_plot_all_fallback <- shiny::renderPlot(
      {
        d <- ora_combined_plot_df()
        if (!is.null(d$error)) {
          return(.ora_msg_plot(d$error))
        }
        combined <- d$combined
        all_lbls <- vapply(study_ids, study_label_for, character(1L))
        present <- unique(as.character(combined$study_label))
        study_label_levels <- all_lbls[all_lbls %in% present]

        tryCatch(
          .ora_faceted_comparison_plot(combined, study_label_levels, d$ob$pathway_level_order),
          error = function(err) {
            .ora_msg_plot(paste("Plot:", conditionMessage(err)))
          }
        )
      }
    )

    selected_point <- shiny::reactive({
      id <- input$ora_plot_all_selected
      if (is.null(id) || !nzchar(as.character(id))) return(NULL)
      parts <- strsplit(as.character(id), "\\|\\|\\|", fixed = FALSE)[[1L]]
      if (length(parts) < 5L) return(NULL)
      list(
        study_id = parts[[1L]],
        study_label = parts[[2L]],
        comparison = parts[[3L]],
        pathway_id = parts[[4L]],
        pathway_desc = parts[[5L]]
      )
    })

    output$ora_heatmap_selection_msg <- shiny::renderText({
      sp <- selected_point()
      if (is.null(sp)) {
        return("Click a dot in the ORA panel to generate the heatplot.")
      }
      oi <- ora_input()
      mode <- if (!is.null(oi$click_target) && identical(oi$click_target, "comparison")) "comparison" else "pathway"
      if (identical(mode, "pathway")) {
        paste0(
          "Selected pathway: ", sp$pathway_desc,
          " (", sp$pathway_id, "). Heatmap uses pathway genes x all study comparisons."
        )
      } else {
        paste0(
          "Selected comparison: ", sp$study_label, "::", sp$comparison,
          ". Heatmap uses pathways x genes with log2FC values."
        )
      }
    })

    ora_heatmap_state <- shiny::reactive({
      sp <- selected_point()
      if (is.null(sp)) {
        return(list(
          plot = .ora_msg_plot("Click a dot in the ORA panel to generate the heatplot."),
          mode = "none",
          n_rows = 0L,
          n_cols = 0L,
          mat = NULL
        ))
      }
      pathway_file <- if (is.null(input$pathway_file)) "" else as.character(input$pathway_file)
      co <- if (is.function(custom_ontology)) custom_ontology() else NULL
      term2gene_custom <- if (!is.null(co) && is.list(co) && isTRUE(co$active) && is.data.frame(co$t2g) && nrow(co$t2g) > 0L) {
        co$t2g
      } else {
        NULL
      }
      pathway_file_heat <- if (!is.null(term2gene_custom)) "" else pathway_file
      if (!nzchar(pathway_file_heat) && is.null(term2gene_custom)) {
        return(list(
          plot = .ora_msg_plot("Select a pathway database or load a custom ontology to enable heatplot."),
          mode = "none",
          n_rows = 0L,
          n_cols = 0L,
          mat = NULL
        ))
      }
      oi <- ora_input()
      mode <- if (!is.null(oi$click_target) && identical(oi$click_target, "comparison")) "comparison" else "pathway"
      if (identical(mode, "pathway")) {
        d <- ora_combined_plot_df()
        if (!is.null(d$error)) {
          return(list(plot = .ora_msg_plot(d$error), mode = "none", n_rows = 0L, n_cols = 0L, mat = NULL))
        }
        pm <- ora_build_pathway_mode_matrix(
          pathway_file_heat,
          sp$pathway_id,
          sp$pathway_desc,
          study_ids,
          study_labels,
          d$combined,
          hide_empty_comparisons = isTRUE(oi$hide_empty_pathway_comparisons),
          term2gene_df = term2gene_custom
        )
        if (!is.null(pm$error)) {
          return(list(plot = .ora_msg_plot(pm$error), mode = "none", n_rows = 0L, n_cols = 0L, mat = NULL))
        }
        return(list(
          plot = ora_heatmap_plot_pathway_faceted(
            pm$mat,
            title = paste0("Pathway heatplot: ", sp$pathway_desc),
            study_label_order = vapply(study_ids, study_label_for, character(1L)),
            sig_mat = pm$sig_mat
          ),
          mode = "pathway",
          n_rows = nrow(pm$mat),
          n_cols = ncol(pm$mat),
          mat = pm$mat
        ))
      }
      d <- ora_combined_plot_df()
      if (!is.null(d$error)) {
        return(list(plot = .ora_msg_plot(d$error), mode = "none", n_rows = 0L, n_cols = 0L, mat = NULL))
      }
      cm <- ora_build_comparison_mode_matrix(
        pathway_file_heat,
        sp$study_id,
        sp$comparison,
        d$combined,
        term2gene_df = term2gene_custom
      )
      if (!is.null(cm$error)) {
        return(list(plot = .ora_msg_plot(cm$error), mode = "none", n_rows = 0L, n_cols = 0L, mat = NULL))
      }
      list(
        plot = ora_heatmap_plot(
          cm$mat,
          title = paste0("Comparison heatplot: ", sp$study_label, "::", sp$comparison),
          x_lab = "Genes",
          y_lab = "Pathways",
          sig_mat = cm$sig_mat
        ),
        mode = "comparison",
        n_rows = nrow(cm$mat),
        n_cols = ncol(cm$mat),
        mat = cm$mat
      )
    })

    output$ora_heatmap_plot_ui <- shiny::renderUI({
      st <- ora_heatmap_state()
      w <- .ora_heatmap_width_px(st$n_cols)
      shiny::plotOutput(
        session$ns("ora_heatmap_plot"),
        width = paste0(w, "px"),
        height = "auto",
        click = session$ns("ora_heatmap_plot_click")
      )
    })

    output$ora_heatmap_plot <- shiny::renderPlot({
      st <- ora_heatmap_state()
      st$plot
    }, height = function() {
      st <- ora_heatmap_state()
      if (identical(st$mode, "pathway")) {
        return(.ora_row_heatmap_height_px(st$n_rows))
      }
      if (identical(st$mode, "comparison")) {
        return(.ora_comparison_heatmap_height_px(st$n_rows))
      }
      420L
    })

    shiny::observeEvent(input$ora_heatmap_plot_click, {
      if (is.null(on_gene_select) || !is.function(on_gene_select)) return()
      st <- ora_heatmap_state()
      if (is.null(st$mat) || !is.matrix(st$mat)) return()
      click <- input$ora_heatmap_plot_click
      if (is.null(click)) return()

      gene_symbol <- NULL
      if (identical(st$mode, "pathway")) {
        rows <- rev(rownames(st$mat))
        if (!is.null(rows) && length(rows) > 0L) {
          iy <- .ora_pick_discrete_index(click$y, length(rows))
          if (!is.na(iy)) gene_symbol <- rows[[iy]]
        }
      } else if (identical(st$mode, "comparison")) {
        cols <- colnames(st$mat)
        if (!is.null(cols) && length(cols) > 0L) {
          ix <- .ora_pick_discrete_index(click$x, length(cols))
          if (!is.na(ix)) gene_symbol <- cols[[ix]]
        }
      }

      if (!is.null(gene_symbol) && nzchar(trimws(as.character(gene_symbol)))) {
        on_gene_select(as.character(gene_symbol))
      }
    })

    shiny::observeEvent(copy_genes_trigger(), {
      st <- ora_heatmap_state()
      if (is.null(st$mat) || !is.matrix(st$mat)) {
        shiny::showNotification("Generate ORA heatmap first, then copy genes.", type = "warning")
        return()
      }
      genes <- if (identical(st$mode, "pathway")) rownames(st$mat) else if (identical(st$mode, "comparison")) colnames(st$mat) else character(0)
      genes <- unique(trimws(as.character(genes)))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      if (length(genes) < 1L) {
        shiny::showNotification("No visible genes available to copy.", type = "warning")
        return()
      }
      txt <- paste(genes, collapse = "\n")
      session$sendCustomMessage("exprs_copy_text", txt)
      shiny::showNotification(paste0("Copied ", length(genes), " gene names."), type = "message")
    }, ignoreInit = TRUE)
  })
}
