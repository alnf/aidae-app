pheno_colors <- function(levels) {
  palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")
  setNames(palette[seq_along(levels) %% length(palette) + 1L], levels)
}

# Defaults match make_heatmap() formal arguments (used when study config has no thresholds).
default_heatmap_thresholds <- function() {
  list(fdr = 0.05, base_mean = 0, log2fc = 1, svalue = 0.005)
}

# Row indices into `res` for genes that pass numeric thresholds and finite z-scores (same rules as make_heatmap).
filter_heatmap_row_index <- function(res, mm, fdr, base_mean, log2fc, svalue) {
  if (is.null(res) || is.null(mm) || nrow(res) == 0L) return(NULL)
  gene_id_col <- if ("ens_gene" %in% colnames(res)) "ens_gene" else "symbol"
  gene_ids <- as.character(res[[gene_id_col]])
  mm <- mm[gene_ids, , drop = FALSE]
  l <- res$padj <= fdr & abs(res$log2FoldChange) >= log2fc
  if ("baseMean" %in% colnames(res)) l <- l & res$baseMean >= base_mean
  if ("svalue" %in% colnames(res)) l <- l & res$svalue <= svalue
  l[is.na(l)] <- FALSE
  if (sum(l) == 0L) return(NULL)
  m <- mm[l, , drop = FALSE]
  row_index <- which(l)
  m_z <- t(scale(t(m)))
  keep_row <- rowSums(!is.finite(m_z)) == 0L
  row_index <- row_index[keep_row]
  if (length(row_index) == 0L) return(NULL)
  row_index
}

make_heatmap <- function(res, mm, fdr = 0.05, base_mean = 0, log2fc = 1, svalue = 0.005, col_annot = NULL,
                         show_row_names = FALSE) {
  if (is.null(res) || is.null(mm) || nrow(res) == 0L) return(NULL)
  row_index <- filter_heatmap_row_index(res, mm, fdr, base_mean, log2fc, svalue)
  if (is.null(row_index)) return(NULL)
  gene_id_col <- if ("ens_gene" %in% colnames(res)) "ens_gene" else "symbol"
  gene_ids <- as.character(res[[gene_id_col]])
  mm <- mm[gene_ids, , drop = FALSE]
  m <- mm[row_index, , drop = FALSE]
  m_z <- t(scale(t(m)))
  if (nrow(m_z) == 0L) return(NULL)
  # Prepare row labels (symbols) without touching matrix rownames
  row_lab <- NULL
  if ("symbol" %in% colnames(res)) {
    row_lab <- res$symbol[row_index]
    print(row_lab)
  }
  n_row <- nrow(m_z)
  n_col <- ncol(m_z)
  message("Heatmap dimensions: n_row=", n_row, " n_col=", n_col)
  top_anno <- NULL
  if (!is.null(col_annot) && "PhenoNames" %in% colnames(col_annot) && nrow(col_annot) == n_col) {
    pheno <- col_annot$PhenoNames
    col_list <- list(PhenoNames = pheno_colors(unique(pheno)))
    top_anno <- HeatmapAnnotation(PhenoNames = pheno, col = col_list, show_legend = TRUE)
  }
  row_km_arg <- if (n_row >= 3L && n_col >= 2L) 2L else NULL
  has_ragg <- requireNamespace("ragg", quietly = TRUE)
  raster_device <- if (has_ragg) "agg_png" else "png"
  ht <- Heatmap(
    m_z, name = "z-score",
    show_row_names = show_row_names, show_column_names = FALSE,
    row_names_side = "left",
    row_labels = row_lab,
    row_km = row_km_arg, show_row_dend = !is.null(row_km_arg),
    column_title = paste0(n_row, " significant genes with FDR < ", fdr),
    top_annotation = top_anno,
    use_raster = TRUE,
    raster_device = raster_device,
    raster_quality = 1
  )
  if ("baseMean" %in% colnames(res)) {
    ht <- ht + Heatmap(
      log10(res$baseMean[row_index] + 1), show_row_names = FALSE, width = unit(5, "mm"),
      name = "log10(baseMean+1)", show_column_names = FALSE,
      use_raster = TRUE
    )
  }
  ht <- ht + Heatmap(
    res$log2FoldChange[row_index], show_row_names = FALSE, width = unit(5, "mm"),
    name = "log2FoldChange", show_column_names = FALSE,
    use_raster = TRUE,
    col = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
  )
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  ht <- draw(ht, merge_legend = TRUE)
  list(ht = ht, row_index = row_index)
}

make_maplot <- function(res, highlight = NULL) {
  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if (!is.null(highlight)) {
    col[highlight] <- "red"
    cex[highlight] <- 1
  }
  x <- if ("baseMean" %in% colnames(res)) res$baseMean else seq_len(nrow(res))
  y <- res$log2FoldChange
  y[y > 2] <- 2
  y[y < -2] <- -2
  col[col == "red" & y < 0] <- "darkgreen"
  par(mar = c(4, 4, 1, 1))
  xlab <- if ("baseMean" %in% colnames(res)) "baseMean" else "rank"
  log_x <- if ("baseMean" %in% colnames(res)) "x" else ""

  suppressWarnings(
    plot(
      x, y, col = col,
      pch = ifelse(res$log2FoldChange > 2 | res$log2FoldChange < -2, 1, 16),
      cex = cex, log = log_x,
      xlab = xlab, ylab = "log2 fold change"
    )
  )
}

make_volcano <- function(res, highlight = NULL) {
  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if (!is.null(highlight)) {
    col[highlight] <- "red"
    cex[highlight] <- 1
  }
  x <- res$log2FoldChange
  y <- -log10(res$padj)
  col[col == "red" & x < 0] <- "darkgreen"
  par(mar = c(4, 4, 1, 1))

  suppressWarnings(
    plot(
      x, y, col = col,
      pch = 16,
      cex = cex,
      xlab = "log2 fold change", ylab = "-log10(FDR)"
    )
  )
}

