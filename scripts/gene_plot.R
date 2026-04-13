# Gene-level boxplots for the Gene tab (one study).
# Requires: ggplot2, ggpubr, ggnewscale (same stack as old/app.R).
# Depends on pheno_colors() from heatmap_utils.R — load heatmap_utils before this file (see app.R).

# ggplot with only a message (no gene / missing packages).
.gene_plot_empty <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = msg,
      size = 4.2
    ) +
    ggplot2::theme_void() +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))
}

# Assign bracket y positions: stack by comparison within each facet (or globally).
.assign_pval_y <- function(df_counts, pv, facet_column) {
  pv$y.position <- NA_real_
  step <- 0.85
  pad <- 0.5
  comps_order <- unique(as.character(pv$comp))
  use_facet <- !is.null(facet_column) && nzchar(facet_column) &&
    facet_column %in% names(df_counts) && facet_column %in% names(pv)

  if (use_facet) {
    regions <- unique(as.character(df_counts[[facet_column]][!is.na(df_counts[[facet_column]])]))
    for (r in regions) {
      sub_r <- as.character(df_counts[[facet_column]]) == r
      maxy <- max(df_counts$Counts[sub_r], na.rm = TRUE)
      y0 <- maxy + pad
      for (cmp in comps_order) {
        w <- which(as.character(pv[[facet_column]]) == r & as.character(pv$comp) == cmp)
        if (length(w) == 0L) next
        pv$y.position[w] <- y0
        y0 <- y0 + step
      }
    }
  } else {
    maxy <- max(df_counts$Counts, na.rm = TRUE)
    y0 <- maxy + pad
    for (cmp in comps_order) {
      w <- which(as.character(pv$comp) == cmp)
      if (length(w) == 0L) next
      pv$y.position[w] <- y0
      y0 <- y0 + step
    }
  }
  pv
}

#' Boxplot of one gene across `PhenoNames`, optional DE brackets from `gdegs` long table.
#'
#' @param facet_column Optional: name of a **metadata** column to facet by; must match `gene_tab_facet` in study config and the extra column in `gdegs` when brackets are used. If NULL or not in `metadata`, no faceting.
#' @param ens_gene Row name in `mm` (Ensembl or matrix row ID).
#' @param mm Expression matrix: genes x samples (same as `load_gene_tab_study()$mm`).
#' @param metadata Samples aligned to `colnames(mm)`; must include `PhenoNames`, `SampleNumber`.
#' @param gene_symbol Optional display name.
#' @param pval_df Slice of global DEGs for this gene (same columns as `gdegs_long.tsv`), or NULL.
#' @param study_title Optional study label for the title.
#' @param mcols Named vector of colours for `PhenoNames` levels; if NULL, uses [pheno_colors()].
plot_gene_study <- function(
    ens_gene,
    mm,
    metadata,
    gene_symbol = NULL,
    pval_df = NULL,
    study_title = NULL,
    mcols = NULL,
    facet_column = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.gene_plot_empty("ggplot2 is required."))
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    return(.gene_plot_empty("Install ggpubr for gene plots: install.packages(\"ggpubr\")"))
  }
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    return(.gene_plot_empty("Install ggnewscale: install.packages(\"ggnewscale\")"))
  }

  if (is.null(mm) || nrow(mm) == 0L || ncol(mm) == 0L) {
    return(.gene_plot_empty("No expression data for this study."))
  }
  if (is.null(metadata) || nrow(metadata) == 0L) {
    return(.gene_plot_empty("No metadata for this study."))
  }
  if (!ens_gene %in% rownames(mm)) {
    return(.gene_plot_empty("Gene not found in this study."))
  }
  if (!all(c("PhenoNames", "SampleNumber") %in% colnames(metadata))) {
    return(.gene_plot_empty("Metadata must contain PhenoNames and SampleNumber."))
  }

  fc <- facet_column
  if (!is.null(fc)) {
    fc <- as.character(fc)
    if (length(fc) != 1L || !nzchar(fc)) {
      fc <- NULL
    } else if (!fc %in% colnames(metadata)) {
      return(.gene_plot_empty(paste0("gene_tab_facet column \"", fc, "\" not found in metadata.")))
    }
  }

  counts <- as.numeric(mm[ens_gene, , drop = TRUE])
  df <- data.frame(
    SampleNumber = colnames(mm),
    PhenoNames = metadata$PhenoNames,
    Counts = counts,
    stringsAsFactors = FALSE
  )
  if (!is.null(fc)) {
    df[[fc]] <- metadata[[fc]]
  }

  pn_order <- unique(as.character(df$PhenoNames))
  df$PhenoNames <- factor(df$PhenoNames, levels = pn_order)

  if (is.null(mcols)) {
    mcols <- pheno_colors(pn_order)
  } else {
    base <- pheno_colors(pn_order)
    nm_ok <- intersect(names(mcols), pn_order)
    if (length(nm_ok) > 0L) {
      base[nm_ok] <- mcols[nm_ok]
    }
    mcols <- base
  }

  has_facet <- !is.null(fc) && fc %in% names(df) && any(!is.na(df[[fc]]))

  ttl <- paste(
    c(
      if (!is.null(study_title) && nzchar(study_title)) study_title else NULL,
      if (!is.null(gene_symbol) && nzchar(gene_symbol)) gene_symbol else ens_gene
    ),
    collapse = " — "
  )

  p <- ggpubr::ggboxplot(
    df,
    x = "PhenoNames",
    y = "Counts",
    color = "PhenoNames",
    add = "jitter",
    title = ttl,
    ggtheme = ggplot2::theme_gray()
  ) +
    ggplot2::scale_color_manual(values = mcols, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.18))) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15, face = "bold"),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 13, colour = "gray15"),
      axis.text.x = ggplot2::element_text(size = 13, colour = "gray15"),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12.5),
      strip.text = ggplot2::element_text(size = 13, face = "bold")
    )

  if (has_facet) {
    p <- p + ggplot2::facet_wrap(as.formula(paste("~", fc)), nrow = 1L)
  }

  pv <- NULL
  if (!is.null(pval_df) && nrow(pval_df) > 0L) {
    pv <- pval_df
    if (!all(c("group1", "group2", "padj", "comp") %in% colnames(pv))) {
      pv <- NULL
    } else {
      pv <- pv[!is.na(pv$padj) & pv$padj <= 0.05, , drop = FALSE]
      levs <- levels(df$PhenoNames)
      pv$group1 <- as.character(pv$group1)
      pv$group2 <- as.character(pv$group2)
      pv <- pv[pv$group1 %in% levs & pv$group2 %in% levs, , drop = FALSE]
      if (nrow(pv) > 0L) {
        comps <- unique(as.character(pv$comp))
        pal <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")
        pv$color <- pal[(match(pv$comp, comps) - 1L) %% length(pal) + 1L]
        pv$padj_lab <- formatC(pv$padj, format = "e", digits = 2)
        if (has_facet && !is.null(fc) && fc %in% colnames(pv)) {
          facet_u <- unique(as.character(df[[fc]][!is.na(df[[fc]])]))
          keep <- !is.na(pv[[fc]]) & as.character(pv[[fc]]) %in% facet_u
          pv <- pv[keep, , drop = FALSE]
        }
        if (nrow(pv) > 0L) {
          use_facet_y <- has_facet && !is.null(fc) && fc %in% colnames(pv) && all(!is.na(pv[[fc]]))
          pv <- .assign_pval_y(df, pv, facet_column = if (use_facet_y) fc else NULL)
        }
      }
    }
  }

  if (!is.null(pv) && nrow(pv) > 0L && !any(is.na(pv$y.position))) {
    p <- p +
      ggnewscale::new_scale_color() +
      ggpubr::stat_pvalue_manual(
        pv,
        label = "padj_lab",
        tip.length = 0.03,
        color = "color",
        label.size = 4,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_color_identity()
  }

  p
}
