# ORA tab: clusterProfiler::enricher, one faceted ggplot (facet = study; x = DEG lists; color = GeneRatio).
# Depends: scripts/ora_cache.R, scripts/heatmap_utils.R, scripts/study_data.R, scripts/pathway_signatures.R;
#   Bioconductor clusterProfiler; ggplot2. Optional: ggh4x (shared y-axis labels on outer margin only).

# Compact dimensions for message-only ggplot.
.ora_msg_plot_px <- function() {
  list(width = 380L, height = 420L)
}

.ora_msg_plot <- function(msg) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }
  msg_txt <- as.character(msg)
  msg_wrapped <- paste(strwrap(msg_txt, width = 48), collapse = "\n")
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.01,
      y = 0.5,
      hjust = 0,
      vjust = 0.5,
      label = msg_wrapped,
      size = 4.6
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
}

# Pixel size for Shiny renderPlot: total width = sum of per-study panel widths (matches facet_grid space = "free_x").
.ora_wrap_strip_text <- function(labels, width = 9L) {
  labs <- as.character(labels)
  vapply(labs, function(lbl) {
    parts <- strwrap(lbl, width = max(8L, as.integer(width)))
    if (length(parts) < 1L) lbl else paste(parts, collapse = "\n")
  }, character(1L))
}

.ora_strip_text_sizes <- function(labels, base_size = 12.5, comp_per_panel = NULL) {
  labs <- as.character(labels)
  labs <- labs[!is.na(labs) & nzchar(labs)]
  if (length(labs) < 1L) return(numeric(0))
  cps <- suppressWarnings(as.integer(comp_per_panel))
  cps <- cps[!is.na(cps) & cps >= 1L]
  cps_default <- if (length(cps) > 0L) stats::median(cps) else 3
  cps_by_label <- if (!is.null(names(comp_per_panel)) && length(comp_per_panel) > 0L) {
    stats::setNames(suppressWarnings(as.integer(comp_per_panel)), names(comp_per_panel))
  } else {
    NULL
  }
  out <- vapply(labs, function(lbl) {
    chars <- nchar(lbl, type = "width")
    n_comp <- cps_default
    if (!is.null(cps_by_label) && lbl %in% names(cps_by_label)) {
      n_comp <- cps_by_label[[lbl]]
      if (is.na(n_comp) || n_comp < 1L) n_comp <- cps_default
    }
    # Per-facet adaptation: long labels and narrow facets get smaller text.
    shrink_label <- max(0, chars - 14) * 0.14
    shrink_narrow <- max(0, 3L - n_comp) * 0.45
    max(9.2, base_size - shrink_label - shrink_narrow)
  }, numeric(1L))
  stats::setNames(out, labs)
}

.ora_strip_text_size <- function(labels, base_size = 12.5, comp_per_panel = NULL) {
  sizes <- .ora_strip_text_sizes(labels, base_size = base_size, comp_per_panel = comp_per_panel)
  if (length(sizes) < 1L) return(base_size)
  # Robust global fallback: avoid over-shrinking all strips because of one narrow facet.
  as.numeric(stats::quantile(sizes, probs = 0.85, na.rm = TRUE, type = 7))
}

.ora_facet_x_expand_add <- function(comp_per_panel) {
  cps <- suppressWarnings(as.integer(comp_per_panel))
  cps <- cps[!is.na(cps) & cps >= 1L]
  if (length(cps) < 1L) return(0.7)
  if (min(cps) <= 2L) return(0.95)
  if (min(cps) == 3L) return(0.8)
  0.65
}

.ora_add_facet_spacers <- function(plot_df, min_slots = 3L) {
  if (is.null(plot_df) || nrow(plot_df) < 1L) return(plot_df)
  if (!("study_label" %in% colnames(plot_df)) || !("comparison" %in% colnames(plot_df))) return(plot_df)
  min_slots <- max(1L, as.integer(min_slots))
  out <- plot_df
  out$comparison <- as.character(out$comparison)
  out$comparison_display <- out$comparison
  out$comparison_slot <- NA_character_
  studies <- unique(as.character(plot_df$study_label))
  for (s in studies) {
    idx <- which(as.character(out$study_label) == s)
    sub <- out[idx, , drop = FALSE]
    cmp_vals <- unique(as.character(sub$comparison))
    n_cmp <- length(cmp_vals)
    slots_total <- max(min_slots, n_cmp)
    slot_keys <- paste0("..slot__", s, "__", seq_len(slots_total))
    slot_idx <- if (n_cmp == 1L) {
      # Keep a single real comparison centered in padded facets.
      as.integer(ceiling(slots_total / 2))
    } else {
      # Spread real comparisons across available slots (e.g., 2 in 3 -> slots 1 and 3).
      as.integer(round(seq(1, slots_total, length.out = n_cmp)))
    }
    cmp_to_slot <- stats::setNames(slot_keys[slot_idx], cmp_vals)
    out$comparison_slot[idx] <- unname(cmp_to_slot[as.character(out$comparison[idx])])
    n_add <- slots_total - n_cmp
    if (n_add < 1L) next
    # One y-level anchor keeps spacers in the facet x scale while producing no visible points.
    y_anchor <- as.character(sub$Description[[1L]])
    if (!nzchar(y_anchor)) y_anchor <- as.character(sub$Description[which(nzchar(as.character(sub$Description)))[1L]])
    if (!nzchar(y_anchor)) next
    spacer_slots <- setdiff(slot_keys, unname(cmp_to_slot))
    for (k in seq_len(length(spacer_slots))) {
      r <- sub[1L, , drop = FALSE]
      r$comparison <- as.character(sub$comparison[[1L]])
      r$comparison_display <- ""
      r$comparison_slot <- spacer_slots[[k]]
      r$Description <- y_anchor
      if ("Count" %in% colnames(r)) r$Count <- NA_real_
      if ("gene_ratio" %in% colnames(r)) r$gene_ratio <- NA_real_
      if ("selection_id" %in% colnames(r)) r$selection_id <- paste0("..spacer__", s, "__", k)
      if ("tooltip" %in% colnames(r)) r$tooltip <- ""
      out <- rbind(out, r)
    }
  }
  na_slots <- is.na(out$comparison_slot) | !nzchar(as.character(out$comparison_slot))
  if (any(na_slots)) {
    out$comparison_slot[na_slots] <- paste0("..slot__", as.character(out$study_label[na_slots]), "__fallback")
  }
  out
}

# Build x-axis labels from comparison slot -> displayed label.
.ora_comparison_axis_labels <- function(x, label_map) {
  xx <- as.character(x)
  out <- unname(label_map[xx])
  out[is.na(out)] <- ""
  out
}

.ora_assay_palette <- function(assay_types) {
  ats <- as.character(assay_types)
  ats <- ats[!is.na(ats) & nzchar(ats)]
  if (length(ats) < 1L) return(character(0))
  ats_u <- unique(ats)
  out <- stats::setNames(rep("#D9D9D9", length(ats_u)), ats_u) # Unknown / unspecified
  key <- toupper(gsub("[^A-Z0-9]", "", ats_u))
  is_rna <- grepl("RNA", key)
  is_ms <- grepl("MS", key) & !is_rna
  out[is_ms] <- "#F4C2D7"
  out[is_rna] <- "#BFE3FF"
  out
}

.ora_study_assay_types <- function(study_ids, study_labels) {
  if (length(study_ids) < 1L) return(character(0))
  out <- character(0)
  for (sid in study_ids) {
    cfg_path <- file.path("data", sid, "config.yaml")
    if (!file.exists(cfg_path)) next
    cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
    if (is.null(cfg)) next
    lbl <- study_labels[[sid]]
    if (is.null(lbl) || !nzchar(as.character(lbl))) lbl <- sid
    at <- cfg$assay_type
    if (is.null(at) || !nzchar(as.character(at))) at <- "Unknown"
    out[[as.character(lbl)]] <- as.character(at)
  }
  out
}

# Pixel size for Shiny renderPlot: total width = sum of per-study panel widths (matches facet_grid space = "free_x").
.ora_faceted_plot_dims <- function(comp_per_panel, max_path, panel_labels = NULL) {
  np <- max(1L, as.integer(max_path))
  if (length(comp_per_panel) < 1L) {
    return(list(width = 520L, height = 340L))
  }
  w_body <- 0L
  for (i in seq_along(comp_per_panel)) {
    nc <- comp_per_panel[[i]]
    nc <- max(1L, as.integer(nc))
    # Concave non-linear growth:
    # - noticeably wider when nc is small (1-3 comparisons)
    # - slower growth for larger nc to avoid over-expanding dense facets
    panel_width <- max(138L, round(84 + 106 * log1p(nc)))
    if (!is.null(panel_labels) && length(panel_labels) >= i) {
      lbl <- as.character(panel_labels[[i]])
      if (!is.na(lbl) && nzchar(lbl)) {
        # Ensure narrow facets can still host long study names in strip labels.
        chars <- nchar(lbl, type = "width")
        wrapped_lines <- max(1L, ceiling(chars / 9))
        # Reserve width for long labels, but avoid exploding total plot width.
        panel_width <- max(panel_width, 74L + 7L * min(34L, chars) + 7L * (wrapped_lines - 1L))
      }
    }
    w_body <- w_body + panel_width
  }
  # Stronger compression when many facets are present.
  n_facets <- length(comp_per_panel)
  compress <- if (n_facets <= 2L) 1 else max(0.64, 1 - 0.085 * (n_facets - 2L))
  w <- min(4200L, max(420L, 130L + round(w_body * compress)))
  # Sublinear height scaling to avoid excessive trailing whitespace for large pathway sets.
  npf <- as.numeric(np)
  h <- 320 + 26 * (npf^0.88) + 7 * log1p(npf)
  h <- min(5200L, max(420L, round(h)))
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
.ora_faceted_comparison_plot <- function(combined_df, study_label_levels, pathway_level_order = NULL, assay_type_by_label = NULL) {
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
  comp_per_panel <- vapply(
    split(as.character(combined_df$comparison), as.character(combined_df$study_label)),
    function(x) length(unique(x)),
    integer(1L)
  )
  plot_df$comparison <- as.character(plot_df$comparison)
  plot_df <- .ora_add_facet_spacers(plot_df, min_slots = 3L)
  cmp_slot_map <- stats::setNames(
    as.character(plot_df$comparison_display),
    as.character(plot_df$comparison_slot)
  )

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
    strip_sizes <- .ora_strip_text_sizes(study_label_levels, base_size = 12.5, comp_per_panel = comp_per_panel)
    assay_vals <- if (!is.null(assay_type_by_label) && length(assay_type_by_label) > 0L) {
      unname(assay_type_by_label[study_label_levels])
    } else {
      rep("Unknown", length(study_label_levels))
    }
    assay_vals[is.na(assay_vals) | !nzchar(assay_vals)] <- "Unknown"
    assay_cols <- .ora_assay_palette(assay_vals)
    strip_fills <- unname(assay_cols[assay_vals])
    strip_fills[is.na(strip_fills)] <- "#e3e3e3"
    strip_obj <- ggh4x::strip_themed(
      text_x = ggh4x::elem_list_text(
        face = rep("bold", length(strip_sizes)),
        size = as.numeric(strip_sizes),
        lineheight = rep(0.95, length(strip_sizes))
      ),
      background_x = ggh4x::elem_list_rect(
        fill = strip_fills,
        colour = rep("grey55", length(strip_fills))
      )
    )
    ggh4x::facet_grid2(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study_label = function(x) .ora_wrap_strip_text(x, width = 9L)),
      strip = strip_obj,
      axes = "margins",
      remove_labels = "none",
      drop = FALSE
    )
  } else {
    ggplot2::facet_grid(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study_label = function(x) .ora_wrap_strip_text(x, width = 9L)),
      drop = FALSE
    )
  }
  strip_size <- .ora_strip_text_size(study_label_levels, base_size = 12.5, comp_per_panel = comp_per_panel)
  x_expand_add <- .ora_facet_x_expand_add(comp_per_panel)

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = comparison_slot,
    y = Description,
    color = gene_ratio,
    size = Count
  )) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(mult = 0, add = x_expand_add),
      labels = function(x) .ora_comparison_axis_labels(x, cmp_slot_map)
    ) +
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
      strip.text = ggplot2::element_text(face = "bold", size = strip_size, lineheight = 0.95),
      panel.spacing.x = grid::unit(12, "pt"),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 11),
      axis.text.y = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 10.5),
      legend.title = ggplot2::element_text(size = 11)
    )
}

.ora_faceted_comparison_plot_girafe <- function(combined_df, study_label_levels, pathway_level_order = NULL, assay_type_by_label = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggiraph", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(combined_df) || nrow(combined_df) < 1L) {
    return(.ora_msg_plot("No pathways to display."))
  }
  plot_df <- combined_df
  comp_per_panel <- vapply(
    split(as.character(combined_df$comparison), as.character(combined_df$study_label)),
    function(x) length(unique(x)),
    integer(1L)
  )
  plot_df$study_label <- factor(
    as.character(plot_df$study_label),
    levels = study_label_levels
  )
  plot_df$comparison <- as.character(plot_df$comparison)
  plot_df <- .ora_add_facet_spacers(plot_df, min_slots = 3L)
  cmp_slot_map <- stats::setNames(
    as.character(plot_df$comparison_display),
    as.character(plot_df$comparison_slot)
  )
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
  p_raw <- if ("pvalue" %in% colnames(plot_df)) {
    as.numeric(plot_df$pvalue)
  } else if ("p_value" %in% colnames(plot_df)) {
    as.numeric(plot_df$p_value)
  } else {
    rep(NA_real_, nrow(plot_df))
  }
  p_adj <- if ("p_adj" %in% colnames(plot_df)) {
    as.numeric(plot_df$p_adj)
  } else if ("p.adjust" %in% colnames(plot_df)) {
    as.numeric(plot_df[["p.adjust"]])
  } else if ("padj" %in% colnames(plot_df)) {
    as.numeric(plot_df$padj)
  } else {
    rep(NA_real_, nrow(plot_df))
  }
  p_raw_lbl <- ifelse(is.na(p_raw), "not available in cache", formatC(p_raw, digits = 3, format = "f"))
  p_adj_lbl <- ifelse(is.na(p_adj), "NA", formatC(p_adj, digits = 3, format = "f"))

  plot_df$tooltip <- paste0(
    "<b>Study:</b> ", as.character(plot_df$study_label),
    "<br><b>Comparison:</b> ", as.character(plot_df$comparison_display),
    "<br><b>Pathway:</b> ", as.character(plot_df$Description),
    "<br><b>Count:</b> ", as.character(plot_df$Count),
    "<br><b>Gene ratio:</b> ", formatC(as.numeric(plot_df$gene_ratio), digits = 3, format = "f"),
    "<br><b>P-value:</b> ", p_raw_lbl,
    "<br><b>Adjusted p-value:</b> ", p_adj_lbl
  )

  facet_obj <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    strip_sizes <- .ora_strip_text_sizes(study_label_levels, base_size = 12.5, comp_per_panel = comp_per_panel)
    assay_vals <- if (!is.null(assay_type_by_label) && length(assay_type_by_label) > 0L) {
      unname(assay_type_by_label[study_label_levels])
    } else {
      rep("Unknown", length(study_label_levels))
    }
    assay_vals[is.na(assay_vals) | !nzchar(assay_vals)] <- "Unknown"
    assay_cols <- .ora_assay_palette(assay_vals)
    strip_fills <- unname(assay_cols[assay_vals])
    strip_fills[is.na(strip_fills)] <- "#e3e3e3"
    strip_obj <- ggh4x::strip_themed(
      text_x = ggh4x::elem_list_text(
        face = rep("bold", length(strip_sizes)),
        size = as.numeric(strip_sizes),
        lineheight = rep(0.95, length(strip_sizes))
      ),
      background_x = ggh4x::elem_list_rect(
        fill = strip_fills,
        colour = rep("grey55", length(strip_fills))
      )
    )
    ggh4x::facet_grid2(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study_label = function(x) .ora_wrap_strip_text(x, width = 9L)),
      strip = strip_obj,
      axes = "margins",
      remove_labels = "none",
      drop = FALSE
    )
  } else {
    ggplot2::facet_grid(
      cols = ggplot2::vars(study_label),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study_label = function(x) .ora_wrap_strip_text(x, width = 9L)),
      drop = FALSE
    )
  }
  strip_size <- .ora_strip_text_size(study_label_levels, base_size = 12.5, comp_per_panel = comp_per_panel)
  x_expand_add <- .ora_facet_x_expand_add(comp_per_panel)

  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = comparison_slot,
    y = Description,
    color = gene_ratio,
    size = Count
  )) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(data_id = selection_id, tooltip = tooltip),
      na.rm = TRUE
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(mult = 0, add = x_expand_add),
      labels = function(x) .ora_comparison_axis_labels(x, cmp_slot_map)
    ) +
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
      strip.text = ggplot2::element_text(face = "bold", size = strip_size, lineheight = 0.95),
      panel.spacing.x = grid::unit(12, "pt"),
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

#' Build a minimal term–gene table from ORA `long_df` rows (e.g. precomputed RDS)
#' when `databases/pathways/<pathway_file>` is absent (thin Connect deploy).
ora_t2g_from_ora_subframe <- function(sub) {
  if (is.null(sub) || !is.data.frame(sub) || nrow(sub) < 1L) {
    return(NULL)
  }
  if (!"geneID" %in% colnames(sub) || !"ID" %in% colnames(sub)) {
    return(NULL)
  }
  ids <- as.character(sub$ID)
  genes_split <- lapply(seq_len(nrow(sub)), function(i) ora_parse_gene_set(sub$geneID[i]))
  if (all(lengths(genes_split) == 0L)) {
    return(NULL)
  }
  term <- rep(ids, lengths(genes_split))
  gene <- unlist(genes_split, use.names = FALSE)
  data.frame(
    term = term,
    gene = toupper(trimws(as.character(gene))),
    stringsAsFactors = FALSE
  )
}

ora_extract_p_adj <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 1L) return(rep(NA_real_, 0L))
  if ("p_adj" %in% colnames(df)) return(as.numeric(df$p_adj))
  if ("p.adjust" %in% colnames(df)) return(as.numeric(df[["p.adjust"]]))
  if ("padj" %in% colnames(df)) return(as.numeric(df$padj))
  rep(NA_real_, nrow(df))
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
    pabs <- file.path("databases", "pathways", pathway_file)
    if (file.exists(pabs)) {
      parse_pathway_file_to_term2gene(pathway_file)
    } else {
      NULL
    }
  }
  if (is.null(t2g) || nrow(t2g) < 1L) return(character(0))
  term <- as.character(t2g$term)
  keep <- term %in% c(as.character(pathway_id), as.character(pathway_desc))
  genes <- toupper(trimws(as.character(t2g$gene[keep])))
  unique(genes[!is.na(genes) & nzchar(genes)])
}

ora_read_deg_logfc <- function(study_id, comparison_label, deg_file = NULL) {
  .pick_num_col <- function(df, candidates) {
    for (nm in candidates) {
      if (nm %in% colnames(df)) return(as.numeric(df[[nm]]))
    }
    rep(NA_real_, nrow(df))
  }
  gdegs <- load_gene_tab_gdegs(study_id)
  if (!is.null(gdegs) && is.data.frame(gdegs) &&
      all(c("symbol", "log2FC") %in% colnames(gdegs))) {
    cmp_col <- if ("joint" %in% colnames(gdegs)) "joint" else if ("comp" %in% colnames(gdegs)) "comp" else NULL
    if (!is.null(cmp_col)) {
      gsub <- gdegs[as.character(gdegs[[cmp_col]]) == as.character(comparison_label), , drop = FALSE]
      if (nrow(gsub) > 0L) {
        p_raw_vec <- .pick_num_col(gsub, c("pvalue", "p_value", "p.value", "pval", "P.Value"))
        p_adj_vec <- .pick_num_col(gsub, c("padj", "p.adjust", "p_adj", "adj_p_val", "FDR", "qvalue", "q.value", "pvalue"))
        out <- data.frame(
          symbol = toupper(trimws(as.character(gsub$symbol))),
          log2FC = as.numeric(gsub$log2FC),
          p_value = p_raw_vec,
          p_adj = p_adj_vec,
          stringsAsFactors = FALSE
        )
        out <- out[!is.na(out$symbol) & nzchar(out$symbol), , drop = FALSE]
        out <- out[!is.na(out$log2FC), , drop = FALSE]
        if (nrow(out) > 0L) {
          # gdegs_long can omit raw p-values; in that case fall back to DEG file
          # (when available) so heatmap tooltips/markers can use p-value.
          has_raw <- any(!is.na(as.numeric(out$p_value)))
          if (has_raw || is.null(deg_file) || !nzchar(as.character(deg_file))) {
            return(out)
          }
        }
      }
    }
  }

  if (is.null(deg_file) || !nzchar(as.character(deg_file))) return(NULL)
  loaded <- load_study_data(study_id, deg_file)
  res <- loaded$res
  if (is.null(res) || nrow(res) < 1L || !all(c("symbol", "log2FoldChange") %in% colnames(res))) {
    return(NULL)
  }
  p_adj_vec <- .pick_num_col(res, c("padj", "p.adjust", "p_adj", "adj_p_val", "FDR", "qvalue", "q.value", "pvalue"))
  p_raw_vec <- .pick_num_col(res, c("pvalue", "p_value", "p.value", "pval", "P.Value"))
  out <- data.frame(
    symbol = toupper(trimws(as.character(res$symbol))),
    log2FC = as.numeric(res$log2FoldChange),
    p_value = p_raw_vec,
    p_adj = p_adj_vec,
    stringsAsFactors = FALSE
  )
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
  t2g_fallback <- NULL
  if ((is.null(term2gene_df) || nrow(term2gene_df) < 1L) &&
        !file.exists(file.path("databases", "pathways", pathway_file)) &&
        !is.null(sub) && nrow(sub) > 0L) {
    t2g_fallback <- ora_t2g_from_ora_subframe(sub)
  }
  t2g_for_terms <- if (!is.null(term2gene_df) && is.data.frame(term2gene_df) && nrow(term2gene_df) > 0L) {
    term2gene_df
  } else {
    t2g_fallback
  }
  if (length(pathway_genes) < 1L && !is.null(sub) && nrow(sub) > 0L) {
    term_genes <- ora_term_genes(pathway_file, pathway_id, pathway_desc, term2gene_df = t2g_for_terms)
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
    pathway_genes <- ora_pathway_genes(pathway_file, pathway_id, pathway_desc, term2gene_df = t2g_for_terms)
  }
  if (length(pathway_genes) < 1L) {
    return(list(error = "No genes found for selected pathway.", mat = NULL))
  }

  vals <- list()
  sig_vals <- list()
  padj_vals <- list()
  pval_vals <- list()
  p_adj_threshold <- 0.05
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
      p_adj_vec <- rep(NA_real_, length(pathway_genes))
      p_raw_vec <- rep(NA_real_, length(pathway_genes))
      if ("p_adj" %in% colnames(tab)) {
        p_adj_vec[ok] <- as.numeric(tab$p_adj[hit[ok]])
      }
      if ("p_value" %in% colnames(tab)) {
        p_raw_vec[ok] <- as.numeric(tab$p_value[hit[ok]])
      }
      padj_vals[[col_name]] <- p_adj_vec
      pval_vals[[col_name]] <- p_raw_vec
      sig_vals[[col_name]] <- !is.na(x) & is.finite(p_adj_vec) & (p_adj_vec < p_adj_threshold)
    }
  }
  if (length(vals) < 1L) {
    return(list(error = "No DEG tables available for selected pathway.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  mat <- do.call(cbind, vals)
  sig_mat <- do.call(cbind, sig_vals)
  padj_mat <- do.call(cbind, padj_vals)
  pval_mat <- do.call(cbind, pval_vals)
  rownames(mat) <- pathway_genes
  rownames(sig_mat) <- pathway_genes
  rownames(padj_mat) <- pathway_genes
  rownames(pval_mat) <- pathway_genes
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
        padj_mat <- padj_mat[, keep_cols, drop = FALSE]
        pval_mat <- pval_mat[, keep_cols, drop = FALSE]
        reorder_names <- keep_names[keep_names %in% colnames(mat)]
        if (length(reorder_names) > 0L) {
          mat <- mat[, reorder_names, drop = FALSE]
          sig_mat <- sig_mat[, reorder_names, drop = FALSE]
          padj_mat <- padj_mat[, reorder_names, drop = FALSE]
          pval_mat <- pval_mat[, reorder_names, drop = FALSE]
        }
      }
    }
  }
  if (ncol(mat) < 1L) {
    return(list(error = "No comparisons remained after pathway-significance filter.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  keep_rows <- rowSums(!is.na(mat)) > 0L
  mat <- mat[keep_rows, , drop = FALSE]
  sig_mat <- sig_mat[keep_rows, , drop = FALSE]
  padj_mat <- padj_mat[keep_rows, , drop = FALSE]
  pval_mat <- pval_mat[keep_rows, , drop = FALSE]
  if (nrow(mat) < 1L) {
    return(list(error = "Selected pathway genes are absent from all DEG tables.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  list(error = NULL, mat = mat, sig_mat = sig_mat, padj_mat = padj_mat, pval_mat = pval_mat)
}

ora_build_comparison_mode_matrix <- function(
    pathway_file,
    selected_sid,
    selected_comparison,
    combined_df,
    term2gene_df = NULL) {
  if (is.null(combined_df) || nrow(combined_df) < 1L) {
    return(list(error = "No ORA rows available for selected comparison.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  sub <- combined_df[
    as.character(combined_df$study_id) == as.character(selected_sid) &
      as.character(combined_df$comparison) == as.character(selected_comparison),
    ,
    drop = FALSE
  ]
  if (nrow(sub) < 1L) {
    return(list(error = "Selected comparison has no pathways in current ORA view.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  lists <- study_deg_lists(selected_sid)
  idx <- which(vapply(lists, function(e) {
    identical(as.character(ora_deg_entry_comparison_label(e)), as.character(selected_comparison))
  }, logical(1L)))
  if (length(idx) != 1L) {
    return(list(error = "Could not map selected comparison to one DEG list.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  tab <- ora_read_deg_logfc(selected_sid, selected_comparison, lists[[idx]]$deg_file)
  if (is.null(tab) || nrow(tab) < 1L) {
    return(list(error = "Selected comparison DEG table is empty.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  pathway_abs <- file.path("databases", "pathways", pathway_file)
  t2g <- if (!is.null(term2gene_df) && is.data.frame(term2gene_df) && nrow(term2gene_df) > 0L) {
    term2gene_df
  } else if (file.exists(pathway_abs)) {
    parse_pathway_file_to_term2gene(pathway_file)
  } else {
    ora_t2g_from_ora_subframe(sub)
  }
  if (is.null(t2g) || nrow(t2g) < 1L) {
    return(list(error = "Pathway database has no term-to-gene mappings.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
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
    return(list(error = "No overlapping genes between selected pathways and DEG table.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
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
    return(list(error = "No genes remained after overlap with selected comparison.", mat = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL))
  }
  p_adj_threshold <- 0.05
  sig_mat <- matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
  rownames(sig_mat) <- rownames(mat)
  colnames(sig_mat) <- colnames(mat)
  padj_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  rownames(padj_mat) <- rownames(mat)
  colnames(padj_mat) <- colnames(mat)
  pval_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  rownames(pval_mat) <- rownames(mat)
  colnames(pval_mat) <- colnames(mat)
  p_adj_by_gene <- stats::setNames(as.numeric(tab$p_adj), as.character(tab$symbol))
  p_raw_by_gene <- if ("p_value" %in% colnames(tab)) {
    stats::setNames(as.numeric(tab$p_value), as.character(tab$symbol))
  } else {
    NULL
  }
  for (j in seq_len(ncol(mat))) {
    gene <- colnames(mat)[[j]]
    pv <- suppressWarnings(as.numeric(p_adj_by_gene[[gene]]))
    pr <- if (!is.null(p_raw_by_gene)) suppressWarnings(as.numeric(p_raw_by_gene[[gene]])) else NA_real_
    padj_mat[, j] <- pv
    pval_mat[, j] <- pr
    sig_mat[, j] <- !is.na(mat[, j]) & is.finite(pv) & (pv < p_adj_threshold)
  }
  list(error = NULL, mat = mat, sig_mat = sig_mat, padj_mat = padj_mat, pval_mat = pval_mat)
}

ora_heatmap_df_with_sig_padj <- function(mat, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL) {
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
  if (!is.null(padj_mat) && all(dim(padj_mat) == dim(mat))) {
    p_df <- as.data.frame(as.table(padj_mat), stringsAsFactors = FALSE)
    colnames(p_df) <- c("row_key", "col_key", "p_adj")
    df <- merge(df, p_df, by = c("row_key", "col_key"), all.x = TRUE, sort = FALSE)
    df$p_adj <- as.numeric(df$p_adj)
  } else {
    df$p_adj <- NA_real_
  }
  if (!is.null(pval_mat) && all(dim(pval_mat) == dim(mat))) {
    p_raw_df <- as.data.frame(as.table(pval_mat), stringsAsFactors = FALSE)
    colnames(p_raw_df) <- c("row_key", "col_key", "p_value")
    df <- merge(df, p_raw_df, by = c("row_key", "col_key"), all.x = TRUE, sort = FALSE)
    df$p_value <- as.numeric(df$p_value)
  } else {
    df$p_value <- NA_real_
  }
  df$marker <- ifelse(
    df$is_sig,
    "*",
    ifelse(!is.na(df$p_value) & is.finite(df$p_value) & (df$p_value <= 0.05), "\u00B7", "")
  )
  df
}

ora_heatmap_plot <- function(mat, title, x_lab, y_lab, sig_mat = NULL, pval_mat = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(.ora_msg_plot("Install ggplot2 to render heatplot."))
  }
  if (is.null(mat) || nrow(mat) < 1L || ncol(mat) < 1L) {
    return(.ora_msg_plot("No matrix values to display."))
  }
  df <- ora_heatmap_df_with_sig_padj(mat, sig_mat = sig_mat, padj_mat = NULL, pval_mat = pval_mat)
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
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 11.5),
      axis.text.y = ggplot2::element_text(size = row_text_size),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right"
    )
  mark_df <- df[nzchar(as.character(df$marker)), , drop = FALSE]
  if (nrow(mark_df) > 0L) {
    p <- p + ggplot2::geom_text(
      data = mark_df,
      ggplot2::aes(x = col_key, y = row_key, label = marker),
      size = 4.2,
      color = "black",
      inherit.aes = FALSE
    )
  }
  p
}

ora_heatmap_plot_girafe <- function(mat, title, x_lab, y_lab, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggiraph", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(mat) || nrow(mat) < 1L || ncol(mat) < 1L) {
    return(.ora_msg_plot("No matrix values to display."))
  }
  df <- ora_heatmap_df_with_sig_padj(mat, sig_mat = sig_mat, padj_mat = padj_mat, pval_mat = pval_mat)
  p_adj_lbl <- ifelse(is.na(df$p_adj), "NA", formatC(df$p_adj, digits = 3, format = "f"))
  logfc_lbl <- ifelse(is.na(df$value), "NA", formatC(df$value, digits = 3, format = "f"))
  df$tooltip <- paste0(
    "<b>Pathway:</b> ", as.character(df$row_key),
    "<br><b>Gene:</b> ", as.character(df$col_key),
    "<br><b>log2FC:</b> ", logfc_lbl,
    "<br><b>P-value:</b> ", ifelse(is.na(df$p_value), "NA", formatC(df$p_value, digits = 3, format = "f")),
    "<br><b>Adjusted p-value:</b> ", p_adj_lbl
  )
  df$data_id <- paste(as.character(df$row_key), as.character(df$col_key), sep = "|||")
  n_rows <- nrow(mat)
  row_text_size <- 18.2 - 2.95 * log1p(max(0, n_rows - 1L) / 24)
  row_text_size <- min(18.2, max(8.6, row_text_size))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_key, y = row_key, fill = value)) +
    ggiraph::geom_tile_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = data_id),
      color = "grey90",
      linewidth = 0.15
    ) +
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
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 11.5),
      axis.text.y = ggplot2::element_text(size = row_text_size),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right"
    )
  mark_df <- df[nzchar(as.character(df$marker)), , drop = FALSE]
  if (nrow(mark_df) > 0L) {
    p <- p + ggplot2::geom_text(
      data = mark_df,
      ggplot2::aes(x = col_key, y = row_key, label = marker),
      size = 4.2,
      color = "black",
      inherit.aes = FALSE
    )
  }
  p
}

ora_heatmap_plot_pathway_faceted <- function(mat, title, study_label_order = NULL, sig_mat = NULL, pval_mat = NULL, assay_type_by_label = NULL) {
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

  df <- ora_heatmap_df_with_sig_padj(mat, sig_mat = sig_mat, padj_mat = NULL, pval_mat = pval_mat)
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
    assay_vals <- if (!is.null(assay_type_by_label) && length(assay_type_by_label) > 0L) {
      unname(assay_type_by_label[study_label_order])
    } else {
      rep("Unknown", length(study_label_order))
    }
    assay_vals[is.na(assay_vals) | !nzchar(assay_vals)] <- "Unknown"
    assay_cols <- .ora_assay_palette(assay_vals)
    strip_fills <- unname(assay_cols[assay_vals])
    strip_fills[is.na(strip_fills)] <- "#e3e3e3"
    strip_obj <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = strip_fills,
        colour = rep("grey55", length(strip_fills))
      )
    )
    ggh4x::facet_grid2(
      cols = ggplot2::vars(study),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study = function(x) .ora_wrap_strip_text(x, width = 12L)),
      strip = strip_obj,
      axes = "margins",
      remove_labels = "none",
      drop = TRUE
    )
  } else {
    ggplot2::facet_grid(
      cols = ggplot2::vars(study),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study = function(x) .ora_wrap_strip_text(x, width = 12L)),
      drop = TRUE
    )
  }
  strip_size <- .ora_strip_text_size(levels(df$study), base_size = 12.5)

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
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      strip.text = ggplot2::element_text(face = "bold", size = strip_size, lineheight = 0.95),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 11.5),
      axis.text.y = ggplot2::element_text(size = row_text_size),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right"
    )
  mark_df <- df[nzchar(as.character(df$marker)), , drop = FALSE]
  if (nrow(mark_df) > 0L) {
    p <- p + ggplot2::geom_text(
      data = mark_df,
      ggplot2::aes(x = comparison, y = row_key, label = marker),
      size = 4.2,
      color = "black",
      inherit.aes = FALSE
    )
  }
  p
}

ora_heatmap_plot_pathway_faceted_girafe <- function(mat, title, study_label_order = NULL, sig_mat = NULL, padj_mat = NULL, pval_mat = NULL, assay_type_by_label = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggiraph", quietly = TRUE)) {
    return(NULL)
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
  df <- ora_heatmap_df_with_sig_padj(mat, sig_mat = sig_mat, padj_mat = padj_mat, pval_mat = pval_mat)
  key_df <- data.frame(
    col_key = cn,
    study = study_part,
    comparison = comp_part,
    stringsAsFactors = FALSE
  )
  df <- merge(df, key_df, by = "col_key", all.x = TRUE, sort = FALSE)
  df$study <- factor(as.character(df$study), levels = study_label_order)
  df$comparison <- factor(as.character(df$comparison), levels = unique(comp_part))
  p_adj_lbl <- ifelse(is.na(df$p_adj), "NA", formatC(df$p_adj, digits = 3, format = "f"))
  logfc_lbl <- ifelse(is.na(df$value), "NA", formatC(df$value, digits = 3, format = "f"))
  df$tooltip <- paste0(
    "<b>Study:</b> ", as.character(df$study),
    "<br><b>Comparison:</b> ", as.character(df$comparison),
    "<br><b>Gene:</b> ", as.character(df$row_key),
    "<br><b>log2FC:</b> ", logfc_lbl,
    "<br><b>P-value:</b> ", ifelse(is.na(df$p_value), "NA", formatC(df$p_value, digits = 3, format = "f")),
    "<br><b>Adjusted p-value:</b> ", p_adj_lbl
  )
  df$data_id <- paste(as.character(df$row_key), as.character(df$comparison), as.character(df$study), sep = "|||")
  n_rows <- nrow(mat)
  row_text_size <- 18.2 - 2.95 * log1p(max(0, n_rows - 1L) / 24)
  row_text_size <- min(18.2, max(8.6, row_text_size))
  facet_obj <- if (requireNamespace("ggh4x", quietly = TRUE)) {
    assay_vals <- if (!is.null(assay_type_by_label) && length(assay_type_by_label) > 0L) {
      unname(assay_type_by_label[study_label_order])
    } else {
      rep("Unknown", length(study_label_order))
    }
    assay_vals[is.na(assay_vals) | !nzchar(assay_vals)] <- "Unknown"
    assay_cols <- .ora_assay_palette(assay_vals)
    strip_fills <- unname(assay_cols[assay_vals])
    strip_fills[is.na(strip_fills)] <- "#e3e3e3"
    strip_obj <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(
        fill = strip_fills,
        colour = rep("grey55", length(strip_fills))
      )
    )
    ggh4x::facet_grid2(
      cols = ggplot2::vars(study),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study = function(x) .ora_wrap_strip_text(x, width = 12L)),
      strip = strip_obj,
      axes = "margins",
      remove_labels = "none",
      drop = TRUE
    )
  } else {
    ggplot2::facet_grid(
      cols = ggplot2::vars(study),
      scales = "free_x",
      space = "free_x",
      labeller = ggplot2::labeller(study = function(x) .ora_wrap_strip_text(x, width = 12L)),
      drop = TRUE
    )
  }
  strip_size <- .ora_strip_text_size(levels(df$study), base_size = 12.5)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = comparison, y = row_key, fill = value)) +
    ggiraph::geom_tile_interactive(
      ggplot2::aes(tooltip = tooltip, data_id = data_id),
      color = "grey90",
      linewidth = 0.15
    ) +
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
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      strip.text = ggplot2::element_text(face = "bold", size = strip_size, lineheight = 0.95),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, size = 11.5),
      axis.text.y = ggplot2::element_text(size = row_text_size),
      legend.title = ggplot2::element_text(size = 13),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right"
    )
  mark_df <- df[nzchar(as.character(df$marker)), , drop = FALSE]
  if (nrow(mark_df) > 0L) {
    p <- p + ggplot2::geom_text(
      data = mark_df,
      ggplot2::aes(x = comparison, y = row_key, label = marker),
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
            shiny::uiOutput(ns("ora_assay_legend_ui")),
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
    assay_type_by_label <- .ora_study_assay_types(study_ids, study_labels)
    output$ora_assay_legend_ui <- shiny::renderUI({
      if (length(assay_type_by_label) < 1L) return(NULL)
      assays <- unique(as.character(assay_type_by_label))
      assays <- assays[!is.na(assays) & nzchar(assays)]
      if (length(assays) < 1L) return(NULL)
      pal <- .ora_assay_palette(assays)
      shiny::tags$div(
        style = "display:flex;align-items:center;flex-wrap:wrap;gap:10px;margin:2px 0 10px 0;",
        shiny::tags$span(style = "font-size:0.92rem;font-weight:600;color:#555;", "Assay type:"),
        lapply(assays, function(a) {
          col <- pal[[a]]
          if (is.null(col) || !nzchar(col)) col <- "#D9D9D9"
          shiny::tags$span(
            style = "display:inline-flex;align-items:center;gap:6px;font-size:0.9rem;color:#333;",
            shiny::tags$span(
              style = paste0(
                "display:inline-block;width:12px;height:12px;border-radius:2px;",
                "background:", col, ";border:1px solid #8d8d8d;"
              )
            ),
            as.character(a)
          )
        })
      )
    })
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
      max_p_value <- suppressWarnings(as.numeric(oi$max_p_value))
      if (is.na(max_p_value) || max_p_value <= 0) max_p_value <- 1
      if (max_p_value > 1) max_p_value <- 1
      n_show <- suppressWarnings(as.integer(oi$show_category))
      if (is.na(n_show) || n_show < 1L) n_show <- 20L
      n_show <- min(100L, max(1L, n_show))

      src <- ora_long_source()
      if (!is.null(src$error)) {
        .ora_debug_log("ora_long_source error: ", as.character(src$error))
        return(list(error = src$error, by_study = NULL, pathway_level_order = NULL))
      }
      ob <- ora_build_plot_payload(
        src$long_by_sid,
        study_ids,
        study_labels,
        min_ol,
        min_ct,
        min_gr,
        n_show,
        max_p_value = max_p_value,
        max_p_adj = max_p_adj
      )
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
        " max_p=", format(max_p_value, digits = 3),
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
      # Defensive bind: studies can contribute slightly different optional columns
      # when caches were produced across app versions.
      all_cols <- character(0)
      for (d in dfs) all_cols <- union(all_cols, colnames(d))
      dfs_aligned <- lapply(dfs, function(d) {
        miss <- setdiff(all_cols, colnames(d))
        if (length(miss) > 0L) {
          for (nm in miss) d[[nm]] <- NA
        }
        d[, all_cols, drop = FALSE]
      })
      combined <- do.call(rbind, dfs_aligned)
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
      panel_labels <- character(0)
      max_path <- 1L
      if (!is.null(d$combined) && nrow(d$combined) > 0L) {
        for (sid in study_ids) {
          sub <- d$combined[as.character(d$combined$study_id) == sid, , drop = FALSE]
          if (nrow(sub) < 1L) next
          comp_per_panel <- c(comp_per_panel, max(3L, length(unique(as.character(sub$comparison)))))
          panel_labels <- c(panel_labels, unique(as.character(sub$study_label))[[1L]])
          max_path <- max(max_path, length(unique(as.character(sub$ID))))
        }
      }
      dims <- if (length(comp_per_panel) > 0L) {
        .ora_faceted_plot_dims(comp_per_panel, max_path, panel_labels = panel_labels)
      } else {
        .ora_msg_plot_px()
      }
      has_real_dotplot <- is.null(d$error) && !is.null(d$combined) && nrow(d$combined) > 0L
      if (requireNamespace("ggiraph", quietly = TRUE) && has_real_dotplot) {
        ggiraph::girafeOutput(session$ns("ora_plot_all"), width = paste0(dims$width, "px"), height = paste0(dims$height, "px"))
      } else {
        shiny::plotOutput(session$ns("ora_plot_all_fallback"), width = paste0(dims$width, "px"), height = paste0(dims$height, "px"))
      }
    })

    if (requireNamespace("ggiraph", quietly = TRUE)) {
      output$ora_plot_all <- ggiraph::renderGirafe({
        d <- ora_combined_plot_df()
        comp_per_panel <- integer(0)
        panel_labels <- character(0)
        max_path <- 1L
        if (!is.null(d$combined) && nrow(d$combined) > 0L) {
          for (sid in study_ids) {
            sub <- d$combined[as.character(d$combined$study_id) == sid, , drop = FALSE]
            if (nrow(sub) < 1L) next
            comp_per_panel <- c(comp_per_panel, max(3L, length(unique(as.character(sub$comparison)))))
            panel_labels <- c(panel_labels, unique(as.character(sub$study_label))[[1L]])
            max_path <- max(max_path, length(unique(as.character(sub$ID))))
          }
        }
        dims <- if (length(comp_per_panel) > 0L) {
          .ora_faceted_plot_dims(comp_per_panel, max_path, panel_labels = panel_labels)
        } else {
          .ora_msg_plot_px()
        }
        if (!is.null(d$error)) {
          gp <- .ora_msg_plot(d$error)
        } else {
          all_lbls <- vapply(study_ids, study_label_for, character(1L))
          present <- unique(as.character(d$combined$study_label))
          study_label_levels <- all_lbls[all_lbls %in% present]
          gp <- tryCatch(
            .ora_faceted_comparison_plot_girafe(
              d$combined,
              study_label_levels,
              d$ob$pathway_level_order,
              assay_type_by_label = assay_type_by_label
            ),
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
          width_svg = max(6, dims$width / 96),
          height_svg = max(4, dims$height / 96),
          options = list(
            ggiraph::opts_selection(
              type = "single",
              only_shiny = TRUE,
              css = "stroke:#39ff14;stroke-width:2.2px;"
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
          .ora_faceted_comparison_plot(
            combined,
            study_label_levels,
            d$ob$pathway_level_order,
            assay_type_by_label = assay_type_by_label
          ),
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
            sig_mat = pm$sig_mat,
            pval_mat = pm$pval_mat,
            assay_type_by_label = assay_type_by_label
          ),
          plot_girafe = ora_heatmap_plot_pathway_faceted_girafe(
            pm$mat,
            title = paste0("Pathway heatplot: ", sp$pathway_desc),
            study_label_order = vapply(study_ids, study_label_for, character(1L)),
            sig_mat = pm$sig_mat,
            padj_mat = pm$padj_mat,
            pval_mat = pm$pval_mat,
            assay_type_by_label = assay_type_by_label
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
          sig_mat = cm$sig_mat,
          pval_mat = cm$pval_mat
        ),
        plot_girafe = ora_heatmap_plot_girafe(
          cm$mat,
          title = paste0("Comparison heatplot: ", sp$study_label, "::", sp$comparison),
          x_lab = "Genes",
          y_lab = "Pathways",
          sig_mat = cm$sig_mat,
          padj_mat = cm$padj_mat,
          pval_mat = cm$pval_mat
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
      h <- if (identical(st$mode, "pathway")) {
        .ora_row_heatmap_height_px(st$n_rows)
      } else if (identical(st$mode, "comparison")) {
        .ora_comparison_heatmap_height_px(st$n_rows)
      } else {
        420L
      }
      has_real_heatmap <- !is.null(st$mat) && is.matrix(st$mat) && nrow(st$mat) > 0L && ncol(st$mat) > 0L
      if (requireNamespace("ggiraph", quietly = TRUE) && has_real_heatmap) {
        ggiraph::girafeOutput(
          session$ns("ora_heatmap_plot_girafe"),
          width = paste0(w, "px"),
          height = paste0(h, "px")
        )
      } else {
        shiny::plotOutput(
          session$ns("ora_heatmap_plot"),
          width = paste0(w, "px"),
          height = "auto",
          click = session$ns("ora_heatmap_plot_click")
        )
      }
    })

    if (requireNamespace("ggiraph", quietly = TRUE)) {
      output$ora_heatmap_plot_girafe <- ggiraph::renderGirafe({
        st <- ora_heatmap_state()
        w <- .ora_heatmap_width_px(st$n_cols)
        h <- if (identical(st$mode, "pathway")) {
          .ora_row_heatmap_height_px(st$n_rows)
        } else if (identical(st$mode, "comparison")) {
          .ora_comparison_heatmap_height_px(st$n_rows)
        } else {
          420L
        }
        gp <- st$plot_girafe
        if (is.null(gp)) {
          gp <- st$plot
        }
        ggiraph::girafe(
          ggobj = gp,
          width_svg = max(6, w / 96),
          height_svg = max(4, h / 96),
          options = list(
            ggiraph::opts_selection(
              type = "single",
              only_shiny = TRUE,
              css = "stroke:#39ff14;stroke-width:2.2px;"
            ),
            ggiraph::opts_sizing(rescale = FALSE),
            ggiraph::opts_hover(css = "stroke:#000;stroke-width:1px;")
          )
        )
      })
    }

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

    .ora_gene_from_girafe_selection <- function(sel_id, mode) {
      if (is.null(sel_id) || !nzchar(as.character(sel_id))) return(NULL)
      parts <- strsplit(as.character(sel_id), "\\|\\|\\|", fixed = FALSE)[[1L]]
      if (identical(mode, "comparison")) {
        # data_id = pathway|||gene
        if (length(parts) < 2L) return(NULL)
        return(as.character(parts[[2L]]))
      }
      if (identical(mode, "pathway")) {
        # data_id = gene|||comparison|||study
        if (length(parts) < 1L) return(NULL)
        return(as.character(parts[[1L]]))
      }
      NULL
    }

    shiny::observeEvent(input$ora_heatmap_plot_girafe_selected, {
      if (is.null(on_gene_select) || !is.function(on_gene_select)) return()
      st <- ora_heatmap_state()
      if (is.null(st$mat) || !is.matrix(st$mat)) return()
      gene_symbol <- .ora_gene_from_girafe_selection(input$ora_heatmap_plot_girafe_selected, st$mode)
      if (!is.null(gene_symbol) && nzchar(trimws(as.character(gene_symbol)))) {
        on_gene_select(as.character(gene_symbol))
      }
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
