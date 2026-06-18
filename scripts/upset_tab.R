# UpSet tab: DEG list intersections (ComplexUpset + ggiraph matrix clicks).
# Depends: ComplexUpset, ggplot2, ggiraph (optional interactivity), yaml; study_data, heatmap_utils, ora_cache.
# Sidebar inputs live in app.R (`upset_*`).

`%||%` <- function(x, y) if (is.null(x)) y else x

# Baked into ggiraph data_id between study_id and deg_file (must not appear in those strings).
.upset_siddeg_data_id_sep <- "@@::@@"

# Patchwork row weight for the first annotation row (intersection-size bar). Matrix + set-size row
# stays at ComplexUpset `height_ratio` (1). SVG + Shiny output height scale by (w + rest) / (n_ann + hr)
# so the matrix keeps the same vertical share as default all-1 weights.
.upset_cu_intersection_bar_row_weight <- 0.55
.upset_cu_girafe_base_height_px <- 720L
.upset_cu_girafe_plot_height_px <- function() {
  w <- .upset_cu_intersection_bar_row_weight
  n_ann <- 1L
  hr <- 1
  round(
    as.numeric(.upset_cu_girafe_base_height_px) *
      (w + max(0L, n_ann - 1L) + hr) /
      (as.numeric(n_ann) + hr)
  )
}

.upset_assay_palette <- function(assay_types) {
  ats <- as.character(assay_types)
  ats <- ats[!is.na(ats) & nzchar(ats)]
  if (length(ats) < 1L) return(character(0))
  ats_u <- unique(ats)
  out <- stats::setNames(rep("#D9D9D9", length(ats_u)), ats_u)
  key <- toupper(gsub("[^A-Z0-9]", "", ats_u))
  is_rna <- grepl("RNA", key)
  is_ms <- grepl("MS", key) & !is_rna
  out[is_ms] <- "#F4C2D7"
  out[is_rna] <- "#BFE3FF"
  out
}

.upset_study_assay_type <- function(study_id) {
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return("Unknown")
  cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
  if (is.null(cfg)) return("Unknown")
  at <- cfg$assay_type
  if (is.null(at) || !nzchar(as.character(at))) "Unknown" else as.character(at)
}

.upset_merge_thr <- function(base, over) {
  if (is.null(over) || !is.list(over)) return(base)
  if (!is.null(over$fdr)) base$fdr <- as.numeric(over$fdr)
  if (!is.null(over$base_mean)) base$base_mean <- as.numeric(over$base_mean)
  if (!is.null(over$log2fc)) base$log2fc <- as.numeric(over$log2fc)
  if (!is.null(over$svalue)) base$svalue <- as.numeric(over$svalue)
  base
}

.upset_build_wide_matrix <- function(tasks, thr_overrides) {
  if (length(tasks) < 1L) {
    return(list(
      ok = FALSE,
      msg = "No DEG lists found for configured studies.",
      wide = NULL,
      set_cols = character(0),
      set_display = character(0),
      set_study = character(0),
      set_assay = character(0),
      set_sid = character(0),
      set_deg_file = character(0)
    ))
  }

  gene_sets <- vector("list", length(tasks))
  set_cols <- character(length(tasks))
  set_display <- character(length(tasks))
  set_study <- character(length(tasks))
  set_assay <- character(length(tasks))
  set_sid <- character(length(tasks))
  set_deg_file <- character(length(tasks))

  for (i in seq_along(tasks)) {
    tk <- tasks[[i]]
    sid <- tk$sid
    slbl <- tk$study_label
    entry <- tk$entry
    deg_rel <- entry$deg_file
    cmp_lbl <- ora_deg_entry_comparison_label(entry)
    scol <- sprintf("S%04d", i)
    set_cols[[i]] <- scol
    set_display[[i]] <- paste0(slbl, " :: ", cmp_lbl)
    set_study[[i]] <- slbl
    set_assay[[i]] <- .upset_study_assay_type(sid)
    set_sid[[i]] <- sid
    set_deg_file[[i]] <- if (is.null(deg_rel) || !nzchar(as.character(deg_rel))) "" else as.character(deg_rel)

    if (is.null(deg_rel) || !nzchar(as.character(deg_rel))) {
      gene_sets[[i]] <- character(0)
      next
    }

    thr <- deg_list_threshold_defaults(sid, deg_rel)
    key <- paste(sid, deg_rel, sep = "|")
    if (!is.null(thr_overrides[[key]])) {
      thr <- .upset_merge_thr(thr, thr_overrides[[key]])
    }

    loaded <- load_study_data(sid, deg_rel)
    res <- loaded$res
    mm <- loaded$mm
    if (is.null(res) || is.null(mm) || nrow(res) < 1L) {
      gene_sets[[i]] <- character(0)
      next
    }
    idx <- filter_heatmap_row_index(
      res, mm,
      thr$fdr, thr$base_mean, thr$log2fc, thr$svalue
    )
    if (is.null(idx) || length(idx) < 1L) {
      gene_sets[[i]] <- character(0)
    } else {
      sy <- toupper(trimws(as.character(res$symbol[idx])))
      gene_sets[[i]] <- unique(sy[!is.na(sy) & nzchar(sy)])
    }
  }

  all_genes <- unique(unlist(gene_sets, use.names = FALSE))
  all_genes <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
  if (length(all_genes) < 1L) {
    return(list(
      ok = FALSE,
      msg = "No genes pass thresholds for any DEG list (after overrides).",
      wide = NULL,
      set_cols = set_cols,
      set_display = set_display,
      set_study = set_study,
      set_assay = set_assay,
      set_sid = set_sid,
      set_deg_file = set_deg_file
    ))
  }

  wide <- data.frame(symbol = all_genes, stringsAsFactors = FALSE)
  for (i in seq_along(tasks)) {
    wide[[set_cols[[i]]]] <- wide$symbol %in% gene_sets[[i]]
  }

  list(
    ok = TRUE,
    msg = "",
    wide = wide,
    set_cols = set_cols,
    set_display = set_display,
    set_study = set_study,
    set_assay = set_assay,
    set_sid = set_sid,
    set_deg_file = set_deg_file
  )
}

#' Max intersection bar height (`comb_size`) among combination columns that include each set.
.upset_max_intersect_involvement <- function(m) {
  sn <- ComplexHeatmap::set_name(m)
  cn <- ComplexHeatmap::comb_name(m)
  sz <- ComplexHeatmap::comb_size(m)
  mx <- stats::setNames(rep(0, length(sn)), sn)
  for (k in seq_along(cn)) {
    pat <- as.character(cn[k])
    ch <- strsplit(pat, "", fixed = TRUE)[[1L]]
    if (length(ch) != length(sn)) next
    ck <- suppressWarnings(as.integer(sz[k]))
    if (!is.finite(ck)) next
    for (i in seq_along(sn)) {
      if (identical(ch[[i]], "1")) {
        si <- sn[[i]]
        mx[[si]] <- max(mx[[si]], ck)
      }
    }
  }
  mx
}

#' Named gene lists + combination matrix (max degree).
#' Uses `mode = "intersect"` in [ComplexHeatmap::make_comb_mat] so each bar is genes in the
#' selected sets with **1** (membership) regardless of other sets — e.g. all three pairwise
#' columns `110`, `101`, `011` exist for three lists. `mode = "distinct"` would use exclusive
#' partitions only, which drops empty slices and can omit pairwise patterns users expect.
#' @param row_sort `"study"` (study block order = config task order, then list label) or
#'   `"intersect_max"` (largest involved overlap bar, then label).
.upset_comb_from_pl <- function(pl, max_deg, row_sort = "study") {
  wide <- pl$wide
  set_cols <- pl$set_cols
  dsp <- make.unique(as.character(pl$set_display))
  lt <- stats::setNames(
    lapply(seq_along(set_cols), function(i) {
      unique(as.character(wide$symbol[wide[[set_cols[[i]]]] %in% TRUE]))
    }),
    dsp
  )
  lt <- lt[lengths(lt) > 0L]
  if (length(lt) < 1L) {
    return(list(ok = FALSE, msg = "All DEG lists are empty after filtering.", m = NULL))
  }
  m <- tryCatch(
    ComplexHeatmap::make_comb_mat(lt, mode = "intersect"),
    error = function(e) list(error = conditionMessage(e))
  )
  if (!inherits(m, "comb_mat")) {
    return(list(ok = FALSE, msg = paste0("make_comb_mat: ", m$error %||% "failed"), m = NULL))
  }
  mk <- ComplexHeatmap::comb_degree(m) <= as.integer(max_deg) & ComplexHeatmap::comb_degree(m) >= 2L
  m <- tryCatch(m[, mk], error = function(e) NULL)
  if (is.null(m)) {
    return(list(ok = FALSE, msg = "No intersections after applying max degree filter.", m = NULL))
  }
  sn <- ComplexHeatmap::set_name(m)
  assay_by <- stats::setNames(as.character(pl$set_assay), dsp)[sn]
  study_by <- stats::setNames(trimws(as.character(pl$set_study)), dsp)[sn]
  lab_by <- stats::setNames(trimws(as.character(pl$set_display)), dsp)[sn]
  assay_by[is.na(assay_by)] <- "Unknown"
  study_by[is.na(study_by)] <- "(unknown)"
  lab_by[is.na(lab_by)] <- sn
  # Study column order in sidebar/config (first appearance), not alphabetical on labels.
  study_levels <- unique(trimws(as.character(pl$set_study)))
  rs <- as.character(row_sort %||% "study")
  if (identical(rs, "set_size")) rs <- "intersect_max"
  if (!rs %in% c("study", "intersect_max")) rs <- "study"
  if (identical(rs, "study")) {
    study_rank <- match(study_by[sn], study_levels)
    study_rank[is.na(study_rank)] <- .Machine$integer.max
    ord_sets <- sn[order(study_rank, lab_by[sn])]
  } else {
    imx <- .upset_max_intersect_involvement(m)
    ord_sets <- sn[order(-unname(imx[sn]), lab_by[sn])]
  }
  # Integer row_order is robust if any set name encoding differs from rownames(m).
  ord_ix <- match(ord_sets, sn)
  if (length(ord_ix) != length(sn) || anyNA(ord_ix)) {
    ord_ix <- seq_along(sn)
  }
  list(
    ok = TRUE,
    msg = "",
    m = m,
    set_order = ord_ix,
    set_names_plot_order = ord_sets,
    row_labels = unname(lab_by[ord_sets]),
    assay_ord = unname(assay_by[ord_sets]),
    study_ord = unname(study_by[ord_sets]),
    list_names = dsp
  )
}

#' Genes in each ComplexUpset intersection (`mode = intersect` / inclusive_intersection).
.upset_cu_genes_by_intersection <- function(ud, wide_symbol) {
  ints <- as.character(ud$plot_intersections_subset)
  if (length(ints) < 1L) {
    return(stats::setNames(list(), character(0)))
  }
  ns_lab <- ud$non_sanitized_labels
  out <- stats::setNames(vector("list", length(ints)), ints)
  cn <- setdiff(colnames(wide_symbol), "symbol")
  for (int in ints) {
    toks <- strsplit(int, "-", fixed = TRUE)[[1L]]
    orig <- unname(ns_lab[toks])
    orig <- orig[!is.na(orig) & nzchar(orig)]
    if (length(orig) < 1L) {
      out[[int]] <- character(0)
      next
    }
    miss <- setdiff(orig, cn)
    if (length(miss) > 0L) {
      out[[int]] <- character(0)
      next
    }
    sub <- wide_symbol[, orig, drop = FALSE]
    ok <- apply(sub, 1L, function(r) all(as.logical(r)))
    sy <- as.character(wide_symbol$symbol[ok])
    sy <- unique(sy[!is.na(sy) & nzchar(sy)])
    out[[int]] <- sy
  }
  out
}

.upset_cu_tooltip_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

#' Plain-text tooltips for ggiraph: one line per DEG list in the intersection, then gene count.
#' Lines follow matrix row order (top → bottom) like the plot: use `set_order_y` from
#' `rev(levels(ud$matrix_frame$group))` (ggplot discrete y: first level at bottom).
#' Uses newlines; pair with `opts_tooltip(css = "white-space: pre-line; ...")`.
.upset_cu_tooltip_by_intersection <- function(ud, genes_by_int, set_order_y = NULL) {
  if (length(genes_by_int) < 1L) return(character(0))
  ns_lab <- ud$non_sanitized_labels
  if (is.null(set_order_y) || length(set_order_y) < 1L) {
    gf <- ud$matrix_frame$group
    set_order_y <- if (is.factor(gf)) {
      rev(levels(gf))
    } else {
      rev(unique(as.character(gf)))
    }
  }
  set_order_y <- as.character(set_order_y)
  # With encode_sets=TRUE, matrix rows use encoded ids ("1","2",…); map to display for ordering.
  set_order_disp <- unname(ns_lab[set_order_y])
  na_hit <- is.na(set_order_disp) | !nzchar(as.character(set_order_disp))
  set_order_disp[na_hit] <- set_order_y[na_hit]
  out <- character(length(genes_by_int))
  names(out) <- names(genes_by_int)
  for (int in names(genes_by_int)) {
    toks <- strsplit(int, "-", fixed = TRUE)[[1L]]
    labels <- unname(ns_lab[toks])
    labels <- labels[!is.na(labels) & nzchar(labels)]
    n <- length(genes_by_int[[int]])
    if (length(labels) < 1L) {
      out[[int]] <- paste0(
        .upset_cu_tooltip_escape(int),
        "\nGenes in intersection: ",
        n
      )
      next
    }
    labels_plot <- set_order_disp[set_order_disp %in% labels]
    miss <- labels[!(labels %in% labels_plot)]
    labels_ord <- c(labels_plot, miss)
    lab_lines <- paste0(.upset_cu_tooltip_escape(labels_ord), collapse = "\n")
    out[[int]] <- paste0(lab_lines, "\nGenes in intersection: ", n)
  }
  out
}

# Strip between intersection bars and matrix. In HeatmapAnnotation(top), the *last*
# item is closest to the heatmap body, so the pad must follow intersection_size.
.upset_top_annotation_gap <- function(m, spacer_mm = 1.5) {
  sor <- attr(m, "param")$set_on_rows
  if (!isTRUE(sor)) {
    return(ComplexHeatmap::upset_top_annotation(m))
  }
  ncomb <- length(ComplexHeatmap::comb_size(m))
  pad_vec <- rep(0, ncomb)
  ComplexHeatmap::HeatmapAnnotation(
    intersection_size = ComplexHeatmap::anno_barplot(
      ComplexHeatmap::comb_size(m),
      border = FALSE,
      gp = grid::gpar(fill = "black"),
      height = grid::unit(2.85, "cm"),
      which = "column",
      # Large `extend` adds most of the visible white band under the bars (y-scale padding).
      extend = 0.03,
      axis_param = list(side = "left")
    ),
    matrix_pad = ComplexHeatmap::anno_barplot(
      pad_vec,
      which = "column",
      border = FALSE,
      axis = FALSE,
      gp = grid::gpar(fill = "#F0F0F0", col = NA),
      height = grid::unit(spacer_mm, "mm"),
      bar_width = 1
    ),
    gap = grid::unit(0, "mm"),
    show_annotation_name = c(TRUE, FALSE),
    annotation_label = c("Intersection\nsize", "")
  )
}

.upset_multiplicity_df <- function(wide, set_cols) {
  if (is.null(wide) || length(set_cols) < 1L) {
    return(data.frame(degree = integer(0), n_genes = integer(0)))
  }
  m <- as.matrix(wide[, set_cols, drop = FALSE])
  d <- as.integer(rowSums(m))
  tb <- table(d)
  if (length(tb) < 1L) {
    return(data.frame(degree = integer(0), n_genes = integer(0)))
  }
  deg <- as.integer(names(tb))
  n <- as.integer(tb)
  keep <- n > 0L
  deg <- deg[keep]
  n <- n[keep]
  if (length(deg) < 1L) {
    return(data.frame(degree = integer(0), n_genes = integer(0)))
  }
  ord <- order(deg)
  data.frame(degree = deg[ord], n_genes = n[ord], stringsAsFactors = FALSE)
}

# Left of matrix: DEG list label + assay + study only (`show_row_names = FALSE`).
# Set-size bars are a separate [ComplexHeatmap::rowAnnotation] on the right (see below).
.upset_left_row_annotation <- function(m, assay_ord, study_ord, row_labels) {
  assays_u <- unique(as.character(assay_ord))
  studies_u <- unique(as.character(study_ord))
  pal_a <- .upset_assay_palette(assays_u)
  pal_s <- pheno_colors(studies_u)
  lab_w <- .upset_row_names_max_width(row_labels)
  ComplexHeatmap::rowAnnotation(
    Label = ComplexHeatmap::anno_text(
      as.character(row_labels),
      gp = grid::gpar(fontsize = 10),
      just = "left"
    ),
    Assay = ComplexHeatmap::anno_simple(assay_ord, col = pal_a),
    Study = ComplexHeatmap::anno_simple(study_ord, col = pal_s),
    annotation_name_side = "bottom",
    show_annotation_name = TRUE,
    gap = grid::unit(1.5, "mm"),
    annotation_width = grid::unit.c(
      lab_w,
      grid::unit(5, "mm"),
      grid::unit(5, "mm")
    )
  )
}

# Right of matrix: per-set sizes only (matches usual UpSet set-size strip placement).
.upset_right_set_size_annotation <- function(m) {
  ComplexHeatmap::rowAnnotation(
    set_size = ComplexHeatmap::anno_barplot(
      ComplexHeatmap::set_size(m),
      border = FALSE,
      gp = grid::gpar(fill = "black")
    ),
    annotation_name_side = "bottom",
    show_annotation_name = TRUE,
    annotation_width = grid::unit(2.2, "cm")
  )
}

# Horizontal legend drawn with the UpSet (assay + study).
.upset_horizontal_legends <- function(assay_ord, study_ord) {
  assays_u <- unique(as.character(assay_ord))
  studies_u <- unique(as.character(study_ord))
  assays_u <- assays_u[!is.na(assays_u) & nzchar(assays_u)]
  studies_u <- studies_u[!is.na(studies_u) & nzchar(studies_u)]
  pal_a <- .upset_assay_palette(assays_u)
  pal_s <- pheno_colors(studies_u)
  legs <- list()
  if (length(assays_u) >= 1L) {
    fill_a <- unname(pal_a[assays_u])
    fill_a[is.na(fill_a)] <- "#D9D9D9"
    legs[[length(legs) + 1L]] <- ComplexHeatmap::Legend(
      title = "Assay",
      at = assays_u,
      legend_gp = grid::gpar(fill = fill_a),
      nrow = 1,
      column_gap = grid::unit(4, "mm")
    )
  }
  if (length(studies_u) >= 1L) {
    fill_s <- unname(pal_s[studies_u])
    fill_s[is.na(fill_s)] <- "#BBBBBB"
    legs[[length(legs) + 1L]] <- ComplexHeatmap::Legend(
      title = "Study",
      at = studies_u,
      legend_gp = grid::gpar(fill = fill_s),
      nrow = 1,
      column_gap = grid::unit(10, "mm")
    )
  }
  if (length(legs) == 0L) return(NULL)
  if (length(legs) == 1L) return(legs[[1L]])
  lg <- legs[[1L]]
  for (k in seq_len(length(legs) - 1L) + 1L) {
    lg <- ComplexHeatmap::packLegend(
      lg,
      legs[[k]],
      direction = "horizontal",
      gap = grid::unit(14, "mm")
    )
  }
  lg
}

# Width for the Label `anno_text` column (Heatmap row names default to 6 cm and clip).
.upset_row_names_max_width <- function(row_labels) {
  rl <- as.character(row_labels %||% character(0))
  mx <- if (length(rl) > 0L) suppressWarnings(max(nchar(rl), na.rm = TRUE)) else 20L
  if (!is.finite(mx) || mx < 1L) mx <- 20L
  cm <- max(8, as.numeric(mx) * 0.12)
  cm <- min(cm, 32)
  grid::unit(cm, "cm")
}

# SVG width for ComplexUpset (patchwork): long y-labels + intersections + set-size strip.
# ggplot uses one discrete x scale (equal column width); total SVG width sums a sublinear
# per-intersection budget from column degree so high-degree columns get more room overall.
.upset_cu_intersection_degrees <- function(intersection_ids) {
  ids <- as.character(intersection_ids %||% character(0))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) < 1L) return(integer(0))
  vapply(ids, function(s) as.integer(length(strsplit(s, "-", fixed = TRUE)[[1L]])), integer(1L), USE.NAMES = FALSE)
}

.upset_cu_intersection_pretty_labels <- function(intersection_ids, ud) {
  ids <- as.character(intersection_ids %||% character(0))
  if (length(ids) < 1L) return(character(0))
  ns_lab <- ud$non_sanitized_labels
  stats::setNames(
    vapply(ids, function(id) {
      toks <- strsplit(id, "-", fixed = TRUE)[[1L]]
      labs <- unname(ns_lab[toks])
      labs <- labs[!is.na(labs) & nzchar(as.character(labs))]
      if (length(labs) < 1L) return(id)
      paste(as.character(labs), collapse = " & ")
    }, character(1L), USE.NAMES = FALSE),
    ids
  )
}

.upset_cu_plot_width_svg <- function(mat, n_intersections, intersection_degrees = NULL) {
  cn <- colnames(mat)
  lw <- if (length(cn) > 0L) suppressWarnings(max(nchar(as.character(cn)), na.rm = TRUE)) else 20L
  if (!is.finite(lw) || lw < 12L) lw <- 12L
  ni <- as.integer(n_intersections %||% max(length(cn), 1L))
  ni <- max(ni, 1L)

  degs <- intersection_degrees
  if (is.null(degs) || length(degs) < 1L) {
    degs <- rep(2L, ni)
  }
  degs <- as.integer(degs)
  if (length(degs) != ni) {
    if (length(degs) > ni) {
      degs <- degs[seq_len(ni)]
    } else if (length(degs) >= 1L) {
      pad <- rep(max(degs, na.rm = TRUE), ni - length(degs))
      pad[is.na(pad) | !is.finite(pad)] <- 2L
      degs <- c(degs, pad)
    } else {
      degs <- rep(2L, ni)
    }
  }
  degs[!is.finite(degs) | degs < 2L] <- 2L
  degs[degs > 48L] <- 48L

  # Per-column inch budget: sqrt(degree-1) grows sublinearly; cap + floor so bars do not
  # visually merge when the SVG is scaled down in the browser. Slightly generous vs an older
  # baseline so many columns stay readable at typical browser widths.
  d1 <- pmax(0, as.numeric(degs - 1L))
  per_col <- 0.108 + 0.09 * sqrt(d1)
  per_col <- pmin(per_col, 0.33)
  per_col <- pmax(per_col, 0.158)
  w_cols <- sum(per_col)

  w_in <- 1.02 + w_cols + length(cn) * 0.10 + 1.32
  w_in <- max(w_in, 10)
  min(w_in, 120)
}

# Max n*(n-1)/2 pairwise columns before upset_data(intersections='all') is skipped (observed only).
.upset_cu_max_pair_grid_columns <- 550L

.upset_cu_intersections_setting <- function(max_degree, n_sets) {
  md <- suppressWarnings(as.integer(max_degree))
  ns <- suppressWarnings(as.integer(n_sets))
  if (length(md) != 1L || is.na(md) || md < 2L) md <- 2L
  if (length(ns) != 1L || is.na(ns) || ns < 2L) return("observed")
  if (md != 2L) return("observed")
  n_pair <- (as.numeric(ns) * (as.numeric(ns) - 1)) / 2
  if (!is.finite(n_pair) || n_pair > as.numeric(.upset_cu_max_pair_grid_columns)) {
    return("observed")
  }
  "all"
}

# ComplexUpset::upset() ends with plot_layout(heights = c(rep(1, n_ann), height_ratio)): every
# annotation row (intersection-size bar first) has weight 1. To shorten only the bar row without
# stretching the dot matrix, use weight < 1 on that row and scale SVG height so the matrix row
# keeps the same share of total height as with all weights 1 (see upset() deparse: heights = ...).
.upset_cu_shrink_bar_layout <- function(
    p,
    height_ratio = 1,
    bar_row_weight = .upset_cu_intersection_bar_row_weight,
    base_height_svg = 7.6,
    n_annotation_rows = 1L
) {
  if (!inherits(p, "patchwork") || !requireNamespace("patchwork", quietly = TRUE)) {
    return(list(plot = p, height_svg = base_height_svg))
  }
  n_ann <- suppressWarnings(as.integer(n_annotation_rows[[1L]]))
  if (!is.finite(n_ann) || n_ann < 1L) n_ann <- 1L
  bar_w <- as.numeric(bar_row_weight)
  if (!is.finite(bar_w) || bar_w <= 0 || bar_w > 1) bar_w <- 1
  hr <- as.numeric(height_ratio)
  if (!is.finite(hr) || hr <= 0) hr <- 1
  t_old <- as.numeric(n_ann) + hr
  t_new <- bar_w + as.numeric(max(0L, n_ann - 1L)) + hr
  h_svg <- base_height_svg * t_new / t_old
  heights <- c(bar_w, rep(1, max(0L, n_ann - 1L)), hr)
  list(
    plot = p + patchwork::plot_layout(heights = heights),
    height_svg = h_svg
  )
}

#' Minimal theme tweaks only. Never set `axis.text.x` on `intersections_matrix`: ComplexUpset
#' keeps that axis blank; forcing labels there paints long intersection strings across the
#' matrix and destroys the patchwork layout. The intersection-size bar panel keeps discrete x
#' axis text blank; bar heights encode size (no count labels on bars).
.upset_cu_upset_themes <- function() {
  m_side <- ggplot2::margin(t = 1, r = 2.5, b = 2, l = 2.5, unit = "mm")
  ComplexUpset::upset_modify_themes(list(
    intersections_matrix = ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin = m_side
    ),
    overall_sizes = ggplot2::theme(plot.margin = m_side),
    `Intersection size` = ggplot2::theme(
      plot.margin = ggplot2::margin(t = 2, r = 2, b = 0, l = 2, unit = "mm"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  ))
}

#' Row stripes coloured by study display label; `label_colors` is named vector label -> hex.
.upset_cu_stripes <- function(mat, study_ord, label_colors) {
  if (is.null(mat) || ncol(mat) < 1L) return(ComplexUpset::upset_stripes())
  sets <- colnames(mat)
  so <- trimws(as.character(study_ord))
  if (length(so) != length(sets)) return(ComplexUpset::upset_stripes())
  sd <- data.frame(set = sets, study = so, stringsAsFactors = FALSE)
  studies_u <- unique(sd$study)
  pal_named <- stats::setNames(rep("#e8eaef", length(studies_u)), studies_u)
  lc <- label_colors
  if (is.null(lc)) lc <- character(0)
  for (st in studies_u) {
    hit <- lc[[st]]
    if (!is.null(hit) && !is.na(hit) && nzchar(as.character(hit))) {
      pal_named[[st]] <- as.character(hit)[[1L]]
    }
  }
  ComplexUpset::upset_stripes(
    mapping = ggplot2::aes(color = .data$study),
    data = sd,
    colors = pal_named,
    geom = ggplot2::geom_segment(linewidth = 5)
  )
}

#' Small HTML legend for assay + study colours (matches stripes / prior ComplexHeatmap legend).
.upset_cu_color_key_ui <- function(assay_ord, study_ord, study_color_map = NULL) {
  assays_u <- unique(trimws(as.character(assay_ord %||% character(0))))
  studies_u <- unique(trimws(as.character(study_ord %||% character(0))))
  assays_u <- assays_u[!is.na(assays_u) & nzchar(assays_u)]
  studies_u <- studies_u[!is.na(studies_u) & nzchar(studies_u)]
  pal_a <- .upset_assay_palette(assays_u)
  pal_s <- pheno_colors(studies_u)
  if (!is.null(study_color_map) && length(study_color_map) > 0L) {
    for (nm in studies_u) {
      hit <- study_color_map[[nm]]
      if (!is.null(hit) && !is.na(hit) && nzchar(as.character(hit))) {
        pal_s[[nm]] <- as.character(hit)[[1L]]
      }
    }
  }
  chip_row <- function(title, labels, pal) {
    chips <- lapply(labels, function(nm) {
      col <- pal[[as.character(nm)]] %||% "#cccccc"
      if (is.na(col) || !nzchar(col)) col <- "#cccccc"
      shiny::tags$span(
        style = "display:inline-flex;align-items:center;margin:2px 12px 2px 0;font-size:12px;",
        shiny::tags$span(
          style = sprintf(
            "display:inline-block;width:14px;height:14px;border-radius:2px;background:%s;margin-right:6px;border:1px solid #bbb;",
            col
          )
        ),
        shiny::tags$span(as.character(nm))
      )
    })
    shiny::tags$div(
      style = "margin-bottom:8px;",
      shiny::tags$strong(title),
      shiny::tags$div(style = "display:flex;flex-wrap:wrap;margin-top:4px;", chips)
    )
  }
  shiny::tagList(
    if (length(assays_u) > 0L) chip_row("Assay", assays_u, pal_a),
    if (length(studies_u) > 0L) chip_row("Study", studies_u, pal_s)
  )
}

upsetTabUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(shiny::HTML(
      ".upset-plot-wrap img { max-width: none !important; width: auto !important; height: auto !important; }
      .upset-mult-wrap img { max-width: 100%; width: auto; height: auto; }"
    )),
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::uiOutput(ns("status_msg")),
        shiny::uiOutput(ns("color_key")),
        shiny::tags$div(
          style = "max-width: 380px;",
          shiny::selectInput(
            ns("row_sort"),
            label = "Row (set) order:",
            choices = c(
              "By study (then list label)" = "study",
              "By largest intersection overlap (desc.)" = "intersect_max"
            ),
            selected = "study",
            width = "100%"
          )
        ),
        shiny::tags$div(
          class = "upset-plot-wrap",
          style = "overflow-x: auto; width: 100%; max-width: 100%; min-width: 0;",
          shiny::tags$div(
            style = "display: inline-block; vertical-align: top; max-width: none;",
            ggiraph::girafeOutput(ns("upset_plot"), height = paste0(.upset_cu_girafe_plot_height_px(), "px"))
          )
        ),
        shiny::tags$hr(),
        shiny::tags$div(
          class = "upset-mult-wrap",
          style = "display: inline-block; max-width: 520px;",
          shiny::plotOutput(ns("plot_multiplicity"), height = "220px")
        ),
        shiny::tags$hr(),
        shiny::fluidRow(
          shiny::column(
            width = 6,
            shiny::selectInput(ns("int_pick"), "Intersection (copy genes)", choices = NULL, width = "100%")
          ),
          shiny::column(
            width = 3,
            shiny::tags$div(
              style = "margin-top: 24px;",
              shiny::actionButton(ns("int_copy"), "Copy genes", class = "btn-primary")
            )
          )
        ),
      )
    )
  )
}

#' @param cfg Reactive: `list(max_degree = int, min_intersection = num, refresh = int, row_sort = chr)` invalidates rebuilds.
#' @param thr_overrides_parent `reactiveValues()` from parent; updated on Refresh before `cfg()$refresh` increments.
#' @param on_dot_click Optional `function(study_id, deg_file, genes_chr)` when user clicks
#'   an active dot in the intersection matrix (opens DEGs tab with that list + intersection genes).
#' @param deg_by_study reactive or static named list: study_id -> visible deg_file paths
upsetTabServer <- function(id, study_ids, study_labels, cfg, thr_overrides_parent, on_dot_click = NULL, deg_by_study = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    app_scope <- function() {
      list(
        study_ids = if (is.function(study_ids)) study_ids() else study_ids,
        study_labels = if (is.function(study_labels)) study_labels() else study_labels,
        deg_by_study = if (is.function(deg_by_study)) deg_by_study() else deg_by_study
      )
    }

    thr_overrides <- thr_overrides_parent
    last_payload <- shiny::reactiveVal(NULL)

    upset_payload <- shiny::reactive({
      sc <- app_scope()
      cfg()$max_degree
      cfg()$min_intersection
      cfg()$refresh
      cfg()$row_sort
      ov <- shiny::reactiveValuesToList(thr_overrides)
      tasks <- ora_build_deg_tasks(sc$study_ids, sc$study_labels, deg_filter = NULL, deg_by_study = sc$deg_by_study)
      pl <- .upset_build_wide_matrix(tasks, ov)
      if (!isTRUE(pl$ok)) {
        return(pl)
      }
      max_deg <- suppressWarnings(as.integer(cfg()$max_degree))
      if (length(max_deg) != 1L || is.na(max_deg) || max_deg < 2L) {
        max_deg <- length(pl$set_cols)
      }
      max_deg <- min(max_deg, length(pl$set_cols))
      max_deg <- max(2L, max_deg)

      min_sz <- suppressWarnings(as.numeric(cfg()$min_intersection))
      if (length(min_sz) != 1L || is.na(min_sz) || min_sz < 1) min_sz <- 1
      min_sz <- min(min_sz, 1e7)

      ch <- .upset_comb_from_pl(pl, max_deg, row_sort = cfg()$row_sort %||% "study")
      if (!isTRUE(ch$ok)) {
        pl$ok <- FALSE
        pl$msg <- ch$msg
        return(pl)
      }
      dsp <- make.unique(as.character(pl$set_display))
      mat0 <- pl$wide[, pl$set_cols, drop = FALSE]
      colnames(mat0) <- dsp
      sn_order <- ch$set_names_plot_order
      mat <- mat0[, sn_order, drop = FALSE]
      if (ncol(mat) < 2L) {
        pl$ok <- FALSE
        pl$msg <- "Need at least two non-empty DEG lists to draw an UpSet."
        return(pl)
      }
      if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
        pl$ok <- FALSE
        pl$msg <- "Install CRAN package ComplexUpset for the UpSet tab."
        return(pl)
      }
      int_sets <- .upset_cu_intersections_setting(max_deg, ncol(mat))
      ud <- tryCatch(
        ComplexUpset::upset_data(
          mat,
          intersect = colnames(mat),
          mode = "intersect",
          min_size = min_sz,
          min_degree = 2L,
          max_degree = max_deg,
          intersections = int_sets,
          sort_sets = FALSE,
          sort_intersections = "descending",
          sort_intersections_by = "cardinality",
          encode_sets = TRUE,
          group_by = "degree"
        ),
        error = function(e) list(error = conditionMessage(e))
      )
      if (!is.null(ud$error)) {
        pl$ok <- FALSE
        pl$msg <- paste0("ComplexUpset: ", ud$error)
        return(pl)
      }
      if (length(ud$plot_intersections_subset) < 1L) {
        pl$ok <- FALSE
        pl$msg <- "No intersections after applying max degree filter."
        return(pl)
      }
      wide_sym <- cbind(symbol = pl$wide$symbol, mat, stringsAsFactors = FALSE)
      genes_by_int <- .upset_cu_genes_by_intersection(ud, wide_sym)
      gf <- ud$matrix_frame$group
      plot_y_order <- if (is.factor(gf)) rev(levels(gf)) else rev(unique(as.character(gf)))
      pl$m_comb <- ch$m
      pl$row_labels <- ch$row_labels
      pl$assay_ord <- ch$assay_ord
      pl$study_ord <- ch$study_ord
      pl$set_order_ch <- ch$set_order
      pl$max_degree <- max_deg
      pl$min_intersection <- min_sz
      scm <- study_color_map_for_labels(sc$study_ids, sc$study_labels)
      int_ids <- as.character(ud$plot_intersections_subset %||% character(0))
      int_deg <- .upset_cu_intersection_degrees(int_ids)
      int_labs <- .upset_cu_intersection_pretty_labels(int_ids, ud)
      # CRITICAL: ComplexUpset(encode_sets=TRUE) does NOT assign encoded ids in colnames(mat)
      # order. The truth table is ud$non_sanitized_labels: names = encoded id strings, values =
      # display names (= keys of disp_to_sid / disp_to_deg).
      disp_to_sid_v <- stats::setNames(pl$set_sid, dsp)
      disp_to_deg_v <- stats::setNames(pl$set_deg_file, dsp)
      ns_lab <- ud$non_sanitized_labels
      enc_ids <- as.character(names(ns_lab) %||% character(0))
      sid_by_enc <- stats::setNames(
        vapply(enc_ids, function(enc) {
          cnm <- unname(ns_lab[[enc]])
          if (is.null(cnm) || is.na(cnm) || !nzchar(cnm)) return(NA_character_)
          v <- disp_to_sid_v[[cnm]]
          if (is.null(v)) NA_character_ else as.character(v)
        }, character(1L)),
        enc_ids
      )
      deg_by_enc <- stats::setNames(
        vapply(enc_ids, function(enc) {
          cnm <- unname(ns_lab[[enc]])
          if (is.null(cnm) || is.na(cnm) || !nzchar(cnm)) return(NA_character_)
          v <- disp_to_deg_v[[cnm]]
          if (is.null(v)) NA_character_ else as.character(v)
        }, character(1L)),
        enc_ids
      )
      pl$cu <- list(
        mat = mat,
        genes_by_int = genes_by_int,
        tooltip_by_int = .upset_cu_tooltip_by_intersection(ud, genes_by_int, plot_y_order),
        disp_to_sid = stats::setNames(pl$set_sid, dsp),
        disp_to_deg = stats::setNames(pl$set_deg_file, dsp),
        sid_by_enc = sid_by_enc,
        deg_by_enc = deg_by_enc,
        n_intersections = length(int_ids),
        intersection_degrees = int_deg,
        intersection_label = int_labs,
        intersections_setting = int_sets,
        study_color_map = scm
      )
      pl
    })

    output$status_msg <- shiny::renderUI({
      pl <- upset_payload()
      if (isTRUE(pl$ok)) {
        nset <- ncol(pl$cu$mat)
        nint <- pl$cu$n_intersections %||% 0L
        int_set <- pl$cu$intersections_setting %||% "observed"
        min_i <- pl$min_intersection %||% 1
        pair_note <- if (identical(int_set, "all")) {
          " Every pair of lists with ≥1 shared gene has a column (inclusive intersect). "
        } else if (identical(pl$max_degree, 2L)) {
          " Pair grid omitted (too many lists for full pair matrix); showing observed intersections only. "
        } else {
          " "
        }
        min_note <- if (is.numeric(min_i) && is.finite(min_i) && min_i > 1) {
          paste0(" Minimum intersection size ≥ ", as.integer(min_i), ". ")
        } else {
          " Empty intersections (0 genes) are never shown. "
        }
        shiny::tags$div(
          class = "text-muted",
          style = "font-size: 0.9rem;",
          paste0(
            "Lists: ", nset,
            ". Columns: ", nint,
            ". Max degree ", pl$max_degree,
            ". Mode: intersect (inclusive), not distinct.",
            min_note,
            pair_note,
            "Study row striping uses ColorBrewer Set3 in main-config study order (swatches below). ",
            "Click a black (active) dot on a row to open that DEG list on the DEGs tab with the intersection gene set. ",
            "Gene multiplicity and Copy genes use the same intersection keys as the plot."
          )
        )
      } else {
        shiny::tags$div(class = "text-warning", pl$msg %||% "No data.")
      }
    })

    output$upset_plot <- ggiraph::renderGirafe({
      pl <- upset_payload()
      mk_empty <- function(msg) {
        gp <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg) +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::xlim(0, 1) +
          ggplot2::ylim(0, 1) +
          ggplot2::theme_void()
        ggiraph::girafe(ggobj = gp, options = list(ggiraph::opts_sizing(rescale = FALSE)))
      }
      if (!requireNamespace("ggiraph", quietly = TRUE)) {
        return(mk_empty("Package ggiraph is required for the UpSet plot."))
      }
      if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
        return(mk_empty("Package ComplexUpset is required for the UpSet plot."))
      }
      if (!isTRUE(pl$ok)) {
        return(mk_empty(pl$msg %||% "No data."))
      }
      mat <- pl$cu$mat
      if (is.null(mat) || ncol(mat) < 2L) {
        return(mk_empty("Need at least two non-empty DEG lists to draw an UpSet."))
      }
      tips <- pl$cu$tooltip_by_int
      if (is.null(tips) || length(tips) < 1L) {
        tips <- character(0)
      }
      # Bake study_id + deg_file into data_id keyed by ComplexUpset's ENCODED id (not by
      # colnames(mat) position — encode_sets=TRUE does NOT preserve column order). Use the
      # sid_by_enc / deg_by_enc maps built in the payload from ud$non_sanitized_labels.
      siddeg_sep <- .upset_siddeg_data_id_sep
      sid_map <- pl$cu$sid_by_enc %||% character(0)
      deg_map <- pl$cu$deg_by_enc %||% character(0)
      enc_keys <- as.character(names(sid_map))
      siddeg_by_enc <- stats::setNames(
        vapply(enc_keys, function(enc) {
          s <- as.character(sid_map[[enc]])
          d <- as.character(deg_map[[enc]])
          if (is.na(s)) s <- ""
          if (is.na(d)) d <- ""
          paste(s, d, sep = siddeg_sep)
        }, character(1L)),
        enc_keys
      )
      mx <- ComplexUpset::intersection_matrix(
        geom = ggiraph::geom_point_interactive(
          ggplot2::aes(
            data_id = paste0(
              as.character(.data$intersection),
              "\x01",
              unname(siddeg_by_enc[as.character(.data$group)]),
              "\x01",
              ifelse(as.logical(.data$value), "1", "0")
            ),
            tooltip = ifelse(
              as.logical(.data$value),
              ifelse(
                is.na(tips[as.character(.data$intersection)]),
                paste0(
                  .upset_cu_tooltip_escape(as.character(.data$intersection)),
                  "\nGenes in intersection: ?"
                ),
                tips[as.character(.data$intersection)]
              ),
              NA_character_
            )
          ),
          size = 3
        )
      )
      # Tighter intersection columns only; do not add coord_fixed or y-scale expand here —
      # those apply only to the matrix and break vertical alignment with the set-size panel.
      mx <- mx + ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.012, add = 0.02))
      int_sets <- pl$cu$intersections_setting %||% "observed"
      # Custom bar annotation: must set mode to inclusive_intersection (same as upset mode "intersect").
      # counts = FALSE: no numeric labels on bars (sizes still in tooltips / status / picker).
      base_ann <- list(
        `Intersection size` = ComplexUpset::intersection_size(
          counts = FALSE,
          mode = "inclusive_intersection"
        )
      )
      p <- ComplexUpset::upset(
        mat,
        intersect = colnames(mat),
        name = "Intersection",
        mode = "intersect",
        min_size = pl$min_intersection %||% 1,
        min_degree = 2,
        max_degree = pl$max_degree,
        intersections = int_sets,
        sort_sets = FALSE,
        sort_intersections = "descending",
        sort_intersections_by = "cardinality",
        encode_sets = TRUE,
        base_annotations = base_ann,
        matrix = mx,
        stripes = .upset_cu_stripes(mat, pl$study_ord, pl$cu$study_color_map),
        height_ratio = 1,
        width_ratio = 0.18,
        themes = .upset_cu_upset_themes()
      )
      # One base annotation ("Intersection size"); shrink that row only + scale SVG height so
      # the matrix row matches default row density (ComplexUpset uses rep(1, n_ann) for bars).
      lay <- .upset_cu_shrink_bar_layout(
        p,
        height_ratio = 1,
        bar_row_weight = .upset_cu_intersection_bar_row_weight,
        base_height_svg = 7.6,
        n_annotation_rows = 1L
      )
      wsvg <- .upset_cu_plot_width_svg(mat, pl$cu$n_intersections, pl$cu$intersection_degrees)
      ggiraph::girafe(
        ggobj = lay$plot,
        width_svg = wsvg,
        height_svg = lay$height_svg,
        options = list(
          ggiraph::opts_selection(
            type = "single",
            only_shiny = TRUE,
            css = "stroke:#39ff14;stroke-width:2.2px;"
          ),
          ggiraph::opts_hover(css = "stroke:#000;stroke-width:1.2px;"),
          ggiraph::opts_tooltip(
            opacity = 0.9,
            css = paste0(
              "padding:6px 10px;background:black;color:white;border-radius:2px;",
              "max-width:380px;white-space:pre-line;word-wrap:break-word;",
              "text-align:left;line-height:1.45;font-size:15px;"
            )
          ),
          ggiraph::opts_sizing(rescale = FALSE)
        )
      )
    })

    shiny::observe({
      pl <- upset_payload()
      if (!isTRUE(pl$ok) || is.null(pl$cu) || is.null(pl$cu$genes_by_int)) {
        last_payload(NULL)
        shiny::updateSelectInput(session, "int_pick", choices = NULL, selected = character(0))
        return()
      }
      genes_by_int <- pl$cu$genes_by_int
      ints <- names(genes_by_int)
      last_payload(list(genes_by_int = genes_by_int))
      ilab <- pl$cu$intersection_label
      pref <- as.character(ilab[ints])
      na_pref <- is.na(pref) | !nzchar(pref)
      pref[na_pref] <- ints[na_pref]
      labs <- paste0(pref, " (n=", lengths(genes_by_int), ")")
      shiny::updateSelectInput(
        session,
        "int_pick",
        choices = stats::setNames(ints, labs),
        selected = if (length(ints) > 0L) ints[[1L]] else NULL
      )
    })

    shiny::observeEvent(input$int_copy, {
      pl <- last_payload()
      if (is.null(pl)) return()
      int <- input$int_pick
      if (is.null(int) || !nzchar(int)) return()
      genes <- pl$genes_by_int[[int]]
      if (length(genes) < 1L) {
        shiny::showNotification("No genes in this intersection.", type = "warning")
        return()
      }
      session$sendCustomMessage("exprs_copy_text", paste(genes, collapse = "\n"))
      shiny::showNotification(paste(length(genes), "gene symbols copied."), type = "message")
    }, ignoreInit = TRUE)

    output$color_key <- shiny::renderUI({
      pl <- upset_payload()
      if (!isTRUE(pl$ok)) return(NULL)
      inner <- .upset_cu_color_key_ui(pl$assay_ord, pl$study_ord, pl$cu$study_color_map)
      shiny::tags$div(style = "margin: 0 0 12px 0;", inner)
    })

    shiny::observeEvent(input$upset_plot_selected, {
      if (!is.function(on_dot_click)) return()
      sel <- input$upset_plot_selected
      if (is.null(sel) || (is.character(sel) && !nzchar(sel))) return()
      if (is.character(sel) && length(sel) > 1L) sel <- sel[[1L]]
      parts <- strsplit(as.character(sel), "\x01", fixed = TRUE)[[1L]]
      if (length(parts) != 3L || !identical(parts[[3L]], "1")) return()
      int <- parts[[1L]]
      mid <- parts[[2L]]
      pl <- shiny::isolate(upset_payload())
      if (!isTRUE(pl$ok) || is.null(pl$cu)) return()
      genes <- pl$cu$genes_by_int[[int]]
      if (is.null(genes)) genes <- character(0)
      sid <- ""
      deg <- ""
      siddeg_sep <- .upset_siddeg_data_id_sep
      if (nzchar(mid) && grepl(siddeg_sep, mid, fixed = TRUE)) {
        bits <- strsplit(mid, siddeg_sep, fixed = TRUE)[[1L]]
        if (length(bits) >= 2L) {
          sid <- trimws(as.character(bits[[1L]]))
          deg <- trimws(as.character(bits[[2L]]))
        }
      } else {
        # Legacy data_id with only the encoded group id; resolve via sid_by_enc.
        grp <- mid
        sid <- as.character(pl$cu$sid_by_enc[[as.character(grp)]]) %||% ""
        deg <- as.character(pl$cu$deg_by_enc[[as.character(grp)]]) %||% ""
      }
      if (!nzchar(sid) || !nzchar(deg)) {
        shiny::showNotification("Could not resolve study / DEG list for that dot.", type = "warning")
        return()
      }
      on_dot_click(sid, deg, genes)
    }, ignoreInit = TRUE, ignoreNULL = TRUE)

    output$plot_multiplicity <- shiny::renderPlot(
      {
        pl <- upset_payload()
        shiny::req(isTRUE(pl$ok))
        df <- .upset_multiplicity_df(pl$wide, pl$set_cols)
        if (nrow(df) < 1L || !requireNamespace("ggplot2", quietly = TRUE)) {
          return(invisible(NULL))
        }
        lvl <- sort(unique(df$degree))
        df$x <- factor(df$degree, levels = lvl)
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$n_genes)) +
          ggplot2::geom_col(fill = "#5c7cfa", width = 0.78) +
          ggplot2::scale_x_discrete(
            name = "Number of DEG lists",
            expand = ggplot2::expansion(add = 0.55)
          ) +
          ggplot2::scale_y_continuous("Genes", expand = ggplot2::expansion(mult = c(0, 0.06))) +
          ggplot2::labs(title = "Gene multiplicity") +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::theme(
            aspect.ratio = 0.58,
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank()
          )
      },
      height = 220,
      width = 480,
      res = 100
    )

  })
}
