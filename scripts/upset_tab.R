# UpSet tab: DEG list intersections (ComplexHeatmap::UpSet + row annotations).
# Depends: ComplexHeatmap, grid, circlize, yaml, ggplot2; study_data, heatmap_utils, ora_cache.
# Sidebar inputs live in app.R (`upset_*`).

`%||%` <- function(x, y) if (is.null(x)) y else x

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
      set_sid = character(0)
    ))
  }

  gene_sets <- vector("list", length(tasks))
  set_cols <- character(length(tasks))
  set_display <- character(length(tasks))
  set_study <- character(length(tasks))
  set_assay <- character(length(tasks))
  set_sid <- character(length(tasks))

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
      set_sid = set_sid
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
    set_sid = set_sid
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
    row_labels = unname(lab_by[ord_sets]),
    assay_ord = unname(assay_by[ord_sets]),
    study_ord = unname(study_by[ord_sets]),
    list_names = dsp
  )
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

# Device width (px) for renderPlot: wide enough for long row labels + many intersections +
# left annotations; outer div scrolls horizontally (img must not shrink — see UI CSS).
.upset_plot_width_px <- function(m, row_labels) {
  ni <- length(ComplexHeatmap::comb_name(m))
  rl <- as.character(row_labels %||% character(0))
  lw <- if (length(rl) > 0L) suppressWarnings(max(nchar(rl), na.rm = TRUE)) else 20L
  if (!is.finite(lw) || lw < 12L) lw <- 12L
  lbl_px <- min(ceiling(as.numeric(lw) * 7.2), 980L)
  int_px <- max(ni, 1L) * 42L
  left_ann_px <- 200L
  right_sz_px <- 220L
  margins_px <- 140L
  w <- lbl_px + int_px + left_ann_px + right_sz_px + margins_px
  w <- max(w, 960L)
  w <- min(w, 9000L)
  as.integer(w)
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
            shiny::plotOutput(ns("upset_plot"), height = "600px")
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

#' @param cfg Reactive: `list(max_degree = int, refresh = int, row_sort = chr)` invalidates rebuilds.
#' @param thr_overrides_parent `reactiveValues()` from parent; updated on Refresh before `cfg()$refresh` increments.
upsetTabServer <- function(id, study_ids, study_labels, cfg, thr_overrides_parent) {
  shiny::moduleServer(id, function(input, output, session) {
    thr_overrides <- thr_overrides_parent
    last_payload <- shiny::reactiveVal(NULL)

    upset_payload <- shiny::reactive({
      cfg()$max_degree
      cfg()$refresh
      cfg()$row_sort
      ov <- shiny::reactiveValuesToList(thr_overrides)
      tasks <- ora_build_deg_tasks(study_ids, study_labels, deg_filter = NULL)
      pl <- .upset_build_wide_matrix(tasks, ov)
      if (!isTRUE(pl$ok)) {
        return(pl)
      }
      max_deg <- suppressWarnings(as.integer(cfg()$max_degree))
      if (length(max_deg) != 1L || is.na(max_deg) || max_deg < 2L) {
        max_deg <- min(2L, length(pl$set_cols))
      }
      max_deg <- min(max_deg, length(pl$set_cols))

      ch <- .upset_comb_from_pl(pl, max_deg, row_sort = cfg()$row_sort %||% "study")
      if (!isTRUE(ch$ok)) {
        pl$ok <- FALSE
        pl$msg <- ch$msg
        return(pl)
      }
      pl$m_comb <- ch$m
      pl$row_labels <- ch$row_labels
      pl$assay_ord <- ch$assay_ord
      pl$study_ord <- ch$study_ord
      pl$set_order_ch <- ch$set_order
      pl$max_degree <- max_deg
      pl
    })

    output$status_msg <- shiny::renderUI({
      pl <- upset_payload()
      if (isTRUE(pl$ok)) {
        shiny::tags$div(
          class = "text-muted",
          style = "font-size: 0.9rem;",
          paste0(
            "Lists: ", length(ComplexHeatmap::set_name(pl$m_comb)),
            ". Max degree: ", pl$max_degree,
            ". Intersect mode; overlaps only (degree ≥ 2). Label + assay + study on the left, set-size on the right; assay/study legend is drawn above the matrix in the figure. See ComplexHeatmap UpSet book."
          )
        )
      } else {
        shiny::tags$div(class = "text-warning", pl$msg %||% "No data.")
      }
    })

    output$upset_plot <- shiny::renderPlot(
      {
        pl <- upset_payload()
        if (!isTRUE(pl$ok)) {
          grid::grid.newpage()
          grid::grid.text(pl$msg %||% "No data.")
          return(invisible(NULL))
        }
        m <- pl$m_comb
        if (length(ComplexHeatmap::set_name(m)) < 2L) {
          grid::grid.newpage()
          grid::grid.text("Need at least two non-empty DEG lists to draw an UpSet.")
          return(invisible(NULL))
        }
        la <- .upset_left_row_annotation(m, pl$assay_ord, pl$study_ord, pl$row_labels)
        ra <- .upset_right_set_size_annotation(m)
        rs <- cfg()$row_sort %||% "study"
        if (identical(rs, "set_size")) rs <- "intersect_max"
        if (!rs %in% c("study", "intersect_max")) rs <- "study"
        so <- pl$set_order_ch
        sor <- isTRUE(attr(m, "param")$set_on_rows)
        if (identical(rs, "intersect_max")) {
          csz <- ComplexHeatmap::comb_size(m)
          cdeg <- ComplexHeatmap::comb_degree(m)
          cnm <- ComplexHeatmap::comb_name(m)
          comb_order <- order(-csz, -cdeg, cnm)
        } else {
          # Match ComplexHeatmap::UpSet default: columns follow degree/pattern, not bar height.
          m_ord <- if (sor) m[so, , drop = FALSE] else m[, so, drop = FALSE]
          comb_order <- ComplexHeatmap::order.comb_mat(m_ord, decreasing = TRUE)
        }
        grid::grid.newpage()
        ht <- ComplexHeatmap::UpSet(
          m,
          set_order = so,
          comb_order = comb_order,
          top_annotation = .upset_top_annotation_gap(m, spacer_mm = 1.5),
          show_row_names = FALSE,
          left_annotation = la,
          right_annotation = ra,
          gap = grid::unit(1, "mm")
        )
        leg <- .upset_horizontal_legends(pl$assay_ord, pl$study_ord)
        if (!is.null(leg)) {
          # Top legend: ComplexHeatmap centers `heatmap_legend_list`; draw it top-left in post_fun.
          pad <- grid::unit(c(6, 14, 18, 10), "mm")
          ComplexHeatmap::draw(
            ht,
            newpage = FALSE,
            show_heatmap_legend = FALSE,
            padding = pad,
            post_fun = function(obj) {
              grid::pushViewport(grid::viewport(
                x = grid::unit(2, "mm"),
                y = grid::unit(1, "npc") - grid::unit(2, "mm"),
                width = grid::unit(0.65, "npc"),
                height = grid::unit(12, "mm"),
                just = c("left", "top")
              ))
              ComplexHeatmap::draw(
                leg,
                x = grid::unit(0, "npc"),
                y = grid::unit(1, "npc"),
                just = c("left", "top")
              )
              grid::upViewport()
            }
          )
        } else {
          ComplexHeatmap::draw(ht, newpage = FALSE, padding = grid::unit(c(6, 14, 6, 10), "mm"))
        }
      },
      height = 600,
      width = function() {
        pl <- upset_payload()
        if (!isTRUE(pl$ok) || is.null(pl$m_comb)) return(960L)
        if (length(ComplexHeatmap::set_name(pl$m_comb)) < 2L) return(960L)
        .upset_plot_width_px(pl$m_comb, pl$row_labels)
      },
      res = 100
    )

    shiny::observe({
      pl <- upset_payload()
      if (!isTRUE(pl$ok)) {
        last_payload(NULL)
        shiny::updateSelectInput(session, "int_pick", choices = NULL, selected = character(0))
        return()
      }
      m <- pl$m_comb
      ints <- ComplexHeatmap::comb_name(m)
      genes_by_int <- stats::setNames(
        lapply(ints, function(nm) ComplexHeatmap::extract_comb(m, nm)),
        ints
      )
      last_payload(list(genes_by_int = genes_by_int))
      labs <- paste0(ints, " (n=", lengths(genes_by_int), ")")
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
