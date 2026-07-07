# Panels tab: category-annotated per-study expression heatmaps.
# Depends: ComplexHeatmap, study_data, panel_heatmap_utils, pathway_signatures.

`%||%` <- function(x, y) if (is.null(x)) y else x

# Plot width (px): from built heatmap object when available.
.panel_heatmap_width_px <- function(ht, n_categories = 1L, show_row_names = TRUE) {
  panel_heatmap_plot_width_px(ht, show_row_names, n_categories)
}

.panel_plot_height_px <- function(panels) {
  if (is.null(panels) || length(panels) < 1L) return(420L)
  hs <- vapply(
    panels,
    function(p) {
      if (!is.null(p$error) || is.null(p$ht)) return(0L)
      panel_heatmap_plot_height_px(p$ht)
    },
    integer(1L)
  )
  hs <- hs[hs > 0L]
  if (length(hs) < 1L) return(420L)
  as.integer(min(6800L, max(280L, sum(hs))))
}

.panel_plot_width_px <- function(panels) {
  if (is.null(panels) || length(panels) < 1L) return(620L)
  ws <- vapply(
    panels,
    function(p) {
      if (!is.null(p$error) || is.null(p$ht)) return(0L)
      .panel_heatmap_width_px(
        p$ht,
        p$n_categories %||% 1L,
        isTRUE(p$show_row_names)
      )
    },
    integer(1L)
  )
  ws <- ws[ws > 0L]
  if (length(ws) < 1L) return(620L)
  max(ws)
}

.panel_default_thr <- function() {
  if (exists("default_heatmap_thresholds", mode = "function")) {
    default_heatmap_thresholds()
  } else {
    list(fdr = 0.05, log2fc = 1)
  }
}

.panel_display_options <- function(input) {
  thr <- .panel_default_thr()
  fdr <- suppressWarnings(as.numeric(input$fdr))
  if (is.na(fdr)) fdr <- thr$fdr
  log2fc <- suppressWarnings(as.numeric(input$log2fc))
  if (is.na(log2fc)) log2fc <- thr$log2fc
  list(
    deduplicate_genes = isTRUE(input$deduplicate_genes),
    cluster_rows = isTRUE(input$cluster_rows),
    cluster_columns = isTRUE(input$cluster_columns),
    show_row_names = isTRUE(input$show_rownames),
    significant_only = isTRUE(input$significant_only),
    fdr = fdr,
    log2fc = log2fc
  )
}

.panel_heatmap_ncol <- function(ht) {
  if (is.null(ht)) return(NA_integer_)
  suppressWarnings(as.integer(ncol(ht@matrix)))
}

.panel_columns_alignable <- function(ht_list) {
  ht_list <- ht_list[!vapply(ht_list, is.null, logical(1))]
  if (length(ht_list) < 2L) return(TRUE)
  ns <- vapply(ht_list, .panel_heatmap_ncol, integer(1L))
  ns <- ns[!is.na(ns)]
  length(ns) < 2L || length(unique(ns)) == 1L
}

.panel_draw_stacked_viewports <- function(ht_list, pad) {
  heights_mm <- vapply(
    ht_list,
    function(ht) {
      as.numeric(grid::convertHeight(
        ComplexHeatmap::ht_size(ht)$height,
        "mm",
        valueOnly = TRUE
      ))
    },
    numeric(1L)
  )
  total_mm <- sum(heights_mm)
  y_top <- 1
  n <- length(ht_list)
  grid::grid.newpage()
  for (i in seq_len(n)) {
    h_frac <- heights_mm[i] / total_mm
    grid::pushViewport(grid::viewport(y = y_top, height = h_frac, just = "top"))
    ComplexHeatmap::draw(
      ht_list[[i]],
      newpage = FALSE,
      merge_legend = TRUE,
      show_heatmap_legend = (i == n),
      padding = pad
    )
    grid::popViewport()
    y_top <- y_top - h_frac
  }
  invisible(ht_list)
}

.panel_draw_combined <- function(ht_list) {
  ht_list <- ht_list[!vapply(ht_list, is.null, logical(1))]
  if (length(ht_list) < 1L) return(invisible(NULL))
  pad <- grid::unit(c(2, 6, 2, 2), "mm")
  if (length(ht_list) == 1L) {
    ComplexHeatmap::draw(ht_list[[1L]], merge_legend = TRUE, padding = pad)
    return(invisible(ht_list[[1L]]))
  }
  if (.panel_columns_alignable(ht_list)) {
    combined <- tryCatch({
      cmb <- ht_list[[1L]]
      for (i in seq_along(ht_list)[-1L]) {
        cmb <- cmb %v% ht_list[[i]]
      }
      cmb
    }, error = function(e) NULL)
    if (!is.null(combined)) {
      ComplexHeatmap::draw(combined, merge_legend = TRUE, padding = pad)
      return(invisible(combined))
    }
  }
  .panel_draw_stacked_viewports(ht_list, pad)
}

panelTabUI <- function(id) {
  ns <- shiny::NS(id)
  thr <- .panel_default_thr()
  shiny::tagList(
    shiny::tags$style(shiny::HTML(sprintf(
      "#%s .panel-toolbar .form-group { margin-bottom: 0; }
       #%s .panel-toolbar .checkbox { margin-top: 0; margin-bottom: 0; min-height: 0; }
       #%s .panel-toolbar .control-label { margin-bottom: 2px; font-size: 0.82rem; font-weight: 600; }
       #%s .panel-toolbar .shiny-input-container { width: auto; max-width: 100%%; }",
      ns("toolbar"), ns("toolbar"), ns("toolbar"), ns("toolbar")
    ))),
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::uiOutput(ns("status_ui")),
        shiny::tags$div(
          id = ns("toolbar"),
          class = "panel-toolbar",
          style = paste(
            "margin: 8px 0 12px 0; padding: 8px 12px;",
            "background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 4px;"
          ),
          shiny::tags$div(
            style = "display: flex; flex-wrap: wrap; align-items: center; gap: 8px 16px;",
            shiny::tags$div(
              style = "display: flex; flex-wrap: wrap; align-items: center; gap: 4px 14px;",
              shiny::checkboxInput(ns("deduplicate_genes"), "Deduplicate", value = FALSE),
              shiny::checkboxInput(ns("cluster_rows"), "Cluster rows", value = TRUE),
              shiny::checkboxInput(ns("cluster_columns"), "Cluster columns", value = TRUE),
              shiny::checkboxInput(ns("show_rownames"), "Gene names", value = TRUE),
              shiny::checkboxInput(ns("significant_only"), "Significant only", value = FALSE)
            ),
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("significant_only")),
              shiny::tags$div(
                style = "display: flex; flex-wrap: wrap; align-items: flex-end; gap: 8px 12px;",
                shiny::selectInput(
                  ns("fdr"),
                  label = "FDR",
                  choices = c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1, "0.5" = 0.5),
                  selected = thr$fdr,
                  width = "96px"
                ),
                shiny::numericInput(
                  ns("log2fc"),
                  label = "|log2FC|",
                  value = thr$log2fc,
                  min = 0,
                  step = 0.1,
                  width = "88px"
                )
              )
            ),
            shiny::actionButton(
              ns("update"),
              "Update",
              class = "btn-sm btn-primary",
              style = "margin: 0;"
            )
          )
        ),
        shiny::uiOutput(ns("heatmap_plot_ui"))
      )
    )
  )
}

panelTabServer <- function(
    id,
    study_id,
    custom_ontology,
    panels_input,
    generate_trigger = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    panel_results <- shiny::reactiveVal(NULL)

    build_one_panel <- function(panel, t2g, pi, disp, study) {
      matrix_rn <- rownames(panel$mm)
      resolved <- resolve_symbols_in_study_matrix(
        study, unique(t2g$gene), matrix_ids = matrix_rn
      )
      row_df <- prepare_panel_heatmap_rows(t2g, resolved, deduplicate_genes = isTRUE(disp$deduplicate_genes))
      missing_syms <- setdiff(unique(toupper(t2g$gene)), resolved$symbol[resolved$found])
      matrix_title <- if (!is.null(panel$title) && nzchar(panel$title)) {
        as.character(panel$title)
      } else {
        as.character(panel$counts_file)
      }
      if (is.null(row_df) || nrow(row_df) < 1L) {
        return(list(
          error = "No ontology genes matched this expression matrix.",
          title = matrix_title,
          counts_file = panel$counts_file,
          missing = missing_syms,
          found = 0L,
          requested = length(unique(t2g$gene))
        ))
      }
      sample_filter <- pi$sample_filter
      filter_col <- NULL
      filter_val <- NULL
      if (is.list(sample_filter) && !is.null(sample_filter$col) && !is.null(sample_filter$value)) {
        filter_col <- as.character(sample_filter$col)
        filter_val <- as.character(sample_filter$value)
      }
      sig_n_before <- nrow(row_df)
      sig_n_after <- sig_n_before
      sig_msg <- NULL
      if (isTRUE(disp$significant_only)) {
        sig_res <- filter_panel_row_df_significant(
          row_df,
          study_id = study,
          panel = panel,
          sample_filter_col = filter_col,
          sample_filter_value = filter_val,
          fdr = disp$fdr,
          log2fc = disp$log2fc
        )
        sig_msg <- sig_res$sig_msg
        sig_n_before <- sig_res$n_before
        sig_n_after <- sig_res$n_after
        row_df <- sig_res$rows
        if (!is.null(sig_msg) && nzchar(sig_msg)) {
          return(list(
            error = sig_msg,
            title = matrix_title,
            counts_file = panel$counts_file,
            missing = missing_syms,
            found = 0L,
            requested = length(unique(t2g$gene))
          ))
        }
        if (is.null(row_df) || nrow(row_df) < 1L) {
          return(list(
            error = paste0(
              "No significant ontology genes at FDR \u2264 ", disp$fdr,
              ", |log2FC| \u2265 ", disp$log2fc, "."
            ),
            title = matrix_title,
            counts_file = panel$counts_file,
            missing = missing_syms,
            found = 0L,
            requested = length(unique(t2g$gene)),
            sig_n_before = sig_n_before,
            sig_n_after = 0L
          ))
        }
      }
      pheno_col <- "PhenoNames"
      pheno_order <- pi$pheno_order
      if (nzchar(filter_val %||% "")) {
        matrix_title <- paste0(matrix_title, " (", filter_val, ")")
      }
      pd <- prepare_panel_heatmap_data(
        mm = panel$mm,
        metadata = panel$metadata,
        row_df = row_df,
        pheno_col = pheno_col,
        pheno_order = pheno_order,
        sample_filter_col = filter_col,
        sample_filter_value = filter_val,
        is_count_like = isTRUE(panel$is_count_like),
        cluster_rows = isTRUE(disp$cluster_rows),
        cluster_columns = isTRUE(disp$cluster_columns),
        show_row_names = isTRUE(disp$show_row_names),
        matrix_title = matrix_title
      )
      if (!is.null(pd$error)) {
        return(list(
          error = pd$error,
          title = matrix_title,
          counts_file = panel$counts_file,
          missing = missing_syms,
          found = nrow(row_df),
          requested = length(unique(t2g$gene))
        ))
      }
      ht <- panel_heatmap_build(pd)
      if (is.null(ht)) {
        return(list(
          error = "Could not build panel heatmap.",
          title = matrix_title,
          counts_file = panel$counts_file,
          missing = missing_syms,
          found = nrow(row_df),
          requested = length(unique(t2g$gene))
        ))
      }
      list(
        error = NULL,
        ht = ht,
        title = matrix_title,
        counts_file = panel$counts_file,
        missing = missing_syms,
        found = pd$n_genes,
        requested = length(unique(t2g$gene)),
        n_categories = length(pd$categories),
        n_samples = pd$n_samples,
        n_pheno_groups = pd$n_pheno_groups,
        show_row_names = isTRUE(disp$show_row_names),
        plot_width_px = panel_heatmap_plot_width_px(
          ht,
          isTRUE(disp$show_row_names),
          length(pd$categories)
        ),
        sig_n_before = if (isTRUE(disp$significant_only)) sig_n_before else NA_integer_,
        sig_n_after = if (isTRUE(disp$significant_only)) sig_n_after else NA_integer_,
        key = panel$key
      )
    }

    run_generate <- function() {
      panel_results(list(loading = TRUE))
      shiny::withProgress(
        message = "Building panel heatmap…",
        value = 0,
        min = 0,
        max = 1,
        {
          study <- study_id()
          if (is.null(study) || !nzchar(study)) {
            panel_results(list(error = "Select a study in the sidebar."))
            return(invisible(NULL))
          }
          co <- if (is.function(custom_ontology)) custom_ontology() else NULL
          if (is.null(co) || !is.list(co) || !isTRUE(co$active) || is.null(co$t2g) || nrow(co$t2g) < 1L) {
            panel_results(list(error = "Load a gene ontology (sidebar) before generating the panel heatmap."))
            return(invisible(NULL))
          }
          t2g <- co$t2g
          pi <- if (is.function(panels_input)) panels_input() else list()
          disp <- .panel_display_options(input)
          shiny::incProgress(0.15, detail = "Loading expression matrix")
          panel_set <- load_gene_tab_study_panels(study)
          if (is.null(panel_set) || length(panel_set$panels) < 1L) {
            panel_results(list(error = "Could not load expression data for this study."))
            return(invisible(NULL))
          }
          panels <- panel_set$panels
          matrix_key <- if (!is.null(pi$matrix_key) && nzchar(as.character(pi$matrix_key)) && !identical(pi$matrix_key, "__all__")) {
            as.character(pi$matrix_key)
          } else {
            NULL
          }
          if (!is.null(matrix_key)) {
            panels <- panels[vapply(panels, function(p) identical(p$key, matrix_key), logical(1L))]
          }
          if (length(panels) < 1L) {
            panel_results(list(error = "No expression matrix matched the current filter."))
            return(invisible(NULL))
          }
          shiny::incProgress(0.25, detail = "Matching genes")
          out <- lapply(panels, function(p) build_one_panel(p, t2g, pi, disp, study))
          shiny::incProgress(0.9, detail = "Finishing heatmap layout")
          ok_ht <- vapply(out, function(p) is.null(p$error) && !is.null(p$ht), logical(1L))
          ht_list <- lapply(out[ok_ht], function(p) p$ht)
          panel_results(list(
            error = NULL,
            panels = out,
            study = study,
            source = co$source %||% "",
            significant_only = isTRUE(disp$significant_only),
            fdr = disp$fdr,
            log2fc = disp$log2fc,
            column_aligned = .panel_columns_alignable(ht_list)
          ))
        }
      )
      invisible(NULL)
    }

    shiny::observeEvent(input$update, {
      run_generate()
    }, ignoreInit = TRUE)

    shiny::observeEvent(
      {
        pi <- if (is.function(panels_input)) panels_input() else NULL
        if (is.null(pi)) 0L else pi$generate %||% 0L
      },
      {
        gen <- if (is.function(panels_input)) panels_input()$generate else 0L
        if (is.null(gen) || gen < 1L) return()
        run_generate()
      },
      ignoreInit = TRUE
    )

    if (!is.null(generate_trigger) && is.function(generate_trigger)) {
      shiny::observeEvent(generate_trigger(), {
        trig <- generate_trigger()
        if (is.null(trig) || trig < 1L) return()
        run_generate()
      }, ignoreInit = TRUE)
    }

    output$status_ui <- shiny::renderUI({
      res <- panel_results()
      if (is.null(res)) {
        return(shiny::tags$div(
          class = "text-muted",
          style = "margin-bottom: 12px;",
          "Choose a study, load an ontology, and click ",
          shiny::tags$strong("Generate panel heatmap"),
          " in the sidebar."
        ))
      }
      if (isTRUE(res$loading)) {
        return(shiny::tags$div(
          class = "text-muted",
          style = "margin-bottom: 12px;",
          "Generating panel heatmap… (may take a few seconds on large matrices)"
        ))
      }
      if (!is.null(res$error) && nzchar(res$error)) {
        return(shiny::tags$div(class = "text-danger", style = "margin-bottom: 12px;", res$error))
      }
      src <- if (nzchar(res$source %||% "")) paste0(" (", res$source, ")") else ""
      lines <- lapply(res$panels, function(p) {
        title <- if (!is.null(p$title) && nzchar(p$title)) p$title else p$counts_file
        miss_n <- length(p$missing %||% character(0))
        if (!is.null(p$error)) {
          return(shiny::tags$li(paste0(title, ": ", p$error)))
        }
        parts <- character(0)
        if (!is.null(p$sig_n_after) && !is.na(p$sig_n_after)) {
          parts <- c(parts, paste0(p$sig_n_after, " significant of ", p$sig_n_before))
        }
        if (miss_n > 0L) {
          parts <- c(parts, paste0(miss_n, " symbols not in matrix"))
        }
        if (!is.null(p$n_samples) && !is.na(p$n_samples)) {
          parts <- c(parts, paste0(p$n_samples, " samples"))
        }
        shiny::tags$li(paste0(
          title, ": ", p$found, " genes plotted",
          if (length(parts) > 0L) paste0(" (", paste(parts, collapse = "; "), ")") else ""
        ))
      })
      note <- NULL
      if (isFALSE(res$column_aligned %||% TRUE)) {
        note <- shiny::tags$p(
          class = "text-muted",
          style = "margin: 6px 0 0 0; font-size: 0.9rem;",
          "This study has separate expression matrices per region with different sample counts. ",
          "They are shown as stacked heatmaps (not column-aligned). ",
          "Pick one matrix in the sidebar to view a single panel."
        )
      }
      shiny::tags$div(
        style = "margin-bottom: 12px;",
        shiny::tags$strong("Panel heatmaps"),
        src,
        shiny::tags$ul(style = "margin: 6px 0 0 18px;", lines),
        note
      )
    })

    output$heatmap_plot_ui <- shiny::renderUI({
      res <- panel_results()
      if (is.null(res) || isTRUE(res$loading) || !is.null(res$error) || is.null(res$panels)) {
        return(NULL)
      }
      has_ht <- any(vapply(res$panels, function(p) is.null(p$error) && !is.null(p$ht), logical(1L)))
      if (!has_ht) return(NULL)
      w_px <- .panel_plot_width_px(res$panels)
      h_px <- .panel_plot_height_px(res$panels)
      shiny::tags$div(
        style = "overflow-x: auto; width: 100%;",
        shiny::plotOutput(
          session$ns("panel_heatmaps"),
          width = paste0(w_px, "px"),
          height = paste0(h_px, "px")
        )
      )
    })

    output$panel_heatmaps <- shiny::renderPlot(
      {
        res <- panel_results()
        shiny::req(res, is.null(res$error), !isTRUE(res$loading), !is.null(res$panels))
        pds <- lapply(res$panels, function(p) if (is.null(p$error)) p$ht else NULL)
        shiny::req(any(!vapply(pds, is.null, logical(1))))
        .panel_draw_combined(pds)
      },
      res = 96
    )
  })
}
