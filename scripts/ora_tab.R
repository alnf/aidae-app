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

oraTabUI <- function(id, study_ids, study_labels, ora_file_choices, pathway_default = "") {
  ns <- shiny::NS(id)
  sel <- if (nzchar(as.character(pathway_default)) && as.character(pathway_default) %in% unname(ora_file_choices)) {
    as.character(pathway_default)
  } else if (length(ora_file_choices)) {
    unname(ora_file_choices)[[1L]]
  } else {
    NULL
  }
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::tags$div(
          class = "ora-pathway-db-select",
          style = "max-width: 260px;",
          shiny::selectInput(
            ns("pathway_file"),
            "Pathway database:",
            choices = ora_file_choices,
            selected = sel,
            width = "100%"
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
        shiny::tags$div(
          class = "ora-plot-wrap",
          style = "overflow-x: auto; width: 100%; max-width: 100%; min-width: 0;",
          shiny::tags$div(
            style = "display: inline-block; vertical-align: top; max-width: none;",
            shiny::plotOutput(ns("ora_plot_all"), width = "auto", height = "auto")
          )
        )
      )
    )
  )
}

#' @param ora_input reactive: list(min_overlap, min_count, min_gene_ratio, show_category)
oraTabServer <- function(id, study_ids, study_labels, ora_input) {
  shiny::moduleServer(id, function(input, output, session) {
    # Heavy: cache load or enricher — depends on pathway file only (not sidebar thresholds).
    ora_long_source <- shiny::reactive({
      shiny::req(!is.null(input$pathway_file))
      shiny::req(nzchar(as.character(input$pathway_file)))

      pathway_file <- as.character(input$pathway_file)

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
      n_show <- suppressWarnings(as.integer(oi$show_category))
      if (is.na(n_show) || n_show < 1L) n_show <- 20L
      n_show <- min(100L, max(1L, n_show))

      src <- ora_long_source()
      if (!is.null(src$error)) {
        return(list(error = src$error, by_study = NULL, pathway_level_order = NULL))
      }
      ora_build_plot_payload(src$long_by_sid, study_ids, study_labels, min_ol, min_ct, min_gr, n_show)
    })

    study_label_for <- function(sid) {
      lbl <- study_labels[[sid]]
      if (is.null(lbl) || !nzchar(as.character(lbl))) as.character(sid) else as.character(lbl)
    }

    output$ora_plot_all <- shiny::renderPlot(
      {
        ob <- ora_by_study()
        if (!is.null(ob$error)) {
          return(.ora_msg_plot(ob$error))
        }

        dfs <- list()
        comp_per_panel <- integer(0)
        max_path <- 1L
        first_msg <- NULL
        for (sid in study_ids) {
          b <- ob$by_study[[sid]]
          if (is.null(b)) next
          if (!is.null(b$msg)) {
            if (is.null(first_msg)) first_msg <- b$msg
            next
          }
          if (!is.null(b$plot_df) && nrow(b$plot_df) > 0L) {
            dfs[[length(dfs) + 1L]] <- b$plot_df
            comp_per_panel <- c(comp_per_panel, length(unique(b$plot_df$comparison)))
            max_path <- max(max_path, length(unique(as.character(b$plot_df$ID))))
          }
        }
        if (!is.null(ob$pathway_level_order) && length(ob$pathway_level_order) > 0L) {
          max_path <- max(max_path, length(ob$pathway_level_order))
        }

        if (length(dfs) == 0L) {
          return(.ora_msg_plot(if (!is.null(first_msg)) first_msg else "No pathways to display."))
        }

        combined <- do.call(rbind, dfs)
        if (!is.null(ob$pathway_level_order) && length(ob$pathway_level_order) > 0L) {
          combined$Description <- factor(
            as.character(combined$Description),
            levels = ob$pathway_level_order
          )
        }
        all_lbls <- vapply(study_ids, study_label_for, character(1L))
        present <- unique(as.character(combined$study_label))
        study_label_levels <- all_lbls[all_lbls %in% present]

        tryCatch(
          .ora_faceted_comparison_plot(combined, study_label_levels, ob$pathway_level_order),
          error = function(err) {
            .ora_msg_plot(paste("Plot:", conditionMessage(err)))
          }
        )
      },
      width = function() {
        ob <- ora_by_study()
        if (!is.null(ob$error)) return(.ora_msg_plot_px()$width)
        comp_per_panel <- integer(0)
        max_path <- 1L
        for (sid in study_ids) {
          b <- ob$by_study[[sid]]
          if (is.null(b) || !is.null(b$msg) || is.null(b$plot_df) || nrow(b$plot_df) < 1L) next
          comp_per_panel <- c(comp_per_panel, length(unique(b$plot_df$comparison)))
          max_path <- max(max_path, length(unique(as.character(b$plot_df$ID))))
        }
        if (!is.null(ob$pathway_level_order) && length(ob$pathway_level_order) > 0L) {
          max_path <- max(max_path, length(ob$pathway_level_order))
        }
        if (length(comp_per_panel) == 0L) return(.ora_msg_plot_px()$width)
        .ora_faceted_plot_dims(comp_per_panel, max_path)$width
      },
      height = function() {
        ob <- ora_by_study()
        if (!is.null(ob$error)) return(.ora_msg_plot_px()$height)
        comp_per_panel <- integer(0)
        max_path <- 1L
        for (sid in study_ids) {
          b <- ob$by_study[[sid]]
          if (is.null(b) || !is.null(b$msg) || is.null(b$plot_df) || nrow(b$plot_df) < 1L) next
          comp_per_panel <- c(comp_per_panel, length(unique(b$plot_df$comparison)))
          max_path <- max(max_path, length(unique(as.character(b$plot_df$ID))))
        }
        if (!is.null(ob$pathway_level_order) && length(ob$pathway_level_order) > 0L) {
          max_path <- max(max_path, length(ob$pathway_level_order))
        }
        if (length(comp_per_panel) == 0L) return(.ora_msg_plot_px()$height)
        .ora_faceted_plot_dims(comp_per_panel, max_path)$height
      }
    )
  })
}
