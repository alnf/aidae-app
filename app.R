# app.R - Dashboard for DESeq2 results visualization
# Licensed under GPL-3.0. See LICENSE in the repository.
#
# Inspired by the Shiny app for visualizing DESeq2 results by Zuguang Gu
# (InteractiveComplexHeatmap vignette, MIT License):
# https://github.com/jokergoo/InteractiveComplexHeatmap/blob/master/vignettes/deseq2_app.Rmd

library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(circlize)
library(yaml)

source("scripts/study_data.R")
source("scripts/heatmap_utils.R")
source("scripts/perf_utils.R")

# Main app config (title, studies list; fallback if file missing)
main_config_path <- "config.yaml"
main_config <- if (file.exists(main_config_path)) {
  yaml::read_yaml(main_config_path)
} else {
  list(title = "Gene expression dashboard", studies = character(0))
}
if (is.null(main_config$title)) main_config$title <- "Gene expression dashboard"
if (is.null(main_config$studies)) main_config$studies <- character(0)

# Study list: from main config or scan data/ for subdirs with config.yaml
study_ids <- if (length(main_config$studies) > 0L) {
  main_config$studies
} else {
  data_dirs <- list.dirs("data", full.names = TRUE, recursive = FALSE)
  keep <- vapply(data_dirs, function(d) file.exists(file.path(d, "config.yaml")), logical(1L))
  basename(data_dirs[keep])
}

# Build study dropdown choices: value = study id (sent to server), name = display label
# Shiny selectInput: names(choices) = label shown, values(choices) = input$study value
study_labels <- character(length(study_ids))
for (i in seq_along(study_ids)) {
  id <- study_ids[i]
  cfg_path <- file.path("data", id, "config.yaml")
  if (file.exists(cfg_path)) {
    sc <- yaml::read_yaml(cfg_path)
    study_labels[i] <- if (!is.null(sc$name)) sc$name else id
  } else {
    study_labels[i] <- id
  }
}
study_choices <- setNames(study_ids, study_labels)
if (length(study_choices) == 0L) study_choices <- c("(no studies)" = "")

 

# Default study = first; default DEG list = first list of first study.
default_study <- if (length(study_ids) > 0L) study_ids[1L] else ""
first_study_lists <- study_deg_lists(default_study)
default_deg_choices <- if (length(first_study_lists) > 0L) {
  setNames(vapply(first_study_lists, function(x) x$deg_file, character(1L)),
           vapply(first_study_lists, function(x) x$label, character(1L)))
} else {
  c("(no DEG lists)" = "")
}
default_deg <- if (length(first_study_lists) > 0L) first_study_lists[[1L]]$deg_file else ""

# Brush action is defined in server and uses current study data from reactiveValues.
# It updates the MA-plot, volcano plot and result table for the selected genes.
library(DT)
library(GetoptLong)

# The dashboard body contains three columns:
# 1. the original heatmap
# 2. the sub-heatmap and the default output
# 3. the self-defined output
library(shiny)
library(bs4Dash)
library(shinymanager)

creds <- read.table("data/creds.txt", sep="\t", header = T)
credentials <- data.frame(
  user     = c(creds$user),
  password = c(creds$password),
  start    = c("2025-09-11"),
  expire   = c(NA),
  admin    = c(FALSE),
  stringsAsFactors = FALSE,
  is_hashed_password = TRUE
)


body <- dashboardBody(
  tags$style(HTML("
    .main-header .navbar-nav.ml-auto { display: none !important; }
    /* Heatmap toolbar: icons are Font Awesome (.fa) inside .nav-tabs */
    [id$='_heatmap_control'] .nav-tabs .fa,
    [id$='_heatmap_control'] .nav-tabs .fas,
    [id$='_heatmap_control'] .nav-tabs .far,
    [id$='_heatmap_control'] .nav-tabs .fab {
      font-size: 14px !important;
      -webkit-font-smoothing: antialiased;
    }
    [id$='_heatmap_control'] .nav-tabs > li {
      margin-left: 16px;
    }
    [id$='_heatmap_control'] .nav-tabs > li:first-child {
      margin-left: 0;
    }
    /* Align sidebar description with other inputs: remove indent from shiny-html-output */
    .main-sidebar .shiny-html-output {
      padding-left: 0;
      margin-left: 0;
    }
  ")),
  tabItems(
    tabItem(
      tabName = "degs",
      fluidRow(
        column(width = 4,
          box(title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
            originalHeatmapOutput("ht", height = 800, containment = TRUE)
          )
        ),
        column(width = 4,
          box(title = uiOutput("res_table_title"), width = NULL, solidHeader = TRUE, status = "primary",
            DTOutput("res_table")
          )
        ),
        column(width = 4,
          id = "column3",
          box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
            subHeatmapOutput("ht", title = NULL, containment = TRUE, height = 800)
          ),
          box(title = "MA-plot", width = NULL, solidHeader = TRUE, status = "primary",
            plotOutput("ma_plot")
          ),
          box(title = "Volcano plot", width = NULL, solidHeader = TRUE, status = "primary",
            plotOutput("volcano_plot")
          )
        ),
        tags$style("
          .content-wrapper, .right-side {
            overflow-x: auto;
          }
          .content {
            min-width:1500px;
          }
        ")
      )
    ),
    tabItem(
      tabName = "gene",
      box(title = "Gene", width = 12, solidHeader = TRUE, status = "secondary",
        p("Gene-level content to be added.")
      )
    )
  )
)

# Side bar: study selector, then DEG list, then cutoffs for significant genes.
ui <- secure_app(dashboardPage(
  title = main_config$title,
  fullscreen = FALSE,
  dark = FALSE,
  help = FALSE,
  header = dashboardHeader(
    title = main_config$title,
    navbarMenu(
      id = "navtabs",
      navbarTab(tabName = "degs", text = "DEGs"),
      navbarTab(tabName = "gene", text = "Gene")
    )
  ),
  sidebar = dashboardSidebar(
    minified = FALSE,
    selectInput("study", label = "Study", choices = study_choices, selected = default_study),
    selectInput("deg_list", label = "DEG list", choices = default_deg_choices, selected = default_deg),
    uiOutput("deg_description"),
    selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1), selected = 0.05),
    uiOutput("svalue_ui"),
    uiOutput("base_mean_ui"),
    numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 1),
    checkboxInput("show_rownames", label = "Show row names on heatmap", value = FALSE),
    actionButton("filter", label = "Generate heatmap"),
    br(),
    actionButton("select_genes_btn", label = "Select genes"),
    actionButton("clear_genes_btn", label = "Clear genes"),
    checkboxInput("lock_gene_list", label = "Apply selected genes across studies", value = TRUE)
  ),
  controlbar = dashboardControlbar(disable = TRUE),
  body = body
))

server <- function(input, output, session) {
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  rv <- reactiveValues(
    current_res = NULL,
    current_mm = NULL,
    current_col_annot = NULL,
    row_index = NULL,
    selected_rows = NULL,
    custom_genes = NULL,
    custom_genes_study = NULL
  )

  # Dynamic title for result table: all genes by default, selected genes after sub-heatmap selection
  output$res_table_title <- renderText({
    if (!is.null(rv$selected_rows) && length(rv$selected_rows) > 0) {
      "Result table of the selected genes"
    } else {
      "Result table of all genes"
    }
  })

  # Helper: apply current custom gene list (rv$custom_genes) to the loaded data
  # If apply_thresholds = FALSE, show all selected genes without numeric filtering.
  # If apply_thresholds = TRUE, apply current FDR/baseMean/log2FC/svalue thresholds
  # within the selected subset.
  apply_custom_gene_list <- function(apply_thresholds = TRUE) {
    genes_vec <- rv$custom_genes
    if (is.null(genes_vec)) return()
    genes_vec <- unique(trimws(genes_vec))
    genes_vec <- genes_vec[nzchar(genes_vec)]
    if (length(genes_vec) == 0) return()
    res <- rv$current_res
    mm <- rv$current_mm
    if (is.null(res) || is.null(mm) || !("symbol" %in% colnames(res))) return()
    sel <- tolower(res$symbol) %in% tolower(genes_vec)
    if (!any(sel)) {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("No genes from the list were found in the current DEG table.")
      })
      rv$row_index <- NULL
      rv$selected_rows <- NULL
      return()
    }
    # Subset DEGs and matrix to the selected genes only
    res_sub <- res[sel, , drop = FALSE]
    mm_sub <- mm[res_sub$ens_gene, , drop = FALSE]

    if (isTRUE(apply_thresholds)) {
      # Apply numeric thresholds within the selected genes only
      sval <- if ("svalue" %in% colnames(res_sub)) as.numeric(if (!is.null(input$svalue)) input$svalue else 0.005) else 0.005
      bmean <- if ("baseMean" %in% colnames(res_sub)) as.numeric(if (!is.null(input$base_mean)) input$base_mean else 20) else 0
      message(
        "[perf] make_heatmap (custom gene list, with thresholds): start; study=", input$study,
        " deg_list=", input$deg_list,
        " n_genes_input=", length(genes_vec),
        " genes_input=", paste(genes_vec, collapse = ", "),
        " n_genes_matched=", nrow(res_sub),
        " nrow(res)=", nrow(res),
        " fdr=", input$fdr,
        " base_mean=", bmean,
        " log2fc=", input$log2fc,
        " svalue=", sval
      )
      out <- perf_time(
        "make_heatmap_custom_genes",
        make_heatmap(
          res_sub,
          mm_sub,
          fdr = as.numeric(input$fdr),
          base_mean = bmean,
          log2fc = input$log2fc,
          svalue = sval,
          col_annot = rv$current_col_annot,
          show_row_names = isTRUE(input$show_rownames)
        )
      )
    } else {
      # Show all selected genes without any numeric thresholds
      message(
        "[perf] make_heatmap (custom gene list, no thresholds): start; study=", input$study,
        " deg_list=", input$deg_list,
        " n_genes_input=", length(genes_vec),
        " genes_input=", paste(genes_vec, collapse = ", "),
        " n_genes_matched=", nrow(res_sub),
        " nrow(res)=", nrow(res)
      )
      out <- perf_time(
        "make_heatmap_custom_genes_no_thresholds",
        make_heatmap(
          res_sub,
          mm_sub,
          fdr = 1,
          base_mean = 0,
          log2fc = 0,
          svalue = 1,
          col_annot = rv$current_col_annot,
          show_row_names = isTRUE(input$show_rownames)
        )
      )
    }
    if (is.null(out)) {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("No row exists after processing the custom gene list.")
      })
      rv$row_index <- NULL
      rv$selected_rows <- NULL
      return()
    }
    # Map row indices back to the original DEG table
    selected_idx_full <- which(sel)
    rv$row_index <- selected_idx_full[out$row_index]
    rv$selected_rows <- NULL
    render_res_table(rv$row_index)
    message("[perf] makeInteractiveComplexHeatmap (custom gene list): start; study=", input$study, " deg_list=", input$deg_list)
    perf_time(
      "makeInteractiveComplexHeatmap_custom_genes",
      makeInteractiveComplexHeatmap(
        input, output, session, out$ht, "ht",
        brush_action = brush_action
      )
    )
  }

  # Open modal to input custom gene list
  observeEvent(input$select_genes_btn, {
    showModal(
      modalDialog(
        title = "Select genes",
        textAreaInput(
          "gene_list_input",
          label = "Gene symbols (one per line):",
          value = "",
          rows = 10,
          placeholder = "TP53\nBRCA1\nEGFR"
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("apply_gene_list", "Apply")
        ),
        easyClose = TRUE
      )
    )
  })

  # Helper to render result table for either all or a subset of genes
  render_res_table <- function(selected_idx = NULL) {
    res <- rv$current_res
    if (is.null(res)) return()
    tbl_cols <- c("symbol", "baseMean", "log2FoldChange", "padj")
    tbl_cols <- intersect(tbl_cols, colnames(res))
    num_idx <- which(tbl_cols %in% c("baseMean", "log2FoldChange", "padj"))
    rows <- if (!is.null(selected_idx)) selected_idx else seq_len(nrow(res))
    output[["res_table"]] <- renderDT(
      formatRound(
        datatable(
          res[rows, tbl_cols, drop = FALSE],
          rownames = FALSE,
          options = list(
            pageLength = 20,
            lengthMenu = list(c(20, 50, 100), c("20", "50", "100"))
          )
        ),
        columns = num_idx,
        digits = 3
      )
    )
  }

  # Keep result table in sync with the full heatmap when there is no sub-heatmap selection
  observe({
    if (!is.null(rv$current_res) &&
        !is.null(rv$row_index) &&
        length(rv$row_index) > 0 &&
        (is.null(rv$selected_rows) || length(rv$selected_rows) == 0)) {
      render_res_table(rv$row_index)
    }
  })

  # Apply custom gene list (case-insensitive match on symbol; ignore numeric thresholds)
  observeEvent(input$apply_gene_list, {
    isolate({
      genes_raw <- input$gene_list_input
      removeModal()
      if (is.null(genes_raw)) return()
      genes_vec <- unique(trimws(unlist(strsplit(genes_raw, "[\r\n]+"))))
      genes_vec <- genes_vec[nzchar(genes_vec)]
      if (length(genes_vec) == 0) return()
      rv$custom_genes <- genes_vec
      rv$custom_genes_study <- input$study
      # First display: show all selected genes, without thresholds
      apply_custom_gene_list(apply_thresholds = FALSE)
    })
  })

  # Clear current custom gene list and revert to threshold-based heatmap for the current study/DEG list
  observeEvent(input$clear_genes_btn, {
    rv$custom_genes <- NULL
    rv$custom_genes_study <- NULL
    rv$selected_rows <- NULL
    # Regenerate heatmap using current numeric thresholds
    if (is.null(rv$current_res) || is.null(rv$current_mm)) {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("Select a study and load data.")
      })
      return()
    }
    # Trigger the same logic as in the main heatmap regeneration observer
    sval <- if ("svalue" %in% colnames(rv$current_res)) as.numeric(if (!is.null(input$svalue)) input$svalue else 0.005) else 0.005
    bmean <- if ("baseMean" %in% colnames(rv$current_res)) as.numeric(if (!is.null(input$base_mean)) input$base_mean else 20) else 0
    message(
      "[perf] make_heatmap (after clear genes): start; study=", input$study,
      " deg_list=", input$deg_list,
      " fdr=", input$fdr,
      " base_mean=", bmean,
      " log2fc=", input$log2fc,
      " svalue=", sval
    )
    out <- perf_time(
      "make_heatmap_after_clear_genes",
      make_heatmap(
        rv$current_res, rv$current_mm,
        fdr = as.numeric(input$fdr), base_mean = bmean, log2fc = input$log2fc, svalue = sval,
        col_annot = rv$current_col_annot,
        show_row_names = isTRUE(input$show_rownames)
      )
    )
    if (!is.null(out)) {
      rv$row_index <- out$row_index
      # By default, show all genes (within the filtered/heatmap subset) in the result table
      if (!is.null(rv$row_index)) {
        render_res_table(rv$row_index)
      } else {
        render_res_table(NULL)
      }
      message("[perf] makeInteractiveComplexHeatmap (after clear genes): start; study=", input$study, " deg_list=", input$deg_list)
      perf_time(
        "makeInteractiveComplexHeatmap_after_clear_genes",
        makeInteractiveComplexHeatmap(
          input, output, session, out$ht, "ht",
          brush_action = brush_action
        )
      )
    } else {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  })

  # When study changes, update DEG list dropdown to that study's lists (first selected).
  observeEvent(input$study, {
    rv$current_res <- NULL
    rv$current_mm <- NULL
    rv$row_index <- NULL
    rv$selected_rows <- NULL
    if (is.null(input$study) || input$study == "") return()
    lists <- study_deg_lists(input$study)
    if (length(lists) == 0L) {
      updateSelectInput(session, "deg_list", choices = c("(no DEG lists)" = ""), selected = "")
      return()
    }
    choices <- setNames(vapply(lists, function(x) x$deg_file, character(1L)),
                       vapply(lists, function(x) x$label, character(1L)))
    updateSelectInput(session, "deg_list", choices = choices, selected = lists[[1L]]$deg_file)
  }, ignoreNULL = FALSE)

  # Show svalue cutoff only when the loaded DEG table has a svalue column
  output$svalue_ui <- renderUI({
    res <- rv$current_res
    if (is.null(res) || !("svalue" %in% colnames(res))) return(NULL)
    numericInput("svalue", label = "Cutoff for svalue:", value = 0.005)
  })

  # Show base mean cutoff only when the loaded DEG table has a baseMean column
  output$base_mean_ui <- renderUI({
    res <- rv$current_res
    if (is.null(res) || !("baseMean" %in% colnames(res))) return(NULL)
    numericInput("base_mean", label = "Minimal base mean:", value = 20)
  })

  # Description for the selected DEG list (from config)
  output$deg_description <- renderUI({
    if (is.null(input$study) || input$study == "" || is.null(input$deg_list) || input$deg_list == "") return(NULL)
    lists <- study_deg_lists(input$study)
    idx <- match(input$deg_list, vapply(lists, function(x) x$deg_file, character(1L)))
    desc <- if (!is.na(idx) && !is.null(lists[[idx]]$description)) lists[[idx]]$description else ""
    if (desc == "") return(NULL)
    tags$div(
      class = "form-group shiny-input-container",
      tags$label("Description", class = "control-label"),
      tags$div(desc, class = "text-muted", style = "margin-top: 0.25rem; font-size: 0.9rem;")
    )
  })

  # Load data when study or DEG list selection changes (runs on init so default study+list load).
  observeEvent(list(input$study, input$deg_list), {
    rv$current_res <- NULL
    rv$current_mm <- NULL
    rv$current_col_annot <- NULL
    rv$row_index <- NULL
    rv$selected_rows <- NULL
    if (is.null(input$study) || input$study == "" || is.null(input$deg_list) || input$deg_list == "") return()
    message("[perf] load_study_data: start; study=", input$study, " deg_list=", input$deg_list)
    loaded <- perf_time("load_study_data", load_study_data(input$study, input$deg_list))
    rv$current_res <- loaded$res
    rv$current_mm <- loaded$mm
    rv$current_col_annot <- loaded$col_annot
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  # Brush action uses current study data from rv
  brush_action <- function(df, input, output, session) {
    res <- rv$current_res
    row_index <- rv$row_index
    if (is.null(res) || is.null(row_index)) return()
    idx <- unique(unlist(df$row_index))
    selected <- row_index[idx]
    rv$selected_rows <- selected
    output[["ma_plot"]] <- renderPlot({
      make_maplot(res, selected)
    })
    output[["volcano_plot"]] <- renderPlot({
      make_volcano(res, selected)
    })
    render_res_table(selected)
  }

  # Regenerate heatmap when filter is clicked or study/DEG list changes
  observeEvent(list(input$filter, input$study, input$deg_list), {
    perf_time(
      "regenerate_heatmap_observer",
      {
        # Use custom gene list when available and either locked across studies
        # or when viewing the study where it was defined; otherwise fall back
        # to threshold-based filtering.
        use_custom <- !is.null(rv$custom_genes) &&
          (isTRUE(input$lock_gene_list) ||
             (!is.null(rv$custom_genes_study) && identical(input$study, rv$custom_genes_study)))
        if (use_custom) {
          # On Generate heatmap, re-apply custom genes with thresholds
          apply_custom_gene_list(apply_thresholds = TRUE)
        } else {
          if (is.null(rv$current_res) || is.null(rv$current_mm)) {
            output$ht_heatmap <- renderPlot({
              grid.newpage()
              grid.text("Select a study and load data.")
            })
            return()
          }
          sval <- if ("svalue" %in% colnames(rv$current_res)) as.numeric(if (!is.null(input$svalue)) input$svalue else 0.005) else 0.005
          bmean <- if ("baseMean" %in% colnames(rv$current_res)) as.numeric(if (!is.null(input$base_mean)) input$base_mean else 20) else 0
          message(
            "[perf] make_heatmap: start; study=", input$study,
            " deg_list=", input$deg_list,
            " fdr=", input$fdr,
            " base_mean=", bmean,
            " log2fc=", input$log2fc,
            " svalue=", sval
          )
          out <- perf_time(
            "make_heatmap",
            make_heatmap(
              rv$current_res, rv$current_mm,
              fdr = as.numeric(input$fdr), base_mean = bmean, log2fc = input$log2fc, svalue = sval,
              col_annot = rv$current_col_annot,
              show_row_names = isTRUE(input$show_rownames)
            )
          )
          if (!is.null(out)) {
            rv$row_index <- out$row_index
            rv$selected_rows <- NULL
            # By default, show all genes (within the filtered/heatmap subset) in the result table
            if (!is.null(rv$row_index)) {
              render_res_table(rv$row_index)
            } else {
              render_res_table(NULL)
            }
            message("[perf] makeInteractiveComplexHeatmap: start; study=", input$study, " deg_list=", input$deg_list)
            perf_time(
              "makeInteractiveComplexHeatmap",
              makeInteractiveComplexHeatmap(
                input, output, session, out$ht, "ht",
                brush_action = brush_action
              )
            )
          } else {
            output$ht_heatmap <- renderPlot({
              grid.newpage()
              grid.text("No row exists after filtering.")
            })
          }
        }
      }
    )
  }, ignoreNULL = FALSE)
}

shinyApp(ui, server)