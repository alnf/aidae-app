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

source("scripts/heatmap_utils.R")
source("scripts/result_table_indices.R")
source("scripts/study_data.R")
source("scripts/perf_utils.R")
source("scripts/gene_plot.R")
source("scripts/gene_tab.R")
source("scripts/pathway_signatures.R")
source("scripts/ora_cache.R")
source("scripts/ora_tab.R")

# UI defaults align with make_heatmap() / default_heatmap_thresholds(); study+DEG list can override via config.
default_thr <- default_heatmap_thresholds()

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

pathway_txt_files <- list_ora_pathway_files(study_ids)
# Start empty so ORA does not load until the user picks a pathway database.
ora_pathway_default <- ""
ora_file_choices <- if (length(pathway_txt_files) > 0L) {
  stats::setNames(
    pathway_txt_files,
    gsub("_", " ", tools::file_path_sans_ext(pathway_txt_files), fixed = TRUE)
  )
} else {
  c("(no pathway files in databases/pathways)" = "")
}

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
  tags$head(
    tags$script(HTML("
      Shiny.addCustomMessageHandler('exprs_copy_text', function(text) {
        if (text === null || text === undefined) return;
        var s = String(text);
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(s).catch(function(e) { console.error(e); });
        } else {
          var ta = document.createElement('textarea');
          ta.value = s;
          ta.style.position = 'fixed';
          ta.style.left = '-9999px';
          document.body.appendChild(ta);
          ta.select();
          try { document.execCommand('copy'); } catch (err) { console.error(err); }
          document.body.removeChild(ta);
        }
      });
    "))
  ),
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
    /* withProgress() bar: default Shiny placement is bottom; anchor to top of viewport */
    .shiny-progress {
      position: fixed !important;
      top: 0 !important;
      bottom: auto !important;
      left: 0 !important;
      right: 0 !important;
      width: 100% !important;
      padding: 10px 16px 8px 16px !important;
      margin: 0 !important;
      z-index: 2000 !important;
      background: rgba(255, 255, 255, 0.97) !important;
      box-shadow: 0 1px 4px rgba(0, 0, 0, 0.12) !important;
    }
    .shiny-progress .progress {
      margin-bottom: 4px !important;
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
          box(
            title = uiOutput("res_table_title"),
            width = NULL, solidHeader = TRUE, status = "primary",
            actionButton("copy_table_symbols", label = "Copy gene symbols", class = "btn-sm btn-outline-secondary mb-2"),
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
      box(
        title = "Gene expression by study",
        width = 12, solidHeader = TRUE, status = "secondary",
        geneTabUI("gene", study_ids, stats::setNames(study_labels, study_ids))
      )
    ),
    tabItem(
      tabName = "ora",
      box(
        title = "Overrepresentation analysis (ORA)",
        width = 12, solidHeader = TRUE, status = "secondary",
        oraTabUI("ora", study_ids, stats::setNames(study_labels, study_ids), ora_file_choices, ora_pathway_default)
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
      navbarTab(tabName = "gene", text = "Gene"),
      navbarTab(tabName = "ora", text = "ORA")
    )
  ),
  sidebar = dashboardSidebar(
    minified = FALSE,
    # DEGs tab only: heatmap / table controls.
    shiny::conditionalPanel(
      condition = "input.navtabs == 'degs'",
      selectInput("study", label = "Study", choices = study_choices, selected = default_study),
      selectInput("deg_list", label = "DEG list", choices = default_deg_choices, selected = default_deg),
      uiOutput("deg_description"),
      selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1), selected = default_thr$fdr),
      uiOutput("svalue_ui"),
      uiOutput("base_mean_ui"),
      numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = default_thr$log2fc),
      checkboxInput("show_rownames", label = "Show row names on heatmap", value = FALSE),
      actionButton("filter", label = "Generate heatmap"),
      br(),
      actionButton("select_genes_btn", label = "Select genes"),
      actionButton("clear_genes_btn", label = "Clear genes"),
      checkboxInput("lock_gene_list", label = "Apply selected genes across studies", value = TRUE)
    ),
    shiny::conditionalPanel(
      condition = "input.navtabs == 'ora'",
      numericInput("ora_min_overlap", label = "Minimum pathway size (minGSSize):", value = 10L, min = 1L, step = 1L),
      numericInput("ora_min_count", label = "Minimum overlap (Count):", value = 5L, min = 1L, step = 1L),
      numericInput("ora_min_gene_ratio", label = "Minimum gene ratio:", value = 0.1, min = 0, max = 1, step = 0.01),
      numericInput("ora_show_category", label = "Max pathways to show:", value = 20L, min = 1L, step = 1L),
      radioButtons(
        "ora_heatmap_click_target",
        "Dot click interprets as:",
        choices = c("Pathway (row intent)" = "pathway", "Comparison (column intent)" = "comparison"),
        selected = "pathway"
      ),
      shiny::conditionalPanel(
        condition = "input.navtabs == 'ora' && input.ora_heatmap_click_target == 'pathway'",
        checkboxInput(
          "ora_pathway_hide_empty_comparisons",
          label = "Hide comparisons without significant pathways",
          value = TRUE
        )
      ),
      tags$div(
        class = "text-muted",
        style = "padding: 8px 0; font-size: 0.9rem;",
        "Gene sets use thresholds from config (study-level and per DEG list), or app defaults — not the DEGs sidebar sliders."
      )
    ),
    shiny::conditionalPanel(
      condition = "input.navtabs == 'gene'",
      tags$div(
        class = "text-muted",
        style = "padding: 10px 12px; font-size: 0.9rem;",
        "Study, DEG list, and heatmap filters apply to the DEGs tab. Choose a gene in the main panel."
      )
    )
  ),
  controlbar = dashboardControlbar(disable = TRUE),
  body = body
))

server <- function(input, output, session) {
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  gene_jump_symbol <- shiny::reactiveVal(NULL)

  geneTabServer(
    "gene",
    study_ids,
    stats::setNames(study_labels, study_ids),
    external_symbol = shiny::reactive(gene_jump_symbol())
  )

  rv <- reactiveValues(
    current_res = NULL,
    current_mm = NULL,
    current_col_annot = NULL,
    row_index = NULL,
    selected_rows = NULL,
    custom_genes = NULL,
    custom_genes_study = NULL,
    threshold_defaults = default_heatmap_thresholds(),
    deg_snapshot = list(
      study = default_study,
      deg_list = default_deg,
      fdr = default_thr$fdr,
      log2fc = default_thr$log2fc,
      base_mean = default_thr$base_mean,
      svalue = default_thr$svalue,
      lock_gene_list = TRUE
    )
  )

  shiny::observe({
    shiny::req(input$navtabs == "degs")
    shiny::req(input$study, input$deg_list)
    if (is.null(input$study) || input$study == "" || is.null(input$deg_list) || input$deg_list == "") {
      return()
    }
    d <- default_heatmap_thresholds()
    rv$deg_snapshot <- list(
      study = input$study,
      deg_list = input$deg_list,
      fdr = as.numeric(input$fdr),
      log2fc = input$log2fc,
      base_mean = if (!is.null(input$base_mean)) as.numeric(input$base_mean) else d$base_mean,
      svalue = if (!is.null(input$svalue)) as.numeric(input$svalue) else d$svalue,
      lock_gene_list = isTRUE(input$lock_gene_list)
    )
  })

  ora_input <- shiny::reactive({
    mo <- if (is.null(input$ora_min_overlap)) 10L else input$ora_min_overlap
    mc <- if (is.null(input$ora_min_count)) 5L else input$ora_min_count
    mgr <- if (is.null(input$ora_min_gene_ratio)) 0.1 else as.numeric(input$ora_min_gene_ratio)
    if (is.na(mgr) || mgr < 0) mgr <- 0
    if (mgr > 1) mgr <- 1
    sc <- if (is.null(input$ora_show_category)) 20L else input$ora_show_category
    click_target <- if (is.null(input$ora_heatmap_click_target)) "pathway" else as.character(input$ora_heatmap_click_target)
    if (!click_target %in% c("pathway", "comparison")) click_target <- "pathway"
    hide_empty <- if (is.null(input$ora_pathway_hide_empty_comparisons)) TRUE else isTRUE(input$ora_pathway_hide_empty_comparisons)
    list(
      min_overlap = mo,
      min_count = mc,
      min_gene_ratio = mgr,
      show_category = sc,
      click_target = click_target,
      hide_empty_pathway_comparisons = hide_empty
    )
  })

  oraTabServer(
    "ora",
    study_ids,
    stats::setNames(study_labels, study_ids),
    ora_input,
    on_gene_select = function(symbol) {
      sym <- trimws(as.character(symbol))
      if (!nzchar(sym)) return(invisible(NULL))
      gene_jump_symbol(NULL)
      gene_jump_symbol(sym)
      if (requireNamespace("shinydashboard", quietly = TRUE)) {
        shinydashboard::updateTabItems(session, "navtabs", selected = "gene")
      }
      invisible(NULL)
    }
  )

  # Dynamic title for result table: threshold-filtered genes, or sub-heatmap selection
  output$res_table_title <- renderText({
    if (!is.null(rv$selected_rows) && length(rv$selected_rows) > 0) {
      "Result table of the selected genes"
    } else {
      "Result table of genes passing current thresholds"
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
      sval <- if ("svalue" %in% colnames(res_sub)) as.numeric(if (!is.null(input$svalue)) input$svalue else default_heatmap_thresholds()$svalue) else default_heatmap_thresholds()$svalue
      bmean <- if ("baseMean" %in% colnames(res_sub)) as.numeric(if (!is.null(input$base_mean)) input$base_mean else default_heatmap_thresholds()$base_mean) else 0
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

  # Result table: selected genes from brush, else genes passing current thresholds (and custom list when active)
  render_res_table <- function(selected_idx) {
    res <- rv$current_res
    if (is.null(res)) return()
    tbl_cols <- c("symbol", "baseMean", "log2FoldChange", "padj")
    tbl_cols <- intersect(tbl_cols, colnames(res))
    num_idx <- which(tbl_cols %in% c("baseMean", "log2FoldChange", "padj"))
    rows <- if (length(selected_idx) == 0L) integer(0) else selected_idx
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

  # Same row set as the result table (brush selection or threshold-filtered rows)
  current_res_table_row_indices <- function() {
    res <- rv$current_res
    mm <- rv$current_mm
    if (is.null(res) || is.null(mm)) return(NULL)
    sval <- if ("svalue" %in% colnames(res)) {
      as.numeric(if (!is.null(input$svalue)) input$svalue else default_heatmap_thresholds()$svalue)
    } else {
      default_heatmap_thresholds()$svalue
    }
    bmean <- if ("baseMean" %in% colnames(res)) {
      as.numeric(if (!is.null(input$base_mean)) input$base_mean else default_heatmap_thresholds()$base_mean)
    } else {
      0
    }
    result_table_row_indices_for_study(
      res, mm,
      active_study_id = input$study,
      this_study_id = input$study,
      selected_rows = rv$selected_rows,
      custom_genes = rv$custom_genes,
      custom_genes_study = rv$custom_genes_study,
      lock_gene_list = input$lock_gene_list,
      fdr = as.numeric(input$fdr),
      log2fc = input$log2fc,
      base_mean = bmean,
      svalue = sval
    )
  }

  observe({
    idx <- current_res_table_row_indices()
    if (is.null(idx)) return()
    render_res_table(idx)
  })

  observeEvent(input$copy_table_symbols, {
    idx <- current_res_table_row_indices()
    res <- rv$current_res
    if (is.null(idx) || is.null(res) || !("symbol" %in% colnames(res))) {
      showNotification("No genes to copy.", type = "warning")
      return()
    }
    if (length(idx) == 0L) {
      showNotification("No genes to copy.", type = "warning")
      return()
    }
    syms <- as.character(res$symbol[idx])
    syms <- syms[!is.na(syms) & nzchar(syms)]
    if (length(syms) == 0L) {
      showNotification("No gene symbols to copy.", type = "warning")
      return()
    }
    session$sendCustomMessage("exprs_copy_text", paste(syms, collapse = "\n"))
    showNotification(paste(length(syms), "gene symbol(s) copied to clipboard."), type = "message")
  })

  # Apply custom gene list (case-insensitive match on symbol; heatmap updates on Generate heatmap)
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
      rv$row_index <- NULL
      rv$selected_rows <- NULL
      if (is.null(rv$current_res) || is.null(rv$current_mm)) {
        output$ht_heatmap <- renderPlot({
          grid.newpage()
          grid.text("Select a study and load data.")
        })
      } else {
        output$ht_heatmap <- renderPlot({
          grid.newpage()
          grid.text("Gene list saved. Click \"Generate heatmap\" to display the heatmap.")
        })
      }
      showNotification("Gene list saved. Click \"Generate heatmap\" to update the view.", type = "message")
    })
  })

  # Clear current custom gene list; heatmap updates on Generate heatmap
  observeEvent(input$clear_genes_btn, {
    rv$custom_genes <- NULL
    rv$custom_genes_study <- NULL
    rv$selected_rows <- NULL
    rv$row_index <- NULL
    if (is.null(rv$current_res) || is.null(rv$current_mm)) {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("Select a study and load data.")
      })
      return()
    }
    output$ht_heatmap <- renderPlot({
      grid.newpage()
      grid.text("Custom gene list cleared. Click \"Generate heatmap\" to update the view.")
    })
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
    d <- rv$threshold_defaults
    val <- if (!is.null(d$svalue)) d$svalue else default_heatmap_thresholds()$svalue
    numericInput("svalue", label = "Cutoff for svalue:", value = val)
  })

  # Show base mean cutoff only when the loaded DEG table has a baseMean column
  output$base_mean_ui <- renderUI({
    res <- rv$current_res
    if (is.null(res) || !("baseMean" %in% colnames(res))) return(NULL)
    d <- rv$threshold_defaults
    val <- if (!is.null(d$base_mean)) d$base_mean else default_heatmap_thresholds()$base_mean
    numericInput("base_mean", label = "Minimal base mean:", value = val)
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
    d <- deg_list_threshold_defaults(input$study, input$deg_list)
    rv$threshold_defaults <- d
    fdr_choices <- c(0.001, 0.01, 0.05, 0.1)
    fdr_sel <- d$fdr
    if (!fdr_sel %in% fdr_choices) {
      message(
        "[config] fdr=", fdr_sel, " not in UI choices ",
        paste(fdr_choices, collapse = ", "),
        "; using ", default_heatmap_thresholds()$fdr
      )
      fdr_sel <- default_heatmap_thresholds()$fdr
    }
    updateSelectInput(session, "fdr", selected = fdr_sel)
    updateNumericInput(session, "log2fc", value = d$log2fc)
    output$ht_heatmap <- renderPlot({
      grid.newpage()
      grid.text("Click \"Generate heatmap\" to display the heatmap.")
    })
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
  }

  # Regenerate heatmap only when Generate heatmap is clicked
  observeEvent(input$filter, {
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
          sval <- if ("svalue" %in% colnames(rv$current_res)) as.numeric(if (!is.null(input$svalue)) input$svalue else default_heatmap_thresholds()$svalue) else default_heatmap_thresholds()$svalue
          bmean <- if ("baseMean" %in% colnames(rv$current_res)) as.numeric(if (!is.null(input$base_mean)) input$base_mean else default_heatmap_thresholds()$base_mean) else 0
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
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
}

shinyApp(ui, server)