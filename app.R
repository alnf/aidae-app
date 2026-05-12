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
source("scripts/upset_tab.R")

# UI defaults align with make_heatmap() / default_heatmap_thresholds(); study+DEG list can override via config.
default_thr <- default_heatmap_thresholds()

# Main app config (title, studies list; fallback if file missing).
# On Posit Connect, set EXPRS_MAIN_CONFIG to a deploy YAML (e.g. deploy/apps/heart.yaml).
main_config_path_env <- Sys.getenv("EXPRS_MAIN_CONFIG", unset = "")

# Resolve EXPRS_MAIN_CONFIG more defensively. On Posit Connect, working dirs and
# relative paths can differ; this makes sure we still find `deploy/apps/*.yaml`.
resolve_main_config_path <- function(p) {
  if (is.null(p) || !nzchar(p)) return("config.yaml")
  candidates <- unique(c(
    p,
    # If user gives "heart.yaml" / "heart.yml"
    file.path("deploy", "apps", basename(p)),
    # If user gives just "heart"
    if (!grepl("\\.ya?ml$", p, ignore.case = TRUE)) paste0("deploy/apps/", p, ".yaml") else character(0),
    if (!grepl("\\.ya?ml$", p, ignore.case = TRUE)) paste0("deploy/apps/", p, ".yml") else character(0)
  ))
  ok <- candidates[vapply(candidates, file.exists, logical(1L))]
  if (length(ok) < 1L) "config.yaml" else ok[[1L]]
}

main_config_path <- resolve_main_config_path(main_config_path_env)
message(
  "[exprs-app] EXPRS_MAIN_CONFIG env=(", main_config_path_env %||% "", ") resolved=(", main_config_path, ")",
  " exists=", file.exists(main_config_path),
  " cwd=", getwd()
)
repo_root_app <- getwd()
main_config <- if (file.exists(main_config_path)) {
  mc <- read_main_yaml_merged(main_config_path, repo_root = repo_root_app)
  if (length(mc) < 1L) {
    list(title = "Gene expression dashboard", studies = character(0))
  } else {
    mc
  }
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

message("[exprs-app] main_config title=", main_config$title %||% "", " studies=", paste(study_ids, collapse = ","))

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
if (length(study_choices) == 0L) {
  study_choices <- c("(no studies)" = "")
} else {
  study_choices <- c("Select study..." = "", study_choices)
}

 

# Default study = none selected; default DEG list empty until study is chosen.
default_study <- ""
first_study_lists <- study_deg_lists(default_study)
default_deg_choices <- if (length(first_study_lists) > 0L) {
  setNames(vapply(first_study_lists, function(x) x$deg_file, character(1L)),
           vapply(first_study_lists, function(x) x$label, character(1L)))
} else {
  c("(no DEG lists)" = "")
}
default_deg <- if (length(first_study_lists) > 0L) first_study_lists[[1L]]$deg_file else ""

default_upset_study <- if (length(study_ids) > 0L) study_ids[[1L]] else ""
upset_lists0 <- study_deg_lists(default_upset_study)
default_upset_deg <- if (length(upset_lists0) > 0L) upset_lists0[[1L]]$deg_file else ""
default_upset_deg_choices <- if (length(upset_lists0) > 0L) {
  setNames(
    vapply(upset_lists0, function(x) x$deg_file, character(1L)),
    vapply(upset_lists0, function(x) x$label, character(1L))
  )
} else {
  c("(no DEG lists)" = "")
}
upset_study_choices <- if (length(study_ids) > 0L) {
  stats::setNames(study_ids, study_labels)
} else {
  c("(no studies)" = "")
}
n_upset_deg_tasks <- if (length(study_ids) > 0L) {
  length(ora_build_deg_tasks(study_ids, stats::setNames(study_labels, study_ids)))
} else {
  1L
}

pathway_txt_files <- list_ora_pathway_files(
  study_ids,
  main_config_path = main_config_path,
  main_cfg = main_config
)
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
    /* UpSet sidebar section bottom padding only (button widths use inline styles on the control) */
    .sidebar .upset-sidebar-section,
    .main-sidebar .upset-sidebar-section,
    .left-side .upset-sidebar-section {
      padding-bottom: 2rem;
    }
    .sidebar .upset-threshold-actions .btn,
    .main-sidebar .upset-threshold-actions .btn,
    .left-side .upset-threshold-actions .btn {
      width: auto !important;
      max-width: 100% !important;
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
      tabName = "info",
      fluidRow(
        box(
          title = "Study overview",
          width = 12, solidHeader = TRUE, status = "secondary",
          uiOutput("info_study_description"),
          tags$hr(),
          tags$div(
            style = "max-width: 980px; margin: 0 auto;",
            uiOutput("info_pca_ui")
          )
        )
      ),
      fluidRow(
        box(
          title = "Study metadata",
          width = 12, solidHeader = TRUE, status = "secondary",
          uiOutput("info_metadata_table_ui")
        )
      )
    ),
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
    ),
    tabItem(
      tabName = "upset",
      box(
        title = "DEG list intersections (UpSet)",
        width = 12, solidHeader = TRUE, status = "secondary",
        upsetTabUI("upset")
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
      navbarTab(tabName = "info", text = "Info"),
      navbarTab(tabName = "degs", text = "DEGs"),
      navbarTab(tabName = "gene", text = "Gene"),
      navbarTab(tabName = "ora", text = "ORA"),
      navbarTab(tabName = "upset", text = "UpSet")
    )
  ),
  sidebar = dashboardSidebar(
    minified = FALSE,
    shiny::conditionalPanel(
      condition = "input.navtabs == 'info' || input.navtabs == 'degs' || input.navtabs == 'gene'",
      selectInput("study", label = "Study", choices = study_choices, selected = default_study)
    ),
    shiny::conditionalPanel(
      condition = "input.navtabs == 'info'",
      uiOutput("info_color_by_ui")
    ),
    # DEGs tab only: heatmap / table controls.
    shiny::conditionalPanel(
      condition = "input.navtabs == 'degs'",
      selectInput("deg_list", label = "DEG list", choices = default_deg_choices, selected = default_deg),
      uiOutput("deg_description"),
      selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1, "0.5" = 0.5), selected = default_thr$fdr),
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
      shiny::tags$div(
        style = "font-size: 0.88rem; margin-bottom: 6px;",
        shiny::tags$strong("Custom ontology (on-the-fly ORA)")
      ),
      shiny::fileInput(
        "ora_custom_ontology_xlsx",
        label = shiny::tags$span("Excel: symbol + category", style = "font-weight: normal;"),
        buttonLabel = "Choose .xlsx…",
        accept = c(".xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
        width = "100%"
      ),
      shiny::tags$div(style = "margin-top: 2px;"),
      shiny::actionButton(
        "ora_load_custom_ontology",
        "Load custom ontology",
        class = "btn-sm btn-primary"
      ),
      shiny::tags$div(style = "height: 4px;"),
      shiny::actionButton(
        "ora_clear_custom_ontology",
        "Clear custom",
        class = "btn-sm btn-default"
      ),
      shiny::tags$div(
        class = "text-muted",
        style = "font-size: 0.78rem; margin: 6px 0 10px 0;",
        "Long format: col 1 = symbol, col 2 = category (one gene per row; gene sets = all rows per category). ",
        "Unlike dropdown ",
        shiny::tags$code(".txt"),
        " files (one pathway per line). Gene symbols are uppercased on load so matching is case-insensitive. Requires ",
        shiny::tags$code("readxl"),
        "."
      ),
      shiny::tags$hr(),
      numericInput("ora_min_overlap", label = "Minimum pathway size (minGSSize):", value = 10L, min = 1L, step = 1L),
      numericInput("ora_min_count", label = "Minimum overlap (Count):", value = 5L, min = 1L, step = 1L),
      numericInput("ora_min_gene_ratio", label = "Minimum gene ratio:", value = 0.1, min = 0, max = 1, step = 0.01),
      selectInput(
        "ora_fdr_cutoff",
        label = "Maximum adjusted p-value (FDR):",
        choices = c("0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1, "0.5" = 0.5, "1" = 1),
        selected = 1
      ),
      selectInput(
        "ora_pvalue_cutoff",
        label = "Maximum p-value:",
        choices = c("0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1, "0.5" = 0.5, "1" = 1),
        selected = 1
      ),
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
      actionButton("ora_copy_visible_genes", label = "Copy gene symbols"),
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
        "Choose a gene in the main panel. Study selection also controls the Info tab."
      )
    ),
    shiny::conditionalPanel(
      condition = "input.navtabs == 'upset'",
      shiny::tags$div(
        class = "upset-sidebar-section",
        shiny::tags$div(
          class = "text-muted",
          style = "font-size: 0.82rem; margin-bottom: 8px;",
          "UpSet uses ",
          shiny::tags$strong("all"),
          " DEG lists from every study in the main config. Thresholds here apply ",
          shiny::tags$strong("only to the UpSet tab"),
          " (not the DEGs tab sliders). Pick a study and list to edit values, then use the buttons below."
        ),
        shiny::numericInput(
          "upset_max_degree",
          label = "Max lists per intersection (degree):",
          value = min(2L, max(2L, as.integer(n_upset_deg_tasks))),
          min = 2L,
          max = max(8L, as.integer(n_upset_deg_tasks)),
          step = 1L
        ),
        shiny::tags$div(
          class = "text-muted",
          style = "font-size: 0.78rem; margin: -6px 0 10px 0;",
          "Mode is always inclusive ",
          shiny::tags$strong("intersect"),
          " (not distinct). With max degree 2, every pair of DEG lists gets a column when the pair count is below the app cap; raise degree to include triples etc. (then only observed intersections are listed)."
        ),
        selectInput(
          "upset_min_intersection",
          label = "Minimum intersection size (genes):",
          choices = c(
            "No limit" = 0,
            "1" = 1,
            "2" = 2,
            "5" = 5,
            "10" = 10,
            "25" = 25,
            "50" = 50,
            "100" = 100,
            "200" = 200,
            "500" = 500,
            "1,000" = 1000
          ),
          selected = 0,
          width = "100%"
        ),
        selectInput("upset_study", label = "Study (threshold target)", choices = upset_study_choices, selected = default_upset_study),
        selectInput("upset_deg_list", label = "DEG list (threshold target)", choices = default_upset_deg_choices, selected = default_upset_deg),
        uiOutput("upset_deg_description"),
        selectInput("upset_fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1, "0.5" = 0.5), selected = default_thr$fdr),
        uiOutput("upset_svalue_sidebar"),
        uiOutput("upset_base_mean_sidebar"),
        numericInput("upset_log2fc", label = "Minimal abs(log2 fold change):", value = default_thr$log2fc),
        shiny::tags$div(
          class = "form-group upset-threshold-actions",
          style = paste0(
            "width:270px;max-width:100%;box-sizing:border-box;",
            "display:flex;flex-direction:column;gap:6px;",
            "margin-top:2rem;margin-bottom:1.5rem;margin-right:15px;"
          ),
          shiny::tags$div(
            style = "min-width:0;width:100%;max-width:100%;box-sizing:border-box;overflow:hidden;",
            shiny::actionButton(
              "upset_refresh",
              label = "Refresh thresholds for selected list",
              class = "btn-primary",
              style = "margin:6px 15px 6px 15px;box-sizing:border-box;white-space:normal;overflow-wrap:anywhere;word-break:break-word;"
            )
          ),
          shiny::tags$div(
            style = "min-width:0;width:100%;max-width:100%;box-sizing:border-box;overflow:hidden;",
            shiny::actionButton(
              "upset_thr_apply_all",
              label = "Apply current thresholds to all lists",
              class = "btn-secondary",
              style = "margin:6px 15px 6px 15px;box-sizing:border-box;white-space:normal;overflow-wrap:anywhere;word-break:break-word;"
            )
          ),
          shiny::tags$div(
            style = "min-width:0;width:100%;max-width:100%;box-sizing:border-box;overflow:hidden;",
            shiny::actionButton(
              "upset_thr_reset_all",
              label = "Reset all lists to default thresholds",
              class = "btn-outline-secondary",
              style = "margin:6px 15px 6px 15px;box-sizing:border-box;white-space:normal;overflow-wrap:anywhere;word-break:break-word;"
            )
          )
        )
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

  session$onFlushed(function() {
    if (requireNamespace("shinydashboard", quietly = TRUE)) {
      shinydashboard::updateTabItems(session, "navtabs", selected = "info")
    }
  }, once = TRUE)

  gene_jump_symbol <- shiny::reactiveVal(NULL)

  geneTabServer(
    "gene",
    study_ids,
    stats::setNames(study_labels, study_ids),
    external_symbol = shiny::reactive(gene_jump_symbol())
  )

  upset_refresh_n <- shiny::reactiveVal(0L)
  upset_thr_overrides <- shiny::reactiveValues()

  shiny::observeEvent(input$upset_refresh, {
    sid <- input$upset_study %||% ""
    deg <- input$upset_deg_list %||% ""
    if (!nzchar(sid) || !nzchar(deg)) {
      shiny::showNotification("Select a study and DEG list before refreshing.", type = "warning")
      return()
    }
    th <- upset_thr_sidebar()
    key <- paste(sid, deg, sep = "|")
    upset_thr_overrides[[key]] <- list(
      fdr = as.numeric(th$fdr),
      log2fc = as.numeric(th$log2fc),
      base_mean = as.numeric(th$base_mean %||% 0),
      svalue = as.numeric(th$svalue %||% default_heatmap_thresholds()$svalue)
    )
    upset_refresh_n(upset_refresh_n() + 1L)
    shiny::showNotification("Threshold override saved for this DEG list. Rebuilding UpSet.", type = "message")
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$upset_thr_apply_all, {
    tasks <- ora_build_deg_tasks(study_ids, stats::setNames(study_labels, study_ids), deg_filter = NULL)
    if (length(tasks) < 1L) {
      shiny::showNotification("No DEG lists found in config.", type = "warning")
      return()
    }
    th <- upset_thr_sidebar()
    d0 <- default_heatmap_thresholds()
    n <- 0L
    for (tk in tasks) {
      deg_rel <- tk$entry$deg_file
      if (is.null(deg_rel) || !nzchar(as.character(deg_rel))) next
      key <- paste(tk$sid, deg_rel, sep = "|")
      upset_thr_overrides[[key]] <- list(
        fdr = as.numeric(th$fdr %||% d0$fdr),
        log2fc = as.numeric(th$log2fc %||% d0$log2fc),
        base_mean = as.numeric(th$base_mean %||% d0$base_mean),
        svalue = as.numeric(th$svalue %||% d0$svalue)
      )
      n <- n + 1L
    }
    upset_refresh_n(upset_refresh_n() + 1L)
    shiny::showNotification(paste("Applied current sidebar thresholds to", n, "DEG list(s). UpSet only."), type = "message")
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$upset_thr_reset_all, {
    ov <- shiny::reactiveValuesToList(upset_thr_overrides)
    for (nm in names(ov)) {
      upset_thr_overrides[[nm]] <- NULL
    }
    upset_refresh_n(upset_refresh_n() + 1L)
    shiny::showNotification("Cleared all UpSet threshold overrides; lists use YAML defaults.", type = "message")
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$upset_study, {
    if (is.null(input$upset_study) || !nzchar(input$upset_study)) return()
    lists <- study_deg_lists(input$upset_study)
    if (length(lists) < 1L) {
      shiny::updateSelectInput(session, "upset_deg_list", choices = c("(no DEG lists)" = ""), selected = "")
      return()
    }
    ch <- stats::setNames(
      vapply(lists, function(x) x$deg_file, character(1L)),
      vapply(lists, function(x) x$label, character(1L))
    )
    shiny::updateSelectInput(session, "upset_deg_list", choices = ch, selected = lists[[1L]]$deg_file)
  }, ignoreNULL = FALSE)

  shiny::observeEvent(list(input$upset_study, input$upset_deg_list), {
    if (is.null(input$upset_study) || !nzchar(input$upset_study)) return()
    if (is.null(input$upset_deg_list) || !nzchar(input$upset_deg_list)) return()
    d <- deg_list_threshold_defaults(input$upset_study, input$upset_deg_list)
    fc <- c(0.001, 0.01, 0.05, 0.1, 0.5)
    fdr_sel <- if (d$fdr %in% fc) d$fdr else 0.05
    shiny::updateSelectInput(session, "upset_fdr", selected = fdr_sel)
    shiny::updateNumericInput(session, "upset_log2fc", value = d$log2fc)
    if (!is.null(input$upset_base_mean)) {
      shiny::updateNumericInput(session, "upset_base_mean", value = d$base_mean)
    }
    if (!is.null(input$upset_svalue)) {
      shiny::updateNumericInput(session, "upset_svalue", value = d$svalue)
    }
  }, ignoreInit = FALSE)

  output$upset_deg_description <- shiny::renderUI({
    if (is.null(input$upset_study) || !nzchar(input$upset_study)) return(NULL)
    if (is.null(input$upset_deg_list) || !nzchar(input$upset_deg_list)) return(NULL)
    lists <- study_deg_lists(input$upset_study)
    idx <- match(input$upset_deg_list, vapply(lists, function(x) x$deg_file, character(1L)))
    desc <- if (!is.na(idx) && !is.null(lists[[idx]]$description)) lists[[idx]]$description else ""
    if (is.null(desc) || !nzchar(trimws(as.character(desc)))) return(NULL)
    tags$div(
      class = "form-group shiny-input-container",
      tags$label("Description", class = "control-label"),
      tags$div(desc, class = "text-muted", style = "margin-top: 0.25rem; font-size: 0.9rem;")
    )
  })

  output$upset_svalue_sidebar <- shiny::renderUI({
    shiny::req(input$upset_study, input$upset_deg_list)
    shiny::req(nzchar(input$upset_study), nzchar(input$upset_deg_list))
    loaded <- load_study_data(input$upset_study, input$upset_deg_list)
    res <- loaded$res
    if (is.null(res) || !("svalue" %in% colnames(res))) return(NULL)
    d <- deg_list_threshold_defaults(input$upset_study, input$upset_deg_list)
    val <- if (!is.null(d$svalue)) d$svalue else default_heatmap_thresholds()$svalue
    numericInput("upset_svalue", label = "Cutoff for svalue:", value = val)
  })

  output$upset_base_mean_sidebar <- shiny::renderUI({
    shiny::req(input$upset_study, input$upset_deg_list)
    shiny::req(nzchar(input$upset_study), nzchar(input$upset_deg_list))
    loaded <- load_study_data(input$upset_study, input$upset_deg_list)
    res <- loaded$res
    if (is.null(res) || !("baseMean" %in% colnames(res))) return(NULL)
    d <- deg_list_threshold_defaults(input$upset_study, input$upset_deg_list)
    val <- if (!is.null(d$base_mean)) d$base_mean else default_heatmap_thresholds()$base_mean
    numericInput("upset_base_mean", label = "Minimal base mean:", value = val)
  })

  upset_tab_cfg <- shiny::reactive({
    md <- suppressWarnings(as.numeric(input$upset_max_degree %||% 2))
    if (length(md) != 1L || is.na(md) || md < 2L) md <- 2
    rs <- input[["upset-row_sort"]] %||% "study"
    if (!is.character(rs) || length(rs) != 1L || !nzchar(rs)) rs <- "study"
    if (identical(rs, "set_size")) rs <- "intersect_max"
    if (!rs %in% c("study", "intersect_max")) rs <- "study"
    mis <- suppressWarnings(as.numeric(input$upset_min_intersection %||% 0))
    if (length(mis) != 1L || is.na(mis) || mis < 0) mis <- 0
    mis <- min(mis, 1e7)
    list(
      max_degree = md,
      min_intersection = mis,
      refresh = upset_refresh_n(),
      row_sort = rs
    )
  })

  upset_thr_sidebar <- shiny::reactive({
    sid <- input$upset_study %||% ""
    deg <- input$upset_deg_list %||% ""
    d0 <- default_heatmap_thresholds()
    if (!nzchar(sid) || !nzchar(deg)) {
      return(list(
        fdr = d0$fdr,
        log2fc = d0$log2fc,
        base_mean = d0$base_mean,
        svalue = d0$svalue
      ))
    }
    d <- deg_list_threshold_defaults(sid, deg)
    list(
      fdr = as.numeric(input$upset_fdr %||% d$fdr),
      log2fc = as.numeric(input$upset_log2fc %||% d$log2fc),
      base_mean = as.numeric(input$upset_base_mean %||% d$base_mean),
      svalue = as.numeric(input$upset_svalue %||% d$svalue)
    )
  })

  upsetTabServer(
    "upset",
    study_ids,
    stats::setNames(study_labels, study_ids),
    upset_tab_cfg,
    upset_thr_overrides,
    on_dot_click = function(sid, deg, genes) {
      if (!requireNamespace("shinydashboard", quietly = TRUE)) {
        return(invisible(NULL))
      }
      sid <- as.character(sid %||% "")[[1L]]
      deg <- as.character(deg %||% "")[[1L]]
      if (!nzchar(sid) || !nzchar(deg)) return(invisible(NULL))
      genes <- unique(trimws(as.character(genes)))
      genes <- genes[nzchar(genes) & !is.na(genes)]

      lists <- study_deg_lists(sid)
      if (length(lists) < 1L) return(invisible(NULL))
      choices <- stats::setNames(
        vapply(lists, function(x) x$deg_file, character(1L)),
        vapply(lists, function(x) x$label, character(1L))
      )
      deg_vals <- unname(choices)
      if (!deg %in% deg_vals) return(invisible(NULL))

      # pending_deg_list is only consumed in observeEvent(input$study). If the DEGs tab is
      # already on this study, updateSelectInput(study) is a no-op and the list never switched —
      # update deg_list here. Cross-study navigation keeps using pending + study observer.
      cur_study <- as.character(shiny::isolate(input$study %||% ""))[[1L]]
      same_study <- nzchar(cur_study) && identical(cur_study, sid)

      rv$apply_custom_after_load <- length(genes) > 0L
      rv$custom_genes <- if (length(genes) > 0L) genes else NULL
      rv$custom_genes_study <- sid

      shinydashboard::updateTabItems(session, "navtabs", selected = "degs")
      shiny::updateCheckboxInput(session, "lock_gene_list", value = TRUE)

      if (same_study) {
        rv$pending_deg_list <- NULL
        shiny::updateSelectInput(session, "deg_list", choices = choices, selected = deg)
      } else {
        rv$pending_deg_list <- deg
        shiny::updateSelectInput(session, "study", selected = sid)
      }

      # If study+deg were already correct, neither observer re-runs and apply_custom_after_load
      # would never be consumed — flush once and apply when data are ready.
      if (isTRUE(rv$apply_custom_after_load)) {
        session$onFlushed(function() {
          if (!isTRUE(shiny::isolate(rv$apply_custom_after_load))) return()
          if (!identical(as.character(shiny::isolate(input$study) %||% ""), sid)) return()
          if (!identical(as.character(shiny::isolate(input$deg_list) %||% ""), deg)) return()
          if (is.null(shiny::isolate(rv$current_res)) || is.null(shiny::isolate(rv$current_mm))) {
            return()
          }
          rv$apply_custom_after_load <- FALSE
          apply_custom_gene_list(apply_thresholds = TRUE)
        }, once = TRUE)
      }

      invisible(NULL)
    }
  )

  ora_custom_ontology <- shiny::reactiveVal(list(active = FALSE, t2g = NULL))

  shiny::observeEvent(input$ora_load_custom_ontology, {
    shiny::req(input$ora_custom_ontology_xlsx)
    pr <- parse_ontology_xlsx_for_ora(input$ora_custom_ontology_xlsx$datapath)
    if (!isTRUE(pr$ok)) {
      shiny::showNotification(
        if (!is.null(pr$error) && nzchar(as.character(pr$error))) pr$error else "Failed to read custom ontology.",
        type = "error"
      )
      return()
    }
    ora_custom_ontology(list(active = TRUE, t2g = pr$t2g))
    shiny::updateSelectInput(session, "ora-pathway_file", choices = c(" " = ""), selected = "")
    n_cat <- length(unique(pr$t2g$term))
    message(
      "[ORA debug] custom ontology loaded from ",
      as.character(input$ora_custom_ontology_xlsx$name),
      " categories=", n_cat,
      " pairs=", nrow(pr$t2g),
      " unique_genes=", length(unique(as.character(pr$t2g$gene))),
      if (isTRUE(pr$swapped)) " swapped_columns=TRUE" else " swapped_columns=FALSE"
    )
    msg <- paste0("Custom ontology loaded: ", n_cat, " categories, ", nrow(pr$t2g), " gene–category pairs.")
    if (isTRUE(pr$swapped)) {
      msg <- paste0(msg, " (Detected unnamed columns looked swapped; used column 2 as symbols, column 1 as categories.)")
    }
    shiny::showNotification(msg, type = "message", duration = 10)
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$ora_clear_custom_ontology, {
    if (is.null(input$ora_clear_custom_ontology) || input$ora_clear_custom_ontology < 1) return()
    message("[ORA debug] clear custom ontology clicked: n=", input$ora_clear_custom_ontology)
    ora_custom_ontology(list(active = FALSE, t2g = NULL))
    ch <- ora_file_choices
    if (!("" %in% unname(ch))) {
      ch <- c("Select a pathway database..." = "", ch)
    }
    shiny::updateSelectInput(session, "ora-pathway_file", choices = ch, selected = "")
    shiny::showNotification("Custom ontology cleared.", type = "message")
  }, ignoreInit = TRUE)

  rv <- reactiveValues(
    current_res = NULL,
    current_mm = NULL,
    current_col_annot = NULL,
    row_index = NULL,
    res_table_row_index = integer(0),
    selected_rows = NULL,
    pending_deg_list = NULL,
    apply_custom_after_load = FALSE,
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

  info_panel_set <- shiny::reactive({
    if (is.null(input$study) || input$study == "") return(NULL)
    load_gene_tab_study_panels(input$study)
  })

  info_data <- shiny::reactive({
    ps <- info_panel_set()
    if (is.null(ps) || length(ps$panels) < 1L) return(NULL)
    p1 <- ps$panels[[1L]]
    list(mm = p1$mm, metadata = p1$metadata, is_count_like = p1$is_count_like)
  })

  output$info_color_by_ui <- renderUI({
    ps <- info_panel_set()
    if (is.null(ps) || length(ps$panels) < 1L) {
      return(
        tags$div(
          class = "text-muted",
          style = "padding: 6px 0;",
          "No metadata available for selected study."
        )
      )
    }
    dat <- list(metadata = ps$panels[[1L]]$metadata)
    if (is.null(dat$metadata) || nrow(dat$metadata) < 1L) {
      return(
        tags$div(
          class = "text-muted",
          style = "padding: 6px 0;",
          "No metadata available for selected study."
        )
      )
    }
    meta <- dat$metadata
    col_choices <- colnames(meta)
    default_col <- if ("PhenoNames" %in% col_choices) "PhenoNames" else col_choices[[1L]]
    selectInput(
      "info_color_by",
      "Color PCA by metadata column",
      choices = col_choices,
      selected = default_col
    )
  })

  output$info_study_description <- renderUI({
    if (is.null(input$study) || input$study == "") return(NULL)
    cfg_path <- file.path("data", input$study, "config.yaml")
    if (!file.exists(cfg_path)) return(tags$div(class = "text-muted", "No study config found."))
    cfg <- yaml::read_yaml(cfg_path)
    desc <- cfg$description
    if (is.null(desc) || !nzchar(trimws(as.character(desc)))) {
      return(tags$div(class = "text-muted", "No study description provided in config."))
    }
    tags$div(
      style = "font-size: 0.95rem; line-height: 1.5;",
      as.character(desc)
    )
  })

  output$info_pca_ui <- renderUI({
    ps <- info_panel_set()
    if (is.null(ps) || length(ps$panels) < 1L) {
      return(tags$div(class = "text-muted", "No data available for PCA."))
    }
    use_girafe <- requireNamespace("ggiraph", quietly = TRUE)
    out <- list()
    for (i in seq_along(ps$panels)) {
      panel <- ps$panels[[i]]
      p_title <- if (is.null(panel$title) || !nzchar(as.character(panel$title))) {
        paste0("Panel ", i)
      } else {
        as.character(panel$title)
      }
      out[[length(out) + 1L]] <- tags$div(
        style = "font-weight: 600; margin: 6px 0 4px 0;",
        p_title
      )
      out_id <- paste0("info_pca_", i)
      out[[length(out) + 1L]] <- if (use_girafe) {
        ggiraph::girafeOutput(out_id, width = "100%", height = "460px")
      } else {
        plotOutput(out_id, height = 460, width = "100%")
      }
    }
    tagList(out)
  })

  info_pca_data_one <- function(mm, meta, color_by, study_label) {
    if (ncol(mm) < 2L) {
      return(list(ok = FALSE, msg = "At least 2 samples are required for PCA."))
    }
    if (!color_by %in% colnames(meta)) {
      return(list(ok = FALSE, msg = "Selected metadata column is not available."))
    }
    mm_num <- suppressWarnings(apply(mm, 2, as.numeric))
    rownames(mm_num) <- rownames(mm)
    finite_vals <- as.vector(mm_num)
    finite_vals <- finite_vals[is.finite(finite_vals)]
    if (length(finite_vals) < 1L) {
      return(list(ok = FALSE, msg = "No numeric values available for PCA."))
    }
    # Heuristic: integer-like matrix is treated as raw counts and log-transformed.
    is_count_like <- all(abs(finite_vals - round(finite_vals)) < 1e-8)
    mm_for_pca <- if (is_count_like) log2(mm_num + 0.5) else mm_num
    # Remove genes with non-finite or zero variance; prcomp(scale.=TRUE) errors on constants.
    keep_rows <- apply(mm_for_pca, 1, function(x) {
      x <- as.numeric(x)
      x <- x[is.finite(x)]
      length(x) >= 2L && stats::sd(x) > 0
    })
    mm_for_pca <- mm_for_pca[keep_rows, , drop = FALSE]
    if (nrow(mm_for_pca) < 2L) {
      return(list(ok = FALSE, msg = "Not enough variable genes for PCA after filtering constant rows."))
    }
    pca <- stats::prcomp(t(mm_for_pca), center = TRUE, scale. = TRUE)
    var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
    pca_df <- data.frame(
      SampleNumber = rownames(pca$x),
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      ColorBy = as.character(meta[[color_by]]),
      stringsAsFactors = FALSE
    )
    pca_df$ColorBy[is.na(pca_df$ColorBy) | !nzchar(pca_df$ColorBy)] <- "(missing)"
    list(
      ok = TRUE,
      pca_df = pca_df,
      study_label = study_label,
      var_expl = var_expl
    )
  }

  info_pca_base_plot <- function(pd, color_by) {
    if (!isTRUE(pd$ok)) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0, y = 0, label = pd$msg, size = 4.2) +
          ggplot2::theme_void() +
          ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))
      )
    }
    ggplot2::ggplot(pd$pca_df, ggplot2::aes(x = PC1, y = PC2, color = ColorBy)) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::labs(
        title = paste0("PCA - ", pd$study_label),
        x = sprintf("PC1 (%.1f%%)", 100 * pd$var_expl[1]),
        y = sprintf("PC2 (%.1f%%)", 100 * pd$var_expl[2]),
        color = color_by
      ) +
      ggplot2::theme_minimal(base_size = 13)
  }

  observe({
    ps <- info_panel_set()
    shiny::req(ps, ps$panels, input$info_color_by)
    use_girafe <- requireNamespace("ggiraph", quietly = TRUE)
    sid_lbl <- study_labels[match(input$study, study_ids)]
    if (is.na(sid_lbl) || !nzchar(sid_lbl)) sid_lbl <- input$study
    for (i in seq_along(ps$panels)) {
      local({
        idx <- i
        panel <- ps$panels[[idx]]
        panel_title <- if (is.null(panel$title) || !nzchar(as.character(panel$title))) {
          sid_lbl
        } else {
          paste0(sid_lbl, " — ", as.character(panel$title))
        }
        pd <- info_pca_data_one(panel$mm, panel$metadata, input$info_color_by, panel_title)
        out_id <- paste0("info_pca_", idx)
        if (use_girafe) {
          output[[out_id]] <- ggiraph::renderGirafe({
            if (!isTRUE(pd$ok)) {
              return(ggiraph::girafe(ggobj = info_pca_base_plot(pd, input$info_color_by)))
            }
            p <- ggplot2::ggplot(pd$pca_df, ggplot2::aes(x = PC1, y = PC2, color = ColorBy)) +
              ggiraph::geom_point_interactive(
                ggplot2::aes(
                  tooltip = paste0("Sample: ", SampleNumber, "\n", input$info_color_by, ": ", ColorBy),
                  data_id = SampleNumber
                ),
                size = 3, alpha = 0.85
              ) +
              ggplot2::labs(
                title = paste0("PCA - ", pd$study_label),
                x = sprintf("PC1 (%.1f%%)", 100 * pd$var_expl[1]),
                y = sprintf("PC2 (%.1f%%)", 100 * pd$var_expl[2]),
                color = input$info_color_by
              ) +
              ggplot2::theme_minimal(base_size = 13)
            ggiraph::girafe(
              ggobj = p,
              width_svg = 10,
              height_svg = 4.8
            )
          })
        } else {
          output[[out_id]] <- renderPlot({
            info_pca_base_plot(pd, input$info_color_by)
          })
        }
      })
    }
  })

  output$info_metadata_table_ui <- renderUI({
    dat <- info_data()
    if (is.null(dat) || is.null(dat$metadata)) {
      return(tags$div(class = "text-muted", "Select a study to view metadata."))
    }
    n_cols <- ncol(dat$metadata)
    if (isTRUE(n_cols > 8L)) {
      return(
        tags$div(
          style = "overflow-x: auto;",
          DTOutput("info_metadata_table")
        )
      )
    }
    tags$div(
      style = "max-width: 980px; margin: 0 auto;",
      DTOutput("info_metadata_table")
    )
  })

  output$info_metadata_table <- renderDT({
    dat <- info_data()
    req(dat, dat$metadata)
    meta <- dat$metadata
    is_wide <- isTRUE(ncol(meta) > 8L)
    DT::datatable(
      meta,
      rownames = FALSE,
      options = list(
        pageLength = 20,
        lengthMenu = list(c(20, 50, 100), c("20", "50", "100")),
        scrollX = is_wide,
        autoWidth = TRUE
      )
    )
  })

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
    fdr_cut <- if (is.null(input$ora_fdr_cutoff)) 1 else as.numeric(input$ora_fdr_cutoff)
    if (is.na(fdr_cut) || fdr_cut <= 0) fdr_cut <- 1
    if (fdr_cut > 1) fdr_cut <- 1
    p_cut <- if (is.null(input$ora_pvalue_cutoff)) 1 else as.numeric(input$ora_pvalue_cutoff)
    if (is.na(p_cut) || p_cut <= 0) p_cut <- 1
    if (p_cut > 1) p_cut <- 1
    sc <- if (is.null(input$ora_show_category)) 20L else input$ora_show_category
    click_target <- if (is.null(input$ora_heatmap_click_target)) "pathway" else as.character(input$ora_heatmap_click_target)
    if (!click_target %in% c("pathway", "comparison")) click_target <- "pathway"
    hide_empty <- if (is.null(input$ora_pathway_hide_empty_comparisons)) TRUE else isTRUE(input$ora_pathway_hide_empty_comparisons)
    list(
      min_overlap = mo,
      min_count = mc,
      min_gene_ratio = mgr,
      max_p_value = p_cut,
      max_p_adj = fdr_cut,
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
    copy_genes_trigger = shiny::reactive(input$ora_copy_visible_genes),
    custom_ontology = shiny::reactive(ora_custom_ontology()),
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
    gene_id_col <- if ("ens_gene" %in% colnames(res_sub)) "ens_gene" else "symbol"
    mm_sub <- mm[as.character(res_sub[[gene_id_col]]), , drop = FALSE]

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
    rv$res_table_row_index <- rows
    output[["res_table"]] <- renderDT(
      formatRound(
        datatable(
          res[rows, tbl_cols, drop = FALSE],
          rownames = FALSE,
          selection = list(mode = "single", target = "row"),
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

  observeEvent(input$res_table_rows_selected, {
    sel <- input$res_table_rows_selected
    if (is.null(sel) || length(sel) < 1L) return()
    selected_row <- as.integer(sel[[1L]])
    idx_map <- rv$res_table_row_index
    res <- rv$current_res
    if (is.null(res) || length(idx_map) < selected_row) return()
    row_idx <- idx_map[[selected_row]]
    if (is.na(row_idx) || row_idx < 1L || row_idx > nrow(res)) return()
    if (!("symbol" %in% colnames(res))) return()
    sym <- trimws(as.character(res$symbol[[row_idx]]))
    if (!nzchar(sym)) return()
    gene_jump_symbol(NULL)
    gene_jump_symbol(sym)
    if (requireNamespace("shinydashboard", quietly = TRUE)) {
      shinydashboard::updateTabItems(session, "navtabs", selected = "gene")
    }
  }, ignoreInit = TRUE)

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
    pending <- rv$pending_deg_list
    rv$pending_deg_list <- NULL
    deg_vals <- unname(choices)
    sel <- if (!is.null(pending) && nzchar(as.character(pending)) && as.character(pending) %in% deg_vals) {
      as.character(pending)
    } else {
      lists[[1L]]$deg_file
    }
    updateSelectInput(session, "deg_list", choices = choices, selected = sel)
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
    if (is.null(input$study) || input$study == "" || is.null(input$deg_list) || input$deg_list == "") return()
    lists_chk <- study_deg_lists(input$study)
    if (length(lists_chk) < 1L) return()
    deg_ok <- vapply(lists_chk, function(x) x$deg_file, character(1L))
    if (!input$deg_list %in% deg_ok) return()
    rv$current_res <- NULL
    rv$current_mm <- NULL
    rv$current_col_annot <- NULL
    rv$row_index <- NULL
    rv$selected_rows <- NULL
    message("[perf] load_study_data: start; study=", input$study, " deg_list=", input$deg_list)
    loaded <- perf_time("load_study_data", load_study_data(input$study, input$deg_list))
    rv$current_res <- loaded$res
    rv$current_mm <- loaded$mm
    rv$current_col_annot <- loaded$col_annot
    d <- deg_list_threshold_defaults(input$study, input$deg_list)
    rv$threshold_defaults <- d
    fdr_choices <- c(0.001, 0.01, 0.05, 0.1, 0.5)
    fdr_sel <- d$fdr
    if (!fdr_sel %in% fdr_choices) {
      message(
        "[config] fdr=", fdr_sel, " not in UI choices ",
        paste(fdr_choices, collapse = ", "),
        "; using ", default_heatmap_thresholds()$fdr
      )
      fdr_sel <- default_heatmap_thresholds()$fdr
    }
    updateSelectInput(
      session,
      "fdr",
      choices = c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05, "0.1" = 0.1, "0.5" = 0.5),
      selected = fdr_sel
    )
    updateNumericInput(session, "log2fc", value = d$log2fc)
    output$ht_heatmap <- renderPlot({
      grid.newpage()
      grid.text("Click \"Generate heatmap\" to display the heatmap.")
    })
    if (isTRUE(rv$apply_custom_after_load)) {
      rv$apply_custom_after_load <- FALSE
      shiny::isolate(apply_custom_gene_list(apply_thresholds = TRUE))
    }
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