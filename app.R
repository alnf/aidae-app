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

# Return list of {label, deg_file} for a study (supports deg_lists or legacy deg_file).
study_deg_lists <- function(study_id) {
  if (is.null(study_id) || study_id == "") return(list())
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(list())
  cfg <- yaml::read_yaml(cfg_path)
  if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    return(cfg$deg_lists)
  }
  if (!is.null(cfg$deg_file)) {
    return(list(list(label = if (!is.null(cfg$name)) cfg$name else study_id, deg_file = cfg$deg_file)))
  }
  list()
}

# Return SampleNumber vector for samples in the comparison identified by label (joint in
# comparison_file), or NULL if no filtering (missing config/files, or no matching row).
samples_for_comparison <- function(study_id, label) {
  message("samples_for_comparison entered: study_id=", study_id, " label=", label)
  if (is.null(study_id) || study_id == "" || is.null(label) || label == "") {
    message("samples_for_comparison: return NULL (null/empty study_id or label)")
    return(NULL)
  }
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    message("samples_for_comparison: return NULL (config file not found)")
    return(NULL)
  }
  cfg <- yaml::read_yaml(cfg_path)
  if (is.null(cfg$comparison_file) || is.null(cfg$metadata_file)) {
    message("samples_for_comparison: return NULL (comparison_file or metadata_file missing in config)")
    return(NULL)
  }
  comp_path <- file.path("data", study_id, cfg$comparison_file)
  meta_path <- file.path("data", study_id, cfg$metadata_file)
  if (!file.exists(comp_path) || !file.exists(meta_path)) {
    message("samples_for_comparison: return NULL (file not found: comp=", file.exists(comp_path), " meta=", file.exists(meta_path), " meta_path=", meta_path, ")")
    return(NULL)
  }
  comp <- read.table(comp_path, sep = "\t", header = TRUE, check.names = FALSE)
  if (!all(c("joint", "group1", "group2") %in% colnames(comp))) {
    message("samples_for_comparison: return NULL (comparison file missing joint/group1/group2)")
    return(NULL)
  }
  idx <- which(comp$joint == label)
  if (length(idx) == 0L) {
    message("samples_for_comparison: return NULL (no row with joint == label)")
    return(NULL)
  }
  group1 <- comp$group1[idx[1L]]
  group2 <- comp$group2[idx[1L]]
  meta <- read.table(meta_path, sep = "\t", header = TRUE, check.names = FALSE)
  if (!all(c("PhenoNames", "SampleNumber") %in% colnames(meta))) {
    message("samples_for_comparison: return NULL (metadata missing PhenoNames or SampleNumber)")
    return(NULL)
  }
  keep <- meta$PhenoNames %in% c(group1, group2)
  # Optional extra filter(s) from this DEG list: metadata_filter is column -> value(s).
  deg_entry <- if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    idx <- match(label, vapply(cfg$deg_lists, function(x) x$label, character(1L)))
    if (!is.na(idx)) cfg$deg_lists[[idx]] else NULL
  } else if (!is.null(cfg$deg_file) && (if (!is.null(cfg$name)) cfg$name else study_id) == label) {
    list(metadata_filter = NULL)
  } else NULL
  mf <- if (!is.null(deg_entry) && !is.null(deg_entry$metadata_filter)) deg_entry$metadata_filter else NULL
  if (!is.null(mf) && is.list(mf) && length(mf) > 0L) {
    for (col in names(mf)) {
      if (col %in% colnames(meta)) {
        vals <- as.character(unlist(mf[[col]]))
        keep <- keep & (meta[[col]] %in% vals)
      }
    }
  }
  out <- meta$SampleNumber[keep]
  message("Samples for heatmap (SampleNumber): ", paste(out, collapse = ", "))
  out
}

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

# Make the heatmap for differentially expressed genes under certain cutoffs.
# Returns list(ht = ..., row_index = which(l)) or NULL; res and mm must be set.
make_heatmap <- function(res, mm, fdr = 0.01, base_mean = 0, log2fc = 1, svalue = 0.005) {
  if (is.null(res) || is.null(mm) || nrow(res) == 0L) return(NULL)
  mm <- mm[res$ens_gene, , drop = FALSE]
  l <- res$padj <= fdr & res$baseMean >= base_mean & abs(res$log2FoldChange) >= log2fc
  if ("svalue" %in% colnames(res)) l <- l & res$svalue <= svalue
  l[is.na(l)] <- FALSE
  if (sum(l) == 0L) return(NULL)
  m <- mm[    l, , drop = FALSE]
  row_index <- which(l)
  # Z-score by row; drop rows with any NA/NaN/Inf (e.g. constant rows) so kmeans gets valid data
  m_z <- t(scale(t(m)))
  keep_row <- rowSums(!is.finite(m_z)) == 0L
  m_z <- m_z[keep_row, , drop = FALSE]
  row_index <- row_index[keep_row]
  if (nrow(m_z) == 0L) return(NULL)
  n_row <- nrow(m_z)
  n_col <- ncol(m_z)
  message("Heatmap dimensions: n_row=", n_row, " n_col=", n_col)
  # row_km can fail with few rows/columns (e.g. kmeans); only cluster when enough data
  row_km_arg <- if (n_row >= 3L && n_col >= 2L) 2L else NULL
  ht <- Heatmap(m_z, name = "z-score",
      show_row_names = FALSE, show_column_names = FALSE,
      row_km = row_km_arg, show_row_dend = !is.null(row_km_arg),
      column_title = paste0(n_row, " significant genes with FDR < ", fdr)) +
      Heatmap(log10(res$baseMean[row_index] + 1), show_row_names = FALSE, width = unit(5, "mm"),
          name = "log10(baseMean+1)", show_column_names = FALSE) +
      Heatmap(res$log2FoldChange[row_index], show_row_names = FALSE, width = unit(5, "mm"),
          name = "log2FoldChange", show_column_names = FALSE,
          col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  ht <- draw(ht, merge_legend = TRUE)
  list(ht = ht, row_index = row_index)
}

# make the MA-plot with some genes highlighted
make_maplot <- function(res, highlight = NULL) {
  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if(!is.null(highlight)) {
    col[highlight] <- "red"
    cex[highlight] <- 1
  }
  x <- res$baseMean
  y <- res$log2FoldChange
  y[y > 2] <- 2
  y[y < -2] <- -2
  col[col == "red" & y < 0] <- "darkgreen"
  par(mar = c(4, 4, 1, 1))

  suppressWarnings(
    plot(x, y, col = col, 
      pch = ifelse(res$log2FoldChange > 2 | res$log2FoldChange < -2, 1, 16), 
      cex = cex, log = "x",
      xlab = "baseMean", ylab = "log2 fold change")
  )
}

# make the volcano plot with some genes highlited
make_volcano <- function(res, highlight = NULL) {
  col <- rep("#00000020", nrow(res))
  cex <- rep(0.5, nrow(res))
  names(col) <- rownames(res)
  names(cex) <- rownames(res)
  if(!is.null(highlight)) {
    col[highlight] <- "red"
    cex[highlight] <- 1
  }
  x <- res$log2FoldChange
  y <- -log10(res$padj)
  col[col == "red" & x < 0] <- "darkgreen"
  par(mar = c(4, 4, 1, 1))

  suppressWarnings(
    plot(x, y, col = col, 
      pch = 16, 
      cex = cex,
      xlab = "log2 fold change", ylab = "-log10(FDR)")
  )
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
          id = "column2",
          box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
            subHeatmapOutput("ht", title = NULL, containment = TRUE)
          ),
          box(title = "Result table of the selected genes", width = NULL, solidHeader = TRUE, status = "primary",
            DTOutput("res_table")
          )
        ),
        column(width = 4,
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
    selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05), selected = 0.05),
    numericInput("svalue", label = "Cutoff for svalue:", value = 0.005),
    numericInput("base_mean", label = "Minimal base mean:", value = 20),
    numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 1),
    actionButton("filter", label = "Generate heatmap")
  ),
  controlbar = dashboardControlbar(disable = TRUE),
  body = body
))

# Load study data (res, mm) from data/<study_id>/ using study config and selected deg_file.
load_study_data <- function(study_id, deg_file) {
  if (is.null(study_id) || study_id == "" || is.null(deg_file) || deg_file == "") return(list(res = NULL, mm = NULL))
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(list(res = NULL, mm = NULL))
  cfg <- yaml::read_yaml(cfg_path)
  deg_path <- file.path("data", study_id, deg_file)
  counts_path <- file.path("data", study_id, cfg$counts_file)
  if (!file.exists(deg_path) || !file.exists(counts_path)) return(list(res = NULL, mm = NULL))
  res <- read.table(deg_path, sep = "\t", header = TRUE, check.names = FALSE)
  mm <- read.table(counts_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  # Resolve label for selected DEG list and optionally filter mm to comparison samples.
  lists <- if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) cfg$deg_lists else
    if (!is.null(cfg$deg_file)) list(list(label = if (!is.null(cfg$name)) cfg$name else study_id, deg_file = cfg$deg_file)) else list()
  idx <- match(deg_file, vapply(lists, function(x) x$deg_file, character(1L)))
  label <- if (!is.na(idx) && !is.null(lists[[idx]]$label)) lists[[idx]]$label else NULL
  print("TEST")
  sample_names <- if (!is.null(label)) samples_for_comparison(study_id, label) else NULL
  if (length(sample_names) > 0L) {
    keep <- intersect(sample_names, colnames(mm))
    if (length(keep) > 0L) mm <- mm[, keep, drop = FALSE]
  }
  required <- c("padj", "baseMean", "log2FoldChange", "symbol")
  if (!all(required %in% colnames(res)) || nrow(res) != nrow(mm)) return(list(res = NULL, mm = NULL))
  list(res = res, mm = mm)
}

server <- function(input, output, session) {
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  rv <- reactiveValues(current_res = NULL, current_mm = NULL, row_index = NULL)

  # When study changes, update DEG list dropdown to that study's lists (first selected).
  observeEvent(input$study, {
    rv$current_res <- NULL
    rv$current_mm <- NULL
    rv$row_index <- NULL
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
    rv$row_index <- NULL
    if (is.null(input$study) || input$study == "" || is.null(input$deg_list) || input$deg_list == "") return()
    loaded <- load_study_data(input$study, input$deg_list)
    rv$current_res <- loaded$res
    rv$current_mm <- loaded$mm
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  # Brush action uses current study data from rv
  brush_action <- function(df, input, output, session) {
    res <- rv$current_res
    row_index <- rv$row_index
    if (is.null(res) || is.null(row_index)) return()
    idx <- unique(unlist(df$row_index))
    selected <- row_index[idx]
    output[["ma_plot"]] <- renderPlot({
      make_maplot(res, selected)
    })
    output[["volcano_plot"]] <- renderPlot({
      make_volcano(res, selected)
    })
    output[["res_table"]] <- renderDT(
      formatRound(datatable(res[selected, c("symbol", "baseMean", "log2FoldChange", "padj")], rownames = FALSE), columns = 2:4, digits = 3)
    )
  }

  # Regenerate heatmap when filter is clicked or study/DEG list changes
  observeEvent(list(input$filter, input$study, input$deg_list), {
    if (is.null(rv$current_res) || is.null(rv$current_mm)) {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("Select a study and load data.")
      })
      return()
    }
    out <- make_heatmap(rv$current_res, rv$current_mm,
      fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc, svalue = input$svalue)
    if (!is.null(out)) {
      rv$row_index <- out$row_index
      makeInteractiveComplexHeatmap(input, output, session, out$ht, "ht",
        brush_action = brush_action)
    } else {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui, server)