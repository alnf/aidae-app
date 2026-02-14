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

# Make the heatmap for differentially expressed genes under certain cutoffs.
# Returns list(ht = ..., row_index = which(l)) or NULL; res and mm must be set.
make_heatmap <- function(res, mm, fdr = 0.01, base_mean = 0, log2fc = 1) {
  if (is.null(res) || is.null(mm) || nrow(res) == 0L) return(NULL)
  l <- res$padj <= fdr & res$baseMean >= base_mean & abs(res$log2FoldChange) >= log2fc
  l[is.na(l)] <- FALSE
  if (sum(l) == 0L) return(NULL)
  m <- mm[l, , drop = FALSE]
  row_index <- which(l)
  ht <- Heatmap(t(scale(t(m))), name = "z-score",
      show_row_names = FALSE, show_column_names = FALSE, row_km = 2,
      column_title = paste0(sum(l), " significant genes with FDR < ", fdr),
      show_row_dend = FALSE) +
      Heatmap(log10(res$baseMean[l] + 1), show_row_names = FALSE, width = unit(5, "mm"),
          name = "log10(baseMean+1)", show_column_names = FALSE) +
      Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
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

# Side bar: study selector then cutoffs for significant genes.
# selected must be the choice value (study id), not the label
default_study <- if (length(study_ids) == 1L) study_ids[1L] else ""
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
    selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05)),
    numericInput("base_mean", label = "Minimal base mean:", value = 0),
    numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 1),
    actionButton("filter", label = "Generate heatmap")
  ),
  controlbar = dashboardControlbar(disable = TRUE),
  body = body
))

# Load study data (res, mm) from data/<study_id>/ using study config.
load_study_data <- function(study_id) {
  if (is.null(study_id) || study_id == "") return(list(res = NULL, mm = NULL))
  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) return(list(res = NULL, mm = NULL))
  cfg <- yaml::read_yaml(cfg_path)
  deg_path <- file.path("data", study_id, cfg$deg_file)
  counts_path <- file.path("data", study_id, cfg$counts_file)
  if (!file.exists(deg_path) || !file.exists(counts_path)) return(list(res = NULL, mm = NULL))
  res <- read.table(deg_path, sep = "\t", header = TRUE, check.names = FALSE)
  mm <- read.table(counts_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  required <- c("padj", "baseMean", "log2FoldChange", "symbol")
  if (!all(required %in% colnames(res)) || nrow(res) != nrow(mm)) return(list(res = NULL, mm = NULL))
  list(res = res, mm = mm)
}

server <- function(input, output, session) {
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  rv <- reactiveValues(current_res = NULL, current_mm = NULL, row_index = NULL)

  # Load data when study selection changes
  observeEvent(input$study, {
    rv$current_res <- NULL
    rv$current_mm <- NULL
    rv$row_index <- NULL
    if (is.null(input$study) || input$study == "") return()
    loaded <- load_study_data(input$study)
    rv$current_res <- loaded$res
    rv$current_mm <- loaded$mm
  }, ignoreNULL = FALSE)

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

  # Regenerate heatmap when filter is clicked or study changes
  observeEvent(list(input$filter, input$study), {
    if (is.null(rv$current_res) || is.null(rv$current_mm)) {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("Select a study and load data.")
      })
      return()
    }
    out <- make_heatmap(rv$current_res, rv$current_mm,
      fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc)
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