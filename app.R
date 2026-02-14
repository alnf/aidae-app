# app.R - Dashboard for DESeq2 results visualization
# Licensed under GPL-3.0. See LICENSE in the repository.
#
# Inspired by the Shiny app for visualizing DESeq2 results by Zuguang Gu
# (InteractiveComplexHeatmap vignette, MIT License):
# https://github.com/jokergoo/InteractiveComplexHeatmap/blob/master/vignettes/deseq2_app.Rmd

# First we perform DESeq2 analysis on the airway dataset.
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
#remotes::install_github("datastorm-open/shinymanager", ref = "de42a23")


res <- read.table("data/res_annot_D0.tsv", sep="\t", header = T, check.names=F)
mm <- read.table("data/res_annot_D0_counts.tsv", sep="\t", header = T, row.names = 1, check.names=F)
print(colnames(mm))
print(colnames(res))

library(ComplexHeatmap)
library(circlize)

env <- new.env()

# Make the heatmap for differentially expressed genes under certain cutoffs.
make_heatmap <- function(fdr = 0.01, base_mean = 0, log2fc = 1) {
	l = res$padj <= fdr & res$baseMean >= base_mean & abs(res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE

	if(sum(l) == 0) return(NULL)

	m = mm[l, ]

  env$row_index <- which(l)

  ht <- Heatmap(t(scale(t(m))), name = "z-score",
      #top_annotation = HeatmapAnnotation(
      #    dex = colData(dds)$dex,
      #    sizeFactor = anno_points(colData(dds)$sizeFactor)
      #),
      show_row_names = FALSE, show_column_names = FALSE, row_km = 2,
      column_title = paste0(sum(l), " significant genes with FDR < ", fdr),
      show_row_dend = FALSE) + 
      Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
          name = "log10(baseMean+1)", show_column_names = FALSE) +
      Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
          name = "log2FoldChange", show_column_names = FALSE,
          col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)          
  ht <- draw(ht, merge_legend = TRUE)
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

# A self-defined action to respond brush event. It updates the MA-plot, the volcano plot
# and a table which contains DESeq2 results for the selected genes.
library(DT)
library(GetoptLong)
brush_action <- function(df, input, output, session) {

  row_index <- unique(unlist(df$row_index))
  selected <- env$row_index[row_index]

  output[["ma_plot"]] <- renderPlot({
    make_maplot(res, selected)
  })

  output[["volcano_plot"]] <- renderPlot({
    make_volcano(res, selected)
  })

  output[["res_table"]] <- renderDT(
    formatRound(datatable(res[selected, c("symbol", "baseMean", "log2FoldChange", "padj")], rownames = F), columns = 2:4, digits = 3)
  )

}

# The dashboard body contains three columns:
# 1. the original heatmap
# 2. the sub-heatmap and the default output
# 3. the self-defined output
library(shiny)
library(shinydashboard)
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
)

# Side bar contains settings for certain cutoffs to select significant genes.
ui <- secure_app(dashboardPage(
  dashboardHeader(title = "DESeq2 results"),
  dashboardSidebar(
    selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05)),
    numericInput("base_mean", label = "Minimal base mean:", value = 0),
    numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 1),
    actionButton("filter", label = "Generate heatmap")
  ),
  body
))

# makeInteractiveComplexHeatmap() is put inside observeEvent() so that changes on the cutoffs can regenerate the heatmap.
server <- function(input, output, session) {
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  observeEvent(input$filter, {
    ht <- make_heatmap(fdr = as.numeric(input$fdr), base_mean = input$base_mean, log2fc = input$log2fc)
    if(!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
        brush_action = brush_action)
    } else {
      # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap, thus, it is ht_heatmap here.
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui, server)