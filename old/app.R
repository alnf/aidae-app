#library(remotes)
#remotes::install_github('datastorm-open/shinymanager')
#library(conflicted)
library(tidyverse)
library(shiny)
#library(CEMiTool)
library(shinymanager)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggnewscale)
library(ggh4x)
library(ggiraph)
library(grid)
library(tidyr)
#library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(stringr)
library(circlize)
library(igraph)
library(tidygraph)
library(ggraph)
library(rbioapi)
library(STRINGdb)
library(cowplot)
library(patchwork)
#library(mygene)
#library(rentrez)
source("UI/geneTab.R")
source("UI/pathwayTab.R")
source("UI/wgcnaTab.R")
source("UI/networkTab.R")
source("scripts/utils.R")

creds <- read.table("data/creds.txt", sep="\t", header = T)
#labs <- get_labels("en")
#labs[["Login"]] <- "Login"
#do.call(set_labels, c(list("en"), labs))

# https://github.com/datastorm-open/shinymanager/issues/195
credentials <- data.frame(
  user     = c(creds$user),
  password = c(creds$password),
  start    = c("2025-01-01"),
  expire   = c(NA),
  admin    = c(FALSE),
  stringsAsFactors = FALSE,
  is_hashed_password = TRUE
)

ui <- secure_app(
  fluidPage(
    titlePanel("Rat heart RNA-seq Shiny App"),
    selectInput("view_selector", "Select View:", choices = c("PE", "norma"), selected = "PE"),
    tabsetPanel(
      tabPanel("Gene", geneTabUI("gene")),
      tabPanel("WGCNA", wgcnaTabUI("wgcna")),
      tabPanel("Pathway", pathwayTabUI("pathway")),
      tabPanel("Network", networkTabUI("network"))
    )
  )
)

path <- "data/"
readData(path)

server <- function(input, output, session) {
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  view <- reactive({input$view_selector})
  groups <- reactive({
    if (view() == "PE") {
      c("WTpreg", "PEpreg", "WTpost", "PEpost")
    } else if (view() == "norma") {
      c("np", "WTpreg", "WTpost")
    }
  })

  shared_values <- reactiveValues(
    ontology = NULL,
    selected_row = NULL
  )
  
  geneTabServer("gene", t2g, exprs, mdata, mcols, gdf, groups)
  wgcnaTabServer("wgcna", mcols, exprs, mdata, groups, view)
  pathwayTabServer("pathway", ontologies, t2gh, gdf, path, groups, view, shared_values)
  networkTabServer("network", ontologies, path, shared_values, view)
}

shinyApp(ui, server)
