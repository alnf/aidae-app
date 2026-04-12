# UI module for WGCNA tab (initially empty)
wgcnaTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns("wgcna_heatmap"), height = "860px", width = "1000")
  )
}

# Server module for WGCNA tab (initially empty)
wgcnaTabServer <- function(id, mcols, exprs, mdata, groups, view) {
  moduleServer(id, function(input, output, session) {
    output$wgcna_heatmap <- renderPlot({
      # Assuming plotWGCNA() returns a ComplexHeatmap object
      ht <- plotWGCNA(mcols, exprs, mdata, groups(), view())
      
      # Properly render the heatmap
      ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
    }, res = 100)
    
  })    
}