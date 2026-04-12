networkTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    div(style = "width:90%;", 
        plotOutput(ns("network"), width = "80%", height = "100%"))
  )
}

# Server module for WGCNA tab
networkTabServer <- function(id, ontologies, path, shared_values, view) {
  moduleServer(id, function(input, output, session) {
    
    output$network <- renderPlot({
      par(mar = c(2, 2, 2, 2))
      plotNetwork(ontologies, shared_values$ontology, shared_values$selected_row, path, view())
    }, 
    width = function() {
      session$clientData[[paste0("output_", session$ns("network"), "_width")]]
    },
    height = function() {
      2 * session$clientData[[paste0("output_", session$ns("network"), "_width")]]
    },
    res = 100)
    
  })
}