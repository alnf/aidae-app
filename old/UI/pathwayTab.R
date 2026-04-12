# Define the module UI for Tab 1
pathwayTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Include custom CSS to control heatmap container width
    tags$head(
      tags$style(HTML("
        /* Default to 100% width for narrow screens */
        #heatmap-container {
          width: 100%;
        }
        /* For screens wider than 768px, set the container to a fixed percentage (e.g., 70%) */
        @media (min-width: 768px) {
          #heatmap-container {
            width: 90%;
            margin-left: auto;
            margin-right: auto;          }
        }
        @media (min-width: 992px) {
          #heatmap-container {
            width: 60%;
            margin-left: auto;
            margin-right: auto;          }
          }
        }
        @media (min-width: 1200px) {
          #heatmap-container {
            width: 40%;
            margin-left: auto;
            margin-right: auto;          }
          }
        }
      "))
    ),
    # selectizeInput provides autocomplete suggestions
    selectizeInput(
      ns("ontology"),
      "Select ontology:",
      choices = NULL,
      options = list(
        placeholder = 'Start typing ...'
      )
    ),
    girafeOutput(ns("lollipop"), width="90%"),
    div(id = "heatmap-container", uiOutput(ns("heatmap")))
  )
}

# Define the module server for Tab 1
pathwayTabServer <- function(id, ontologies, t2gh, gdf, path, groups, view, shared_values) {
  moduleServer(id, function(input, output, session) {
    cdata <- session$clientData
    #print(cdata)

    ns <- session$ns
    
    # Update choices using gene names from your data frame.
    updateSelectizeInput(session, "ontology",
                         choices = names(ontologies),
                         selected = character(0),
                         server = TRUE)
    
    # Reactive expression to filter the data based on the selected gene.
    selectedData <- reactive({
      req(input$ontology)
      
      # Check if the gene exists in the data.
      if (input$ontology %in% names(ontologies)) {
        input$ontology
      } else {
        NULL
      }
    })
    
    # Create a reactive expression that returns both the plot and the row count.
    plotData <- reactive({
      cnames <- names(cdata)
      allvalues <- lapply(cnames, function(name) {
       paste(name, cdata[[name]], sep = " = ")
      })
      paste(allvalues, collapse = "\n")
      #print(allvalues)
      
      ontology <- selectedData()
      container_width <- cdata[["output_pathway-lollipop_width"]] 
      out <- plotPathways(ontologies, ontology, t2gh, ns, container_width, path, view())
      list(plot = out[[1]], n = out[[2]])
    })
    
    output$lollipop <- renderGirafe({
      #req(session$clientData$output_pathway-lollipop_width)
      # Convert the container width (in pixels) to inches (assuming 96 dpi)
      #width_in <- session$clientData$output_pathway-lollipop_width / 96
      width_in = 13
      # Convert the desired pixel height to inches (assuming 96 pixels per inch)
      height_in <- ((plotData()$n+5) * 40) / 96
      girafe(
        ggobj = plotData()$plot,
        width_svg = width_in,                # set a fixed width in inches, or compute as needed
        height_svg = height_in, 
        options = list(
          opts_selection(type = "single", only_shiny = FALSE),
          opts_sizing(rescale = TRUE)  # Makes the SVG responsive to container
        )
      )
    })  
    # Reactive value to store selected row click
    selected_row <- reactiveVal(NULL)
    
    # Reset the clicked row whenever ontology selection changes
    observeEvent(input$ontology, {
      selected_row(NULL)  # reset clicked row
      shared_values$ontology <- NULL
      shared_values$selected_row <- NULL
      output$heatmap <- renderUI(NULL)  # clear the second plot
    })
    
    # Listen for clicks on the row names and display a heatmap plot
    observeEvent(input$row_click, {
      selected_row(input$row_click)
      
      req(selectedData(), selected_row())
      
      shared_values$ontology <- selectedData()
      shared_values$selected_row <- selected_row()
      
      heatmapData <- plotPathwayHeatmap(selectedData(), selected_row(), gdf, groups(), path, view())
      n_rows <- heatmapData[["n_rows"]] 

      # Create heatmap plot
      output$heatmap <- renderUI({
        plotOutput(ns("heatmap_plot"), height = paste0(120 * n_rows ^ 0.5, "px"))
      })
      output$heatmap_plot <- renderPlot({
        heatmapData[["ht"]]  
      })
    })

  })
}