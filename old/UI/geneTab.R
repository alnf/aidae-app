# Define the module UI for Tab 1
geneTabUI <- function(id) {
  ns <- NS(id)
  tagList(
    # selectizeInput provides autocomplete suggestions
    selectizeInput(
      ns("symbol"),
      "Enter gene name:",
      choices = NULL,
      options = list(
        placeholder = 'Type gene name ...'
      )
    ),
    # Boxplot
    plotOutput(ns("boxplot"), width = "100%"),
    # Gene information panel
    conditionalPanel(
      condition = "input.symbol != ''",
      ns = ns,
      div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 15px; border-left: 4px solid #3498db;",
        h4("Gene Information", style = "color: #2c3e50; margin-bottom: 15px;"),
        uiOutput(ns("gene_info"))
      )
    )
  )
}

# Define the module server for Tab 1
geneTabServer <- function(id, t2g, exprs, mdata, mcols, gdf, groups) {
  moduleServer(id, function(input, output, session) {
    # Update choices using gene names from your data frame.
    updateSelectizeInput(session, "symbol",
                         choices = unique(t2g$symbol),
                         selected = character(0),
                         server = TRUE)

    # Reactive expression to filter the data based on the selected gene.
    selectedData <- reactive({
      req(input$symbol)
      # Check if the gene exists in the data.
      if (input$symbol %in% t2g$symbol) {
        input$symbol
      } else {
        NULL
      }
    })

    # Render gene information
    output$gene_info <- renderUI({
      req(selectedData())
      gene <- selectedData()

      # Fetch information
      ens_gene = unique(t2g$ens_gene[which(t2g$symbol == gene)])
      if (length(ens_gene)>1) {
        ens_gene = ens_gene[1]
      }
      info <- fetchGeneInfo(ens_gene, t2g, "rat")

      if (!is.null(info) && length(info) > 0) {
        info_html <- createGeneInfoHTML(info, gene)
        HTML(info_html)
      } else {
        HTML("<div style='color: #e74c3c; text-align: center; padding: 20px;'>Unable to fetch gene information</div>")
      }
    })

    # Render the boxplot when data is available.
    output$boxplot <- renderPlot({
      gene <- selectedData()
      plotGene(gene, exprs, groups(), NULL, mdata, mcols, gdf, id_type="symbol")
    }, res=100)
  })
}

# Fetch gene information
fetchGeneInfo <- function(ens_gene, t2g, species) {
  # Initialize info list
  info_list <- list()

  # Get basic info from t2g data
  gene_data <- t2g[t2g$ens_gene == ens_gene, ]
  if (nrow(gene_data) > 0) {
    info_list$basic <- list(
      symbol = gene_data$symbol[1],
      ensembl_id = gene_data$ens_gene[1],
      description = gene_data$description[1]
    )
  }

  # mapping <- mygene::queryMany(
  #   ens_gene,
  #   scopes = "ensembl.gene",
  #   fields = c("entrezgene", "type_of_gene"),
  #   species = species
  # )
  # 
  # entrez_id <- mapping$entrezgene[!is.na(mapping$entrezgene)][1]

  # Then fetch the NCBI Gene summary
  # summary <- tryCatch({
  #   s <- rentrez::entrez_summary(db = "gene", id = entrez_id)
  # }, error = function(e) NA_character_)
  
  # summary <- NA
  # print(summary)
  # if (!is.na(summary)) {
  #   info_list$ncbi <- list(
  #     entrez_id = entrez_id,
  #     name = summary[["name"]],
  #     chromosome = summary[["chromosome"]],
  #     map_location = summary[["maplocation"]],
  #     type_of_gene = mapping$type_of_gene[1],
  #     summary = summary[["summary"]]
  #   )
  # }

  return(info_list)
}


# Helper function to create HTML for gene information
createGeneInfoHTML <- function(info, gene_symbol) {
  html_parts <- c()

  # Basic information section
  if (!is.null(info$basic)) {
    html_parts <- c(html_parts, 
      "<div style='margin-bottom: 20px; padding: 10px; background-color: white; border-radius: 3px;'>",
      "<h5 style='color: #2c3e50; margin-bottom: 10px;'>Basic Information</h5>",
      "<table style='width: 100%; border-collapse: collapse;'>",
      "<tr><td style='padding: 5px; font-weight: bold; width: 120px;'>Symbol:</td><td style='padding: 5px;'>", info$basic$symbol, "</td></tr>",
      "<tr><td style='padding: 5px; font-weight: bold;'>Ensembl ID:</td><td style='padding: 5px;'>", info$basic$ensembl_id, "</td></tr>"
    )

    if (!is.na(info$basic$description) && info$basic$description != "") {
      html_parts <- c(html_parts,
        "<tr><td style='padding: 5px; font-weight: bold;'>Description:</td><td style='padding: 5px;'>", info$basic$description, "</td></tr>"
      )
    }

    html_parts <- c(html_parts, "</table></div>")
  }

  # NCBI detailed information
  if (!is.null(info$ncbi)) {
    ncbi <- info$ncbi
    html_parts <- c(html_parts,
      "<div style='margin-bottom: 20px; padding: 10px; background-color: white; border-radius: 3px;'>",
      "<h5 style='color: #2c3e50; margin-bottom: 10px;'>NCBI Gene Information</h5>",
      "<table style='width: 100%; border-collapse: collapse;'>"
    )

    if (!is.null(ncbi$name)) {
      html_parts <- c(html_parts,
        "<tr><td style='padding: 5px; font-weight: bold; width: 120px;'>Official Name:</td><td style='padding: 5px;'>", ncbi$name, "</td></tr>"
      )
    }

    if (!is.null(ncbi$entrez_id)) {
      html_parts <- c(html_parts,
        "<tr><td style='padding: 5px; font-weight: bold; width: 120px;'>Entrez ID:</td><td style='padding: 5px;'>", ncbi$entrez_id, "</td></tr>"
      )
    }

    if (!is.null(ncbi$chromosome)) {
      html_parts <- c(html_parts,
        "<tr><td style='padding: 5px; font-weight: bold;'>Chromosome:</td><td style='padding: 5px;'>", ncbi$chromosome, "</td></tr>"
      )
    }

    if (!is.null(ncbi$map_location)) {
      html_parts <- c(html_parts,
        "<tr><td style='padding: 5px; font-weight: bold;'>Map Location:</td><td style='padding: 5px;'>", ncbi$map_location, "</td></tr>"
      )
    }

    if (!is.null(ncbi$type_of_gene)) {
      html_parts <- c(html_parts,
        "<tr><td style='padding: 5px; font-weight: bold;'>Gene Type:</td><td style='padding: 5px;'>", ncbi$type_of_gene, "</td></tr>"
      )
    }

    if (!is.na(ncbi$summary) && ncbi$summary != "") {
      html_parts <- c(html_parts,
        "<tr><td style='padding: 5px; font-weight: bold;'>Description:</td><td style='padding: 5px;'>", ncbi$summary, "</td></tr>"
      )
    }

    html_parts <- c(html_parts, "</table></div>")
  }

  # External links section
  html_parts <- c(html_parts,
    "<div style='margin-top: 15px; padding: 10px; background-color: #ecf0f1; border-radius: 3px;'>",
    "<h5 style='color: #2c3e50; margin-bottom: 10px;'>External Links</h5>",
    "<div style='display: flex; gap: 15px; flex-wrap: wrap;'>"
  )

  # NCBI link
  if (!is.null(info$ncbi$entrez_id)) {
    html_parts <- c(html_parts,
      "<a href='https://www.ncbi.nlm.nih.gov/gene/", info$ncbi$entrez_id, "' target='_blank' style='color: #3498db; text-decoration: none; padding: 5px 10px; background-color: white; border-radius: 3px; border: 1px solid #bdc3c7;'>",
      "NCBI Gene</a>"
    )
  }

  # Ensembl link
  if (!is.null(info$basic$ensembl_id)) {
    html_parts <- c(html_parts,
      "<a href='https://ensembl.org/Rattus_norvegicus/Gene/Summary?g=", info$basic$ensembl_id, "' target='_blank' style='color: #3498db; text-decoration: none; padding: 5px 10px; background-color: white; border-radius: 3px; border: 1px solid #bdc3c7;'>",
      "Ensembl</a>"
    )
  }

  # GeneCards link
  if (!is.null(info$basic$symbol)) {
    html_parts <- c(html_parts,
      "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", info$basic$symbol, "' target='_blank' style='color: #3498db; text-decoration: none; padding: 5px 10px; background-color: white; border-radius: 3px; border: 1px solid #bdc3c7;'>",
      "GeneCards</a>"
    )
  }

  # General search links
  html_parts <- c(html_parts,
    "<a href='https://pubmed.ncbi.nlm.nih.gov/?term=", gene_symbol, "+heart' target='_blank' style='color: #3498db; text-decoration: none; padding: 5px 10px; background-color: white; border-radius: 3px; border: 1px solid #bdc3c7;'>",
    "NCBI Search</a>",
    "<a href='https://www.uniprot.org/uniprot/?query=", gene_symbol, "+organism_id:10116' target='_blank' style='color: #3498db; text-decoration: none; padding: 5px 10px; background-color: white; border-radius: 3px; border: 1px solid #bdc3c7;'>",
    "UniProt</a>",
    "</div></div>"
  )

  paste(html_parts, collapse = "")
}