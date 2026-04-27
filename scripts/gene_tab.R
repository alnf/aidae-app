# Gene tab Shiny module: gene search, HTML gene info (from DEG tables), one boxplot per study.
# Depends on study_data.R, gene_plot.R (and heatmap_utils.R via gene_plot).

.gene_tab_msg_plot <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = msg, size = 4.2) +
    ggplot2::theme_void() +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))
}

# Plot height (px) for renderPlot: compact for messages, taller for real boxplots.
.gene_tab_plot_height_px <- function(symbol, study_id) {
  sym <- trimws(symbol)
  if (!nzchar(sym)) return(72L)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(80L)
  if (!requireNamespace("ggpubr", quietly = TRUE) || !requireNamespace("ggnewscale", quietly = TRUE)) {
    return(110L)
  }
  gt <- load_gene_tab_study(study_id)
  if (is.null(gt)) return(100L)
  ens <- resolve_ens_for_symbol_in_study(sym, study_id)
  if (is.null(ens)) return(130L)
  if (!ens %in% rownames(gt$mm)) return(130L)
  360L
}

# Plot width (px): # of PhenoNames groups (x categories) and facet count (gene_tab_facet).
.gene_tab_plot_width_px <- function(symbol, study_id) {
  sym <- trimws(symbol)
  if (!nzchar(sym)) return(360L)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(380L)
  if (!requireNamespace("ggpubr", quietly = TRUE) || !requireNamespace("ggnewscale", quietly = TRUE)) {
    return(420L)
  }
  gt <- load_gene_tab_study(study_id)
  if (is.null(gt)) return(400L)
  ens <- resolve_ens_for_symbol_in_study(sym, study_id)
  if (is.null(ens) || !ens %in% rownames(gt$mm)) return(440L)

  md <- gt$metadata
  pn_u <- unique(as.character(md$PhenoNames))
  pn_u <- pn_u[!is.na(pn_u) & nzchar(pn_u)]
  n_pheno <- length(pn_u)
  if (n_pheno < 1L) n_pheno <- 1L

  cfg <- gt$cfg
  fc_raw <- cfg$gene_tab_facet
  fc <- if (is.null(fc_raw) || !nzchar(as.character(fc_raw))) NULL else as.character(fc_raw)
  has_facet <- !is.null(fc) && fc %in% colnames(md) && any(!is.na(md[[fc]]))
  n_facet <- if (!has_facet) {
    1L
  } else {
    nf <- length(unique(as.character(md[[fc]])[!is.na(md[[fc]])]))
    max(1L, as.integer(nf))
  }

  np <- max(1L, as.integer(n_pheno))
  npf <- as.numeric(np)

  if (n_facet == 1L) {
    panel_w <- 200L + as.integer(round(58 * npf + 40 * sqrt(npf)))
  } else {
    pheno_core <- 100L + as.integer(round(52 * (npf^0.65) + 10 * sqrt(npf)))
    panel_w <- pheno_core + 120L + 48L
  }
  total <- as.integer(n_facet) * panel_w + if (n_facet > 1L) 88L else 64L
  wmax <- if (n_facet > 1L) 1800L else 1400L
  as.integer(min(wmax, max(280L, total)))
}

if (!exists(".gene_tab_deg_cache", inherits = FALSE)) {
  .gene_tab_deg_cache <- new.env(parent = emptyenv())
}

read_deg_table_cached <- function(study_id, deg_rel) {
  key <- paste(study_id, deg_rel, sep = "|")
  if (exists(key, envir = .gene_tab_deg_cache, inherits = FALSE)) {
    return(get(key, envir = .gene_tab_deg_cache, inherits = FALSE))
  }
  p <- file.path("data", study_id, deg_rel)
  if (!file.exists(p)) return(NULL)
  tab <- read.table(p, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  assign(key, tab, envir = .gene_tab_deg_cache)
  tab
}

# One display string per gene: keep first spelling seen, drop later case-only duplicates.
.dedup_symbols_ci <- function(syms) {
  syms <- trimws(as.character(syms))
  syms <- syms[!is.na(syms) & nzchar(syms)]
  if (length(syms) == 0L) return(character(0))
  keys <- tolower(syms)
  keep <- !duplicated(keys)
  out <- syms[keep]
  out[order(tolower(out))]
}

#' Union of `symbol` values from all DEG lists across studies (for selectize).
collect_all_gene_symbols <- function(study_ids) {
  syms <- character(0)
  for (sid in study_ids) {
    lists <- study_deg_lists(sid)
    for (entry in lists) {
      deg_rel <- entry$deg_file
      if (is.null(deg_rel) || deg_rel == "") next
      res <- read_deg_table_cached(sid, deg_rel)
      if (is.null(res) || !"symbol" %in% colnames(res)) next
      syms <- c(syms, as.character(res$symbol))
    }
  }
  .dedup_symbols_ci(syms)
}

resolve_ens_for_symbol_in_study <- function(symbol, study_id) {
  if (is.null(symbol) || !nzchar(trimws(symbol))) return(NULL)
  sym_l <- tolower(trimws(symbol))
  lists <- study_deg_lists(study_id)
  for (entry in lists) {
    deg_rel <- entry$deg_file
    if (is.null(deg_rel) || deg_rel == "") next
    res <- read_deg_table_cached(study_id, deg_rel)
    if (is.null(res) || !all(c("symbol", "ens_gene") %in% colnames(res))) next
    m <- tolower(as.character(res$symbol)) == sym_l
    if (any(m, na.rm = TRUE)) {
      eg <- as.character(res$ens_gene[which(m)[1L]])
      if (nzchar(eg)) return(eg)
    }
  }
  NULL
}

# First matching DEG row across studies (order = study_ids) for gene info HTML.
fetch_gene_info_from_degs <- function(symbol, study_ids) {
  if (is.null(symbol) || !nzchar(trimws(symbol))) return(NULL)
  sym_l <- tolower(trimws(symbol))
  for (sid in study_ids) {
    lists <- study_deg_lists(sid)
    for (entry in lists) {
      deg_rel <- entry$deg_file
      if (is.null(deg_rel) || deg_rel == "") next
      res <- read_deg_table_cached(sid, deg_rel)
      if (is.null(res) || !all(c("symbol", "ens_gene") %in% colnames(res))) next
      m <- tolower(as.character(res$symbol)) == sym_l
      if (!any(m, na.rm = TRUE)) next
      i <- which(m)[1L]
      desc <- if ("description" %in% colnames(res)) as.character(res$description[i]) else NA_character_
      return(list(
        basic = list(
          symbol = as.character(res[[i, "symbol"]]),
          ensembl_id = as.character(res[[i, "ens_gene"]]),
          description = desc
        ),
        study_id = sid
      ))
    }
  }
  NULL
}

createGeneInfoHTML <- function(info, gene_symbol) {
  if (is.null(info$basic)) {
    return("<div style='color:#e74c3c;padding:12px;'>No matching gene in DEG tables.</div>")
  }
  b <- info$basic
  ens <- htmltools::htmlEscape(as.character(b$ensembl_id), attribute = TRUE)
  sym <- htmltools::htmlEscape(as.character(b$symbol), attribute = TRUE)
  desc <- if (!is.na(b$description) && nzchar(as.character(b$description))) {
    paste0(
      "<tr><td style='padding:5px;font-weight:bold;'>Description:</td><td style='padding:5px;'>",
      htmltools::htmlEscape(as.character(b$description)),
      "</td></tr>"
    )
  } else {
    ""
  }
  ensembl_url <- paste0("https://www.ensembl.org/id/", ens)
  gc_url <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", sym)
  paste0(
    "<div style='margin-bottom:16px;padding:12px;background:#fff;border-radius:4px;'>",
    "<h5 style='color:#2c3e50;margin-bottom:10px;'>Gene information</h5>",
    "<table style='width:100%;border-collapse:collapse;'>",
    "<tr><td style='padding:5px;font-weight:bold;width:120px;'>Symbol:</td><td style='padding:5px;'>", sym, "</td></tr>",
    "<tr><td style='padding:5px;font-weight:bold;'>Ensembl:</td><td style='padding:5px;'>", ens, "</td></tr>",
    desc,
    "</table></div>",
    "<div style='padding:10px;background:#ecf0f1;border-radius:4px;'>",
    "<span style='font-weight:bold;margin-right:8px;'>Links</span>",
    "<a href='", ensembl_url, "' target='_blank' rel='noopener' style='color:#3498db;margin-right:12px;'>Ensembl</a>",
    "<a href='", gc_url, "' target='_blank' rel='noopener' style='color:#3498db;margin-right:12px;'>GeneCards</a>",
    "<a href='https://pubmed.ncbi.nlm.nih.gov/?term=", utils::URLencode(gene_symbol), "' target='_blank' rel='noopener' style='color:#3498db;'>PubMed</a>",
    "</div>"
  )
}

geneTabUI <- function(id, study_ids, study_labels) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 12,
        shiny::selectizeInput(
          ns("symbol"),
          "Gene symbol:",
          choices = NULL,
          width = "320px",
          options = list(placeholder = "Type a gene symbol…")
        ),
        shiny::conditionalPanel(
          condition = "input.symbol != ''",
          ns = ns,
          shiny::div(
            style = "background-color:#f8f9fa;padding:15px;border-radius:5px;margin-bottom:20px;border-left:4px solid #3498db;",
            shiny::h4("Gene information", style = "color:#2c3e50;margin-bottom:12px;"),
            shiny::uiOutput(ns("gene_info"))
          )
        ),
        shiny::tags$hr(),
        shiny::tags$div(
          lapply(study_ids, function(sid) {
            lbl <- study_labels[[sid]]
            if (is.null(lbl) || !nzchar(as.character(lbl))) lbl <- sid
            shiny::tags$div(
              class = "mb-3",
              shiny::h4(as.character(lbl), class = "text-primary", style = "margin-bottom:6px;"),
            shiny::tags$div(
              style = "overflow-x: auto; width: 100%;",
              shiny::tags$div(
                style = "display: inline-block;",
                shiny::plotOutput(ns(paste0("gene_plot_", sid)), width = "auto", height = "auto")
              )
            )
            )
          })
        )
      )
    )
  )
}

geneTabServer <- function(id, study_ids, study_labels, external_symbol = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    gene_choices <- shiny::reactiveVal(character(0))

    shiny::observe({
      syms <- collect_all_gene_symbols(study_ids)
      gene_choices(syms)
      shiny::updateSelectizeInput(
        session, "symbol",
        choices = syms,
        selected = character(0),
        server = TRUE
      )
    }) # runs once: no reactive inputs inside

    output$gene_info <- shiny::renderUI({
      sym_raw <- input$symbol
      if (is.null(sym_raw)) sym_raw <- ""
      sym <- trimws(sym_raw)
      if (!nzchar(sym)) return(NULL)
      info <- fetch_gene_info_from_degs(sym, study_ids)
      if (is.null(info)) {
        return(shiny::HTML("<div style='color:#e74c3c;padding:8px;'>No matching gene in any DEG table for this dashboard.</div>"))
      }
      shiny::HTML(createGeneInfoHTML(info, sym))
    })

    shiny::observeEvent(external_symbol(), {
      sym_raw <- external_symbol()
      sym <- trimws(as.character(sym_raw))
      if (!nzchar(sym)) return()
      choices <- gene_choices()
      if (length(choices) < 1L) {
        choices <- collect_all_gene_symbols(study_ids)
        gene_choices(choices)
      }
      key <- tolower(sym)
      hit <- which(tolower(as.character(choices)) == key)
      sym_sel <- if (length(hit) > 0L) as.character(choices[[hit[[1L]]]]) else sym
      shiny::updateSelectizeInput(
        session,
        "symbol",
        choices = choices,
        selected = sym_sel,
        server = TRUE
      )
    }, ignoreInit = TRUE)

    for (sid in study_ids) {
      local({
        study_id <- sid
        out_name <- paste0("gene_plot_", study_id)
        output[[out_name]] <- shiny::renderPlot(
          {
            sym_raw <- input$symbol
            if (is.null(sym_raw)) sym_raw <- ""
            sym <- trimws(sym_raw)
            if (!nzchar(sym)) {
              return(.gene_tab_msg_plot("Enter a gene symbol."))
            }
            gt <- load_gene_tab_study(study_id)
            if (is.null(gt)) {
              return(.gene_tab_msg_plot("Could not load counts/metadata for this study."))
            }
            ens <- resolve_ens_for_symbol_in_study(sym, study_id)
            if (is.null(ens)) {
              return(plot_gene_study(
                ens_gene = "__not_in_study__",
                mm = gt$mm,
                metadata = gt$metadata,
                gene_symbol = sym,
                study_title = study_labels[[study_id]],
                is_count_like = isTRUE(gt$is_count_like)
              ))
            }
            gdegs <- load_gene_tab_gdegs(study_id)
            pv <- if (!is.null(gdegs) && nrow(gdegs) > 0L && "ens_gene" %in% colnames(gdegs)) {
              gdegs[as.character(gdegs$ens_gene) == ens, , drop = FALSE]
            } else {
              NULL
            }
            cfg <- gt$cfg
            fc <- cfg$gene_tab_facet
            if (is.null(fc) || !nzchar(as.character(fc))) {
              fc <- NULL
            } else {
              fc <- as.character(fc)
            }
            plot_gene_study(
              ens_gene = ens,
              mm = gt$mm,
              metadata = gt$metadata,
              gene_symbol = sym,
              pval_df = pv,
              study_title = study_labels[[study_id]],
              facet_column = fc,
              is_count_like = isTRUE(gt$is_count_like)
            )
          },
          height = function() {
            s <- input$symbol
            if (is.null(s)) s <- ""
            .gene_tab_plot_height_px(s, study_id)
          },
          width = function() {
            s <- input$symbol
            if (is.null(s)) s <- ""
            .gene_tab_plot_width_px(s, study_id)
          }
        )
      })
    }
  })
}
