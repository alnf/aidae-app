# Config tab: session visibility toggles for studies and DEG lists.
# Depends: scripts/study_data.R (build_study_deg_catalog, parse_visibility_from_inputs, default_visibility_state).

configTabUI <- function(id, catalog) {
  ns <- shiny::NS(id)
  study_boxes <- lapply(catalog, function(item) {
    sid <- item$study_id
    lists <- item$deg_entries
    deg_choices <- stats::setNames(
      vapply(lists, function(e) as.character(e$deg_file), character(1L)),
      vapply(lists, function(e) {
        lbl <- e$label
        if (is.null(lbl) || !nzchar(as.character(lbl))) {
          lbl <- basename(as.character(e$deg_file))
        }
        as.character(lbl)
      }, character(1L))
    )
    all_deg <- unname(deg_choices)
    bs4Dash::box(
      title = item$study_label,
      width = 12,
      solidHeader = TRUE,
      status = "secondary",
      collapsible = TRUE,
      shiny::checkboxInput(
        ns(paste0("study_", sid)),
        label = "Include study",
        value = TRUE
      ),
      shiny::tags$div(
        style = "margin-left: 1.25rem;",
        shiny::checkboxGroupInput(
          ns(paste0("deg_", sid)),
          label = "DEG lists",
          choices = deg_choices,
          selected = all_deg
        )
      )
    )
  })
  shiny::tagList(
    bs4Dash::box(
      title = "Study and DEG list visibility",
      width = 12,
      solidHeader = TRUE,
      status = "primary",
      shiny::tags$p(
        class = "text-muted",
        "Choose which studies and DEG lists appear across the dashboard (ORA, UpSet, Gene, sidebar selectors). ",
        "Choices apply for this browser session only. ",
        "The main config.yaml still controls which studies are available in the app."
      ),
      shiny::fluidRow(
        shiny::column(
          width = 12,
          shiny::actionButton(ns("select_all"), "Select all", class = "btn-sm btn-outline-primary"),
          shiny::actionButton(ns("deselect_all"), "Deselect all", class = "btn-sm btn-outline-secondary"),
          shiny::actionButton(ns("reset_defaults"), "Reset to defaults", class = "btn-sm btn-outline-secondary")
        )
      ),
      shiny::tags$hr(),
      shiny::textOutput(ns("visibility_summary"))
    ),
    study_boxes
  )
}

configTabServer <- function(id, catalog) {
  shiny::moduleServer(id, function(input, output, session) {
    apply_defaults <- function() {
      d <- default_visibility_state(catalog)
      for (item in catalog) {
        sid <- item$study_id
        all_deg <- vapply(item$deg_entries, function(e) as.character(e$deg_file), character(1L))
        on <- sid %in% d$study_ids
        shiny::updateCheckboxInput(session, paste0("study_", sid), value = on)
        shiny::updateCheckboxGroupInput(
          session,
          paste0("deg_", sid),
          selected = if (on) d$deg_by_study[[sid]] else character(0)
        )
      }
    }

    for (item in catalog) {
      local({
        sid <- item$study_id
        all_deg <- vapply(item$deg_entries, function(e) as.character(e$deg_file), character(1L))
        shiny::observeEvent(input[[paste0("study_", sid)]], {
          on <- isTRUE(input[[paste0("study_", sid)]])
          shiny::updateCheckboxGroupInput(
            session,
            paste0("deg_", sid),
            selected = if (on) all_deg else character(0)
          )
        }, ignoreInit = TRUE)
      })
    }

    shiny::observeEvent(input$select_all, {
      for (item in catalog) {
        sid <- item$study_id
        all_deg <- vapply(item$deg_entries, function(e) as.character(e$deg_file), character(1L))
        shiny::updateCheckboxInput(session, paste0("study_", sid), value = TRUE)
        shiny::updateCheckboxGroupInput(session, paste0("deg_", sid), selected = all_deg)
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$deselect_all, {
      for (item in catalog) {
        sid <- item$study_id
        shiny::updateCheckboxInput(session, paste0("study_", sid), value = FALSE)
        shiny::updateCheckboxGroupInput(session, paste0("deg_", sid), selected = character(0))
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$reset_defaults, apply_defaults(), ignoreInit = TRUE)

    visibility <- shiny::reactive({
      # Module `input` ids are already namespaced (e.g. `study_<sid>`), not `config-study_<sid>`.
      parse_visibility_from_inputs(catalog, input, "")
    })

    output$visibility_summary <- shiny::renderText({
      vis <- visibility()
      n_studies <- length(vis$study_ids)
      n_degs <- sum(vapply(vis$deg_by_study, length, integer(1L)))
      sprintf("Showing %d study/studies and %d DEG list(s).", n_studies, n_degs)
    })

    list(visibility = visibility)
  })
}
