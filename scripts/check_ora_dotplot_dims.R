#!/usr/bin/env Rscript
# Dry-run ORA dot-plot layout metrics for representative scenarios.
# Usage (from repo root): Rscript scripts/check_ora_dotplot_dims.R

source("scripts/ora_tab.R")

mock_combined <- function(
    study_specs,
    n_pathways = 10L,
    pathway_prefix = "PATHWAY",
    long_pathway = FALSE) {
  rows <- list()
  pid <- 0L
  for (spec in study_specs) {
    sid <- spec$study_id
    slbl <- spec$study_label
    for (cmp in spec$comparisons) {
      for (p in seq_len(n_pathways)) {
        pid <- pid + 1L
        desc <- if (long_pathway && p == 1L) {
          paste0(
            "Very long hallmark pathway name that should wrap without squeezing panels ",
            p
          )
        } else {
          paste0(pathway_prefix, "_", p)
        }
        rows[[length(rows) + 1L]] <- data.frame(
          study_id = sid,
          study_label = slbl,
          comparison = cmp,
          ID = paste0("ID_", p),
          Description = desc,
          gene_ratio = runif(1, 0.1, 0.9),
          Count = sample(5:30, 1L),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

print_metrics <- function(label, m) {
  cat(
    sprintf(
      "%-16s width=%4d height=%4d pathways=%2d facets=%d\n",
      paste0(label, ":"),
      m$width,
      m$height,
      m$n_pathways,
      m$n_facets
    )
  )
}

run_case <- function(label, combined, study_ids, study_labels) {
  m <- ora_dotplot_layout_metrics(combined, study_ids, study_labels)
  print_metrics(label, m)
  invisible(m)
}

cat("ORA dot-plot layout metrics (scenario matrix)\n\n")

sparse <- mock_combined(
  list(list(study_id = "s1", study_label = "Study A", comparisons = "cmp_a")),
  n_pathways = 5L
)
m_sparse <- run_case("sparse", sparse, "s1", c(s1 = "Study A"))

activin_like <- mock_combined(
  list(list(
    study_id = "activin",
    study_label = "Activin",
    comparisons = c(
      "Activin_2_control_3h_CTB",
      "Activin_10_control_3h_CTB",
      "Activin_2_A83_control_3h_CTB",
      "Activin_10_A83_control_3h_CTB",
      "Activin_2_control_24h_CTB",
      "Activin_10_control_24h_CTB",
      "Activin_10_control_24h_expl"
    )
  )),
  n_pathways = 20L
)
m_activin <- run_case(
  "activin-like",
  activin_like,
  "activin",
  c(activin = "Activin")
)

multi_study <- mock_combined(
  lapply(seq_len(5L), function(i) {
    list(
      study_id = paste0("s", i),
      study_label = paste("Study", i),
      comparisons = paste0("cmp_", seq_len(3L))
    )
  }),
  n_pathways = 20L
)
m_multi <- run_case(
  "multi-study",
  multi_study,
  paste0("s", 1:5),
  stats::setNames(paste("Study", 1:5), paste0("s", 1:5))
)

long_labels <- mock_combined(
  list(list(
    study_id = "activin",
    study_label = "Activin",
    comparisons = c(
      "Activin_2_control_3h_CTB",
      "Activin_10_A83_control_3h_CTB",
      "Activin_10_control_24h_expl",
      "Activin_2_control_24h_CTB",
      "Activin_10_control_3h_CTB",
      "Activin_2_A83_control_3h_CTB",
      "Activin_10_control_24h_CTB"
    )
  )),
  n_pathways = 15L,
  long_pathway = TRUE
)
m_long <- run_case(
  "long-labels",
  long_labels,
  "activin",
  c(activin = "Activin")
)

cat("\nMonotonicity checks:\n")
stopifnot(m_activin$width > m_sparse$width)
stopifnot(m_activin$height > m_sparse$height)
stopifnot(m_multi$width > m_activin$width)
stopifnot(m_long$width >= m_activin$width)
stopifnot(m_long$height >= m_activin$height)
cat("OK: more DEGs/pathways/longer labels -> wider/taller (not smaller)\n")
