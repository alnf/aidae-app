#!/usr/bin/env Rscript
# Dry-run ORA pipeline (load → filter → render) on a synthetic fixture with shards.
# Usage (from repo root): Rscript scripts/profile_ora_pipeline.R

Sys.setenv(EXPRS_ORA_PROFILE = "1")

source("scripts/perf_utils.R")
source("scripts/heatmap_utils.R")
source("scripts/study_data.R")
source("scripts/ora_cache.R")
source("scripts/ora_tab.R")

mock_long_rows <- function(pathway_files, comparisons, study_id, study_label, n_pathways, n_genes_per_pathway = 40L) {
  rows <- list()
  for (pf in pathway_files) {
    for (cmp in comparisons) {
      for (p in seq_len(n_pathways)) {
        genes <- paste0("GENE", seq_len(n_genes_per_pathway), collapse = "/")
        rows[[length(rows) + 1L]] <- data.frame(
          comparison = cmp,
          ID = paste0(tools::file_path_sans_ext(basename(pf)), "_ID_", p),
          Description = paste0("PATH_", tools::file_path_sans_ext(basename(pf)), "_", p),
          geneID = genes,
          gene_ratio = runif(1, 0.1, 0.9),
          Count = sample(5:30, 1L),
          p_value = runif(1, 1e-6, 0.05),
          p_adj = runif(1, 1e-5, 0.1),
          setSize = sample(50:500, 1L),
          study_label = study_label,
          pathway_file = pf,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

write_profile_fixture <- function(root) {
  pathway_a <- "Ontology_A.txt"
  pathway_b <- "Ontology_B.txt"
  studies <- c("profile_s1", "profile_s2", "profile_s3")
  comparisons <- c("cmp_one", "cmp_two", "cmp_three")

  for (sid in studies) {
    dir.create(file.path(root, "data", sid, "ora"), recursive = TRUE, showWarnings = FALSE)
    cfg <- list(
      name = sid,
      deg_lists = lapply(comparisons, function(cmp) {
        list(label = cmp, deg_file = paste0(cmp, ".tsv"))
      }),
      ora_file = "ora/enrichment.rds"
    )
    yaml::write_yaml(cfg, file.path(root, "data", sid, "config.yaml"))

    long_df <- mock_long_rows(
      c(pathway_a, pathway_b),
      comparisons,
      study_id = sid,
      study_label = paste("Study", sub("profile_s", "", sid)),
      n_pathways = 80L
    )
    obj <- list(
      version = ORA_CACHE_VERSION,
      study_id = sid,
      pathway_files = c(pathway_a, pathway_b),
      created = Sys.time(),
      long_df = long_df
    )
    saveRDS(obj, file.path(root, "data", sid, "ora", "enrichment.rds"))
    ora_sync_shards_from_long_df(sid, long_df, c(pathway_a, pathway_b), created = obj$created)
  }
  list(study_ids = studies, pathway_files = c(pathway_a, pathway_b))
}

profile_render_stage <- function(combined, study_ids, study_labels, pathway_level_order, pathway_file) {
  study_label_levels <- .ora_facet_study_label_levels(combined, study_ids, study_labels)
  t_prep <- ora_profile_timer_start()
  gp <- .ora_faceted_comparison_plot_girafe(
    combined,
    study_label_levels,
    pathway_level_order
  )
  prep_ms <- ora_profile_elapsed_ms(t_prep)
  dims <- ora_dotplot_layout_metrics(combined, study_ids, study_labels, pathway_level_order)
  t_girafe <- ora_profile_timer_start()
  if (requireNamespace("ggiraph", quietly = TRUE) && !is.null(gp)) {
    invisible(ggiraph::girafe(
      ggobj = gp,
      width_svg = max(6, dims$width / 96),
      height_svg = max(4, dims$height / 96),
      options = list(ggiraph::opts_sizing(rescale = FALSE))
    ))
  }
  girafe_ms <- ora_profile_elapsed_ms(t_girafe)
  ora_profile_log(
    "pathway=", pathway_file,
    " render_prep=", round(prep_ms, 1), "ms",
    " girafe=", round(girafe_ms, 1), "ms",
    " render_total=", round(prep_ms + girafe_ms, 1), "ms"
  )
}

run_pipeline <- function(study_ids, study_labels, pathway_file) {
  ora_profile_set_ctx(pathway_file)
  per <- ora_try_load_per_study_caches(study_ids, pathway_file)
  if (!isTRUE(per$ok_all)) {
    stop("Fixture load failed for ", pathway_file, call. = FALSE)
  }
  ob <- ora_build_plot_payload(
    per$long_by_sid,
    study_ids,
    study_labels,
    min_ol = 1L,
    min_ct = 1L,
    min_gene_ratio = 0,
    n_show = 20L
  )
  dfs <- list()
  for (sid in study_ids) {
    b <- ob$by_study[[sid]]
    if (!is.null(b$plot_df) && nrow(b$plot_df) > 0L) {
      df_sid <- b$plot_df
      df_sid$study_id <- sid
      dfs[[length(dfs) + 1L]] <- df_sid
    }
  }
  if (length(dfs) < 1L) {
    stop("No plot rows for ", pathway_file, call. = FALSE)
  }
  combined <- do.call(rbind, dfs)
  if (!is.null(ob$pathway_level_order) && length(ob$pathway_level_order) > 0L) {
    combined$Description <- factor(as.character(combined$Description), levels = ob$pathway_level_order)
  }
  profile_render_stage(combined, study_ids, study_labels, ob$pathway_level_order, pathway_file)
  invisible(TRUE)
}

old_wd <- getwd()
fixture_root <- tempfile("ora_profile_")
dir.create(fixture_root)
on.exit({
  setwd(old_wd)
  unlink(fixture_root, recursive = TRUE)
}, add = TRUE)

meta <- write_profile_fixture(fixture_root)
setwd(fixture_root)

study_ids <- meta$study_ids
study_labels <- stats::setNames(paste("Study", seq_along(study_ids)), study_ids)

cat("ORA pipeline profile (synthetic fixture: 3 studies, 2 ontologies, 80 pathways each)\n\n")

for (run in seq_len(2L)) {
  cat(sprintf("--- Pass %d ---\n", run))
  for (pf in meta$pathway_files) {
    run_pipeline(study_ids, study_labels, pf)
  }
  cat("\n")
}

cat("Done.\n")
