#!/usr/bin/env Rscript
# Precompute ORA: one RDS per study under data/<study_id>/<ora_file from config>.
# Default: bundle every *.txt under databases/pathways/ into a single long_df per study
# (column `pathway_file` + `comparison`). RDS version 2; legacy version-1 single-pathway
# files still load when selecting that pathway.
#
# Usage (from project root):
#   Rscript scripts/precompute_ora.R
#   Rscript scripts/precompute_ora.R MSigDB_Hallmark_2020.txt

args <- commandArgs(trailingOnly = TRUE)

normalize_pathway_arg <- function(s) {
  pathway_file <- s
  if (grepl("^databases/pathways/", s)) {
    pathway_file <- sub("^databases/pathways/", "", s)
  }
  pathway_file
}

if (length(args) >= 1L && args[[1L]] %in% c("-h", "--help")) {
  message(
    "Usage:\n",
    "  Rscript scripts/precompute_ora.R\n",
    "    Precompute ORA for all *.txt files under databases/pathways/.\n",
    "  Rscript scripts/precompute_ora.R <pathway_filename>\n",
    "    Precompute for one file only (under databases/pathways/)."
  )
  quit(status = 0L)
}

root <- Sys.getenv("EXPRS_APP_ROOT", unset = "")
if (nzchar(root)) {
  setwd(root)
}

source("scripts/heatmap_utils.R")
source("scripts/result_table_indices.R")
source("scripts/study_data.R")
source("scripts/perf_utils.R")
source("scripts/pathway_signatures.R")
source("scripts/ora_cache.R")

main_config <- if (file.exists("config.yaml")) yaml::read_yaml("config.yaml") else list(studies = character(0))
study_ids <- if (!is.null(main_config$studies) && length(main_config$studies) > 0L) {
  main_config$studies
} else {
  data_dirs <- list.dirs("data", full.names = TRUE, recursive = FALSE)
  keep <- vapply(data_dirs, function(d) file.exists(file.path(d, "config.yaml")), logical(1L))
  basename(data_dirs[keep])
}

study_labels <- character(length(study_ids))
for (i in seq_along(study_ids)) {
  id <- study_ids[[i]]
  cfg_path <- file.path("data", id, "config.yaml")
  if (file.exists(cfg_path)) {
    sc <- yaml::read_yaml(cfg_path)
    study_labels[i] <- if (!is.null(sc$name)) sc$name else id
  } else {
    study_labels[i] <- id
  }
}
names(study_labels) <- study_ids

if (length(args) < 1L || !nzchar(trimws(paste(args, collapse = "")))) {
  pathway_files_to_run <- list_pathway_txt_files()
} else {
  pathway_files_to_run <- normalize_pathway_arg(args[[1L]])
  pathway_files_to_run <- c(pathway_files_to_run)
}

if (length(pathway_files_to_run) < 1L) {
  stop("No pathway .txt files to process (check databases/pathways/).", call. = FALSE)
}

message(
  "Pathway file(s): ",
  length(pathway_files_to_run),
  if (length(pathway_files_to_run) <= 5L) {
    paste0(" — ", paste(pathway_files_to_run, collapse = ", "))
  } else {
    paste0(" (showing first 5: ", paste(head(pathway_files_to_run, 5L), collapse = ", "), ", …)")
  }
)
message("Studies: ", paste(study_ids, collapse = ", "))

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  stop("Install clusterProfiler: BiocManager::install(\"clusterProfiler\")", call. = FALSE)
}

n_pf <- length(pathway_files_to_run)

for (sid in study_ids) {
  out_path <- study_ora_rds_abs_path(sid)
  if (is.na(out_path) || !nzchar(out_path)) {
    message("[skip] ", sid, " — no config")
    next
  }
  if (length(study_deg_lists(sid)) == 0L) {
    message("[skip] ", sid, " — no deg_lists")
    next
  }

  slbl <- study_labels[[sid]]
  n_lists <- length(study_deg_lists(sid))
  den <- max(1L, n_lists * n_pf)

  parts <- list()
  pathway_files_ok <- character(0)

  for (pi in seq_len(n_pf)) {
    pf <- pathway_files_to_run[[pi]]
    t2g <- parse_pathway_file_to_term2gene(pf)
    if (is.null(t2g) || nrow(t2g) == 0L) {
      message("[skip pathway] ", sid, " | ", pf, " — missing or empty")
      next
    }

    t2g_u <- data.frame(
      term = t2g$term,
      gene = toupper(trimws(as.character(t2g$gene))),
      stringsAsFactors = FALSE
    )

    chunk <- ora_enrichment_long_df_single_study(
      sid,
      slbl,
      t2g_u,
      min_gs_size = 1L,
      max_gs_size = 50000L,
      progress = NULL,
      den = den,
      pathway_rel_file = pf
    )

    if (is.null(chunk)) {
      chunk <- ora_long_df_empty_rows(pf)
    }

    pathway_files_ok <- c(pathway_files_ok, pf)
    parts[[length(parts) + 1L]] <- chunk
  }

  if (length(parts) < 1L) {
    message("[skip] ", sid, " — no usable pathway files")
    next
  }

  long_df <- do.call(rbind, parts)

  obj <- list(
    version = ORA_CACHE_VERSION,
    study_id = sid,
    pathway_files = pathway_files_ok,
    created = Sys.time(),
    long_df = long_df
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, out_path)
  message(
    "[saved] ", sid, " -> ", out_path,
    " (", nrow(long_df), " rows; ", length(pathway_files_ok), " pathway files)"
  )
}

message("Done.")
