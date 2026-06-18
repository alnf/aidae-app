#!/usr/bin/env Rscript
# Split per-study monolithic ora/enrichment.rds into per-ontology shards (no enricher recompute).
#
# Creates data/<study_id>/ora/enrichment/<pathway_basename>.rds from existing long_df rows.
# Backs up enrichment.rds to enrichment.rds.bak (never deletes the monolith).
#
# Usage (from repo root):
#   Rscript scripts/split_ora_enrichment.R
#   Rscript scripts/split_ora_enrichment.R --study <study_id>
#   Rscript scripts/split_ora_enrichment.R --dry-run
#   Rscript scripts/split_ora_enrichment.R --force   # overwrite existing .bak backup

args_all <- commandArgs(trailingOnly = TRUE)

parse_split_cli <- function(args) {
  out <- list(
    study = NA_character_,
    dry_run = FALSE,
    force = FALSE,
    show_help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      out$show_help <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--dry-run")) {
      out$dry_run <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--force")) {
      out$force <- TRUE
      i <- i + 1L
      next
    }
    if (identical(a, "--study")) {
      if (i >= length(args)) {
        stop("--study requires a study id", call. = FALSE)
      }
      out$study <- as.character(args[[i + 1L]])
      i <- i + 2L
      next
    }
    stop("Unknown argument: ", a, call. = FALSE)
  }
  out
}

print_split_help <- function() {
  cat(
    "Split monolithic ORA enrichment.rds into per-ontology shards (no enricher).\n\n",
    "  Rscript scripts/split_ora_enrichment.R [--study <id>] [--dry-run] [--force]\n\n",
    "Shards: data/<study>/ora/enrichment/<pathway_basename>.rds\n",
    "Backup: enrichment.rds.bak (skipped if backup exists unless --force)\n",
    "Monolith enrichment.rds is never removed.\n",
    sep = ""
  )
}

source("scripts/study_data.R")
source("scripts/ora_cache.R")

parsed <- parse_split_cli(args_all)
if (isTRUE(parsed$show_help)) {
  print_split_help()
  quit(status = 0)
}

resolve_study_ids <- function(single_study) {
  if (!is.na(single_study) && nzchar(single_study)) {
    return(as.character(single_study))
  }
  cfg_path <- "config.yaml"
  if (!file.exists(cfg_path)) {
    stop("config.yaml not found; pass --study <id>", call. = FALSE)
  }
  cfg <- yaml::read_yaml(cfg_path)
  ids <- as.character(unlist(cfg$studies))
  ids <- ids[nzchar(ids)]
  if (length(ids) < 1L) {
    stop("No studies in config.yaml", call. = FALSE)
  }
  ids
}

split_one_study <- function(sid, dry_run = FALSE, force = FALSE) {
  monolith <- study_ora_rds_abs_path(sid)
  if (is.na(monolith) || !file.exists(monolith)) {
    message("[skip] ", sid, ": no monolith at ", monolith)
    return(invisible(list(wrote = 0L, skipped = TRUE)))
  }
  obj <- tryCatch(readRDS(monolith), error = function(e) NULL)
  if (is.null(obj) || !is.data.frame(obj$long_df)) {
    stop("Unreadable or empty monolith: ", monolith, call. = FALSE)
  }
  ver <- suppressWarnings(as.integer(obj$version))
  if (length(ver) != 1L || is.na(ver) || ver != ORA_CACHE_VERSION) {
    stop("Monolith is not ORA cache version ", ORA_CACHE_VERSION, ": ", monolith, call. = FALSE)
  }
  if (!"pathway_file" %in% colnames(obj$long_df)) {
    stop("Monolith long_df has no pathway_file column: ", monolith, call. = FALSE)
  }
  pfs <- unique(as.character(obj$long_df$pathway_file))
  pfs <- pfs[nzchar(pfs)]
  if (length(pfs) < 1L) {
    message("[skip] ", sid, ": no pathway_file values in long_df")
    return(invisible(list(wrote = 0L, skipped = TRUE)))
  }

  backup <- paste0(monolith, ".bak")
  if (file.exists(backup) && !isTRUE(force)) {
    message("[backup] ", sid, ": ", basename(backup), " exists (use --force to recopy)")
  } else if (!dry_run) {
    ok <- file.copy(monolith, backup, overwrite = isTRUE(force))
    if (!ok) {
      stop("Failed to copy backup: ", backup, call. = FALSE)
    }
    message("[backup] ", sid, ": ", monolith, " -> ", backup)
  } else {
    message("[dry-run] would backup ", monolith, " -> ", backup)
  }

  cre <- if (!is.null(obj$created)) obj$created else Sys.time()
  wrote <- 0L
  for (pf in pfs) {
    shard <- study_ora_shard_abs_path(sid, pf)
    n_rows <- sum(as.character(obj$long_df$pathway_file) == pf)
    if (dry_run) {
      message("[dry-run] ", sid, " ", pf, " -> ", shard, " (", n_rows, " rows)")
      wrote <- wrote + 1L
      next
    }
    ora_write_study_pathway_shard(sid, pf, obj$long_df, created = cre)
    message("[shard] ", sid, " ", pf, " -> ", shard, " (", n_rows, " rows)")
    wrote <- wrote + 1L
  }
  invisible(list(wrote = wrote, skipped = FALSE))
}

study_ids <- resolve_study_ids(parsed$study)
message(
  "Split ORA enrichment shards for ", length(study_ids), " study/studies",
  if (parsed$dry_run) " (dry-run)" else ""
)

total <- 0L
for (sid in study_ids) {
  r <- split_one_study(sid, dry_run = parsed$dry_run, force = parsed$force)
  if (!is.null(r$wrote)) {
    total <- total + as.integer(r$wrote)
  }
}
message("Done. Shards written (or planned): ", total)
