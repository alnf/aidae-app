# Helpers for Posit Connect bundle file lists (sourced by deploy/write_manifest.R).
# Depends: yaml (same as the app).

ORA_CACHE_VERSION_DEPLOY <- 2L
ORA_CACHE_VERSION_LEGACY_DEPLOY <- 1L

deploy_reserved_keys <- function() {
  c("deploy", "bundle", "extends")
}

default_bundle_options <- function() {
  list(
    omit_databases = TRUE,
    omit_orig_dirs = TRUE,
    omit_counts_tabular_when_rds = TRUE,
    omit_rsconnect_python = TRUE,
    ora_validate = "strict" # strict | warn | skip
  )
}

exprs_deploy_debug <- function() {
  isTRUE(getOption("exprs.deploy.debug", FALSE))
}

exprs_deploy_verbose <- function() {
  isTRUE(getOption("exprs.deploy.verbose", FALSE)) || exprs_deploy_debug()
}

#' Shallow merge: overlay wins for same top-level names.
merge_overlay <- function(base, overlay) {
  if (length(base) < 1L) {
    return(overlay)
  }
  out <- base
  for (nm in names(overlay)) {
    out[[nm]] <- overlay[[nm]]
  }
  out
}

#' Read a deploy app YAML: runtime keys + optional `extends`, `deploy`, `bundle`.
parse_deploy_app <- function(repo_root, app_yaml_path) {
  app_yaml_path <- normalizePath(app_yaml_path, winslash = "/", mustWork = TRUE)
  root_n <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  pref <- paste0(root_n, "/")
  if (!startsWith(app_yaml_path, pref) && !identical(app_yaml_path, root_n)) {
    stop("app_yaml_path must be inside repo_root: ", app_yaml_path, call. = FALSE)
  }
  raw <- yaml::read_yaml(app_yaml_path)
  if (is.null(raw) || !is.list(raw)) {
    stop("Invalid or empty deploy YAML: ", app_yaml_path)
  }
  reserved <- deploy_reserved_keys()
  ext <- raw$extends
  deploy_block <- if (!is.null(raw$deploy)) raw$deploy else list()
  bundle_block <- if (!is.null(raw$bundle)) raw$bundle else list()
  overlay <- raw[setdiff(names(raw), reserved)]
  base <- if (!is.null(ext) && nzchar(as.character(ext))) {
    p <- normalizePath(file.path(repo_root, as.character(ext)), winslash = "/", mustWork = FALSE)
    if (!file.exists(p)) {
      stop("extends file not found: ", p)
    }
    yaml::read_yaml(p)
  } else {
    list()
  }
  main <- merge_overlay(base, overlay)
  yaml_n <- app_yaml_path
  if (startsWith(yaml_n, pref)) {
    rel <- substring(yaml_n, nchar(pref) + 1L)
  } else {
    rel <- basename(yaml_n)
  }
  extends_relpath <- NA_character_
  if (!is.null(ext) && nzchar(trimws(as.character(ext)[[1L]]))) {
    extends_relpath <- trimws(as.character(ext)[[1L]])
  }
  list(
    main = main,
    deploy = deploy_block,
    bundle = utils::modifyList(default_bundle_options(), bundle_block),
    app_yaml_relpath = rel,
    extends_relpath = extends_relpath
  )
}

study_ora_rds_abs_path_deploy <- function(repo_root, study_id) {
  cfg_path <- file.path(repo_root, "data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    return(NA_character_)
  }
  cfg <- yaml::read_yaml(cfg_path)
  rel <- if (!is.null(cfg$ora_file) && nzchar(as.character(cfg$ora_file))) {
    as.character(cfg$ora_file)
  } else {
    "ora/enrichment.rds"
  }
  normalizePath(file.path(repo_root, "data", study_id, rel), winslash = "/", mustWork = FALSE)
}

study_ora_shard_abs_path_deploy <- function(repo_root, study_id, pathway_rel_file) {
  monolith <- study_ora_rds_abs_path_deploy(repo_root, study_id)
  if (is.na(monolith) || !nzchar(monolith)) {
    return(NA_character_)
  }
  base <- tools::file_path_sans_ext(basename(as.character(pathway_rel_file)))
  normalizePath(
    file.path(dirname(monolith), "enrichment", paste0(base, ".rds")),
    winslash = "/",
    mustWork = FALSE
  )
}

study_ora_rds_compatible_deploy <- function(obj, pathway_rel_file, study_id) {
  if (is.null(obj) || !is.list(obj)) {
    return(FALSE)
  }
  if (!identical(as.character(obj$study_id), as.character(study_id))) {
    return(FALSE)
  }
  if (!is.data.frame(obj$long_df)) {
    return(FALSE)
  }
  pf <- as.character(pathway_rel_file)
  ver <- suppressWarnings(as.integer(obj$version))
  if (length(ver) != 1L || is.na(ver)) {
    return(FALSE)
  }
  if (ver == ORA_CACHE_VERSION_LEGACY_DEPLOY) {
    if (is.null(obj$pathway_file)) {
      return(FALSE)
    }
    return(identical(as.character(obj$pathway_file), pf))
  }
  if (ver == ORA_CACHE_VERSION_DEPLOY) {
    pfs <- obj$pathway_files
    if (is.null(pfs) || length(pfs) < 1L) {
      if ("pathway_file" %in% colnames(obj$long_df)) {
        pfs <- unique(as.character(obj$long_df$pathway_file))
      } else {
        return(FALSE)
      }
    } else {
      pfs <- as.character(pfs)
    }
    return(pf %in% pfs)
  }
  FALSE
}

study_needs_ora <- function(repo_root, study_id) {
  cfg_path <- file.path(repo_root, "data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    return(FALSE)
  }
  cfg <- yaml::read_yaml(cfg_path)
  has_lists <- !is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L
  has_legacy <- !is.null(cfg$deg_file) && nzchar(as.character(cfg$deg_file))
  isTRUE(has_lists || has_legacy)
}

#' Resolve pathways_list from merged main config (same rules as pathway_signatures.R).
resolve_pathways_for_main <- function(repo_root, main_cfg) {
  raw <- main_cfg$pathways_list
  if (is.null(raw)) {
    return(character(0))
  }
  vals <- trimws(as.character(unlist(raw)))
  vals <- vals[nzchar(vals)]
  if (length(vals) < 1L) {
    return(character(0))
  }
  if (length(vals) == 1L) {
    p1 <- vals[[1L]]
    full <- normalizePath(file.path(repo_root, p1), winslash = "/", mustWork = FALSE)
    if (file.exists(full)) {
      ext <- tolower(tools::file_ext(p1))
      if (ext %in% c("yaml", "yml")) {
        obj <- yaml::read_yaml(full)
        if (is.list(obj) && !is.null(obj$pathways_list)) {
          v2 <- trimws(as.character(unlist(obj$pathways_list)))
          return(sort(unique(v2[nzchar(v2)])))
        }
        v2 <- trimws(as.character(unlist(obj)))
        return(sort(unique(v2[nzchar(v2)])))
      }
      lines <- readLines(full, warn = FALSE, encoding = "UTF-8")
      lines <- trimws(lines)
      lines <- lines[nzchar(lines)]
      lines <- lines[!startsWith(lines, "#")]
      if (length(lines) < 1L) {
        return(character(0))
      }
      v2 <- trimws(unlist(strsplit(lines, ",", fixed = TRUE)))
      return(sort(unique(v2[nzchar(v2)])))
    }
  }
  sort(unique(vals))
}

collect_study_file_paths <- function(repo_root, study_id) {
  cfg_path <- file.path(repo_root, "data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    stop("Missing study config: ", cfg_path)
  }
  cfg <- yaml::read_yaml(cfg_path)
  rels <- file.path("data", study_id, "config.yaml")
  add_rel <- function(rel) {
    if (is.null(rel) || length(rel) < 1L) {
      return(NULL)
    }
    rel <- as.character(rel)
    rel <- rel[nzchar(trimws(rel))]
    if (length(rel) < 1L) {
      return(NULL)
    }
    file.path("data", study_id, rel)
  }
  rels <- c(rels, add_rel(cfg$counts_file))
  rels <- c(rels, add_rel(cfg$metadata_file))
  rels <- c(rels, add_rel(cfg$comparison_file))
  rels <- c(rels, add_rel(cfg$gdegs_file))
  rels <- c(rels, add_rel(cfg$ora_file))
  monolith_rel <- if (!is.null(cfg$ora_file) && nzchar(as.character(cfg$ora_file))) {
    as.character(cfg$ora_file)
  } else {
    "ora/enrichment.rds"
  }
  shard_dir <- file.path(repo_root, "data", study_id, dirname(monolith_rel), "enrichment")
  if (dir.exists(shard_dir)) {
    shard_files <- list.files(shard_dir, pattern = "\\.rds$", full.names = FALSE)
    for (bn in shard_files) {
      rels <- c(rels, file.path("data", study_id, dirname(monolith_rel), "enrichment", bn))
    }
  }
  rels <- c(rels, add_rel(cfg$deg_file))
  if (!is.null(cfg$deg_lists) && length(cfg$deg_lists) > 0L) {
    for (e in cfg$deg_lists) {
      rels <- c(rels, add_rel(e$deg_file))
      rels <- c(rels, add_rel(e$counts_file))
    }
  }
  unique(rels[!is.na(rels)])
}

path_under_orig <- function(rel_path) {
  grepl("(^|/)orig(/|$)", rel_path, ignore.case = TRUE)
}

is_tabular_counts <- function(rel_path) {
  tolower(tools::file_ext(rel_path)) %in% c("tsv", "txt", "csv")
}

maybe_drop_tabular_counts <- function(repo_root, rel_paths, omit_counts_tabular_when_rds) {
  if (!isTRUE(omit_counts_tabular_when_rds)) {
    return(rel_paths)
  }
  # When config points at counts.tsv but counts.rds exists, we omit the TSV to
  # shrink the bundle — but we must add the .rds / .eds path (not only drop the
  # TSV), otherwise the deployed app has no counts file at all (same basename
  # rules as scripts/study_data.R).
  drop_idx <- logical(length(rel_paths))
  add <- character(0)
  for (i in seq_along(rel_paths)) {
    rel <- rel_paths[[i]]
    if (!is_tabular_counts(rel)) {
      next
    }
    abs_p <- file.path(repo_root, rel)
    if (!file.exists(abs_p)) {
      next
    }
    rds_abs <- sub("\\.(tsv|txt|csv)$", ".rds", abs_p, ignore.case = TRUE)
    eds_abs <- sub("\\.(tsv|txt|csv|rds)$", ".eds", abs_p, ignore.case = TRUE)
    if (identical(eds_abs, abs_p)) {
      eds_abs <- paste0(abs_p, ".eds")
    }
    has_rds <- file.exists(rds_abs)
    has_eds <- file.exists(eds_abs)
    if (!has_rds && !has_eds) {
      next
    }
    drop_idx[[i]] <- TRUE
    if (has_rds) {
      rds_rel <- sub("\\.(tsv|txt|csv)$", ".rds", rel, ignore.case = TRUE)
      add <- c(add, rds_rel)
    }
    if (has_eds) {
      eds_rel <- sub("\\.(tsv|txt|csv|rds)$", ".eds", rel, ignore.case = TRUE)
      if (identical(eds_rel, rel)) {
        eds_rel <- paste0(rel, ".eds")
      }
      add <- c(add, eds_rel)
    }
  }
  unique(c(rel_paths[!drop_idx], add))
}

#' Repo-relative pathway list *files* (yaml/yml/txt with a slash), for thin bundles.
#' Inline basename-only lists do not need a pointer file on disk.
collect_pathway_list_pointer_files <- function(repo_root, main_cfg, study_ids) {
  ptr_one <- function(raw) {
    if (is.null(raw)) {
      return(character(0))
    }
    um <- unlist(raw, use.names = FALSE)
    if (length(um) != 1L) {
      return(character(0))
    }
    p1 <- trimws(as.character(um[[1L]]))
    if (!nzchar(p1)) {
      return(character(0))
    }
    ext <- tolower(tools::file_ext(p1))
    if (!ext %in% c("yaml", "yml", "txt")) {
      return(character(0))
    }
    if (!grepl("/", p1, fixed = TRUE)) {
      return(character(0))
    }
    full <- normalizePath(file.path(repo_root, p1), winslash = "/", mustWork = FALSE)
    if (file.exists(full) && !dir.exists(full)) {
      return(p1)
    }
    character(0)
  }
  out <- if (is.null(main_cfg)) character(0) else ptr_one(main_cfg$pathways_list)
  for (sid in as.character(study_ids)) {
    scp <- file.path(repo_root, "data", sid, "config.yaml")
    if (!file.exists(scp)) {
      next
    }
    sc <- yaml::read_yaml(scp)
    out <- c(out, ptr_one(sc$pathways_list))
  }
  unique(out[nzchar(out)])
}

collect_pathway_database_files <- function(repo_root, main_cfg, study_ids) {
  # Mirror list_ora_pathway_files logic without sourcing the app.
  configured <- resolve_pathways_for_main(repo_root, main_cfg)
  for (sid in as.character(study_ids)) {
    scp <- file.path(repo_root, "data", sid, "config.yaml")
    if (!file.exists(scp)) {
      next
    }
    sc <- yaml::read_yaml(scp)
    raw <- sc$pathways_list
    if (is.null(raw)) {
      next
    }
    v <- trimws(as.character(unlist(raw)))
    v <- v[nzchar(v)]
    if (length(v) == 1L) {
      full <- file.path(repo_root, v[[1L]])
      if (file.exists(full)) {
        ext <- tolower(tools::file_ext(v[[1L]]))
        if (ext %in% c("yaml", "yml")) {
          obj <- yaml::read_yaml(full)
          if (is.list(obj) && !is.null(obj$pathways_list)) {
            v <- trimws(as.character(unlist(obj$pathways_list)))
          } else {
            v <- trimws(as.character(unlist(obj)))
          }
        }
      }
    }
    configured <- c(configured, v[nzchar(v)])
  }
  configured <- sort(unique(configured[nzchar(configured)]))
  if (length(configured) < 1L) {
    return(character(0))
  }
  # List files under databases/pathways
  out <- character(0)
  for (bn in configured) {
    if (grepl("[/\\\\]", bn)) {
      next
    }
    out <- c(out, file.path("databases", "pathways", bn))
  }
  # Include pathways_list pointer file if it is a repo path
  rawm <- main_cfg$pathways_list
  if (!is.null(rawm) && length(as.character(unlist(rawm))) == 1L) {
    p1 <- as.character(unlist(rawm))[[1L]]
    if (nzchar(p1) && file.exists(file.path(repo_root, p1))) {
      out <- c(out, p1)
    }
  }
  unique(out)
}

validate_ora_precompute <- function(
    repo_root,
    study_ids,
    pathway_basenames,
    mode = c("strict", "warn", "skip")) {
  mode <- match.arg(mode)
  if (identical(mode, "skip") || length(pathway_basenames) < 1L) {
    return(invisible(TRUE))
  }
  msgs <- character(0)
  for (sid in study_ids) {
    if (!study_needs_ora(repo_root, sid)) {
      next
    }
    abs_rds <- study_ora_rds_abs_path_deploy(repo_root, sid)
    has_monolith <- !is.na(abs_rds) && file.exists(abs_rds)
    monolith_obj <- NULL
    if (has_monolith) {
      monolith_obj <- tryCatch(readRDS(abs_rds), error = function(e) NULL)
      if (is.null(monolith_obj)) {
        msgs <- c(msgs, sprintf("study %s: unreadable ORA RDS %s", sid, abs_rds))
      }
    }
    for (pf in pathway_basenames) {
      shard_abs <- study_ora_shard_abs_path_deploy(repo_root, sid, pf)
      has_shard <- !is.na(shard_abs) && file.exists(shard_abs)
      ok <- FALSE
      if (has_shard) {
        shard_obj <- tryCatch(readRDS(shard_abs), error = function(e) NULL)
        if (!is.null(shard_obj) && study_ora_rds_compatible_deploy(shard_obj, pf, sid)) {
          ok <- TRUE
        } else {
          msgs <- c(
            msgs,
            sprintf("study %s: shard not compatible with pathway_file=%s (%s)", sid, pf, shard_abs)
          )
        }
      } else if (has_monolith && !is.null(monolith_obj) &&
          study_ora_rds_compatible_deploy(monolith_obj, pf, sid)) {
        ok <- TRUE
      }
      if (!ok && !has_shard) {
        if (!has_monolith) {
          msgs <- c(msgs, sprintf("study %s: missing ORA monolith and shard for %s", sid, pf))
        } else if (is.null(monolith_obj)) {
          NULL
        } else {
          msgs <- c(
            msgs,
            sprintf("study %s: ORA RDS not compatible with pathway_file=%s", sid, pf)
          )
        }
      }
    }
  }
  if (length(msgs) < 1L) {
    return(invisible(TRUE))
  }
  txt <- paste(msgs, collapse = "\n")
  if (identical(mode, "strict")) {
    stop("ORA precompute validation failed:\n", txt, call. = FALSE)
  }
  message("ORA precompute validation (warn):\n", txt)
  invisible(TRUE)
}

collect_bundle_files <- function(repo_root, parsed) {
  bundle <- parsed$bundle
  main <- parsed$main
  study_ids <- main$studies
  if (is.null(study_ids) || length(as.character(unlist(study_ids))) < 1L) {
    stop("Merged main config must define non-empty `studies` for deploy.", call. = FALSE)
  }
  study_ids <- as.character(unlist(study_ids))

  if (exprs_deploy_verbose()) {
    message(
      "[exprs-deploy] bundle: ",
      length(study_ids),
      " studies — ",
      paste(study_ids, collapse = ", ")
    )
    message(
      "[exprs-deploy] bundle: omit_databases=",
      isTRUE(bundle$omit_databases),
      ", omit_orig=",
      isTRUE(bundle$omit_orig_dirs),
      ", omit_counts_tsv_when_rds=",
      isTRUE(bundle$omit_counts_tabular_when_rds),
      ", ora_validate=",
      if (is.null(bundle$ora_validate)) "strict" else as.character(bundle$ora_validate)[1L]
    )
  }

  files <- c("app.R", "LICENSE", "scripts/")
  files <- c(files, parsed$app_yaml_relpath)
  extf <- parsed$extends_relpath
  if (length(extf) == 1L && !is.na(extf) && nzchar(extf)) {
    if (file.exists(file.path(repo_root, extf))) {
      files <- c(files, extf)
    }
  }

  for (sid in study_ids) {
    sp <- collect_study_file_paths(repo_root, sid)
    if (exprs_deploy_verbose()) {
      message(sprintf("[exprs-deploy] bundle: study %s — %d path(s) from config", sid, length(sp)))
    }
    if (exprs_deploy_debug()) {
      message(paste0("  ", sp, collapse = "\n"))
    }
    for (p in sp) {
      abs_p <- normalizePath(file.path(repo_root, p), winslash = "/", mustWork = FALSE)
      if (!file.exists(abs_p)) {
        stop("Bundle reference missing on disk: ", p, call. = FALSE)
      }
      files <- c(files, p)
    }
  }

  if (!isTRUE(bundle$omit_databases)) {
    if (exprs_deploy_verbose()) {
      message("[exprs-deploy] bundle: resolving databases/pathways files …")
    }
    db_files <- collect_pathway_database_files(repo_root, main, study_ids)
    for (p in db_files) {
      abs_p <- file.path(repo_root, p)
      if (file.exists(abs_p)) {
        files <- c(files, p)
      }
    }
  } else if (isTRUE(bundle$omit_databases)) {
    ptrs <- collect_pathway_list_pointer_files(repo_root, main, study_ids)
    for (p in ptrs) {
      if (file.exists(file.path(repo_root, p))) {
        files <- c(files, p)
      }
    }
  }

  if (isTRUE(bundle$omit_rsconnect_python)) {
    # never add rsconnect-python
  }

  files <- unique(files)
  n_before_orig <- length(files)
  if (isTRUE(bundle$omit_orig_dirs)) {
    files <- files[!path_under_orig(files)]
    if (exprs_deploy_verbose() && length(files) < n_before_orig) {
      message(
        sprintf(
          "[exprs-deploy] bundle: dropped %d path(s) under .../orig/",
          n_before_orig - length(files)
        )
      )
    }
  }
  n_before_counts <- length(files)
  files <- maybe_drop_tabular_counts(repo_root, files, bundle$omit_counts_tabular_when_rds)
  if (exprs_deploy_verbose() && length(files) != n_before_counts) {
    message(
      sprintf(
        "[exprs-deploy] bundle: counts tabular ↔ serialized swap changed path count by %+d (see omit_counts_tabular_when_rds)",
        length(files) - n_before_counts
      )
    )
  }

  extra <- bundle$extra_include
  if (!is.null(extra) && length(extra) > 0L) {
    ex <- as.character(unlist(extra))
    ex <- ex[nzchar(trimws(ex))]
    files <- unique(c(files, ex))
  }

  globs <- bundle$extra_omit_globs
  if (!is.null(globs) && length(globs) > 0L) {
    for (g in as.character(unlist(globs))) {
      pat <- utils::glob2rx(g)
      files <- files[!grepl(pat, files, ignore.case = TRUE)]
    }
  }

  # Final existence check
  for (p in files) {
    if (endsWith(p, "/")) {
      d <- file.path(repo_root, sub("/$", "", p))
      if (!dir.exists(d)) {
        stop("Bundle directory missing: ", p, call. = FALSE)
      }
    } else {
      if (!file.exists(file.path(repo_root, p))) {
        stop("Bundle file missing: ", p, call. = FALSE)
      }
    }
  }

  pathways <- resolve_pathways_for_main(repo_root, main)
  if (exprs_deploy_verbose()) {
    message(
      "[exprs-deploy] bundle: ",
      length(pathways),
      " resolved pathway basename(s) for ORA validation"
    )
    if (length(pathways) > 0L && length(pathways) <= 30L) {
      message("  ", paste(pathways, collapse = ", "))
    } else if (length(pathways) > 30L) {
      message("  ", paste(head(pathways, 15L), collapse = ", "), " … (", length(pathways), " total)")
    }
  }
  if (isTRUE(bundle$omit_databases) && length(pathways) < 1L) {
    warning(
      "omit_databases is TRUE but pathways_list resolves empty: ",
      "ORA dropdown may fall back to scanning databases/pathways (empty in thin deploy). ",
      "Set inline pathways_list in the deploy YAML.",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (exprs_deploy_verbose()) {
    message(
      "[exprs-deploy] bundle: ORA precompute check (",
      if (is.null(bundle$ora_validate)) "strict" else as.character(bundle$ora_validate)[1L],
      ") …"
    )
  }
  validate_ora_precompute(
    repo_root,
    study_ids,
    pathways,
    mode = if (is.null(bundle$ora_validate)) "strict" else as.character(bundle$ora_validate)[1L]
  )
  if (exprs_deploy_verbose()) {
    message("[exprs-deploy] bundle: ORA precompute check finished")
  }

  files <- sort(unique(files))
  if (exprs_deploy_debug()) {
    message("[exprs-deploy] bundle: final path list (", length(files), "):")
    message(paste(files, collapse = "\n"))
  }
  files
}

bundle_bytes_and_count <- function(repo_root, files) {
  n <- 0L
  bytes <- 0
  for (p in files) {
    ap <- file.path(repo_root, p)
    if (dir.exists(ap)) {
      inner <- list.files(ap, recursive = TRUE, full.names = TRUE, all.files = FALSE)
      for (f in inner) {
        if (file.exists(f) && !dir.exists(f)) {
          n <- n + 1L
          bytes <- bytes + file.info(f)$size
        }
      }
    } else if (file.exists(ap)) {
      n <- n + 1L
      bytes <- bytes + file.info(ap)$size
    }
  }
  list(n_files = n, bytes = bytes)
}
