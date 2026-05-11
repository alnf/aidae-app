#!/usr/bin/env Rscript
# Build rsconnect manifest from a deploy app YAML (see deploy/apps/*.yaml).
# Usage:
#   Rscript deploy/write_manifest.R deploy/apps/heart.yaml
#   Rscript deploy/write_manifest.R deploy/apps/heart.yaml --dry-run
#   Rscript deploy/write_manifest.R deploy/apps/heart.yaml --verbose   # or -v
#   Rscript deploy/write_manifest.R deploy/apps/heart.yaml --debug     # implies --verbose
# Env: EXPRS_DEPLOY_REPO = repo root (default: current working directory).

`%||%` <- function(x, y) {
  if (is.null(x) || (is.character(x) && !nzchar(x))) {
    y
  } else {
    x
  }
}

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 1L) {
  stop(
    "Usage: Rscript deploy/write_manifest.R <deploy/apps/....yaml> [--dry-run] [--verbose|-v] [--debug]",
    call. = FALSE
  )
}

dry_run <- "--dry-run" %in% argv
verbose_flag <- "--verbose" %in% argv || "-v" %in% argv
debug_flag <- "--debug" %in% argv
argv <- argv[!argv %in% c("--dry-run", "--verbose", "-v", "--debug")]

options(
  exprs.deploy.verbose = isTRUE(verbose_flag || debug_flag),
  exprs.deploy.debug = isTRUE(debug_flag)
)

repo_root <- Sys.getenv("EXPRS_DEPLOY_REPO", unset = "")
if (!nzchar(repo_root)) {
  repo_root <- getwd()
}
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)

app_yaml <- argv[[1L]]
if (!file.exists(app_yaml)) {
  app_yaml <- normalizePath(file.path(repo_root, app_yaml), winslash = "/", mustWork = TRUE)
} else {
  app_yaml <- normalizePath(app_yaml, winslash = "/", mustWork = TRUE)
}

old_wd <- getwd()
setwd(repo_root)
on.exit(setwd(old_wd), add = TRUE)

message("[exprs-deploy] repo: ", repo_root)

source(file.path(repo_root, "deploy/lib/bundle_files.R"), chdir = FALSE)

message("[exprs-deploy] parsing: ", app_yaml)
parsed <- parse_deploy_app(repo_root, app_yaml)
message("[exprs-deploy] deploy app (bundle-relative): ", parsed$app_yaml_relpath)
if (isTRUE(dry_run)) {
  parsed$bundle <- utils::modifyList(parsed$bundle, list(ora_validate = "warn"))
}

t_collect <- proc.time()
message("[exprs-deploy] collecting bundle paths …")
files <- collect_bundle_files(repo_root, parsed)
elapsed_collect <- (proc.time() - t_collect)[["elapsed"]]
message(
  sprintf(
    "[exprs-deploy] collected %d bundle path(s) in %.2fs",
    length(files),
    elapsed_collect
  )
)

if (isTRUE(dry_run)) {
  st <- bundle_bytes_and_count(repo_root, files)
  message("Dry run — paths: ", length(files), " ; file count (recursive): ", st$n_files)
  message("Approx total bytes: ", st$bytes, " (", round(st$bytes / 1024^2, 1), " MiB)")
  if (exprs_deploy_debug()) {
    message("Full path list:")
    message(paste(files, collapse = "\n"))
  } else if (exprs_deploy_verbose()) {
    message("Path list:")
    message(paste(files, collapse = "\n"))
  } else {
    message("Sample paths:")
    message(paste(head(files, 20L), collapse = "\n"))
    if (length(files) > 20L) {
      message("… (use --verbose or --debug for the full list)")
    }
  }
  quit(save = "no", status = 0L)
}

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install rsconnect: install.packages('rsconnect')", call. = FALSE)
}

wm_verbose <- isTRUE(verbose_flag || debug_flag)
wm_hint <- if (wm_verbose) {
  "verbose output enabled"
} else {
  "pass --verbose or --debug for more detail"
}
message(
  "[exprs-deploy] rsconnect::writeManifest() starting ",
  "(dependency capture is often several minutes the first time; ",
  wm_hint,
  ") …"
)

t_wm <- proc.time()
rsconnect::writeManifest(
  appDir = repo_root,
  appFiles = files,
  appPrimaryDoc = "app.R",
  verbose = wm_verbose,
  quiet = FALSE
)
elapsed_wm <- (proc.time() - t_wm)[["elapsed"]]
message(
  sprintf(
    "[exprs-deploy] writeManifest() finished in %.1f min",
    elapsed_wm / 60
  )
)

cache_dir <- file.path(repo_root, "deploy/.cache/last")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dep <- parsed$deploy
ln <- c(
  paste0("RSCONNECT_ACCOUNT=", dep$account %||% ""),
  paste0("RSCONNECT_TITLE=", dep$title %||% ""),
  paste0("RSCONNECT_APP_ID=", dep$app_id %||% "")
)
writeLines(ln, file.path(cache_dir, "rsconnect_deploy.env"))
message("[exprs-deploy] wrote manifest.json and ", file.path(cache_dir, "rsconnect_deploy.env"))
