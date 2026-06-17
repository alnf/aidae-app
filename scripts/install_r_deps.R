#!/usr/bin/env Rscript
# Install R package dependencies for exprs-app-r (see README.md).
#
# Usage:
#   Rscript scripts/install_r_deps.R           # install missing only
#   Rscript scripts/install_r_deps.R --check   # report missing, no install
#   Rscript scripts/install_r_deps.R --full    # also install Suggests deps (slower)
#
# ---------------------------------------------------------------------------
# System packages (Ubuntu/Debian) — install BEFORE running this script:
#
#   sudo apt-get update
#   sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
#     r-base-dev \
#     libcurl4-openssl-dev libxml2-dev libssl-dev libicu-dev \
#     libpng-dev libjpeg-dev libtiff5-dev \
#     libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
#     libcairo2-dev \
#     libuv1-dev \
#     libsodium-dev cmake \
#     libgit2-dev libssh2-1-dev \
#     libgdal-dev libudunits2-dev \
#     pkg-config
#
# Notes:
#   - libcairo2-dev: required for CRAN gdtools -> ggiraph (ORA / UpSet interactivity)
#   - libuv1-dev: required for CRAN fs -> shiny (script sets USE_BUNDLED_LIBUV=1 if absent)
#   - r-base-dev: headers to compile packages from source
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
full_deps <- "--full" %in% args

# --- CRAN -------------------------------------------------------------------
cran_pkgs <- c(
  # Install helper (not loaded by app, but required for bioc_pkgs)
  "BiocManager",

  # app.R — Shiny UI and tables
  "shiny", "bs4Dash", "shinymanager", "shinydashboard",
  "yaml", "DT", "GetoptLong",

  # app.R — heatmap colour scales (ComplexHeatmap itself is Bioconductor)
  "circlize",

  # Gene tab — boxplots (scripts/gene_plot.R, scripts/gene_tab.R)
  "ggplot2", "ggpubr", "ggnewscale",

  # ORA / UpSet tabs — interactive and static plots
  "ggiraph", "ggh4x", "ComplexUpset", "patchwork",

  # ggiraph / InteractiveComplexHeatmap support (often missing on fresh installs)
  "htmlwidgets", "gdtools",
  "knitr", "rmarkdown", "kableExtra",

  # clusterProfiler support on CRAN
  "RSQLite",

  # Heatmap raster backend (scripts/heatmap_utils.R; optional but recommended)
  "ragg",

  # Custom ORA ontology upload (scripts/pathway_signatures.R)
  "readxl",

  # shinymanager password hashing utility (scripts/hash_password.R)
  "scrypt",

  # Posit Connect deploy (deploy/write_manifest.R)
  "rsconnect"
)

# --- Bioconductor -----------------------------------------------------------
bioc_pkgs <- c(
  # app.R — interactive and static heatmaps
  "ComplexHeatmap", "InteractiveComplexHeatmap",

  # ORA stack — explicit deps that failed in cascade when top-level install broke
  "AnnotationDbi", "GO.db", "KEGGREST", "GOSemSim", "DOSE", "ggtree", "enrichplot",

  # ORA tab and precompute (scripts/ora_tab.R, scripts/ora_cache.R, scripts/precompute_ora.R)
  "clusterProfiler"
)

all_pkgs <- unique(c(cran_pkgs, bioc_pkgs))

missing_pkgs <- function(pkgs) {
  ip <- rownames(installed.packages())
  setdiff(pkgs, ip)
}

report <- function(label, pkgs) {
  miss <- missing_pkgs(pkgs)
  have <- setdiff(pkgs, miss)
  message("[install_r_deps] ", label, ": ", length(have), "/", length(pkgs), " installed")
  if (length(miss) > 0L) {
    message("  missing: ", paste(miss, collapse = ", "))
  }
  invisible(miss)
}

message(
  "[install_r_deps] R ", R.version.string,
  " | installed packages: ", nrow(installed.packages()),
  if (check_only) " | check only" else "",
  if (full_deps) " | full deps (Suggests)" else " | deps: Depends + Imports only"
)

miss_cran <- report("CRAN", cran_pkgs)
miss_bioc <- report("Bioconductor", bioc_pkgs)
miss_all <- unique(c(miss_cran, miss_bioc))

if (length(miss_all) < 1L) {
  message("[install_r_deps] All listed packages already installed.")
  quit(status = 0)
}

if (check_only) {
  quit(status = if (length(miss_all) > 0L) 1L else 0L)
}

lib <- .libPaths()[1L]
stale_locks <- Sys.glob(file.path(lib, "00LOCK*"))
if (length(stale_locks) > 0L) {
  message("[install_r_deps] Removing ", length(stale_locks), " stale library lock(s)")
  unlink(stale_locks, recursive = TRUE)
}

ncpus <- max(1L, parallel::detectCores())
options(repos = c(CRAN = "https://cloud.r-project.org"))
# Ncpus = 1: install one package at a time (avoids 00LOCK races when many pkgs build together).
# MAKEFLAGS = -jN: still parallelize compilation within each package.
options(Ncpus = 1L)
Sys.setenv(MAKEFLAGS = paste0("-j", ncpus))
# fs (shiny dependency): bundled libuv if libuv1-dev is not installed.
Sys.setenv(USE_BUNDLED_LIBUV = "1")

dep_mode <- if (full_deps) TRUE else c("Depends", "Imports")

message("[install_r_deps] Compile jobs per package: ", ncpus, " (sequential package installs)")

install_one <- function(pkg, installer = c("cran", "bioc")) {
  installer <- match.arg(installer)
  message("[install_r_deps] Installing ", pkg, " (", installer, ") ...")
  if (installer == "cran") {
    install.packages(pkg, dependencies = dep_mode)
  } else {
    BiocManager::install(pkg, update = FALSE, ask = FALSE, dependencies = dep_mode)
  }
}

if (length(miss_cran) > 0L) {
  for (pkg in miss_cran) {
    if (!pkg %in% missing_pkgs(cran_pkgs)) next
    install_one(pkg, "cran")
  }
}

if (length(miss_bioc) > 0L) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("[install_r_deps] Installing BiocManager")
    install.packages("BiocManager")
  }
  for (pkg in miss_bioc) {
    if (!pkg %in% missing_pkgs(bioc_pkgs)) next
    install_one(pkg, "bioc")
  }
}

message("[install_r_deps] Verification:")
still_missing <- missing_pkgs(all_pkgs)
if (length(still_missing) > 0L) {
  message("  STILL MISSING: ", paste(still_missing, collapse = ", "))
  quit(status = 1L)
}

message("  All listed packages OK")
quit(status = 0)
