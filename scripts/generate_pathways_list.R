#!/usr/bin/env Rscript

# Generate a YAML pathways_list file from databases/pathways/*.txt.
#
# Usage:
#   Rscript scripts/generate_pathways_list.R
#   Rscript scripts/generate_pathways_list.R --out databases/pathways_list.yaml

args <- commandArgs(trailingOnly = TRUE)

out_path <- "databases/pathways_list.yaml"
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (a %in% c("-h", "--help")) {
    message("Usage: Rscript scripts/generate_pathways_list.R [--out <path>]")
    quit(status = 0L)
  }
  if (identical(a, "--out")) {
    if (i >= length(args)) {
      stop("--out requires a value", call. = FALSE)
    }
    out_path <- trimws(args[[i + 1L]])
    i <- i + 2L
    next
  }
  stop("Unknown option: ", a, call. = FALSE)
}

source("scripts/pathway_signatures.R")

pathways <- list_pathway_txt_files("databases/pathways")
if (length(pathways) < 1L) {
  stop("No pathway .txt files found in databases/pathways", call. = FALSE)
}

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
yaml::write_yaml(list(pathways_list = pathways), file = out_path)
message("Saved ", length(pathways), " pathways to ", out_path)
