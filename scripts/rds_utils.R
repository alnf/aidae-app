convert_study_counts_to_rds <- function(study_id, overwrite = FALSE) {
  if (is.null(study_id) || study_id == "") {
    stop("study_id must be a non-empty string")
  }

  cfg_path <- file.path("data", study_id, "config.yaml")
  if (!file.exists(cfg_path)) {
    stop("Config file not found for study_id=", study_id, " at ", cfg_path)
  }

  cfg <- yaml::read_yaml(cfg_path)
  if (is.null(cfg$counts_file)) {
    stop("counts_file is not defined in config.yaml for study_id=", study_id)
  }

  counts_path <- file.path("data", study_id, cfg$counts_file)
  if (!file.exists(counts_path)) {
    stop("Counts file not found for study_id=", study_id, " at ", counts_path)
  }

  counts_rds_path <- sub("\\.(tsv|txt|csv)$", ".rds", counts_path, ignore.case = TRUE)
  if (identical(counts_rds_path, counts_path)) {
    counts_rds_path <- paste0(counts_path, ".rds")
  }

  if (file.exists(counts_rds_path) && !overwrite) {
    stop("RDS file already exists at ", counts_rds_path, ". Set overwrite = TRUE to regenerate it.")
  }

  message("Reading counts table from ", counts_path)
  t <- system.time({
    mm <- read.table(counts_path, sep = "\t", header = TRUE, check.names = FALSE)
  })
  message(sprintf("Finished reading counts table in %.3f s", t[["elapsed"]]))

  saveRDS(mm, counts_rds_path)
  message("Saved counts matrix as RDS to ", counts_rds_path)
  invisible(counts_rds_path)
}

