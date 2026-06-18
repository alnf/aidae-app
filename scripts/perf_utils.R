perf_time <- function(label, expr) {
  t <- system.time(res <- force(expr))
  message(sprintf("[perf] %s: elapsed=%.3f s", label, t[["elapsed"]]))
  res
}

# ORA tab profiling (enable with EXPRS_ORA_PROFILE=1).
.ora_profile_ctx <- new.env(parent = emptyenv())

ora_profile_enabled <- function() {
  identical(Sys.getenv("EXPRS_ORA_PROFILE", unset = ""), "1")
}

ora_profile_set_ctx <- function(pathway_file) {
  if (!ora_profile_enabled()) {
    return(invisible(NULL))
  }
  .ora_profile_ctx$pathway_file <- as.character(pathway_file)
  invisible(NULL)
}

ora_profile_ctx_pathway <- function() {
  pf <- .ora_profile_ctx$pathway_file
  if (is.null(pf) || !nzchar(pf)) "?" else pf
}

ora_profile_timer_start <- function() {
  if (!ora_profile_enabled()) {
    return(NULL)
  }
  proc.time()
}

ora_profile_elapsed_ms <- function(t0) {
  if (is.null(t0)) {
    return(NA_real_)
  }
  as.numeric((proc.time() - t0)["elapsed"]) * 1000
}

ora_profile_log <- function(...) {
  if (!ora_profile_enabled()) {
    return(invisible(NULL))
  }
  message("[ORA profile] ", paste0(..., collapse = ""))
}

