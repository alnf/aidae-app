perf_time <- function(label, expr) {
  t <- system.time(res <- force(expr))
  message(sprintf("[perf] %s: elapsed=%.3f s", label, t[["elapsed"]]))
  res
}

