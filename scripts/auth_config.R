# Shinymanager credentials from merged main YAML (`auth:`) or data/creds.txt fallback.

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_auth_users <- function(auth) {
  if (is.null(auth)) {
    return(NULL)
  }
  users <- if (is.list(auth) && !is.null(auth$users)) auth$users else auth
  if (length(users) < 1L || !is.list(users)) {
    stop("main config `auth` must be a non-empty list of user entries", call. = FALSE)
  }
  lapply(users, function(u) {
    if (!is.list(u)) {
      stop("each auth entry must be a mapping with `user` and `password`", call. = FALSE)
    }
    u_user <- as.character(u$user %||% "")
    u_pass <- as.character(u$password %||% "")
    if (!nzchar(u_user)) {
      stop("auth entry missing `user`", call. = FALSE)
    }
    if (!nzchar(u_pass)) {
      stop("auth entry missing `password` for user ", u_user, call. = FALSE)
    }
    list(
      user = u_user,
      password = u_pass,
      start = as.character(u$start %||% "2025-09-11"),
      expire = u$expire %||% NA,
      admin = isTRUE(u$admin)
    )
  })
}

auth_users_to_credentials <- function(users) {
  data.frame(
    user = vapply(users, `[[`, "", "user"),
    password = vapply(users, `[[`, "", "password"),
    start = vapply(users, `[[`, "", "start"),
    expire = vapply(users, function(x) x$expire, NA),
    admin = vapply(users, `[[`, FALSE, "admin"),
    stringsAsFactors = FALSE,
    is_hashed_password = TRUE
  )
}

read_credentials_tsv <- function(creds_txt_path) {
  if (!file.exists(creds_txt_path)) {
    return(NULL)
  }
  creds <- read.table(
    creds_txt_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  if (nrow(creds) < 1L) {
    return(NULL)
  }
  keep <- !grepl("^#", creds$user)
  creds <- creds[keep, , drop = FALSE]
  if (nrow(creds) < 1L) {
    return(NULL)
  }
  auth_users_to_credentials(lapply(seq_len(nrow(creds)), function(i) {
    list(
      user = creds$user[[i]],
      password = creds$password[[i]],
      start = "2025-09-11",
      expire = NA,
      admin = FALSE
    )
  }))
}

#' Build shinymanager credentials from merged main config, with creds.txt fallback.
credentials_from_main_config <- function(main_config, creds_txt_path = "data/creds.txt") {
  users <- normalize_auth_users(main_config$auth)
  if (!is.null(users)) {
    return(auth_users_to_credentials(users))
  }
  creds <- read_credentials_tsv(creds_txt_path)
  if (!is.null(creds)) {
    message(
      "[exprs-app] auth: using ", creds_txt_path,
      " (", nrow(creds), " user(s); prefer `auth` in main config for deploy apps)"
    )
    return(creds)
  }
  stop(
    "No login credentials configured: add `auth` to the main config YAML ",
    "or provide ", creds_txt_path,
    call. = FALSE
  )
}
