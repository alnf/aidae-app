#!/usr/bin/env Rscript
# Generate a scrypt password hash for shinymanager (`auth:` in config / deploy YAML).
# Usage: Rscript scripts/hash_password.R [password]
# Default password: "test"

args <- commandArgs(trailingOnly = TRUE)
password <- if (length(args) >= 1) args[1] else "test"

if (!requireNamespace("scrypt", quietly = TRUE)) {
  stop("Install the 'scrypt' package: install.packages(\"scrypt\")")
}

hash <- scrypt::hashPassword(password)
cat(hash, "\n")
