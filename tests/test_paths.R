#!/usr/bin/env Rscript
# Work-directory precedence and explicit path overrides.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- dirname(dirname(normalizePath(script)))
source(file.path(repo, "scripts/utils_paths.R"))
expect_error <- function(expression) {
  stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
}

temporary <- tempfile("kc path test ")
dir.create(temporary)
environment_root <- file.path(temporary, "environment")
explicit_root <- file.path(temporary, "explicit run")
Sys.setenv(KC_WORK_DIR = environment_root)
stopifnot(parse_work_directory(character())$root == environment_root)
options <- parse_work_directory(c("stroma", "--work-dir", explicit_root, "--tables-only"))
stopifnot(options$root == explicit_root, identical(options$args, c("stroma", "--tables-only")))
stopifnot(parse_work_directory(paste0("--work-dir=", explicit_root))$root == explicit_root)
expect_error(parse_work_directory("--work-dir"))
expect_error(parse_work_directory("--work-dir="))
expect_error(parse_work_directory(c("--work-dir", "--tables-only")))
expect_error(parse_work_directory(c("--work-dir", "one", "--work-dir", "two")))
stopifnot(resolve_path("data/input.tsv", explicit_root) == file.path(explicit_root, "data/input.tsv"))
stopifnot(resolve_path(explicit_root, environment_root) == explicit_root)
stopifnot(resolve_path("~/input.tsv", explicit_root) == path.expand("~/input.tsv"))

message("Portable R path checks passed.")
