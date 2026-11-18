#!/usr/bin/env Rscript
# Work-directory precedence, explicit path overrides and shared PSI parsing.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- dirname(dirname(normalizePath(script)))
source(file.path(repo, "scripts/utils_paths.R"))
source(file.path(repo, "scripts/utils_psi.R"))
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

expected <- matrix(c(0, 1, 0.5, NA_real_), nrow = 2, byrow = TRUE)
stopifnot(identical(parse_psi(c("0,1", "0.5,NA"), 2), expected))
for (invalid in list("1,2", "0,Inf", "0,word", "0,1,0")) {
  expect_error(parse_psi(invalid, 2))
}
message("Portable R paths and PSI parsing checks passed.")
