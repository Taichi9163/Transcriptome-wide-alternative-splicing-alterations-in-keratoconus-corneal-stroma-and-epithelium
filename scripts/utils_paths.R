# Shared command-line arguments and paths. No working-directory changes occur here.  Priority:
# --work-dir > KC_WORK_DIR > the current working directory.
parse_work_directory <- function(args = commandArgs(trailingOnly = TRUE)) {
  work_dir <- Sys.getenv("KC_WORK_DIR", unset = getwd())
  remaining <- character()
  seen <- FALSE
  i <- 1L
  while (i <= length(args)) {
    argument <- args[i]
    if (argument == "--work-dir" || startsWith(argument, "--work-dir=")) {
      if (seen)
        stop("Specify --work-dir only once.")
      seen <- TRUE
      if (argument == "--work-dir") {
        i <- i + 1L
        if (i > length(args) || startsWith(args[i], "--")) {
          stop("--work-dir requires a directory.")
        }
        work_dir <- args[i]
      } else {
        work_dir <- sub("^--work-dir=", "", argument)
      }
      if (!nzchar(work_dir))
        stop("--work-dir cannot be empty.")
    } else {
      remaining <- c(remaining, argument)
    }
    i <- i + 1L
  }
  if (!nzchar(work_dir))
    work_dir <- getwd()
  work_dir <- path.expand(work_dir)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", work_dir))
    work_dir <- file.path(getwd(), work_dir)
  list(root = normalizePath(work_dir, winslash = "/", mustWork = FALSE), args = remaining)
}

.work_options <- parse_work_directory()
ROOT <- .work_options$root
TABLES_ONLY <- "--tables-only" %in% .work_options$args

script_arguments <- function() .work_options$args

# Resolve explicit input/output arguments against the selected work directory.  Absolute paths and
# paths beginning with ~ remain usable on any machine.
resolve_path <- function(path, root = ROOT) {
  vapply(path, function(value) {
    value <- path.expand(value)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", value))
      value <- file.path(root, value)
    value
  }, character(1), USE.NAMES = FALSE)
}

work_path <- function(...) file.path(ROOT, ...)
