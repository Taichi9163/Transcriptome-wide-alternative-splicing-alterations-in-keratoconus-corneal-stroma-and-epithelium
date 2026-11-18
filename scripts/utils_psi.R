# Parse ordered rMATS PSI vectors without imputing missing values.
parse_psi <- function(values, n_expected, context = "PSI vector") {
  tokens <- strsplit(values, ",", fixed = TRUE)
  if (any(lengths(tokens) != n_expected))
    stop("Unexpected PSI vector length: ", context)
  flat <- trimws(unlist(tokens, use.names = FALSE))
  parsed <- suppressWarnings(as.numeric(flat))
  missing_token <- is.na(flat) | flat %in% c("NA", "NaN", "nan", "")
  if (any(is.na(parsed) & !missing_token))
    stop("Non-numeric PSI token: ", context)
  if (any(is.infinite(parsed)))
    stop("Infinite PSI: ", context)
  if (any(parsed < 0 | parsed > 1, na.rm = TRUE))
    stop("PSI outside [0, 1]: ", context)
  matrix(parsed, ncol = n_expected, byrow = TRUE)
}
