# Exact event identity and directional memberships shared by comparison panels.  Source
# utils_plotting.R first to define ROOT and output helpers.
RMATS_COORDS <- list(
  SE = c(
    "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  ), A3SS = c(
    "chr", "strand", "longExonStart_0base", "longExonEnd",
    "shortES", "shortEE", "flankingES", "flankingEE"
  ), A5SS = c(
    "chr", "strand", "longExonStart_0base",
    "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"
  ), MXE = c(
    "chr", "strand", "1stExonStart_0base",
    "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"
  ),
  RI = c(
    "chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES",
    "downstreamEE"
  )
)
SET_NAMES <- c("Stromal inclusion", "Stromal skipping", "Epithelial inclusion", "Epithelial skipping")

rmats_keys <- function(data, event) {
  columns <- RMATS_COORDS[[event]]
  if (is.null(columns) || !all(columns %in% names(data)))
    stop("Missing coordinates: ", event)
  if (!nrow(data))
    return(character())
  values <- lapply(columns, function(column) {
    value <- trimws(as.character(data[[column]]))
    if (anyNA(value) || any(!nzchar(value)))
      stop("Missing coordinate: ", column)
    if (column == "chr") {
      value <- paste0("chr", sub("^chr", "", value))
    } else if (column == "strand") {
      if (any(!value %in% c("+", "-")))
        stop("Invalid strand")
    } else {
      number <- suppressWarnings(as.numeric(value))
      if (any(!is.finite(number) | number < 0 | number != floor(number)))
        stop("Invalid boundary: ", column)
      value <- format(number, scientific = FALSE, trim = TRUE)
    }
    value
  })
  key <- do.call(paste, c(values, sep = ":"))
  if (anyDuplicated(key))
    stop("Duplicate ", event, " coordinates within one analysis")
  key
}

read_rmats <- function(cohort, event) {
  path <- file.path(ROOT, "data", cohort, "rmats", paste0(event, ".MATS.JC.txt"))
  data <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(data) <- sub("^X([12](st|nd)Exon)", "\\1", names(data))
  required <- c("FDR", "IncLevelDifference", "GeneID", "geneSymbol")
  if (!all(required %in% names(data)))
    stop("Missing rMATS fields: ", path)
  for (column in c("FDR", "IncLevelDifference")) {
    if (!is.numeric(data[[column]]))
      stop("Non-numeric ", column, ": ", path)
  }
  if (any(!is.na(data$FDR) & (!is.finite(data$FDR) | data$FDR < 0 | data$FDR > 1)))
    stop("Invalid FDR")
  if (any(!is.na(data$IncLevelDifference) & (!is.finite(data$IncLevelDifference) | abs(data$IncLevelDifference) >
    1)))
    stop("Invalid delta-PSI")
  data$event_key <- rmats_keys(data, event)
  data
}

is_das <- function(data) {
  is.finite(data$FDR) & data$FDR < 0.05 & is.finite(data$IncLevelDifference) & abs(data$IncLevelDifference) >
    0.1
}

cross_layer_events <- function(event) {
  stroma <- read_rmats("stroma", event)
  epithelium <- read_rmats("epithelium", event)
  sig_s <- is_das(stroma)
  sig_e <- is_das(epithelium)
  keys <- union(stroma$event_key[sig_s], epithelium$event_key[sig_e])
  i <- match(keys, stroma$event_key)
  j <- match(keys, epithelium$event_key)
  data.frame(
    event_key = keys, gene_stroma = stroma$geneSymbol[i], gene_epithelium = epithelium$geneSymbol[j],
    gene_id_stroma = stroma$GeneID[i], gene_id_epithelium = epithelium$GeneID[j], delta_stroma = stroma$IncLevelDifference[i],
    delta_epithelium = epithelium$IncLevelDifference[j], FDR_stroma = stroma$FDR[i], FDR_epithelium = epithelium$FDR[j],
    significant_stroma = keys %in% stroma$event_key[sig_s], significant_epithelium = keys %in% epithelium$event_key[sig_e]
  )
}

four_set_memberships <- function(data) {
  # Bit order is fixed in every output and plot. Opposite directions within a cohort are disjoint;
  # cross-layer overlaps count events, not unique genes.
  flags <- cbind(data$significant_stroma & data$delta_stroma > 0, data$significant_stroma & data$delta_stroma <
    0, data$significant_epithelium & data$delta_epithelium > 0, data$significant_epithelium & data$delta_epithelium <
    0)
  if (anyNA(flags))
    stop("Missing effect for a significant event")
  data$membership <- apply(flags, 1, function(row) paste(as.integer(row), collapse = ""))
  if (any(data$membership == "0000"))
    stop("Venn input must be significant in at least one cohort")
  data
}

four_set_regions <- function(data) {
  masks <- vapply(1:15, function(i) paste(as.integer(intToBits(i))[4:1], collapse = ""), character(1))
  do.call(rbind, lapply(masks, function(mask) {
    members <- data[data$membership == mask, , drop = FALSE]
    genes <- c(members$gene_stroma, members$gene_epithelium)
    genes <- sort(unique(genes[!is.na(genes) & nzchar(genes) & genes != "NA"]))
    data.frame(membership = mask, count = nrow(members), sets = paste(SET_NAMES[strsplit(mask, "")[[1]] ==
      "1"], collapse = "; "), genes = paste(genes, collapse = ", "))
  }))
}

se_support <- function(cohort) {
  raw <- read_rmats(cohort, "SE")
  data <- data.frame(
    event_key = raw$event_key, gene = raw$geneSymbol, gene_id = raw$GeneID, delta = raw$IncLevelDifference,
    FDR = raw$FDR
  )
  data$estimable <- is.finite(data$delta)
  data$supported <- TRUE
  means <- list()
  parse_vector <- function(x) suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
  for (group in 1:2) {
    fields <- c(paste0("IncLevel", group), paste0("IJC_SAMPLE_", group), paste0("SJC_SAMPLE_", group))
    if (!all(fields %in% names(raw)))
      stop("Missing sample PSI or junction counts")
    psi <- lapply(raw[[fields[1]]], parse_vector)
    ijc <- lapply(raw[[fields[2]]], parse_vector)
    sjc <- lapply(raw[[fields[3]]], parse_vector)
    expected <- if (cohort == "epithelium")
      7L else if (group == 1)
      5L else 4L
    if (any(lengths(psi) != expected | lengths(ijc) != expected | lengths(sjc) != expected))
      stop("Wrong sample count")
    finite <- vapply(psi, function(x) sum(is.finite(x)), integer(1))
    support <- vapply(seq_along(psi), function(i) sum(is.finite(psi[[i]]) & is.finite(ijc[[i]] +
      sjc[[i]]) & ijc[[i]] + sjc[[i]] >= 10), integer(1))
    data[[paste0("group", group, "_n_PSI")]] <- finite
    data[[paste0("group", group, "_n_10_reads")]] <- support
    data$estimable <- data$estimable & finite >= 1
    data$supported <- data$supported & support >= ceiling(expected / 2)
    means[[group]] <- vapply(psi, function(x) mean(x[is.finite(x)]), numeric(1))
  }
  if (any(abs(data$delta - (means[[1]] - means[[2]]))[data$estimable] > 0.002))
    stop("Inconsistent PSI direction")
  data$significant <- data$estimable & is_das(raw)
  names(data)[names(data) != "event_key"] <- paste0(names(data)[names(data) != "event_key"], "_", cohort)
  data
}
