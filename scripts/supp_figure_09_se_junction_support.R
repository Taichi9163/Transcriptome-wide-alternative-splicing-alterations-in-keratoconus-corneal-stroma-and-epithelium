# Read-only audit of source rMATS tables; writes new tables, never source data.
# Correlations are descriptive comparisons of independent cohorts, not tests
# of a disease-by-tissue interaction. No additional significance testing.
out <- "results/cross_layer/support"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
coords <- c("chr", "strand", "exonStart_0base", "exonEnd",
            "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
paths <- c(stroma = "data/stroma/rmats/SE.MATS.JC.txt",
           epithelium = "data/epithelium/rmats/SE.MATS.JC.txt")
numbers <- function(x) suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
read_support <- function(path, tissue) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  key <- do.call(paste, c(lapply(x[, coords], function(v) trimws(as.character(v))), sep = ":"))
  stopifnot(!anyDuplicated(key), !anyNA(x[, coords]))
  d <- data.frame(event_key = key, delta = x$IncLevelDifference, FDR = x$FDR)
  d$estimable <- is.finite(d$delta)
  d$half_PSI <- d$all_PSI <- d$half_10_reads <- TRUE
  group_means <- list()
  for (g in 1:2) {
    psi <- lapply(x[[paste0("IncLevel", g)]], numbers)
    ijc <- lapply(x[[paste0("IJC_SAMPLE_", g)]], numbers)
    sjc <- lapply(x[[paste0("SJC_SAMPLE_", g)]], numbers)
    n <- lengths(psi)
    expected <- if (tissue == "epithelium") 7L else if (g == 1) 5L else 4L
    stopifnot(all(n == expected), all(lengths(ijc) == n), all(lengths(sjc) == n))
    nf <- vapply(psi, function(v) sum(is.finite(v)), integer(1))
    n10 <- vapply(seq_along(psi), function(i) sum(is.finite(psi[[i]]) &
      is.finite(ijc[[i]] + sjc[[i]]) & ijc[[i]] + sjc[[i]] >= 10), integer(1))
    d[[paste0("group", g, "_n_PSI")]] <- nf
    d[[paste0("group", g, "_n_10_reads")]] <- n10
    d$estimable <- d$estimable & nf >= 1
    d$half_PSI <- d$half_PSI & nf >= ceiling(n / 2)
    d$all_PSI <- d$all_PSI & nf == n
    d$half_10_reads <- d$half_10_reads & n10 >= ceiling(n / 2)
    group_means[[g]] <- vapply(psi, function(v) mean(v[is.finite(v)]), numeric(1))
  }
  # rMATS reports rounded sample PSI and rounded effects; tolerate rounding.
  delta_error <- abs(d$delta - (group_means[[1]] - group_means[[2]]))
  stopifnot(all(delta_error[d$estimable] <= 0.002))
  d$significant <- d$estimable & is.finite(d$FDR) & d$FDR < 0.05 & abs(d$delta) > 0.1
  names(d)[names(d) != "event_key"] <- paste0(names(d)[names(d) != "event_key"], "_", tissue)
  d
}
s <- read_support(paths[["stroma"]], "stroma")
e <- read_support(paths[["epithelium"]], "epithelium")
paired <- merge(s, e, by = "event_key", sort = FALSE)
paired$paired_estimable <- paired$estimable_stroma & paired$estimable_epithelium
paired$DAS_union <- paired$significant_stroma | paired$significant_epithelium
paired$both_significant <- paired$significant_stroma & paired$significant_epithelium
paired$half_PSI_both <- paired$half_PSI_stroma & paired$half_PSI_epithelium
paired$all_PSI_both <- paired$all_PSI_stroma & paired$all_PSI_epithelium
paired$half_10_reads_both <- paired$half_10_reads_stroma & paired$half_10_reads_epithelium
stats <- list()
rules <- c("paired_estimable", "half_PSI_both", "all_PSI_both", "half_10_reads_both")
for (population in c("all_paired_SE", "paired_SE_DAS_union")) {
  for (rule in rules) {
    keep <- paired$paired_estimable & paired[[rule]]
    if (population == "paired_SE_DAS_union") keep <- keep & paired$DAS_union
    d <- paired[keep, ]
    x <- d$delta_stroma; y <- d$delta_epithelium
    nonzero <- x != 0 & y != 0
    stats[[length(stats) + 1]] <- data.frame(population = population, support_rule = rule,
      n_events = nrow(d), stroma_only_significant = sum(d$significant_stroma & !d$significant_epithelium),
      epithelium_only_significant = sum(d$significant_epithelium & !d$significant_stroma),
      both_significant = sum(d$both_significant), neither_significant = sum(!d$DAS_union),
      same_direction_nonzero = sum(x * y > 0), opposite_direction = sum(x * y < 0),
      measured_zero_in_either = sum(!nonzero),
      same_direction_fraction_nonzero = mean(sign(x[nonzero]) == sign(y[nonzero])),
      pearson_r = cor(x, y), spearman_rho = cor(x, y, method = "spearman"))
  }
}
stats <- do.call(rbind, stats)
write.csv(stats, file.path(out, "paired_SE_support_sensitivity.csv"), row.names = FALSE)
write.csv(paired, file.path(out, "paired_SE_support_audit.csv"), row.names = FALSE, na = "NA")

# Independently reconstruct epithelial significance calls from all raw event tables.
counts <- list()
provenance <- paths
for (event in c("SE", "RI", "A5SS", "A3SS", "MXE")) {
  raw <- file.path("data/epithelium/rmats", paste0(event, ".MATS.JC.txt"))
  saved <- file.path("results/epithelium/das", paste0(event, ".tsv"))
  x <- read.delim(raw, check.names = FALSE)
  old <- read.delim(saved, check.names = FALSE)
  group_ok <- function(v) vapply(lapply(v, numbers), function(z) any(is.finite(z)), logical(1))
  estimable <- is.finite(x$IncLevelDifference) & group_ok(x$IncLevel1) & group_ok(x$IncLevel2)
  fdr <- estimable & is.finite(x$FDR) & x$FDR < 0.05
  keep <- fdr & abs(x$IncLevelDifference) > 0.1
  # Compare original rMATS IDs within each source event class, not across cohorts.
  stopifnot(setequal(as.character(x$ID[keep]), as.character(old$ID)))
  counts[[event]] <- data.frame(event_type = event, total_raw_events = nrow(x),
    estimable_events = sum(estimable), FDR_lt_005 = sum(fdr),
    FDR_lt_005_abs_delta_gt_01 = sum(keep), saved_DAS_count = nrow(old), identical_ID_set = TRUE)
  provenance <- c(provenance, raw, saved)
}
counts <- do.call(rbind, counts)
stopifnot(sum(counts$FDR_lt_005_abs_delta_gt_01) == 805)
write.csv(counts, file.path(out, "epithelium_805_reproduction.csv"), row.names = FALSE)

# Independent fresh-table reconciliation with the revised Figure 4B export.
plot <- read.csv(file.path("results/cross_layer/results", "SE_DAS_union_delta_inc_level_scatter_data.csv"))
expected <- paired[paired$paired_estimable & paired$DAS_union, ]
stopifnot(setequal(plot$event_key, expected$event_key))
idx <- match(plot$event_key, expected$event_key)
stopifnot(all(plot$stroma_delta_inc_for_plot == expected$delta_stroma[idx]),
          all(plot$epithelium_delta_inc_for_plot == expected$delta_epithelium[idx]))
provenance <- unique(unname(provenance))
write.csv(data.frame(path = provenance, md5 = unname(tools::md5sum(provenance))),
          file.path(out, "input_file_checksums.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "audit_R_session.txt"))
print(stats, row.names = FALSE)
print(counts, row.names = FALSE)
cat("Passed: fresh-source effects and event identities match revised Figure 4B data.\n")
