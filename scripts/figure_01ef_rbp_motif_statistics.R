#!/usr/bin/env Rscript
# Reproduce the RBP motif multiplicity analysis from the four saved rMAPS3 tables.  Run with
# Rscript. No packages outside base R or original figure files are needed.  Run settings come from
# each tissue's original log. See README.md for the rMAPS3 code defining the number of tested
# positions and the eight output columns.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "utils_paths.R"))

out <- work_path("results/motif_tests")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cohorts <- c("stroma", "epithelium")
regions <- c(
  "upstreamExon-3prime", "upstreamExonIntron", "upstreamIntron", "targetExon-5prime", "targetExon-3prime",
  "downstreamIntron", "downstreamExonIntron", "downstreamExon-5prime"
)
exonic <- c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE)

read_setting <- function(command, flag) {
  tokens <- strsplit(sub("\\].*$", "", command), "[[:space:]]+")[[1]]
  pos <- match(flag, tokens)
  stopifnot(!is.na(pos))
  value <- as.integer(tokens[pos + 1L])
  stopifnot(is.finite(value), value > 0L)
  value
}

read_family <- function(tissue, direction) {
  src <- work_path("data", tissue, "rmaps")
  # The supplementary workbook records minima from the study's fixed scan.  Validate these settings
  # against the run log when using native rMAPS3 output.
  n_window <- ifelse(exonic, 50L, 250L)
  log_path <- file.path(src, "log.motifMap.txt")
  if (file.exists(log_path)) {
    command <- readLines(log_path, n = 1L)
    stopifnot(read_setting(command, "--step") == 1L, read_setting(command, "--window") == 50L, read_setting(
      command,
      "--exon"
    ) == 50L, read_setting(command, "--intron") == 250L)
  }

  path <- file.path(src, paste0("pVal.", direction, ".vs.bg.RNAmap.txt"))
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(
    identical(names(tab), c("RBP", paste0("smallest_p_in_", regions))), nrow(tab) == 114L,
    !anyDuplicated(tab$RBP)
  )
  mat <- as.matrix(tab[, -1L])
  stopifnot(all(is.finite(mat)), all(mat >= 0 & mat <= 1))

  # as.vector() reads columns, so each region and window count repeats 114 times.
  res <- data.frame(tissue,
    direction = if (direction == "up")
      "inclusion" else "skipping", motif = rep(tab$RBP, length(regions)), region = rep(regions, each = nrow(tab)),
    windows = rep(n_window, each = nrow(tab)), min_window_p = as.vector(mat)
  )
  res$original_BH <- p.adjust(res$min_window_p, "BH")
  res$region_bonferroni_p <- pmin(1, res$windows * res$min_window_p)
  res$window_then_BH <- p.adjust(res$region_bonferroni_p, "BH")
  res$window_then_BY <- p.adjust(res$region_bonferroni_p, "BY")
  res
}

families <- expand.grid(tissue = cohorts, direction = c("up", "dn"), stringsAsFactors = FALSE)
motif_results <- do.call(rbind, Map(read_family, families$tissue, families$direction))
stopifnot(nrow(motif_results) == 3648L, all(motif_results$window_then_BH >= motif_results$original_BH -
  1e-14), all(motif_results$window_then_BY >= motif_results$window_then_BH - 1e-14))

# Map motif labels to the gene symbols used in the expression panels.
display_cutoff <- 5e-05
motif_results$rbp_label <- sub("\\..*$", "", motif_results$motif)
aliases <- c(HuR = "ELAVL1", `SF2-ASF` = "SRSF1")
hit_labels <- unique(motif_results[motif_results$original_BH < display_cutoff, c("tissue", "rbp_label")])
selected <- paste(motif_results$tissue, motif_results$rbp_label) %in% paste(hit_labels$tissue, hit_labels$rbp_label)
matched <- selected & motif_results$rbp_label %in% names(aliases)
motif_results$rbp_label[matched] <- unname(aliases[motif_results$rbp_label[matched]])
motif_results <- motif_results[order(
  motif_results$tissue, motif_results$direction, motif_results$window_then_BH,
  motif_results$motif, motif_results$region
), ]

summary_rows <- list()
for (method in c("original_BH", "window_then_BH", "window_then_BY")) {
  for (cutoff in c(display_cutoff, 0.05)) {
    for (family in split(motif_results, list(motif_results$tissue, motif_results$direction), drop = TRUE)) {
      hits <- family[family[[method]] < cutoff, ]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        tissue = family$tissue[1], direction = family$direction[1],
        method, cutoff, motif_region_pairs = nrow(hits), distinct_motifs = length(unique(hits$motif))
      )
    }
  }
}
write.csv(motif_results, file.path(out, "all_motif_region_results.csv"), row.names = FALSE)
write.csv(motif_results[motif_results$window_then_BH < display_cutoff, ], file.path(out, "hits_at_original_display_cutoff.csv"),
  row.names = FALSE
)
write.csv(do.call(rbind, summary_rows), file.path(out, "summary.csv"), row.names = FALSE)
