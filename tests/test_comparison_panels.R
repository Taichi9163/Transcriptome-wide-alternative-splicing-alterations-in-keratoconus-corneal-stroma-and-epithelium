#!/usr/bin/env Rscript
# Numeric regressions for event identity, strict thresholds and plot memberships.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
repo <- dirname(dirname(normalizePath(sub("^--file=", "", arg))))
source(file.path(repo, "scripts/utils_plotting.R"))
source(file.path(repo, "scripts/utils_splicing_events.R"))
source(file.path(repo, "scripts/utils_venn_plots.R"))
expect_error <- function(expression) {
    stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
}

# Coordinate matching must use all MXE boundaries, chromosome and strand;
# a local numeric ID or changed gene annotation must not change an event key.
one <- as.data.frame(as.list(setNames(c("chr1", "+", "10", "20", "30", "40", "1", "5", "50", "60"), RMATS_COORDS$MXE)), check.names = FALSE)
key <- rmats_keys(one, "MXE")
same <- one; same$chr <- "1"; same$ID <- 999; same$geneSymbol <- "different annotation"
stopifnot(identical(key, rmats_keys(same, "MXE")))
for (column in RMATS_COORDS$MXE) {
    changed <- one
    changed[[column]] <- if (column == "chr") "chr2" else if (column == "strand") "-" else as.character(as.numeric(one[[column]]) + 1)
    stopifnot(key != rmats_keys(changed, "MXE"))
}
expect_error(rmats_keys(rbind(one, one), "MXE"))
bad <- one; bad$upstreamES <- NA; expect_error(rmats_keys(bad, "MXE"))

# Both significance cutoffs are strict, and missing effects are never zero-filled.
stopifnot(identical(is_das(data.frame(FDR = c(.049, .05, .049, NA, .01),
    IncLevelDifference = c(.101, .2, .10, .2, NA))), c(TRUE, FALSE, FALSE, FALSE, FALSE)))

# Eight possible nonempty memberships: four single sets and four cross-layer pairs.
data <- data.frame(event_key = letters[1:8], gene_stroma = paste0("S", 1:8), gene_epithelium = paste0("E", 1:8),
    significant_stroma = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
    significant_epithelium = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    delta_stroma = c(.2, -.2, NA, NA, .2, .2, -.2, -.2),
    delta_epithelium = c(NA, NA, .2, -.2, .2, -.2, .2, -.2))
data <- four_set_memberships(data)
stopifnot(identical(data$membership, c("1000", "0100", "0010", "0001", "1010", "1001", "0110", "0101")))
regions <- four_set_regions(data)
stopifnot(nrow(regions) == 15L, sum(regions$count) == 8L,
          all(regions$count[!regions$membership %in% data$membership] == 0L))
device_before <- grDevices::dev.cur()
plot <- make_four_venn(data, "Synthetic test")
stopifnot(identical(grDevices::dev.cur(), device_before), inherits(plot, "grob"))
# VennDiagram's 15 numeric labels must agree with the independently specified
# eight occupied and seven empty regions; titles and category labels are ignored.
labels <- vapply(plot$children[[2]]$children, function(x) if (inherits(x, "text")) as.character(x$label) else "", character(1))
labels <- as.integer(labels[grepl("^[0-9]+$", labels)])
stopifnot(identical(sort(labels), c(rep(0L, 7), rep(1L, 8))))

# Junction support requires finite PSI AND >=10 reads in the same specimens.
# Three of five cases and two of four controls are required, not a pooled count.
original_root <- ROOT
ROOT <- tempfile("synthetic_support_"); dir.create(file.path(ROOT, "data/stroma/rmats"), recursive = TRUE)
raw <- data.frame(ID = 1:3, GeneID = paste0("G", 1:3), geneSymbol = paste0("G", 1:3),
    chr = "chr1", strand = "+", exonStart_0base = c(100, 200, 300), exonEnd = c(150, 250, 350),
    upstreamES = 1, upstreamEE = 50, downstreamES = 400, downstreamEE = 500,
    FDR = .01, IncLevelDifference = .2,
    IncLevel1 = c(".5,.5,.5,.5,.5", ".5,.5,.5,.5,.5", ".5,.5,NA,.5,.5"), IncLevel2 = ".3,.3,.3,.3",
    IJC_SAMPLE_1 = c("10,10,10,0,0", "10,10,9,0,0", "10,10,10,0,0"), SJC_SAMPLE_1 = "0,0,0,0,0",
    IJC_SAMPLE_2 = "10,10,0,0", SJC_SAMPLE_2 = "0,0,0,0")
write_tsv(raw, file.path(ROOT, "data/stroma/rmats/SE.MATS.JC.txt"))
support <- se_support("stroma")
stopifnot(identical(support$supported_stroma, c(TRUE, FALSE, FALSE)),
          identical(support$group1_n_10_reads_stroma, c(3L, 2L, 2L)))
ROOT <- original_root
cat("Passed: full event keys, strict cutoffs, all directional memberships, Venn labels, no graphics device, and per-group support rules.\n")
