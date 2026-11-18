#!/usr/bin/env Rscript
# Sample-level inclusion for the three stromal RT-PCR targets in Figure 2D-F.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "utils_paths.R"))

events <- read.delim(work_path("results/stroma/das/SE.tsv"), check.names = FALSE)
targets <- c(ZEB2 = 144424795L, ROCK2 = 11222082L, TPM2 = 35684731L)
rows <- lapply(names(targets), function(gene) {
  target_event <- subset(events, geneSymbol == gene & exonStart_0base == targets[[gene]])
  stopifnot(nrow(target_event) == 1L)
  psi <- function(v) as.numeric(strsplit(v, ",", fixed = TRUE)[[1]])
  cases <- psi(target_event$IncLevel1)
  controls <- psi(target_event$IncLevel2)
  stopifnot(length(cases) == 5L, length(controls) == 4L)
  data.frame(gene,
    chr = target_event$chr, exonStart_0base = target_event$exonStart_0base, exonEnd = target_event$exonEnd,
    group = rep(c("Control", "Keratoconus"), c(4, 5)), position = c(1:4, 1:5), PSI = c(
      controls,
      cases
    ), delta = target_event$IncLevelDifference, FDR = target_event$FDR
  )
})
write.table(do.call(rbind, rows), work_path("results/stroma/selected_exons.tsv"),
  sep = "\t", quote = FALSE,
  row.names = FALSE
)

# Export the same complete rMATS rows for the three manuscript sashimi panels.
selected <- do.call(rbind, lapply(names(targets), function(gene) {
  subset(events, geneSymbol == gene & exonStart_0base == targets[[gene]])
}))
# Restore the duplicate ID header expected by native rMATS consumers.
names(selected)[names(selected) == "ID.1"] <- "ID"
write.table(selected, work_path("results/stroma/sashimi_events.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
