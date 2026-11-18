# Sample-level inclusion for the three stromal RT-PCR targets in Figure 2D-F.
x <- read.delim("results/stroma/das/SE.tsv", check.names = FALSE)
targets <- c(ZEB2 = 144424795L, ROCK2 = 11222082L, TPM2 = 35684731L)
rows <- lapply(names(targets), function(gene) {
  r <- subset(x, geneSymbol == gene & exonStart_0base == targets[[gene]])
  stopifnot(nrow(r) == 1L)
  psi <- function(v) as.numeric(strsplit(v, ",", fixed = TRUE)[[1]])
  kc <- psi(r$IncLevel1)
  controls <- psi(r$IncLevel2)
  stopifnot(length(kc) == 5L, length(controls) == 4L)
  data.frame(gene, chr = r$chr, exonStart_0base = r$exonStart_0base, exonEnd = r$exonEnd,
             group = rep(c("Control", "Keratoconus"), c(4, 5)), position = c(1:4, 1:5),
             PSI = c(controls, kc), delta = r$IncLevelDifference, FDR = r$FDR)
})
write.table(do.call(rbind, rows), "results/stroma/selected_exons.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Export the same complete rMATS rows for the three manuscript sashimi panels.
selected <- do.call(rbind, lapply(names(targets), function(gene) {
  subset(x, geneSymbol == gene & exonStart_0base == targets[[gene]])
}))
# Restore the duplicate ID header expected by native rMATS consumers.
names(selected)[names(selected) == "ID.1"] <- "ID"
write.table(selected, "results/stroma/sashimi_events.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
