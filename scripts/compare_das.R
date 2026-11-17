# Coordinate-matched events for all five classes, stratified by direction.
out <- "results/cross_layer"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (event in c("SE", "A5SS", "A3SS", "MXE", "RI")) {
  read_calls <- function(tissue) {
    tab <- read.delim(file.path("results", tissue, "das", paste0(event, ".tsv")), check.names = FALSE)
    coords <- names(tab)[match("chr", names(tab)):(match("ID.1", names(tab)) - 1L)]
    stopifnot(length(coords) >= 8L, !anyNA(tab[coords]), !anyDuplicated(tab[coords]))
    ans <- tab[coords]
    ans[[tissue]] <- ifelse(tab$IncLevelDifference > 0, "inclusion", "skipping")
    ans
  }
  stroma <- read_calls("stroma")
  epithelium <- read_calls("epithelium")
  merged <- merge(stroma, epithelium, all = TRUE, by = intersect(names(stroma), names(epithelium)))
  write.table(merged, file.path(out, paste0(event, "_union.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  print(event)
  print(table(merged$stroma, merged$epithelium, useNA = "ifany"))
}
