#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
repo <- dirname(dirname(normalizePath(sub("^--file=", "", arg))))
source(file.path(repo, "tests/no_graphics.R"))
source(file.path(repo, "scripts/Utils_Reactome_enrichment.R"))
expect_error <- function(expression) {
    failed <- inherits(try(force(expression), silent = TRUE), "try-error")
    stopifnot(failed)
}
expected <- list(
    stroma = c("R-HSA-983168", "R-HSA-983169", "R-HSA-1280218", "R-HSA-5696399", "R-HSA-69239",
               "R-HSA-69306", "R-HSA-390522", "R-HSA-69242", "R-HSA-2467813", "R-HSA-73893"),
    epithelium = c("R-HSA-111465", "R-HSA-75153", "R-HSA-5362517", "R-HSA-1169092", "R-HSA-8873719",
                   "R-HSA-157579", "R-HSA-211976", "R-HSA-5653656", "R-HSA-5658442", "R-HSA-204174"))
for (cohort in names(expected)) {
    result <- reactome_top10(file.path(repo, "inputs", cohort, "nease/Reactome_Python.tsv.gz"))
    stopifnot(result$summary$tested_pathways == 1823L,
              identical(result$top$`Pathway ID`, expected[[cohort]]))
    layer <- ggplot2::ggplot_build(reactome_plot(result$top, cohort))$data[[1]]
    stopifnot(nrow(layer) == 10L, all(is.finite(layer$size)), length(unique(layer$y)) == 10L)
    if (cohort == "epithelium") stopifnot(result$top$`Splicing-affected genes`[8] == 14L)
}
# Synthetic boundary/tie/duplicate-name cases are independent of the study data.
frame <- data.frame("Pathway ID" = sprintf("R-TEST-%02d", 10:1), "Pathway name" = "same name",
    Source = "Reactome", "Spliced genes (number of interactions affecting the pathway)" = "A (2), A (3), B (1), ",
    p_value = .05, "adj p_value" = .05, "Nease score" = 1, check.names = FALSE)
path <- tempfile(fileext = ".tsv")
write_fixture <- function(x) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
write_fixture(frame)
top <- reactome_top10(path)$top
stopifnot(identical(top$`Pathway ID`, sort(frame$`Pathway ID`)), all(top$`Splicing-affected genes` == 2L))
layer <- ggplot2::ggplot_build(reactome_plot(top, "stroma"))$data[[1]]
stopifnot(length(unique(layer$y)) == 10L)
bad <- frame; bad$`Pathway ID`[2] <- bad$`Pathway ID`[1]; write_fixture(bad)
expect_error(reactome_top10(path))
bad <- frame; bad$`adj p_value`[1] <- .02; write_fixture(bad)
expect_error(reactome_top10(path))
bad <- frame; bad$Source[1] <- "KEGG"; write_fixture(bad)
expect_error(reactome_top10(path))
bad <- frame; bad$`Nease score`[1] <- NA; write_fixture(bad)
expect_error(reactome_top10(path))
stopifnot(grDevices::dev.cur() == 1L)
message("Reactome study and edge-case checks passed without graphics devices.")
