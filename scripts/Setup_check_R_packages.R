#!/usr/bin/env Rscript
# Check the R packages needed for the manuscript panels and record versions.
packages <- c("ggplot2", "ggforce", "ragg", "gplots", "RColorBrewer", "dplyr", "tidyr", "purrr",
    "patchwork")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
source(file.path(dirname(normalizePath(sub("^--file=", "", arg))), "Utils_plotting.R"))
dir.create(file.path(ROOT, "results/provenance"), recursive = TRUE, showWarnings = FALSE)
versions <- data.frame(Package = packages, Version = vapply(packages, function(p) as.character(packageVersion(p)),
    character(1)))
write_tsv(versions, file.path(ROOT, "results/provenance/package_versions.tsv"))
writeLines(capture.output(sessionInfo()), file.path(ROOT, "results/provenance/sessionInfo.txt"))
print(versions, row.names = FALSE)
