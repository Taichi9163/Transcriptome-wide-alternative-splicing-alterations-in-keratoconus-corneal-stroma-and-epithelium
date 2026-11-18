#!/usr/bin/env Rscript
# Stromal Reactome top ten, adapted from the recovered enrichment R script.
# --tables-only validates numeric data and plot layers without creating images.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "Utils_plotting.R"))
source(file.path(script_dir, "Utils_Reactome_enrichment.R"))
run_reactome_panel("stroma", "Figure4D", "figure4")
