#!/usr/bin/env Rscript
# Epithelial Reactome top ten, adapted from the recovered enrichment R script.  --tables-only
# validates numeric data and plot layers without creating images.  Paths: --work-dir DIR, then
# KC_WORK_DIR, then the current directory.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "utils_plotting.R"))
source(file.path(script_dir, "utils_reactome_enrichment.R"))
run_reactome_panel("epithelium", "SuppFigure10B", "suppfigure10")
