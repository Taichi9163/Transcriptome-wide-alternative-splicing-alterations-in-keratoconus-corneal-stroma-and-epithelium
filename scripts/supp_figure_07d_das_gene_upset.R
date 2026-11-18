#!/usr/bin/env Rscript
# SuppFigure7D: gene intersections among the five significant DAS event classes.
# KC_WORK_DIR selects the work directory. --tables-only avoids creating images.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "utils_plotting.R"))
source(file.path(script_dir, "utils_das_gene_upset.R"))
run_das_gene_upset("epithelium")
