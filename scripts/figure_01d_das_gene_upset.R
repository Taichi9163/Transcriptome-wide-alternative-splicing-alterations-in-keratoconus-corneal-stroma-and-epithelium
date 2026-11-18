#!/usr/bin/env Rscript
# Figure1D: gene intersections among the five significant DAS event classes.  KC_WORK_DIR selects
# the work directory. --tables-only avoids creating images.  Paths: --work-dir DIR, then
# KC_WORK_DIR, then the current directory.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "utils_plotting.R"))
source(file.path(script_dir, "utils_das_gene_upset.R"))
run_das_gene_upset("stroma")
