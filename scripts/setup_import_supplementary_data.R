#!/usr/bin/env Rscript
# Convert one published supplementary workbook to the input layout used by these scripts.  Create
# data/ in a fresh working directory or keep the original rMATS inputs elsewhere.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "utils_paths.R"))

library(readxl)
args <- script_arguments()
if (length(args) != 2L || !args[1] %in% c("stroma", "epithelium", "motifs")) {
  stop("Usage: Rscript scripts/setup_import_supplementary_data.R stroma|epithelium|motifs workbook.xlsx")
}
write_input <- function(x, tissue, module, filename) {
  dest <- work_path("data", tissue, module)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dest, filename)
  if (file.exists(path))
    stop("Input already exists: ", path)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
if (args[1] == "motifs") {
  for (tissue in c("Stroma", "Epithelium")) for (direction in c("up", "down")) {
    x <- read_excel(resolve_path(args[2]), sheet = paste(tissue, direction), skip = 3, .name_repair = "minimal")
    cols <- c("RBP", grep("^smallest_p_in_", names(x), value = TRUE))
    stopifnot(length(cols) == 9L, nrow(x) == 114L)
    label <- if (direction == "down")
      "dn" else "up"
    write_input(x[cols], tolower(tissue), "rmaps", paste0("pVal.", label, ".vs.bg.RNAmap.txt"))
  }
} else {
  events <- c("SE", "A5SS", "A3SS", "MXE", "RI")
  cohorts <- if (args[1] == "stroma")
    c("stroma", "stroma_with_kc6") else "epithelium"
  for (tissue in cohorts) for (event in events) {
    sheet <- if (tissue == "stroma_with_kc6")
      paste(event, "(all samples)") else event
    x <- read_excel(resolve_path(args[2]), sheet = sheet, skip = 1, col_types = "text", .name_repair = "minimal")
    nm <- sub(" \\(.*$", "", names(x))
    nm[nm == "GeneSymbol"] <- "geneSymbol"
    nm[nm == "P-value"] <- "PValue"
    stopifnot(sum(nm == "ID") == 2L, all(c("GeneID", "geneSymbol", "IncLevelDifference", "FDR") %in%
      nm))
    names(x) <- nm
    write_input(x, tissue, "rmats", paste0(event, ".MATS.JC.txt"))
  }
}
