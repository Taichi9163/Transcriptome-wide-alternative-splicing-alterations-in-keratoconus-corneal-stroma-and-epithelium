#!/usr/bin/env Rscript
# Export significant rMATS events, gene memberships and event-class counts.
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "utils_paths.R"))

read_das_events <- function(path) {
  events <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  # rMATS repeats the ID header after the coordinate columns.
  names(events) <- make.unique(names(events))
  required <- c("ID", "GeneID", "geneSymbol", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference")
  missing_columns <- setdiff(required, names(events))
  if (length(missing_columns)) {
    stop("Missing rMATS columns in ", path, ": ", paste(missing_columns, collapse = ", "))
  }

  # Preserve the strict manuscript thresholds and the original PSI availability rule.
  selected <- is.finite(events$FDR) & events$FDR < 0.05 &
    is.finite(events$IncLevelDifference) & abs(events$IncLevelDifference) > 0.10 &
    !is.na(events$IncLevel1) & !is.na(events$IncLevel2)
  events[selected, , drop = FALSE]
}

save_result <- function(table, output_dir, filename) {
  write.table(table, file.path(output_dir, filename),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

args <- script_arguments()
cohorts <- c("stroma", "epithelium", "stroma_with_kc6")
if (length(args) != 1L || !args %in% cohorts) {
  stop("Choose one cohort: ", paste(cohorts, collapse = ", "))
}
cohort <- args[1]
event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")
input_paths <- work_path("data", cohort, "rmats", paste0(event_types, ".MATS.JC.txt"))
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) stop("Missing rMATS inputs: ", paste(missing_inputs, collapse = ", "))

significant_events <- lapply(input_paths, read_das_events)
genes_by_type <- lapply(significant_events, function(events) unique(events$GeneID))
all_genes <- unique(unlist(genes_by_type))
event_counts <- vapply(significant_events, nrow, integer(1))
summary_table <- data.frame(
  event_class = c(event_types, "All"),
  events = c(event_counts, sum(event_counts)),
  genes = c(lengths(genes_by_type), length(all_genes))
)
membership <- data.frame(gene_id = all_genes)
for (index in seq_along(event_types)) {
  membership[[event_types[index]]] <- all_genes %in% genes_by_type[[index]]
}

output_dir <- work_path("results", cohort, "das")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (index in seq_along(event_types)) {
  save_result(significant_events[[index]], output_dir, paste0(event_types[index], ".tsv"))
}
save_result(membership, output_dir, "gene_membership.tsv")
save_result(summary_table, output_dir, "summary.tsv")
print(summary_table, row.names = FALSE)
