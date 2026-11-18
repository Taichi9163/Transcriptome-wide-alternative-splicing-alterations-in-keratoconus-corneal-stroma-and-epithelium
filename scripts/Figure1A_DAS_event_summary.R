# Extract differential alternative splicing events from rMATS junction-only output.
args <- commandArgs(trailingOnly = TRUE)
cohorts <- c("stroma", "epithelium", "stroma_with_kc6")
if (length(args) != 1L || !args %in% cohorts) stop("Choose one cohort: ", paste(cohorts, collapse = ", "))
out <- file.path("results", args, "das")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
classes <- c("SE", "A5SS", "A3SS", "MXE", "RI")
input <- file.path("data", args, "rmats", paste0(classes, ".MATS.JC.txt"))
stopifnot(all(file.exists(input)))
hits <- lapply(input, function(path) {
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  # rMATS repeats the ID header after the coordinate columns.
  names(tab) <- make.unique(names(tab))
  stopifnot(all(c("ID", "geneSymbol", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference") %in% names(tab)))
  keep <- is.finite(tab$FDR) & tab$FDR < 0.05 & is.finite(tab$IncLevelDifference) &
    abs(tab$IncLevelDifference) > 0.10 & !is.na(tab$IncLevel1) & !is.na(tab$IncLevel2)
  tab[keep, , drop = FALSE]
})
gene_ids <- lapply(hits, function(x) unique(x$GeneID))
all_genes <- unique(unlist(gene_ids))
summary <- data.frame(event_class = c(classes, "All"),
  events = c(vapply(hits, nrow, integer(1)), sum(vapply(hits, nrow, integer(1)))),
  genes = c(lengths(gene_ids), length(all_genes)))
save_table <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
for (i in seq_along(classes)) save_table(hits[[i]], paste0(classes[i], ".tsv"))
membership <- data.frame(gene_id = all_genes)
for (i in seq_along(classes)) membership[[classes[i]]] <- all_genes %in% gene_ids[[i]]
save_table(membership, "gene_membership.tsv")
save_table(summary, "summary.tsv")
print(summary, row.names = FALSE)
