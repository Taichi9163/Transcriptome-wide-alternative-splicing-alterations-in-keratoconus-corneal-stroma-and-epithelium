#!/usr/bin/env Rscript
# Candidate RBP expression, using gene identifiers from Ensembl annotation.  Paths: --work-dir DIR,
# then KC_WORK_DIR, then the current directory.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "utils_paths.R"))

args <- script_arguments()
if (length(args) != 2L || !args[1] %in% c("stroma", "epithelium")) {
  stop("Usage: Rscript scripts/supp_figure_03ab_rbp_expression.R stroma|epithelium gene_symbols.tsv")
}
root <- work_path("results", args[1])
map <- unique(read.delim(resolve_path(args[2]))[c("gene_id", "gene_symbol")])
motif <- read.csv(work_path("results/motif_tests/all_motif_region_results.csv"))
selected <- subset(motif, tissue == args[1] & original_BH < 5e-05)
symbols <- sort(unique(selected$rbp_label))
tpm <- read.delim(file.path(root, "gene/tpm.tsv"), check.names = FALSE)
de <- read.delim(file.path(root, "gene/wald.tsv"))
map <- map[map$gene_symbol %in% symbols & map$gene_id %in% tpm$feature_id, ]
if (anyDuplicated(map$gene_symbol) || !setequal(symbols, map$gene_symbol)) {
  stop("Supply one expressed Ensembl gene ID per candidate symbol from the study annotation.")
}
map <- map[match(symbols, map$gene_symbol), ]
abundance <- as.matrix(tpm[match(map$gene_id, tpm$feature_id), -1, drop = FALSE])
rownames(abundance) <- symbols
z <- t(scale(t(log10(1 + abundance))))
map$padj <- de$padj[match(map$gene_id, de$feature_id)]
out <- file.path(root, "rbp")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(cbind(map, z), file.path(out, "rbp_expression.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
