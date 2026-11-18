# Candidate RBP expression, using gene identifiers from Ensembl annotation.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[1] %in% c("stroma", "epithelium")) {
  stop("Usage: Rscript scripts/SuppFigure3AB_RBP_expression.R stroma|epithelium gene_symbols.tsv")
}
root <- file.path("results", args[1])
map <- unique(read.delim(args[2])[c("gene_id", "gene_symbol")])
motif <- read.csv("results/motif_tests/all_motif_region_results.csv")
selected <- subset(motif, tissue == args[1] & original_BH < 5e-5)
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
