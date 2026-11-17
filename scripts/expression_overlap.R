# Figure 3A-C. Gene membership and the most strongly changed event per DAS gene.
root <- "results/stroma"
classes <- c("SE", "A5SS", "A3SS", "MXE", "RI")
events <- do.call(rbind, lapply(classes, function(k) {
  x <- read.delim(file.path(root, "das", paste0(k, ".tsv")), check.names = FALSE)
  data.frame(gene_id = x$GeneID, gene_symbol = x$geneSymbol, event_class = k, event_id = x$ID,
             delta = x$IncLevelDifference)
}))
events <- events[order(-abs(events$delta), events$event_class, events$event_id), ]
representative <- events[!duplicated(events$gene_id), ]
gene <- read.delim(file.path(root, "gene/wald.tsv"))
iso <- read.delim(file.path(root, "isoform/wald.tsv"))
de_iso_genes <- unique(iso$gene_id[iso$DE])
all_genes <- unique(c(events$gene_id, gene$feature_id[gene$DE], de_iso_genes))
sets <- data.frame(gene_id = all_genes, DAS = all_genes %in% events$gene_id,
                   DEG = all_genes %in% gene$feature_id[gene$DE],
                   DE_isoform = all_genes %in% de_iso_genes)
out <- file.path(root, "overlap")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(sets, file.path(out, "gene_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
points <- merge(representative, gene, by.x = "gene_id", by.y = "feature_id", all.x = TRUE)
write.table(points, file.path(out, "gene_scatter.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(colSums(sets[-1L]))
print(c(DAS_DEG = sum(sets$DAS & sets$DEG), DAS_DE_isoform = sum(sets$DAS & sets$DE_isoform)))
