#!/usr/bin/env Rscript
# Prepare manuscript tables with strict FDR < 0.05 and absolute delta PSI > 0.1.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
source(file.path(dirname(normalizePath(sub("^--file=", "", arg))), "utils_plotting.R"))
out <- file.path(ROOT, "results/stroma/manuscript")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ROOT, "results/stroma/figure3"), recursive = TRUE, showWarnings = FALSE)
types <- c("SE", "RI", "A5SS", "A3SS", "MXE")
normalize_gene <- function(x) sub("[.][0-9]+$", "", trimws(as.character(x)))
events <- lapply(types, function(t) {
    d <- read.delim(file.path(ROOT, "data/stroma/rmats", paste0(t, ".MATS.JC.txt")), check.names = TRUE,
        stringsAsFactors = FALSE)
    d <- d[!is.na(d$FDR) & !is.na(d$IncLevelDifference) & d$FDR < 0.05 & abs(d$IncLevelDifference) >
        0.1 & !is.na(d$IncLevel1) & !is.na(d$IncLevel2), , drop = FALSE]
    d$GeneID <- normalize_gene(d$GeneID)
    d$EventType <- t
    p <- file.path(out, "das", paste0(t, "_JC"))
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    write_tsv(d, file.path(p, paste0("rMATS_", t, "_JC_combined.txt")))
    d
})
cols <- unique(unlist(lapply(events, names)))
all <- do.call(rbind, lapply(events, function(d) {
    for (n in setdiff(cols, names(d))) d[[n]] <- NA
    d[, cols]
}))
stopifnot(nrow(all) == 1073L, length(unique(all$GeneID)) == 807L)
write_tsv(all, file.path(out, "source_DAS_events.tsv"))
psi <- t(vapply(strsplit(paste(all$IncLevel1, all$IncLevel2, sep = ","), ","), as.numeric, numeric(9)))
rownames(psi) <- paste(all$EventType, all$ID, sep = ":")
colnames(psi) <- c(paste0("Keratoconus", 1:5), paste0("Control", 1:4))
psi <- psi[complete.cases(psi), , drop = FALSE]
write_tsv(data.frame(EventKey = rownames(psi), psi, check.names = FALSE), file.path(out, "psi_matrix.tsv"))
summary <- do.call(rbind, lapply(types, function(t) {
    x <- all[all$EventType == t, ]
    data.frame(EventType = t, Events = nrow(x), Genes = length(unique(x$GeneID)))
}))
write_tsv(summary, file.path(out, "das_summary.tsv"))
g <- read_tsv(file.path(ROOT, "data/stroma/expression/gene_deseq2.tsv"))
g$GeneID <- normalize_gene(rownames(g))
g <- g[!duplicated(g$GeneID), ]
sig <- !is.na(g$padj) & is.finite(g$log2FoldChange) & g$padj < 0.05 & abs(g$log2FoldChange) > 1
deg <- g$GeneID[sig]
das <- unique(all$GeneID)
u <- sort(union(deg, das))
m <- data.frame(GeneID = u, DEG = u %in% deg, DAS = u %in% das)
m$GeneSymbol <- all$geneSymbol[match(m$GeneID, all$GeneID)]
stopifnot(sum(m$DEG) == 414L, sum(m$DAS) == 807L, sum(m$DEG & m$DAS) == 8L)
write_tsv(m, file.path(ROOT, "results/stroma/figure3/DEG_DAS_gene_membership.tsv"))
i <- read_tsv(file.path(ROOT, "data/stroma/expression/transcript_deseq2.tsv"))
i$GeneID <- normalize_gene(i$GeneID)
expected <- !is.na(i$padj) & is.finite(i$log2FoldChange) & i$padj < 0.05 & abs(i$log2FoldChange) >
    1
stopifnot(all(expected == (tolower(as.character(i$DE_isoform)) == "true")))
i$DE_isoform <- expected
ide <- unique(i$GeneID[i$DE_isoform])
u <- sort(union(ide, das))
m <- data.frame(GeneID = u, DE_isoform = u %in% ide, DAS = u %in% das)
m$GeneSymbol <- all$geneSymbol[match(m$GeneID, all$GeneID)]
stopifnot(sum(m$DE_isoform) == 410L, sum(m$DE_isoform & m$DAS) == 55L)
write_tsv(m, file.path(out, "isoform_gene_membership.tsv"))
