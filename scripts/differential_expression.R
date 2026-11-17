# RSEM expected counts and effective lengths, analyzed separately per cohort.
suppressPackageStartupMessages({ library(DESeq2); library(tximport) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[1] %in% c("gene", "isoform")) {
  stop("Usage: Rscript scripts/differential_expression.R gene|isoform samples.tsv output_dir")
}
isoform <- args[1] == "isoform"
samples <- read.delim(args[2], stringsAsFactors = FALSE)
stopifnot(all(c("sample", "group", "path") %in% names(samples)), !anyDuplicated(samples$sample),
          !anyNA(samples), setequal(samples$group, c("Control", "Keratoconus")),
          all(table(samples$group) >= 2))
files <- setNames(samples$path, samples$sample)
stopifnot(all(file.exists(files)), !anyDuplicated(normalizePath(files)))
txi <- tximport(files, type = "rsem", txIn = isoform, txOut = isoform, countsFromAbundance = "no")
txi$length[txi$length == 0] <- 1
samples$condition <- factor(samples$group, levels = c("Control", "Keratoconus"))
rownames(samples) <- samples$sample
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)
if (isoform) dds <- dds[rowSums(counts(dds)) >= 10, ]
fit <- DESeq(dds)
# The gene workflow used the DESeq2 default independent-filtering target.
# The isoform workflow set its target to 0.05; both use padj < 0.05 for calls.
filter_alpha <- if (isoform) 0.05 else 0.10
res <- as.data.frame(results(fit, contrast = c("condition", "Keratoconus", "Control"), alpha = filter_alpha))
res$feature_id <- rownames(res)
if (isoform) {
  refs <- read.delim(files[[1]], stringsAsFactors = FALSE)
  res$gene_id <- refs$gene_id[match(res$feature_id, refs$transcript_id)]
  stopifnot(!anyNA(res$gene_id))
}
res$DE <- is.finite(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
write.table(res, file.path(args[3], "wald.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(feature_id = rownames(txi$abundance), txi$abundance, check.names = FALSE),
            file.path(args[3], "tpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(args[3], "session_info.txt"))
message(sum(res$DE), " features meet the differential-expression thresholds.")
