#!/usr/bin/env Rscript
# Figure 1B-C: scaled PCA; Euclidean/Ward.D2 clustering on raw PSI with row Z-score colors.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(gplots)
  library(RColorBrewer)
})
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
source(file.path(dirname(normalizePath(sub("^--file=", "", arg))), "utils_plotting.R"))
out <- file.path(ROOT, "results/stroma/figure1")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- read_tsv(file.path(ROOT, "results/stroma/manuscript/psi_matrix.tsv"))
x <- as.matrix(d[, -1])
rownames(x) <- d$EventKey
stopifnot(nrow(x) == 1073L, ncol(x) == 9L, all(is.finite(x)), all(apply(x, 1, sd) > 0))
pca <- prcomp(t(x), scale. = TRUE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
stopifnot(abs(variance[1] - 46.316) < 0.001, abs(variance[2] - 16.153) < 0.001)
points <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
points$Group <- ifelse(grepl("^Keratoconus", points$Sample), "Keratoconus", "Control")
write_tsv(points, file.path(out, "Figure1B_PCA_scores.tsv"))
write_tsv(data.frame(PC = seq_along(variance), Variance_percent = variance), file.path(out, "Figure1B_variance.tsv"))
p <- ggplot(points, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(
    Control = 16,
    Keratoconus = 16
  )) +
  scale_color_manual(values = c(Control = "black", Keratoconus = "#FFA500")) +
  labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2])) +
  scale_x_reverse(limits = c(
    45,
    -45
  ), breaks = seq(-40, 40, 20)) +
  scale_y_continuous(limits = c(-40, 40), breaks = seq(
    -40, 40,
    20
  )) +
  theme_linedraw(base_size = 15, base_family = FONT) +
  theme(aspect.ratio = 1, panel.grid = element_line(
    linewidth = 1,
    colour = "black", linetype = "dashed"
  ), panel.border = element_rect(
    colour = "black", fill = NA,
    linewidth = 1
  ))
save_panel(p, file.path(out, "Figure1B"), 6.5, 6)
rows <- hclust(dist(x, method = "euclidean"), method = "ward.D2")
columns <- hclust(dist(t(x), method = "euclidean"), method = "ward.D2")
column_tree <- reorder(as.dendrogram(columns), wts = ifelse(grepl("^Control", colnames(x)), 0, 1), agglo.FUN = mean)
write_tsv(data.frame(Sample = colnames(x)[order.dendrogram(column_tree)]), file.path(out, "Figure1C_sample_order.tsv"))
heat <- function() {
  par(family = FONT)
  heatmap.2(x,
    col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), scale = "row", dendrogram = "both",
    key = TRUE, density.info = "none", trace = "none", Rowv = as.dendrogram(rows), Colv = column_tree,
    ColSideColors = ifelse(grepl("^Keratoconus", colnames(x)), "#FFA500", "black"), labRow = FALSE,
    cexCol = 1, margins = c(8, 5), key.title = "Row Z-score", main = NULL
  )
}
if (!TABLES_ONLY) {
  pdf_device(file.path(out, "Figure1C.pdf"), 10, 8)
  heat()
  dev.off()
  ragg::agg_png(file.path(out, "Figure1C.png"), 10, 8, units = "in", res = 600)
  heat()
  dev.off()
  ragg::agg_tiff(file.path(out, "Figure1C.tiff"), 10, 8, units = "in", res = 600, compression = "lzw")
  heat()
  dev.off()
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
