#!/usr/bin/env Rscript
# PCA and Ward.D2 clustering of significant, PSI-complete events in one cohort.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gplots)
  library(RColorBrewer)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "utils_plotting.R"))
args <- setdiff(script_arguments(), "--tables-only")
if (length(args) != 1 || !args %in% c("stroma", "epithelium")) {
  stop("Usage: Rscript scripts/supp_figure_07bc_psi_pca_heatmap.R stroma|epithelium [--tables-only]")
}
cohort <- args[1]
output <- file.path(ROOT, "results", cohort, "figures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
event_types <- c("SE", "RI", "A5SS", "A3SS", "MXE")
matrices <- lapply(event_types, function(event_type) {
  events <- read_tsv(file.path(ROOT, "results", cohort, "das", paste0(event_type, ".tsv")))
  cases <- strsplit(events$IncLevel1, ",", fixed = TRUE)
  controls <- strsplit(events$IncLevel2, ",", fixed = TRUE)
  stopifnot(length(unique(lengths(cases))) == 1, length(unique(lengths(controls))) == 1)
  n_cases <- lengths(cases)[1]
  n_controls <- lengths(controls)[1]
  psi <- cbind(matrix(as.numeric(unlist(cases)), ncol = n_cases, byrow = TRUE), matrix(as.numeric(unlist(controls)),
    ncol = n_controls, byrow = TRUE
  ))
  rownames(psi) <- paste(event_type, events$ID, sep = ":")
  colnames(psi) <- c(paste0("Keratoconus", seq_len(n_cases)), paste0("Control", seq_len(n_controls)))
  psi
})
x <- do.call(rbind, matrices)
x <- x[complete.cases(x), , drop = FALSE]
stopifnot(all(is.finite(x)), all(apply(x, 1, sd) > 0))
case_color <- if (cohort == "stroma") "#FFA500" else "#CC79A7"
point_shape <- if (cohort == "stroma") 16 else 17
sample_group <- ifelse(grepl("^Control", colnames(x)), "Control", "Keratoconus")
row_tree <- as.dendrogram(hclust(dist(x), method = "ward.D2"))
column_tree <- as.dendrogram(hclust(dist(t(x)), method = "ward.D2"))
# Rotate only the top branch to put the group with more controls on the left.
if (sum(grepl("^Control", labels(column_tree[[2]]))) > sum(grepl("^Control", labels(column_tree[[1]])))) {
  tree_attributes <- attributes(column_tree)
  column_tree <- column_tree[c(2, 1)]
  attributes(column_tree) <- tree_attributes
}
pca <- prcomp(t(x), scale. = TRUE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
scores <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = sample_group)
write_tsv(scores, file.path(output, "das_pca_scores.tsv"))
write_tsv(data.frame(PC = seq_along(variance), Variance_percent = variance), file.path(output, "das_pca_variance.tsv"))
write_tsv(data.frame(Sample = labels(column_tree)), file.path(output, "das_heatmap_sample_order.tsv"))
plot <- ggplot(scores, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 5) +
  scale_color_manual(values = c(
    Control = "black",
    Keratoconus = case_color
  )) +
  scale_shape_manual(values = c(Control = point_shape, Keratoconus = point_shape)) +
  labs(x = sprintf("PC1 (%.3g%%)", variance[1]), y = sprintf("PC2 (%.3g%%)", variance[2])) +
  xlim(
    40,
    -40
  ) +
  ylim(-40, 40) +
  theme_linedraw(base_size = 15) +
  theme(aspect.ratio = 1, panel.border = element_rect(
    fill = NA,
    linewidth = 1
  ), panel.grid = element_line(linewidth = 1, colour = "black", linetype = "dashed"))
invisible(ggplot_build(plot))
if (!TABLES_ONLY) {
  pdf_device(file.path(output, "das_heatmap.pdf"), 10, 8)
  heatmap.2(x,
    col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), scale = "row", dendrogram = "both",
    key = TRUE, keysize = 1.5, density.info = "none", trace = "none", Rowv = row_tree, Colv = column_tree,
    ColSideColors = ifelse(sample_group == "Control", "black", case_color), cexRow = 0.4, cexCol = 1,
    margins = c(7, 6), main = ""
  )
  dev.off()
  ggsave(file.path(output, "das_pca.pdf"), plot, width = 6, height = 6, device = pdf_device)
}
