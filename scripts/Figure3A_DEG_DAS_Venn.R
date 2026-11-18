#!/usr/bin/env Rscript
# Figure 3A: DEG/DAS gene overlap; circle area is schematic, counts are exact.
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggforce)
    library(ragg)
})

args0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)

source(file.path(dirname(normalizePath(sub("^--file=", "", args0[1]))), "Utils_plotting.R"))

output_dir <- file.path(ROOT, "results/stroma/figure3")

membership_file <- file.path(output_dir, "DEG_DAS_gene_membership.tsv")

membership <- read.delim(membership_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE)

n_deg <- sum(membership$DEG)

n_das <- sum(membership$DAS)

n_overlap <- sum(membership$DEG & membership$DAS)

n_deg_only <- sum(membership$DEG & !membership$DAS)

n_das_only <- sum(!membership$DEG & membership$DAS)

overlap_genes <- sort(unique(membership$GeneSymbol[membership$DEG & membership$DAS]))

overlap_genes <- overlap_genes[!is.na(overlap_genes) & nzchar(overlap_genes)]

stopifnot(n_deg == 414L, n_das == 807L, n_overlap == 8L, n_deg_only == 406L, n_das_only == 799L,
    length(overlap_genes) == 8L)

deg_color <- "#00A6B2"

das_color <- "#FFD21A"

overlap_color <- "#FF5A5F"

circle_radius <- 1.05

circle_df <- data.frame(set = factor(c("DEG", "DAS"), levels = c("DEG", "DAS")), x0 = c(0.65, -0.65),
    y0 = c(0, 0), r = c(circle_radius, circle_radius))

lens_x <- seq(-0.4, 0.4, length.out = 401)

lens_upper <- pmin(sqrt(pmax(0, circle_radius^2 - (lens_x + 0.65)^2)), sqrt(pmax(0, circle_radius^2 -
    (lens_x - 0.65)^2)))

overlap_df <- data.frame(x = c(lens_x, rev(lens_x)), y = c(lens_upper, -rev(lens_upper)))

overlap_gene_label_left <- paste(overlap_genes[1:4], collapse = "\n")

overlap_gene_label_right <- paste(overlap_genes[5:8], collapse = "\n")

p <- ggplot() + annotate("text", x = c(-1.08, 0, 1.08), y = c(0, 0.15, 0), label = c(n_das_only,
    n_overlap, n_deg_only), family = FONT, size = 5) + annotate("text", x = c(-0.8, 0.8), y = 1.18,
    label = c("DAS genes", "DEGs"), family = FONT, size = 4.5) + geom_circle(data = circle_df, aes(x0 = x0,
    y0 = y0, r = r, fill = set, color = set), alpha = 0.2, linewidth = 0.48, show.legend = FALSE) +
    geom_polygon(data = overlap_df, aes(x = x, y = y), inherit.aes = FALSE, fill = overlap_color,
        color = NA, alpha = 0.55) + annotate("segment", x = 0.05, xend = 0.05, y = -0.18, yend = -1.18,
    linewidth = 0.42, color = "#333333") + annotate("point", x = 0.05, y = -0.18, size = 2, color = "#111111") +
    annotate("rect", xmin = -0.95, xmax = 1.05, ymin = -2.02, ymax = -1.18, fill = "white", color = "#555555",
        linewidth = 0.32) + annotate("text", x = -0.43, y = -1.6, label = overlap_gene_label_left,
    family = FONT, fontface = "italic", size = 3.25, lineheight = 1, color = "#111111") + annotate("text",
    x = 0.53, y = -1.6, label = overlap_gene_label_right, family = FONT, fontface = "italic", size = 3.25,
    lineheight = 1, color = "#111111") + scale_fill_manual(values = c(DEG = deg_color, DAS = das_color)) +
    scale_color_manual(values = c(DEG = deg_color, DAS = das_color)) + coord_fixed(xlim = c(-1.82,
    1.82), ylim = c(-2.12, 1.4), clip = "off") + theme_void(base_family = FONT) + theme(plot.background = element_rect(fill = "white",
    color = NA), plot.margin = margin(4, 4, 4, 4))

save_panel(p, file.path(output_dir, "Figure3A"), 4.7, 4.3)

caption <- paste0("Venn diagram showing the overlap between differentially expressed genes ", "(DEGs; adjusted P < 0.05 and |log2 fold change| > 1) and genes with ",
    "differential alternative splicing (DAS genes; FDR < 0.05 and |delta PSI| ", "> 0.1) in keratoconus stroma. Of 414 DEGs and 807 DAS genes, eight genes ",
    "were shared between the two sets.")

writeLines(caption, file.path(output_dir, "DEG_DAS_venn_figure_legend.txt"))

message("Analysis tables and plot specification validated.")
