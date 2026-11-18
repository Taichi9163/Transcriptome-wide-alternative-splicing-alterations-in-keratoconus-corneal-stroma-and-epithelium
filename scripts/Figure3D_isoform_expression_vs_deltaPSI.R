#!/usr/bin/env Rscript
# Figure 3D: original rMATS delta PSI for both inclusion and skipping-compatible isoforms.
suppressPackageStartupMessages({
    library(ggplot2)
    library(ragg)
})

args0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)

source(file.path(dirname(normalizePath(sub("^--file=", "", args0[1]))), "Utils_plotting.R"))

root <- ROOT

input <- file.path(ROOT, "results/stroma/manuscript/isoform_scatter_displayed_pairs.tsv")

d <- read.delim(input, quote = "", check.names = FALSE)

d$DE_isoform <- tolower(as.character(d$DE_isoform)) == "true"

stopifnot(nrow(d) == 4047, sum(d$DE_isoform) == 42, length(unique(d$GeneID[d$DE_isoform])) == 29,
    !anyDuplicated(paste(d$EventKey, d$TranscriptID)), all(is.finite(d$log2FoldChange)), all(is.finite(d$DeltaPSI)),
    all(d$DAS_FDR < 0.05), all(abs(d$DeltaPSI) > 0.1), setequal(d$Form, c("Inclusion", "Skipping")))

d$Category <- factor(ifelse(d$DE_isoform, "Event-matched DE isoforms", "Not DE isoforms"), levels = c("Not DE isoforms",
    "Event-matched DE isoforms"))

d$Compatibility <- factor(paste0(d$Form, "-compatible"), levels = c("Inclusion-compatible", "Skipping-compatible"))

span <- max(10, ceiling(max(abs(d$log2FoldChange))/10) * 10)

stopifnot(span == 30)

write.table(d, file.path(root, "results/stroma/figure3", "event_transcript_pairs_plotted.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

summary <- do.call(rbind, lapply(c("All", "Inclusion", "Skipping"), function(f) {
    x <- if (f == "All")
        d
    else d[d$Form == f, ]
    data.frame(Form = f, Pairs = nrow(x), Events = length(unique(x$EventKey)), Transcripts = length(unique(x$TranscriptID)),
        Genes = length(unique(x$GeneID)), DE_pairs = sum(x$DE_isoform), DE_transcripts = length(unique(x$TranscriptID[x$DE_isoform])),
        DE_genes = length(unique(x$GeneID[x$DE_isoform])), OutsideX = sum(abs(x$log2FoldChange) >
            span))
}))

write.table(summary, file.path(root, "results/stroma/figure3", "plot_counts.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE)

base_plot <- function(dat) {
    p <- ggplot(dat, aes(x = log2FoldChange, y = DeltaPSI, color = Category)) + geom_hline(yintercept = 0,
        color = "grey65", linewidth = 0.3) + geom_vline(xintercept = 0, color = "grey65", linewidth = 0.3) +
        geom_hline(yintercept = c(-0.1, 0.1), color = "grey65", linetype = "dashed", linewidth = 0.3) +
        geom_vline(xintercept = c(-1, 1), color = "grey65", linetype = "dashed", linewidth = 0.3)
    p <- p + geom_point(aes(shape = Compatibility), size = 1.35, alpha = 0.65) + geom_point(data = dat[dat$DE_isoform,
        ], aes(shape = Compatibility), size = 2, alpha = 0.95) + scale_shape_manual(values = c(`Inclusion-compatible` = 16,
        `Skipping-compatible` = 17), name = NULL)
    p + scale_color_manual(values = c(`Not DE isoforms` = "#FFD21A", `Event-matched DE isoforms` = "#00C878"),
        drop = FALSE) + scale_x_continuous(breaks = seq(-span, span, length.out = 5)) + scale_y_continuous(breaks = seq(-1,
        1, 0.5), labels = function(z) sprintf("%.1f", z)) + coord_cartesian(xlim = c(-span, span),
        ylim = c(-1, 1), expand = TRUE) + labs(x = "Matched isoform expression (log2 fold change)",
        y = "Alternative splicing (ΔPSI; KC − control)", color = NULL) + guides(color = guide_legend(ncol = 1,
        order = 1, override.aes = list(alpha = 1, size = 2.5, shape = 16)), shape = guide_legend(ncol = 1,
        order = 2, override.aes = list(alpha = 1, size = 2.5, color = "black"))) + theme_classic(base_family = FONT,
        base_size = 12) + theme(aspect.ratio = 1, axis.title = element_text(size = 13), axis.text = element_text(color = "black",
        size = 12.6), legend.position = "top", legend.justification = "right", legend.box = "vertical",
        legend.text = element_text(size = 12), strip.background = element_blank(), strip.text = element_text(size = 14),
        panel.spacing = grid::unit(16, "pt"), plot.margin = margin(8, 12, 8, 8), plot.title = element_text(size = 14,
            hjust = 0.5))
}

dl <- d

groups <- c("Not DE isoforms (inclusion)", "Not DE isoforms (skipping)", "Event-matched DE isoforms (inclusion)",
    "Event-matched DE isoforms (skipping)")

dl$Category <- factor(paste0(ifelse(dl$DE_isoform, "Event-matched DE isoforms", "Not DE isoforms"),
    " (", tolower(dl$Form), ")"), levels = groups)

p <- suppressMessages(base_plot(dl) + scale_color_manual(values = setNames(c("#FFD21A", "#FFD21A",
    "#00C878", "#00C878"), groups), drop = FALSE) + guides(shape = "none", color = guide_legend(ncol = 1,
    override.aes = list(alpha = 1, size = 2.5, shape = c(16, 17, 16, 17)))) + theme(legend.position = "top",
    legend.justification = "right", legend.box.just = "right", legend.text = element_text(size = 10.5),
    legend.box.spacing = grid::unit(12, "pt")))

b <- ggplot_build(p)

stopifnot(identical(b$data[[5]]$y, d$DeltaPSI), identical(b$data[[6]]$y, d$DeltaPSI[d$DE_isoform]),
    all(b$data[[5]]$shape[d$Form == "Inclusion"] == 16), all(b$data[[5]]$shape[d$Form == "Skipping"] ==
        17), all(b$data[[5]]$colour[d$DE_isoform] == "#00C878"), all(b$data[[5]]$colour[!d$DE_isoform] ==
        "#FFD21A"))

save_panel(p, file.path(ROOT, "results/stroma/figure3/Figure3D"), 6.1, 7.2)

print(summary, row.names = FALSE)
