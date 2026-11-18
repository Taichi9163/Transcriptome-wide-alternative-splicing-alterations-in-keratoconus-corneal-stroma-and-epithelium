#!/usr/bin/env Rscript
# Figure 3B: largest absolute delta PSI per gene; ties use smaller FDR then source order.
suppressPackageStartupMessages({
    library(ggplot2)
    library(ragg)
})

args0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)

source(file.path(dirname(normalizePath(sub("^--file=", "", args0[1]))), "common.R"))

deg_file <- file.path(ROOT, "data/stroma/expression/gene_deseq2.tsv")

das_base <- file.path(ROOT, "results/stroma/manuscript/das")

output_dir <- file.path(ROOT, "results/stroma/figure3")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

event_types <- c("SE", "RI", "A5SS", "A3SS", "MXE")

deg_padj_cutoff <- 0.05

deg_abs_log2fc_cutoff <- 1

das_fdr_cutoff <- 0.05

das_abs_dpsi_cutoff <- 0.1

normalize_ensembl <- function(x) {
    sub("[.][0-9]+$", "", trimws(as.character(x)))
}

read_tsv_preserve <- function(path) {
    read.delim(path, header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

deg <- read_tsv_preserve(deg_file)

deg$GeneID <- normalize_ensembl(rownames(deg))

stopifnot(all(c("GeneID", "log2FoldChange", "padj") %in% colnames(deg)))

deg <- deg[!duplicated(deg$GeneID), , drop = FALSE]

deg$DEG <- !is.na(deg$padj) & is.finite(deg$log2FoldChange) & deg$padj < deg_padj_cutoff & abs(deg$log2FoldChange) >
    deg_abs_log2fc_cutoff

deg_genes <- unique(deg$GeneID[deg$DEG])

das_list <- lapply(event_types, function(event_type) {
    path <- file.path(das_base, paste0(event_type, "_JC"), paste0("rMATS_", event_type, "_JC_combined.txt"))
    if (!file.exists(path)) {
        stop("Missing DAS input: ", path)
    }
    x <- read_tsv_preserve(path)
    required <- c("GeneID", "geneSymbol", "ID", "FDR", "IncLevelDifference")
    if (!all(required %in% colnames(x))) {
        stop("Required columns are missing from: ", path)
    }
    x$GeneID <- normalize_ensembl(x$GeneID)
    x$EventType <- event_type
    x <- x[!is.na(x$FDR) & !is.na(x$IncLevelDifference) & x$FDR < das_fdr_cutoff & abs(x$IncLevelDifference) >
        das_abs_dpsi_cutoff & nzchar(x$GeneID) & x$GeneID != "NA", c("GeneID", "geneSymbol", "EventType",
        "ID", "FDR", "IncLevelDifference"), drop = FALSE]
    x
})

das_events <- do.call(rbind, das_list)

rownames(das_events) <- NULL

das_genes <- unique(das_events$GeneID)

das_representative <- das_events[order(das_events$GeneID, -abs(das_events$IncLevelDifference), das_events$FDR),
    , drop = FALSE]

das_representative <- das_representative[!duplicated(das_representative$GeneID), , drop = FALSE]

union_genes <- sort(union(deg_genes, das_genes))

plot_data <- data.frame(GeneID = union_genes, stringsAsFactors = FALSE)

deg_match <- match(plot_data$GeneID, deg$GeneID)

das_match <- match(plot_data$GeneID, das_representative$GeneID)

plot_data$GeneSymbol <- das_representative$geneSymbol[das_match]

plot_data$log2FoldChange <- deg$log2FoldChange[deg_match]

plot_data$DEG_padj <- deg$padj[deg_match]

plot_data$DeltaPSI <- das_representative$IncLevelDifference[das_match]

plot_data$DAS_FDR <- das_representative$FDR[das_match]

plot_data$DAS_EventType <- das_representative$EventType[das_match]

plot_data$DAS_EventID <- das_representative$ID[das_match]

plot_data$DEG <- plot_data$GeneID %in% deg_genes

plot_data$DAS <- plot_data$GeneID %in% das_genes

plot_data$Category <- ifelse(plot_data$DEG & plot_data$DAS, "DAS & DEG", ifelse(plot_data$DAS, "DAS",
    "DEG"))

plot_data$Category <- factor(plot_data$Category, levels = c("DAS", "DEG", "DAS & DEG"))

plot_complete <- plot_data[is.finite(plot_data$log2FoldChange) & is.finite(plot_data$DeltaPSI),
    , drop = FALSE]

plot_display <- plot_complete[plot_complete$Category != "DEG", , drop = FALSE]

stopifnot(length(deg_genes) == 414L, length(das_genes) == 807L, sum(plot_display$Category == "DAS") ==
    797L, sum(plot_display$Category == "DAS & DEG") == 8L)

write.table(plot_data, file.path(output_dir, "DEG_DAS_log2FC_deltaPSI_plot_data.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE, na = "")

summary_table <- data.frame(Metric = c("DEG genes", "DAS genes", "DEG and DAS genes", "Union genes",
    "Genes with finite coordinates", "DAS-only genes plotted", "DEG-only genes plotted", "DAS and DEG genes plotted",
    "DAS genes without matched log2FC"), Value = c(length(deg_genes), length(das_genes), length(intersect(deg_genes,
    das_genes)), length(union_genes), nrow(plot_complete), sum(plot_complete$Category == "DAS"),
    sum(plot_complete$Category == "DEG"), sum(plot_complete$Category == "DAS & DEG"), sum(plot_data$DAS &
        !is.finite(plot_data$log2FoldChange))), stringsAsFactors = FALSE)

write.table(summary_table, file.path(output_dir, "DEG_DAS_log2FC_deltaPSI_summary.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)

category_colors <- c(DAS = "#FFD21A", DEG = "#00A6B2", `DAS & DEG` = "#FF5A5F")

p <- ggplot(plot_display, aes(x = log2FoldChange, y = DeltaPSI, color = Category)) + geom_hline(yintercept = 0,
    linewidth = 0.35, color = "#8A8A8A") + geom_vline(xintercept = 0, linewidth = 0.35, color = "#8A8A8A") +
    geom_hline(yintercept = c(-das_abs_dpsi_cutoff, das_abs_dpsi_cutoff), linewidth = 0.32, linetype = "dashed",
        color = "#A6A6A6") + geom_vline(xintercept = c(-deg_abs_log2fc_cutoff, deg_abs_log2fc_cutoff),
    linewidth = 0.32, linetype = "dashed", color = "#A6A6A6") + geom_point(size = 1.75, alpha = 0.86,
    stroke = 0) + geom_point(data = plot_display[plot_display$Category == "DAS & DEG", , drop = FALSE],
    size = 2.35, alpha = 1, stroke = 0) + scale_color_manual(values = category_colors, breaks = c("DAS",
    "DAS & DEG"), name = NULL) + scale_x_continuous(breaks = seq(-10, 10, by = 5), expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = seq(-1, 1, by = 0.5), expand = expansion(mult = 0)) + coord_cartesian(xlim = c(-10,
    10), ylim = c(-1, 1), expand = FALSE) + labs(x = "Gene expression (log2 fold change)", y = "Alternative splicing (ΔPSI)") +
    guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = c(2.2, 2.6),
        alpha = c(0.85, 1)))) + theme_classic(base_family = FONT, base_size = 10) + theme(aspect.ratio = 1,
    axis.title = element_text(size = 12, color = "#111111"), axis.title.x = element_text(margin = margin(t = 7)),
    axis.title.y = element_text(margin = margin(r = 7)), axis.text = element_text(size = 10.5, color = "#111111"),
    axis.text.x = element_text(margin = margin(t = 4)), axis.text.y = element_text(margin = margin(r = 4)),
    axis.line = element_line(linewidth = 0.45, color = "#111111"), axis.ticks = element_line(linewidth = 0.4,
        color = "#111111"), axis.ticks.length = unit(2, "mm"), legend.position = "top", legend.justification = "right",
    legend.direction = "vertical", legend.text = element_text(size = 9.5, color = "#111111"), legend.key.width = unit(4.5,
        "mm"), legend.spacing.x = unit(1.2, "mm"), legend.spacing.y = unit(0.8, "mm"), plot.margin = margin(5,
        7, 8, 8), plot.background = element_rect(fill = "white", color = NA))

save_panel(p, file.path(output_dir, "Figure3B"), 5, 5.5)

caption <- paste0("Relationship between gene expression change and alternative splicing in ", "keratoconus stroma. Each point represents one DAS gene with a finite gene-level fold change. The DEG ",
    "and DAS sets were evaluated separately. DEG status was determined from DESeq2 (adjusted P < 0.05 and ",
    "|log2 fold change| > 1), and DAS status from rMATS (FDR < 0.05 and |delta ", "PSI| > 0.1). For genes with multiple significant DAS events, the event ",
    "with the largest absolute delta PSI was used. DEG-only genes were omitted ", "from the displayed panel. Colors denote DAS-only and DAS-and-DEG genes.")

writeLines(caption, file.path(output_dir, "DEG_DAS_log2FC_deltaPSI_figure_legend.txt"))

message("Analysis tables and plot specification validated.")
