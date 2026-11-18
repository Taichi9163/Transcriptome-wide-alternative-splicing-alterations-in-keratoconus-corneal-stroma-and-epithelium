# Shared numeric preparation and plotting for Figure 4D / Supplementary Figure 10B.
# Input: corrected NEASE Reactome_Python.tsv, not the older Reactome.csv export.

reactome_top10 <- function(path) {
    frame <- read.delim(path, sep = "\t", quote = "", check.names = FALSE,
                        stringsAsFactors = FALSE)
    gene_column <- "Spliced genes (number of interactions affecting the pathway)"
    required <- c("Pathway ID", "Pathway name", "Source", gene_column,
                  "p_value", "adj p_value", "Nease score")
    if (!all(required %in% names(frame))) stop("Missing required NEASE Reactome columns.")
    if (anyNA(frame$`Pathway ID`) || any(!nzchar(frame$`Pathway ID`)) ||
        anyDuplicated(frame$`Pathway ID`)) stop("Pathway IDs must be present and unique.")
    if (anyNA(frame$Source) || any(frame$Source != "Reactome")) {
        stop("Use the complete Reactome-only pathway export.")
    }
    if (anyNA(frame$`Pathway name`) || any(!nzchar(frame$`Pathway name`))) {
        stop("Pathway names must be present.")
    }
    for (column in c("p_value", "adj p_value", "Nease score")) {
        if (!is.numeric(frame[[column]]) || any(!is.finite(frame[[column]]))) {
            stop("Non-finite or nonnumeric values in ", column)
        }
    }
    if (any(frame$p_value < 0 | frame$p_value > 1) ||
        any(frame$`adj p_value` < 0 | frame$`adj p_value` > 1) || any(frame$`Nease score` < 0)) {
        stop("Invalid probability or NEASE score.")
    }
    # Validate, but retain the supplied adjusted values for ranking/color data.
    bh_difference <- max(abs(p.adjust(frame$p_value, method = "BH") - frame$`adj p_value`))
    if (bh_difference > 1e-10) stop("Adjusted P values do not match BH over this Reactome table.")
    eligible <- frame[frame$`adj p_value` <= .05, required, drop = FALSE]
    eligible <- eligible[order(-eligible$`Nease score`, eligible$`Pathway ID`, method = "radix"), ]
    if (nrow(eligible) < 10L) stop("Fewer than ten significant Reactome pathways.")
    top <- head(eligible, 10L)
    count_genes <- function(value) {
        if (is.na(value)) return(0L)
        entries <- strsplit(value, ",", fixed = TRUE)[[1]]
        entries <- trimws(gsub("\\s*\\([^)]*\\)", "", entries, perl = TRUE))
        length(unique(entries[nzchar(entries)]))
    }
    top$`Splicing-affected genes` <- vapply(top[[gene_column]], count_genes, integer(1))
    if (any(top$`Splicing-affected genes` < 1L) || any(top$`adj p_value` <= 0)) {
        stop("Top pathways require affected genes and positive adjusted P values.")
    }
    top$`-log10 adjusted P` <- -log10(top$`adj p_value`)
    top <- cbind(Rank = seq_len(nrow(top)), top)
    rownames(top) <- NULL
    list(top = top, summary = data.frame(
        tested_pathways = nrow(frame), significant_pathways = nrow(eligible),
        displayed_pathways = nrow(top), adjusted_p_cutoff = .05,
        max_bh_absolute_difference = bh_difference))
}

reactome_plot <- function(top, cohort, font = "sans") {
    if (!cohort %in% c("stroma", "epithelium")) stop("Unknown cohort.")
    epithelial <- cohort == "epithelium"
    top$pathway_axis <- factor(top$`Pathway ID`, levels = rev(top$`Pathway ID`))
    labels <- top$`Pathway name`
    if (!epithelial) labels <- vapply(labels, function(x) paste(strwrap(x, width = 38), collapse = "\n"), "")
    names(labels) <- top$`Pathway ID`
    # The old epithelial script capped gene counts at seven, hiding the corrected
    # 14-gene pathway. Keep every selected point using a data-derived scale.
    gene_breaks <- if (epithelial) sort(unique(c(min(top$`Splicing-affected genes`),
                                               7L, max(top$`Splicing-affected genes`)))) else scales::breaks_pretty(n = 3)
    p <- ggplot2::ggplot(top, ggplot2::aes(x = `Nease score`, y = pathway_axis)) +
        ggplot2::geom_point(ggplot2::aes(size = `Splicing-affected genes`, colour = `-log10 adjusted P`)) +
        ggplot2::scale_y_discrete(labels = labels) +
        ggplot2::scale_x_continuous(expand = if (epithelial) ggplot2::expansion(mult = c(.10, .18))
                                  else ggplot2::expansion(add = c(5, 5))) +
        ggplot2::scale_colour_gradient(low = "white", high = if (epithelial) "#CC79A7" else "#FFA500",
            limits = c(0, max(top$`-log10 adjusted P`)), oob = scales::squish,
            name = expression(-Log[10]("P-adjust"))) +
        ggplot2::scale_size_continuous(name = "Splicing-affected\ngenes",
            range = if (epithelial) c(4, 10) else c(4.5, 10), breaks = gene_breaks) +
        ggplot2::labs(x = "NEASE score", y = NULL) +
        ggplot2::theme_bw(base_size = if (epithelial) 11 else 15, base_family = font) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(colour = "black", size = if (epithelial) 10 else 13),
            axis.text.y = ggplot2::element_text(colour = "black", size = if (epithelial) 10 else 14, lineheight = .95),
            axis.title.x = ggplot2::element_text(colour = "black", size = if (epithelial) 10 else 15,
                                                margin = ggplot2::margin(t = 7)),
            legend.position = "right", legend.box = "vertical",
            legend.box.just = "top", legend.spacing.y = grid::unit(6, "mm"),
            legend.title = ggplot2::element_text(size = if (epithelial) 10 else 12.5),
            legend.text = ggplot2::element_text(colour = "black", size = if (epithelial) 9 else 11.5),
            panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .7),
            panel.grid.major = ggplot2::element_line(colour = "black", linetype = "dotted", linewidth = .2),
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(8, 10, 8, 8)) +
        ggplot2::guides(size = ggplot2::guide_legend(order = 1), colour = ggplot2::guide_colourbar(
            order = 2, barwidth = grid::unit(3.5, "mm"), barheight = grid::unit(34, "mm")))
    if (epithelial) p <- p + ggplot2::theme(aspect.ratio = 1.6)
    p
}

run_reactome_panel <- function(cohort, panel, output_folder) {
    args <- commandArgs(trailingOnly = TRUE)
    path <- file.path(ROOT, "data", cohort, "nease", "Reactome_Python.tsv")
    i <- 1L
    while (i <= length(args)) {
        if (args[i] == "--tables-only") {
            i <- i + 1L
        } else if (args[i] == "--pathways" && i < length(args)) {
            path <- path.expand(args[i + 1L])
            if (!grepl("^/", path)) path <- file.path(ROOT, path)
            i <- i + 2L
        } else stop("Usage: [--tables-only] [--pathways PATH]")
    }
    if (!file.exists(path)) stop("Missing pathway table: ", path, "; run setup_prepare_inputs.py first.")
    result <- reactome_top10(path)
    plot <- reactome_plot(result$top, cohort, FONT)
    # Check all ten points before any output. Building layers does not render.
    layer <- ggplot2::ggplot_build(plot)$data[[1]]
    if (nrow(layer) != 10L || any(!is.finite(layer$x)) || any(!is.finite(layer$y)) ||
        any(!is.finite(layer$size)) || anyNA(layer$colour)) stop("A selected pathway is missing from the plot layers.")
    output <- file.path(ROOT, "results", cohort, output_folder)
    write_tsv(result$top, file.path(output, "Reactome_top10_for_figure.tsv"))
    write_tsv(result$summary, file.path(output, "Reactome_summary.tsv"))
    write_tsv(data.frame(input = normalizePath(path), md5 = unname(tools::md5sum(path))),
              file.path(output, "Reactome_input_provenance.tsv"))
    save_panel(plot, file.path(output, paste0(panel, "_Reactome_enrichment")), 7.8,
               if (cohort == "epithelium") 6.5 else 5.2)
    message(panel, ": verified ", result$summary$tested_pathways, " pathways; exported top ten",
            if (TABLES_ONLY) " (tables only; no figures)." else " and figures.")
    invisible(result)
}
