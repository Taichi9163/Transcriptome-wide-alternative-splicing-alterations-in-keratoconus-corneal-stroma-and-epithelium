#!/usr/bin/env Rscript
# Supplementary Figure 3A,B: expression of candidate motif-associated RBP genes.
# Adapted from the recovered stromal and epithelial heatmap scripts.
# Usage: Rscript script.R [stroma|epithelium] [--tables-only]
# Inputs and outputs are relative to KC_WORK_DIR (the current directory by default).

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
source(file.path(dirname(normalizePath(sub("^--file=", "", arg))), "Utils_plotting.R"))
args <- setdiff(commandArgs(trailingOnly = TRUE), "--tables-only")
if (length(args) > 1L || (length(args) == 1L && !args %in% c("stroma", "epithelium"))) {
    stop("Usage: Rscript script.R [stroma|epithelium] [--tables-only]")
}
cohorts <- if (length(args)) args else c("stroma", "epithelium")

read_candidates <- function(cohort) {
    aliases <- c("HuR" = "ELAVL1", "SF2-ASF" = "SRSF1")
    read_direction <- function(direction) {
        path <- file.path(ROOT, "data", cohort, "rmaps",
                          paste0("pVal.", direction, ".vs.bg.RNAmap.txt"))
        dat <- read_tsv(path)
        if (!identical(dim(dat), c(114L, 9L)) || names(dat)[1] != "RBP") {
            stop("Expected 114 motifs and eight regional P values: ", path)
        }
        p <- as.matrix(dat[, -1, drop = FALSE])
        storage.mode(p) <- "double"
        if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid motif P values: ", path)
        # Correct all 912 tests together, separately within each tissue/direction.
        q <- matrix(p.adjust(as.vector(p), method = "BH"), nrow = nrow(p),
                    ncol = ncol(p), dimnames = dimnames(p))
        keep <- apply(q < 5e-5, 1, any)
        rbp_names <- sub("\\..*$", "", dat$RBP[keep])
        symbols <- ifelse(rbp_names %in% names(aliases), unname(aliases[rbp_names]), rbp_names)
        data.frame(gene_symbol = symbols, rmaps_name = rbp_names, motif = dat$RBP[keep],
                   direction = direction, min_BH_FDR = apply(q[keep, , drop = FALSE], 1, min),
                   stringsAsFactors = FALSE)
    }
    hits <- rbind(read_direction("up"), read_direction("dn"))
    metadata <- do.call(rbind, lapply(split(hits, hits$gene_symbol), function(x) {
        directions <- sort(unique(x$direction))
        data.frame(gene_symbol = x$gene_symbol[1],
                   rmaps_name = paste(sort(unique(x$rmaps_name)), collapse = ";"),
                   rmaps_direction = if (all(c("up", "dn") %in% directions)) "both" else directions[1],
                   motif_count = nrow(x), motifs = paste(sort(unique(x$motif)), collapse = ";"),
                   min_BH_FDR = min(x$min_BH_FDR), stringsAsFactors = FALSE)
    }))
    rownames(metadata) <- NULL
    metadata[order(metadata$gene_symbol), , drop = FALSE]
}

run_heatmap <- function(cohort) {
    panel <- if (cohort == "stroma") "SuppFigure3A" else "SuppFigure3B"
    control_count <- if (cohort == "stroma") 4L else 7L
    case_count <- if (cohort == "stroma") 5L else 7L
    expected_genes <- if (cohort == "stroma") 19L else 23L
    sample_labels <- c(paste0("Control ", seq_len(control_count)),
                       paste0("Keratoconus ", seq_len(case_count)))
    metadata <- read_candidates(cohort)
    abundance <- read_tsv(file.path(ROOT, "data", cohort, "rbp", "tpm.tsv"))
    if (!setequal(names(abundance), c("gene_id", "gene_symbol", sample_labels)) ||
        anyDuplicated(names(abundance)) || anyDuplicated(abundance$gene_symbol) ||
        anyDuplicated(abundance$gene_id)) {
        stop("Expected unique gene IDs/symbols and the study sample columns in rbp/tpm.tsv")
    }
    index <- match(metadata$gene_symbol, abundance$gene_symbol)
    if (nrow(metadata) != expected_genes || anyNA(index)) {
        stop("Candidate selection or expression coverage differs from the manuscript: ", cohort)
    }
    metadata$gene_id <- abundance$gene_id[index]
    if (anyNA(metadata$gene_id) || any(metadata$gene_id == "")) stop("Missing Ensembl IDs")
    tpm <- as.matrix(abundance[index, sample_labels, drop = FALSE])
    storage.mode(tpm) <- "double"
    rownames(tpm) <- metadata$gene_symbol
    if (any(!is.finite(tpm)) || any(tpm < 0)) stop("TPM must be finite and nonnegative")

    # Preserve the recovered scripts' transform. Row Z-scores are invariant to
    # the logarithm base, so this also represents row-scaled log10(TPM + 1).
    log2_tpm <- log2(tpm + 1)
    row_zscore <- t(scale(t(log2_tpm)))
    if (any(!is.finite(row_zscore))) stop("Zero-variance or non-finite RBP expression row")
    row_dendrogram <- as.dendrogram(hclust(dist(row_zscore), method = "ward.D2"))
    out_dir <- file.path(ROOT, "results", cohort, "suppfigure3")
    write_tsv(metadata, file.path(out_dir, paste0(panel, "_metadata.tsv")))
    for (name in c("TPM", "log2TPM", "row_zscore")) {
        matrix <- switch(name, TPM = tpm, log2TPM = log2_tpm, row_zscore = row_zscore)
        write_tsv(data.frame(gene_symbol = rownames(matrix), matrix, check.names = FALSE),
                  file.path(out_dir, paste0(panel, "_", name, ".tsv")))
    }
    write_tsv(data.frame(dendrogram_order = seq_len(nrow(tpm)),
                        gene_symbol = labels(row_dendrogram)),
              file.path(out_dir, paste0(panel, "_row_order.tsv")))
    if (TABLES_ONLY) {
        message(panel, ": checked ", nrow(tpm), " genes and ", ncol(tpm), " samples; no images created.")
        return(invisible(NULL))
    }
    if (!requireNamespace("gplots", quietly = TRUE)) stop("Install gplots to draw the heatmap")

    sample_group <- c(rep("Control", control_count), rep("Keratoconus", case_count))
    group_colors <- c(Control = "#222222",
                      Keratoconus = if (cohort == "stroma") "#E69F00" else "#CC79A7")
    heatmap_colors <- grDevices::colorRampPalette(
        c("#2166AC", "#67A9CF", "#D9D9D9", "#EF8A62", "#B2182B"))(256)
    draw_heatmap <- function() {
        gplots::heatmap.2(
            row_zscore, Rowv = row_dendrogram, Colv = FALSE, dendrogram = "row",
            scale = "none", col = heatmap_colors, breaks = seq(-2.5, 2.5, length.out = 257),
            symkey = FALSE, symbreaks = FALSE, density.info = "none", trace = "none",
            key = TRUE, key.title = "", key.xlab = "Row Z-score",
            key.par = list(cex = 1.22, mar = c(4.35, 1, 1, 0.5), mgp = c(2, 0.8, 0)),
            key.xtickfun = function() list(at = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                          labels = c("-2", "-1", "0", "1", "2"), gap.axis = -1),
            keysize = 1.05, ColSideColors = unname(group_colors[sample_group]),
            colsep = control_count, sepcolor = "#666666", sepwidth = c(0.035, 0.035),
            margins = c(9.5, 10), cexRow = 1.65, cexCol = 1.35,
            srtCol = 45, adjCol = c(1, 1), main = NULL,
            lhei = c(1.45, 5.05), lwid = c(1.5, 5))
        legend("topright", legend = names(group_colors), fill = unname(group_colors),
               border = NA, bty = "n", cex = 1.35, inset = c(-0.01, -0.08), xpd = NA)
    }
    width <- if (cohort == "stroma") 2600 else 3000
    height <- if (cohort == "stroma") 2100 else 2200
    # Match the recovered figure dimensions; the shared helper supplies the font.
    png(file.path(out_dir, paste0(panel, ".png")), width = width, height = height,
        res = 300, bg = "white")
    tryCatch({par(family = FONT); draw_heatmap()}, finally = dev.off())
    pdf_device(file.path(out_dir, paste0(panel, ".pdf")), width / 300, height / 300)
    tryCatch(draw_heatmap(), finally = dev.off())
    message(panel, ": wrote PNG and PDF heatmaps.")
}

for (cohort in cohorts) run_heatmap(cohort)
