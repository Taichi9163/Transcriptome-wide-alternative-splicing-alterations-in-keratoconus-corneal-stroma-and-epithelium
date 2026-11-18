# Shared paths and output helpers for the manuscript panels.  KC_WORK_DIR contains data/ and
# results/; it defaults to the current directory.
ROOT <- normalizePath(Sys.getenv("KC_WORK_DIR", getwd()), mustWork = TRUE)
FONT <- Sys.getenv("FIGURE_FONT", "sans")
TABLES_ONLY <- "--tables-only" %in% commandArgs(trailingOnly = TRUE)

read_tsv <- function(path) {
    read.delim(path, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

pdf_device <- function(filename, width, height, ...) {
    if (capabilities("cairo")) {
        grDevices::cairo_pdf(filename, width = width, height = height, family = FONT, ...)
    } else {
        grDevices::pdf(filename, width = width, height = height, family = "Helvetica", ...)
    }
}

save_panel <- function(plot, stem, width, height) {
    # Build layers without opening a graphics device in tables-only mode.
    if (TABLES_ONLY) {
        if (inherits(plot, "patchwork")) {
            for (component in as.list(plot)) invisible(ggplot2::ggplot_build(component))
        } else {
            invisible(ggplot2::ggplot_build(plot))
        }
        return(invisible(NULL))
    }
    dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = pdf_device)
    ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 600, device = ragg::agg_png,
        bg = "white")
    ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = width, height = height, dpi = 600, device = ragg::agg_tiff,
        compression = "lzw", bg = "white")
}
