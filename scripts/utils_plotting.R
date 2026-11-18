# Shared paths and output helpers for the manuscript panels.
source(file.path(dirname(normalizePath(sys.frame(1)$ofile)), "utils_paths.R"))
FONT <- Sys.getenv("FIGURE_FONT", "sans")

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
  ggplot2::ggsave(paste0(stem, ".png"), plot,
    width = width, height = height, dpi = 600, device = ragg::agg_png,
    bg = "white"
  )
  ggplot2::ggsave(paste0(stem, ".tiff"), plot,
    width = width, height = height, dpi = 600, device = ragg::agg_tiff,
    compression = "lzw", bg = "white"
  )
}
