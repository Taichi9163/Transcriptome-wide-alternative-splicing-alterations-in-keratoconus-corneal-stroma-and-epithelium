# Portable Venn reconstruction from event memberships. Areas are not to scale.  VennDiagram
# constructs grid objects with ind=FALSE; drawing is opt-in.
check_venn_package <- function() {
  if (!requireNamespace("VennDiagram", quietly = TRUE))
    stop("Install VennDiagram: install.packages('VennDiagram')")
}

make_four_venn <- function(data, title, show_genes = TRUE) {
  check_venn_package()
  regions <- four_set_regions(data)
  sets <- lapply(1:4, function(i) data$event_key[substr(data$membership, i, i) == "1"])
  arguments <- setNames(lapply(sets, length), paste0("area", 1:4))
  for (n in 2:4) {
    for (indices in combn(1:4, n, simplify = FALSE)) {
      arguments[[paste0("n", paste(indices, collapse = ""))]] <- length(Reduce(intersect, sets[indices]))
    }
  }
  venn <- do.call(VennDiagram::draw.quad.venn, c(arguments, list(
    category = c(
      "Stromal\ninclusion",
      "Stromal\nskipping", "Epithelial\ninclusion", "Epithelial\nskipping"
    ), fill = c(
      "#F4B183", "#ED7D31",
      "#C9B3DD", "#8F5A8F"
    ), alpha = rep(0.35, 4), col = c("#D98F27", "#D98F27", "#8F5A8F", "#8F5A8F"),
    lwd = 1.3, cex = 1.2, cat.cex = 0.85, fontfamily = FONT, cat.fontfamily = FONT, ind = FALSE
  )))
  shared <- regions[regions$membership %in% c("1010", "0101", "1001", "0110") & regions$count > 0, ]
  footer <- "Areas are not proportional to event counts."
  if (show_genes && nrow(shared)) {
    lines <- paste0(shared$sets, ": ", shared$genes)
    footer <- paste(c(
      vapply(lines, function(x) paste(strwrap(x, 95), collapse = "\n"), character(1)),
      footer
    ), collapse = "\n")
  }
  grid::grobTree(grid::textGrob(title, y = 0.97, gp = grid::gpar(
    fontface = "bold", fontsize = 14,
    fontfamily = FONT
  )), grid::grobTree(children = venn, vp = grid::viewport(
    y = 0.59, width = 0.93,
    height = 0.73
  )), grid::textGrob(footer, x = 0.04, y = 0.02, just = c("left", "bottom"), gp = grid::gpar(
    fontsize = 8,
    fontfamily = FONT, lineheight = 1.3
  )))
}

save_venn <- function(plot, stem, width = 8, height = 6) {
  if (TABLES_ONLY)
    return(invisible(plot))
  if (!requireNamespace("ragg", quietly = TRUE))
    stop("Install ragg to generate PNGs")
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  pdf_device(paste0(stem, ".pdf"), width, height)
  grid::grid.newpage()
  grid::grid.draw(plot)
  grDevices::dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 600)
  grid::grid.newpage()
  grid::grid.draw(plot)
  grDevices::dev.off()
}
