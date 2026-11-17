# rMAPS3 motif-enrichment figures for corneal stroma and epithelium
#
# Colour roles:
#   red / blue       = up / down direction
#   orange / pink    = stroma / epithelium (label and target exon only)
#   colour intensity = -log10(BH-FDR)

library(ggplot2)
library(reshape2)
library(patchwork)

fdr_cutoff <- 5e-5
score_cap <- 8

figure_dir <- "results/motif_figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

tissue_config <- list(
  stroma = list(
    dataset = "stroma",
    label = "STROMA",
    accent = "#E69F00"
  ),
  epithelium = list(
    dataset = "epithelium",
    label = "EPITHELIUM",
    accent = "#CC79A7"
  )
)

direction_config <- list(
  up = list(
    title = "Up-regulated splicing events",
    low = "#E6A09A",
    high = "#9F1D20"
  ),
  dn = list(
    title = "Down-regulated splicing events",
    low = "#9AB9D0",
    high = "#174A73"
  )
)

window_size_labels <- c(
  "50 bp", "250 bp", "250 bp", "50 bp",
  "50 bp", "250 bp", "250 bp", "50 bp"
)
cutoff_label <- format(fdr_cutoff, scientific = TRUE, trim = TRUE)

read_adjusted_matrix <- function(dataset, direction) {
  input_file <- file.path(
    "data",
    dataset,
    "rmaps",
    paste0("pVal.", direction, ".vs.bg.RNAmap.txt")
  )
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  input <- read.delim(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  motif <- input[[1]]
  p_matrix <- as.matrix(input[-1])
  storage.mode(p_matrix) <- "double"
  rownames(p_matrix) <- motif

  if (anyNA(p_matrix) || any(p_matrix < 0 | p_matrix > 1)) {
    stop("All p-values must be non-missing values between 0 and 1.")
  }

  q_matrix <- matrix(
    p.adjust(as.vector(p_matrix), method = "BH"),
    nrow = nrow(p_matrix),
    ncol = ncol(p_matrix),
    dimnames = dimnames(p_matrix)
  )
  keep <- apply(q_matrix < fdr_cutoff, 1, any)
  if (!any(keep)) {
    stop("No motif passed FDR < ", fdr_cutoff, " for ", dataset, " ", direction)
  }

  q_filtered <- q_matrix[keep, , drop = FALSE]
  q_filtered[order(apply(q_filtered, 1, min)), , drop = FALSE]
}

make_schematic <- function(tissue) {
  tissue_style <- tissue_config[[tissue]]

  exon_data <- data.frame(
    xmin = c(0.5, 3.5, 7.5),
    xmax = c(1.5, 5.5, 8.5),
    ymin = 0.22,
    ymax = 0.53,
    fill = c("white", tissue_style$accent, "white")
  )

  ggplot() +
    geom_segment(
      aes(x = 1.5, xend = 7.5, y = 0.375, yend = 0.375),
      linewidth = 0.55,
      color = "black"
    ) +
    geom_rect(
      data = exon_data,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      color = "black",
      linewidth = 0.65
    ) +
    scale_fill_identity() +
    geom_segment(
      data = data.frame(x = seq(0.5, 8.5, by = 1)),
      aes(x = x, xend = x, y = 0.22, yend = -0.02),
      linetype = "22",
      linewidth = 0.45,
      color = "#333333"
    ) +
    geom_segment(
      aes(x = 2.34, xend = 2.55, y = 0.27, yend = 0.48),
      linewidth = 0.7,
      color = "black"
    ) +
    geom_segment(
      aes(x = 2.48, xend = 2.69, y = 0.27, yend = 0.48),
      linewidth = 0.7,
      color = "black"
    ) +
    geom_segment(
      aes(x = 6.34, xend = 6.55, y = 0.27, yend = 0.48),
      linewidth = 0.7,
      color = "black"
    ) +
    geom_segment(
      aes(x = 6.48, xend = 6.69, y = 0.27, yend = 0.48),
      linewidth = 0.7,
      color = "black"
    ) +
    annotate(
      "text",
      x = seq_len(8),
      y = 0.70,
      label = window_size_labels,
      family = "Arial",
      size = 3.05,
      color = "black"
    ) +
    scale_x_continuous(limits = c(0.5, 8.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.03, 0.84), expand = c(0, 0)) +
    theme_void(base_family = "Arial") +
    theme(
      plot.margin = margin(t = 2, r = 5, b = 0, l = 5)
    )
}

make_heatmap <- function(direction, q_filtered) {
  direction_style <- direction_config[[direction]]
  motif_order <- rownames(q_filtered)
  score_matrix <- pmin(
    -log10(pmax(q_filtered, .Machine$double.xmin)),
    score_cap
  )
  score_matrix[q_filtered >= fdr_cutoff] <- NA_real_

  plot_data <- melt(
    score_matrix,
    varnames = c("RBP_motif", "Region"),
    value.name = "score"
  )
  plot_data$Region_index <- match(
    as.character(plot_data$Region),
    colnames(score_matrix)
  )
  plot_data$RBP_motif <- factor(
    plot_data$RBP_motif,
    levels = rev(motif_order)
  )

  ggplot(
    plot_data,
    aes(x = Region_index, y = RBP_motif, fill = score)
  ) +
    geom_tile(color = "#606060", linewidth = 0.32) +
    scale_x_continuous(
      limits = c(0.5, 8.5),
      breaks = seq_len(8),
      expand = c(0, 0)
    ) +
    scale_fill_gradientn(
      colours = c(direction_style$low, direction_style$high),
      limits = c(-log10(fdr_cutoff), score_cap),
      breaks = c(4.5, 6, 8),
      na.value = "white",
      name = expression(-plain(Log)[10]~" ("~italic(P)~" adjust)")
    ) +
    guides(
      fill = guide_colorbar(
        barheight = grid::unit(28, "mm"),
        barwidth = grid::unit(4, "mm"),
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = "Arial", base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(
        face = "plain",
        size = 12,
        color = "black",
        margin = margin(r = 5)
      ),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 0, r = 5, b = 8, l = 5)
    )
}

make_panel <- function(tissue, direction, q_filtered, heatmap_height) {
  schematic <- make_schematic(tissue)
  heatmap <- make_heatmap(direction, q_filtered)
  panel <- schematic / heatmap +
    plot_layout(heights = c(0.95, heatmap_height))
  list(plot = panel, n_motifs = nrow(q_filtered))
}

make_tissue_figure <- function(tissue) {
  dataset <- tissue_config[[tissue]]$dataset
  up_q <- read_adjusted_matrix(dataset, "up")
  dn_q <- read_adjusted_matrix(dataset, "dn")
  common_heatmap_height <- max(4.2, 0.42 * max(nrow(up_q), nrow(dn_q)))
  up_panel <- make_panel(tissue, "up", up_q, common_heatmap_height)
  dn_panel <- make_panel(tissue, "dn", dn_q, common_heatmap_height)

  combined_plot <- (wrap_elements(full = up_panel$plot) |
    wrap_elements(full = dn_panel$plot)) +
    plot_layout(widths = c(1, 1)) &
    theme(plot.background = element_rect(fill = "white", color = NA))

  output_dir <- file.path(figure_dir, tissue)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_stem <- file.path(
    output_dir,
    paste0("rMAPS_splicing_map_", tissue, "_horizontal_BH_FDR", cutoff_label)
  )

  # Same landscape page size for both tissues to facilitate figure assembly.
  ggsave(
    paste0(output_stem, ".png"),
    combined_plot,
    width = 18,
    height = 6.2,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    paste0(output_stem, ".pdf"),
    combined_plot,
    width = 18,
    height = 6.2,
    device = cairo_pdf,
    bg = "white"
  )

  message(
    tissue_config[[tissue]]$label,
    ": UP ", up_panel$n_motifs,
    " motifs; DOWN ", dn_panel$n_motifs, " motifs"
  )
  message("Saved: ", output_stem, ".png/.pdf")
}

invisible(lapply(c("stroma", "epithelium"), make_tissue_figure))
