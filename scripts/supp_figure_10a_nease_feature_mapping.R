#!/usr/bin/env Rscript
# SuppFigure10A: NEASE event coverage and overlapping protein-feature classes.  Recovered layout:
# plot_epithelium_with_axis_labels.R.  Counts and labels are calculated from SE_events.tsv;
# --tables-only creates no images.  Paths: --work-dir DIR, then KC_WORK_DIR, then the current
# directory.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "utils_plotting.R"))
source(file.path(script_dir, "utils_nease_feature_summary.R"))
summary <- read_nease_feature_summary("epithelium", "SuppFigure10A")
total <- summary$total
mapped <- summary$mapped
counts <- summary$counts
pairwise <- summary$pairwise
assignments <- summary$assignments
centers <- cumsum(counts) - counts / 2

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = FONT, xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(-2.8, 10.08), ylim = c(-0.85, 5.92), asp = 1)
  x0 <- 0.2
  x1 <- 9.8
  xt <- function(v) x0 + (x1 - x0) * v / total
  xb <- function(v) x0 + (x1 - x0) * v / assignments
  pink <- "#CA78AA"
  grey <- "#919BA4"

  rect(x0, 3.48, xt(mapped), 4.77, col = pink, border = "white", lwd = 1.2)
  rect(xt(mapped), 3.48, x1, 4.77, col = grey, border = "white", lwd = 1.2)
  text(mean(c(x0, xt(mapped))), 4.125, sprintf(
    "Affecting protein\nfeatures\nn = %d (%.1f%%)", mapped,
    100 * mapped / total
  ), cex = 2)
  text(mean(c(xt(mapped), x1)), 4.125, sprintf("No mapped\nprotein feature\nn = %d (%.1f%%)", total -
    mapped, 100 * (total - mapped) / total), cex = 2)
  polygon(c(x0, xt(mapped), x1, x0), c(3.39, 3.39, 1.65, 1.65), col = "#F9F0F6", border = pink, lwd = 0.8)

  starts <- c(0, head(cumsum(counts), -1))
  ends <- cumsum(counts)
  colors <- c(pink, "#DDA0C5", "#EBBFDB")
  for (i in 1:3) {
    rect(xb(starts[i]), 0.36, xb(ends[i]), 1.56, col = colors[i], border = "white", lwd = 1.1)
  }

  text(xb(centers[1]), 1.29, "Domains", cex = 2)
  text(xb(centers[1]), 0.94, sprintf("n = %d (%.1f%%)", counts[1], 100 * counts[1] / assignments), cex = 2)
  text(xb(centers[1]), 0.56, sprintf(
    "Shared with: Residues (%d), Linear motifs (%d)", pairwise[1],
    pairwise[2]
  ), cex = 2)

  # Keep the full category name, n and percentage inside the original narrow segment, using the
  # same font family and separate lines.
  text(xb(centers[2]), 1.29, "Residues", cex = 1.35)
  text(xb(centers[2]), 0.94, sprintf("n = %d", counts[2]), cex = 1.7)
  text(xb(centers[2]), 0.56, sprintf("(%.1f%%)", 100 * counts[2] / assignments), cex = 1.5)
  text(4.9, 2.16, sprintf("Shared with Domains (%d)", pairwise[1]), cex = 2)
  arrows(6.9, 1.98, xb(centers[2]), 1.49, length = 0.08, angle = 25, col = "#666666", lwd = 1)

  text(8.5, 3.22, "Linear motifs", cex = 2)
  text(8.5, 2.85, sprintf("n = %d (%.1f%%)", counts[3], 100 * counts[3] / assignments), cex = 2)
  text(9.8, 2.56, sprintf("Shared: Domains (%d)", pairwise[2]), adj = 1, cex = 1.9)
  segments(9.4, 2.16, xb(centers[3]), 1.58, col = "#777777", lwd = 0.8)

  segments(x0, 4.97, x1, 4.97, lwd = 0.7)
  for (v in c(0, 100, 200, 300, 400)) {
    segments(xt(v), 4.97, xt(v), 5.02, lwd = 0.7)
    text(xt(v), 5.26, as.character(v), cex = 1.5)
  }
  segments(x0, 0.19, x1, 0.19, lwd = 0.7)
  for (v in c(0, 50, 100)) {
    segments(xb(v), 0.19, xb(v), 0.14, lwd = 0.7)
    text(xb(v), -0.1, as.character(v), cex = 1.5)
  }
  text(mean(c(x0, x1)), 5.69, "Number of DAS events", cex = 1.65)
  text(mean(c(x0, x1)), -0.57, "Number of event–feature class assignments", cex = 1.65)
  text(-0.18, 4.125, "DAS events", adj = 1, cex = 1.55, font = 2)
  text(-0.18, 0.96, "Affected protein\nfeatures, expanded", adj = 1, cex = 1.45, font = 2)
}

save_nease_feature_panel(draw, summary$out_dir, "SuppFigure10A", 2576, 1354)
