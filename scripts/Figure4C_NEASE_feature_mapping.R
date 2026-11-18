#!/usr/bin/env Rscript
# Figure4C: NEASE event coverage and overlapping protein-feature classes.
# Recovered layout: plot_overlap_extra_large_text.R.
# Counts and labels are calculated from SE_events.tsv; --tables-only creates no images.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "Utils_plotting.R"))
source(file.path(script_dir, "Utils_NEASE_feature_summary.R"))
summary <- read_nease_feature_summary("stroma", "Figure4C")
total <- summary$total
mapped <- summary$mapped
counts <- summary$counts
pairwise <- summary$pairwise
assignments <- summary$assignments
centers <- cumsum(counts) - counts / 2

draw <- function() {
  par(mar = c(0, 0, 0, 0), family = FONT, xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 10), ylim = c(0, 5.15), asp = 1)
  x0 <- .20
  x1 <- 9.80
  xt <- function(v) x0 + (x1 - x0) * v / total
  xb <- function(v) x0 + (x1 - x0) * v / assignments
  gold <- "#E3A000"
  grey <- "#929DA6"

  # Preserve the original bar lengths, heights, colors and expansion outline.
  rect(x0, 3.48, xt(mapped), 4.77, col = gold, border = "white", lwd = 1.2)
  rect(xt(mapped), 3.48, x1, 4.77, col = grey, border = "white", lwd = 1.2)
  text(mean(c(x0, xt(mapped))), 4.125,
       sprintf("Affecting protein\nfeatures\nn = %d (%.1f%%)", mapped, 100 * mapped / total), cex = 2.0)
  text(mean(c(xt(mapped), x1)), 4.125,
       sprintf("No mapped\nprotein feature\nn = %d (%.1f%%)", total - mapped, 100 * (total - mapped) / total), cex = 2.0)
  polygon(c(x0, xt(mapped), x1, x0), c(3.39, 3.39, 1.65, 1.65),
          col = "#FFF8E6", border = gold, lwd = .8)


  starts <- c(0, head(cumsum(counts), -1))
  ends <- cumsum(counts)
  colors <- c(gold, "#F6C248", "#FFE6A6")
  for (i in 1:3) {
    rect(xb(starts[i]), .36, xb(ends[i]), 1.56,
         col = colors[i], border = "white", lwd = 1.1)
  }

  text(xb(centers[1]), 1.29, "Domains", cex = 2.0)
  text(xb(centers[1]), .94, sprintf("n = %d (%.1f%%)", counts[1], 100 * counts[1] / assignments), cex = 2.0)
  text(xb(centers[1]), .56,
       sprintf("Shared with: Residues (%d), Linear motifs (%d)", pairwise[1], pairwise[2]), cex = 2.0)

  # The narrow segment retains its title and count; the shared counts use
  # the existing expansion area, connected to the residue segment.
  text(xb(centers[2]), 1.29, "Residues", cex = 2.0)
  text(xb(centers[2]), .94, sprintf("n = %d", counts[2]), cex = 2.0)
  text(xb(centers[2]), .56, sprintf("(%.1f%%)", 100 * counts[2] / assignments), cex = 2.0)
  text(4.90, 2.30, sprintf("Shared with Domains (%d)", pairwise[1]), cex = 2.0)
  text(4.90, 1.92, sprintf("and Linear motifs (%d)", pairwise[3]), cex = 2.0)
  arrows(6.80, 1.90, xb(centers[2]), 1.49, length = .08, angle = 25,
         col = "#666666", lwd = 1.0)

  text(8.50, 3.30, "Linear motifs", cex = 2.0)
  text(8.50, 2.97, sprintf("n = %d (%.1f%%)", counts[3], 100 * counts[3] / assignments), cex = 2.0)
  text(9.80, 2.65, sprintf("Shared: Domains (%d)", pairwise[2]), adj = 1, cex = 1.90)
  text(9.80, 2.33, sprintf("Residues (%d)", pairwise[3]), adj = 1, cex = 1.90)
  segments(9.40, 2.06, xb(centers[3]), 1.58, col = "#777777", lwd = .8)

  segments(x0, 4.97, x1, 4.97, lwd = .7)
  for (v in c(0, 200, 400, 600)) segments(xt(v), 4.97, xt(v), 5.02, lwd = .7)
  segments(x0, .19, x1, .19, lwd = .7)
  for (v in c(0, 100, 200, 300)) segments(xb(v), .19, xb(v), .14, lwd = .7)
}

save_nease_feature_panel(draw, summary$out_dir, "Figure4C", 2000, 1030)
