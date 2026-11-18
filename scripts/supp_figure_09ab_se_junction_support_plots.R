#!/usr/bin/env Rscript
# Supplementary Figure 9A,B: coverage-filtered cross-layer SE comparison.  The scatter
# specification is recovered from the 15 October 2026 revision; the Venn is a portable
# reconstruction calculated from exact event memberships.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "utils_plotting.R"))
source(file.path(script_dir, "utils_splicing_events.R"))
source(file.path(script_dir, "utils_venn_plots.R"))
if (any(!script_arguments() %in% "--tables-only")) stop("Unknown argument")
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(ROOT, "results/cross_layer/suppfigure9")
paired <- merge(se_support("stroma"), se_support("epithelium"), by = "event_key", sort = FALSE)
# Require finite PSI and >=10 total junction reads in at least half of EACH case/control group of
# EACH cohort (3/5, 2/4, 4/7 and 4/7 samples).
keep <- paired$estimable_stroma & paired$estimable_epithelium & paired$supported_stroma & paired$supported_epithelium
d <- paired[keep, , drop = FALSE]
stopifnot(nrow(d) == 15343L, all(is.finite(d$delta_stroma)), all(is.finite(d$delta_epithelium)))
lev <- c("Stroma only (n = 427)", "Epithelium only (n = 56)", "Both significant (n = 8)", "Neither significant (n = 14,852)")
d$category <- ifelse(d$significant_stroma & d$significant_epithelium, lev[3], ifelse(d$significant_stroma,
  lev[1], ifelse(d$significant_epithelium, lev[2], lev[4])
))
d$category <- factor(d$category, levels = lev)
stopifnot(identical(as.integer(table(d$category)), c(427L, 56L, 8L, 14852L)))
r <- cor(d$delta_stroma, d$delta_epithelium)
rho <- cor(d$delta_stroma, d$delta_epithelium, method = "spearman")
both <- d$significant_stroma & d$significant_epithelium
summary <- data.frame(Qualifying_events = nrow(d), Stroma_only = sum(d$category == lev[1]), Epithelium_only = sum(d$category ==
  lev[2]), Both = sum(both), Neither = sum(d$category == lev[4]), Shared_same_direction = sum(d$delta_stroma[both] *
  d$delta_epithelium[both] > 0), Shared_opposite_direction = sum(d$delta_stroma[both] * d$delta_epithelium[both] <
  0), Pearson_r = r, Spearman_rho = rho)
stopifnot(summary$Shared_same_direction == 5L, summary$Shared_opposite_direction == 3L)
write_tsv(summary, file.path(out, "SE_support_summary.tsv"))
write_tsv(d, file.path(out, "SE_support_filtered_events.tsv"))
venn_data <- four_set_memberships(d[d$significant_stroma | d$significant_epithelium, , drop = FALSE])
write_tsv(venn_data, file.path(out, "SE_significant_event_memberships.tsv"))
write_tsv(four_set_regions(venn_data), file.path(out, "SE_exclusive_regions.tsv"))
venn <- make_four_venn(venn_data, "SE events meeting the junction-support criterion")
save_venn(venn, file.path(out, "SuppFigure9A_SE_support_Venn"))

colors <- setNames(c("#D98F27", "#8F5A8F", "#009E73", "#AEB4BA"), lev)
shapes <- setNames(c(15, 17, 18, 16), lev)
sizes <- setNames(c(1.45, 1.6, 2.4, 0.7), lev)
alphas <- setNames(c(1, 1, 1, 0.28), lev)
# Render nonsignificant points first so significant events remain visible.
d <- d[order(match(d$category, c(lev[4], lev[1:3]))), ]
p <- ggplot(d, aes(delta_stroma, delta_epithelium)) +
  geom_hline(yintercept = 0, color = "#9EA4AA", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#9EA4AA", linewidth = 0.4) +
  geom_hline(yintercept = c(
    -0.1,
    0.1
  ), color = "#C7CBD0", linewidth = 0.35, linetype = "dashed") +
  geom_vline(xintercept = c(
    -0.1,
    0.1
  ), color = "#C7CBD0", linewidth = 0.35, linetype = "dashed") +
  geom_abline(
    intercept = 0, slope = 1,
    color = "#4C5157", linewidth = 0.4, linetype = "dotted"
  ) +
  geom_point(aes(
    color = category, shape = category,
    size = category, alpha = category
  )) +
  scale_color_manual(values = colors, name = NULL, drop = FALSE) +
  scale_shape_manual(values = shapes, name = NULL, drop = FALSE) +
  scale_size_manual(
    values = sizes,
    guide = "none"
  ) +
  scale_alpha_manual(values = alphas, guide = "none") +
  guides(color = guide_legend(
    ncol = 2,
    byrow = TRUE, override.aes = list(alpha = 1, size = c(2.3, 2.5, 3, 1.8))
  ), shape = guide_legend(
    ncol = 2,
    byrow = TRUE
  )) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), labels = sprintf(
    "%.1f",
    seq(-1, 1, 0.5)
  ), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(
    -1,
    1, 0.5
  ), labels = sprintf("%.1f", seq(-1, 1, 0.5)), expand = expansion(mult = 0)) +
  coord_equal(clip = "off") +
  labs(x = "Stroma IncLevelDifference", y = "Epithelium IncLevelDifference") +
  theme_minimal(
    base_family = FONT,
    base_size = 14
  ) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major = element_line(
      color = "#ECEEF0",
      linewidth = 0.35
    ), axis.line = element_line(color = "#555B61", linewidth = 0.35), axis.ticks = element_line(
      color = "#555B61",
      linewidth = 0.3
    ), axis.text = element_text(color = "#33383D", size = 10.125), axis.title = element_text(
      color = "#20252A",
      size = 10.935
    ), legend.position = "top", legend.justification = "center", legend.box.just = "center",
    legend.location = "plot", legend.text = element_text(size = 10.125), legend.key.width = grid::unit(
      0.32,
      "cm"
    ), legend.key.height = grid::unit(0.4, "cm"), legend.spacing.x = grid::unit(0.06, "cm"),
    plot.margin = margin(4, 5, 4, 4)
  )

# Build layers in numeric-only mode; do not open a graphics device.
if (TABLES_ONLY) {
  invisible(ggplot2::ggplot_build(p))
} else {
  ggsave(file.path(out, "SuppFigure9B_SE_support_scatter.png"), p,
    width = 4.8, height = 4.1, dpi = 300,
    device = ragg::agg_png, bg = "white"
  )
  ggsave(file.path(out, "SuppFigure9B_SE_support_scatter.pdf"), p,
    width = 4.8, height = 4.1, device = pdf_device,
    bg = "white"
  )
}
print(summary, row.names = FALSE)
