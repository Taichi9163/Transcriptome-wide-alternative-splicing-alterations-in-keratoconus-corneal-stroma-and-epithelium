#!/usr/bin/env Rscript
# Figure 3C: structural event-matched DE genes nested within the gene-level intersection.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
  library(ragg)
})

args0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)

source(file.path(dirname(normalizePath(sub("^--file=", "", args0[1]))), "utils_plotting.R"))

out <- file.path(ROOT, "results/stroma/figure3")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read_tsv(file.path(ROOT, "results/stroma/manuscript/isoform_scatter_displayed_pairs.tsv"))

d$DE_isoform <- tolower(as.character(d$DE_isoform)) == "true"

m <- read_tsv(file.path(ROOT, "results/stroma/manuscript/isoform_gene_membership.tsv"))

iso <- read_tsv(file.path(ROOT, "data/stroma/expression/transcript_deseq2.tsv"))

matched <- unique(d$GeneID[d$DE_isoform])

das <- m$GeneID[m$DAS]

de <- m$GeneID[m$DE_isoform]

shared <- intersect(das, de)

stopifnot(length(das) == 807, length(de) == 410, length(shared) == 55, length(matched) == 29, all(matched %in%
  shared), nrow(d) == 4047, sum(d$DE_isoform) == 42, !anyDuplicated(paste(d$EventKey, d$TranscriptID)))

counts <- c(DAS_only = length(setdiff(das, de)), DE_gene_only = length(setdiff(de, das)), Shared_event_unconfirmed = length(setdiff(
  shared,
  matched
)), Event_matched = length(matched))

stopifnot(identical(unname(counts), c(752L, 355L, 26L, 29L)))

matched_iso <- unique(d$TranscriptID[d$DE_isoform])

top <- iso[iso$TranscriptID %in% matched_iso, ]

top$abs_stat <- abs(top$stat)

top <- top[order(-top$abs_stat, top$padj, top$TranscriptID), ]

top <- top[!duplicated(top$GeneID), ]

top$GeneSymbol <- m$GeneSymbol[match(top$GeneID, m$GeneID)]

top <- head(top, 10)

stopifnot(nrow(top) == 10, all(top$GeneID %in% matched))

write.table(top, file.path(out, "event_matched_top10_by_abs_STAT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

m$EventMatchedDE <- m$GeneID %in% matched

m$Category <- ifelse(m$EventMatchedDE, "Event-matched overlap", ifelse(m$DAS & m$DE_isoform, "Gene-level overlap; event correspondence unconfirmed",
  ifelse(m$DAS, "DAS only", "DE isoform gene only")
))

write.table(m, file.path(out, "Figure3C_gene_membership.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

labels <- ifelse(is.na(top$GeneSymbol) | !nzchar(top$GeneSymbol), top$GeneID, top$GeneSymbol)

circles <- data.frame(set = c("DAS", "DE"), x0 = c(-0.65, 0.65), y0 = 0, r = 1.05)

lx <- seq(-0.4, 0.4, length.out = 601)

lu <- pmin(sqrt(pmax(0, 1.05^2 - (lx + 0.65)^2)), sqrt(pmax(0, 1.05^2 - (lx - 0.65)^2)))

lens <- data.frame(x = c(lx, rev(lx)), y = c(lu, -rev(lu)))

inner <- data.frame(x0 = 0, y0 = -0.22, r = 0.3)

theta <- seq(0, 2 * pi, length.out = 1001)

stopifnot(
  all((inner$x0 + inner$r * cos(theta) + 0.65)^2 + (inner$y0 + inner$r * sin(theta))^2 < 1.05^2),
  all((inner$x0 + inner$r * cos(theta) - 0.65)^2 + (inner$y0 + inner$r * sin(theta))^2 < 1.05^2)
)

pC <- ggplot() +
  geom_circle(
    data = circles, aes(x0 = x0, y0 = y0, r = r, fill = set, color = set), alpha = 0.2,
    linewidth = 0.5, show.legend = FALSE
  ) +
  geom_polygon(data = lens, aes(x, y), fill = "#E5E5E5", color = NA) +
  geom_circle(
    data = inner, aes(x0 = x0, y0 = y0, r = r), fill = "#00C878", color = "#008E55", alpha = 0.65,
    linewidth = 0.5
  ) +
  annotate("text",
    x = c(-1.08, 1.08, 0, 0), y = c(0, 0, 0.48, -0.22), label = as.character(counts),
    family = FONT, size = 5
  ) +
  annotate("text",
    x = -0.93, y = 1.25, label = "DAS genes", family = FONT,
    size = 4.8
  ) +
  annotate("text",
    x = 0.93, y = 1.25, label = "Genes with\nDE isoforms", family = FONT,
    size = 4.8, lineheight = 1
  ) +
  annotate("point", x = -1.19, y = 1.74, color = "#00C878", size = 3.8) +
  annotate("text",
    x = -1.03, y = 1.74, label = "Event-matched DE isoform genes", hjust = 0, family = FONT,
    size = 4.1
  ) +
  annotate("point", x = -1.19, y = 1.48, color = "#BDBDBD", size = 3.8) +
  annotate("text",
    x = -1.03, y = 1.48, label = "Other gene-level overlap", hjust = 0, family = FONT, size = 4.1
  ) +
  annotate("point", x = 0, y = -0.47, color = "#111111", size = 1.8) +
  annotate("segment",
    x = 0, xend = 0,
    y = -0.47, yend = -1.16, color = "#444444", linewidth = 0.4
  ) +
  annotate("rect",
    xmin = -1.1, xmax = 1.1,
    ymin = -2.08, ymax = -1.16, fill = "white", color = "#555555", linewidth = 0.32
  ) +
  annotate("text",
    x = -0.55, y = -1.62, label = paste(labels[1:5], collapse = "\n"), family = FONT, fontface = "italic",
    size = 3.5, lineheight = 1
  ) +
  annotate("text",
    x = 0.55, y = -1.62, label = paste(labels[6:10], collapse = "\n"),
    family = FONT, fontface = "italic", size = 3.5, lineheight = 1
  ) +
  annotate("text",
    x = -1.88, y = 1.82,
    label = "C", hjust = 0, family = FONT, size = 7, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    DAS = "#FFD21A",
    DE = "#1565FF"
  )) +
  scale_color_manual(values = c(DAS = "#FFD21A", DE = "#1565FF")) +
  coord_fixed(xlim = c(
    -1.92,
    1.92
  ), ylim = c(-2.18, 1.96), clip = "off") +
  theme_void(base_family = FONT) +
  theme(plot.background = element_rect(
    fill = "white",
    color = NA
  ), plot.margin = margin(5, 5, 5, 5))

save_panel(pC, file.path(out, "Figure3C"), 6.2, 6.65)

print(counts)
