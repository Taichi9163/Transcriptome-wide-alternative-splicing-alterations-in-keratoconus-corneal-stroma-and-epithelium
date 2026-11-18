#!/usr/bin/env Rscript
# Scatter plot for the union of stromal and epithelial SE DAS events.  Counterpart values are
# retrieved by exact genomic coordinates.  Missing/unquantifiable events remain NA and are reported
# separately.  Significance categories describe independent dataset calls, not tissue specificity.


#
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "utils_paths.R"))

args <- script_arguments()
if (any(!args %in% c("--tables-only"))) stop("Unknown argument")
tables_only <- "--tables-only" %in% args
if (!tables_only) suppressPackageStartupMessages(library(ggplot2))

epithelium_raw_path <- work_path("data/epithelium/rmats/SE.MATS.JC.txt")
stroma_raw_path <- work_path("data/stroma/rmats/SE.MATS.JC.txt")
epithelium_das_path <- work_path("results/epithelium/das/SE.tsv")
stroma_das_path <- work_path("results/stroma/das/SE.tsv")
output_root <- work_path("results/cross_layer")
figures_dir <- file.path(output_root, "figures")
results_dir <- file.path(output_root, "results")
if (!tables_only) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

coordinate_columns <- c(
  "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES",
  "downstreamEE"
)
value_columns <- c("GeneID", "geneSymbol", "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference")

read_se_events <- function(path, label, retain_values = TRUE) {
  dat <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- coordinate_columns
  if (retain_values)
    required <- c(required, value_columns)
  missing_columns <- setdiff(required, colnames(dat))
  if (length(missing_columns) > 0) {
    stop(label, " is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  # Convert columns separately: apply(data.frame, ...) pads numeric columns when coercing them to a
  # character matrix, breaking external key lookups.
  dat$event_key <- do.call(paste, c(lapply(dat[, coordinate_columns, drop = FALSE], function(x) trimws(as.character(x))),
    sep = ":"
  ))
  if (anyDuplicated(dat$event_key)) {
    stop(label, " contains duplicated SE event keys.")
  }

  keep <- "event_key"
  if (retain_values) {
    for (group in c("IncLevel1", "IncLevel2")) {
      vals <- strsplit(ifelse(is.na(dat[[group]]), "NA", dat[[group]]), ",", fixed = TRUE)
      dat[[paste0(group, "_n_finite")]] <- vapply(vals, function(x) {
        sum(is.finite(suppressWarnings(as.numeric(x))))
      }, integer(1))
      dat[[paste0(group, "_n_samples")]] <- lengths(vals)
    }
    keep <- c(
      keep, value_columns, "IncLevel1_n_finite", "IncLevel2_n_finite", "IncLevel1_n_samples",
      "IncLevel2_n_samples"
    )
  }
  dat[, keep, drop = FALSE]
}

suffix_columns <- function(dat, suffix) {
  idx <- colnames(dat) != "event_key"
  colnames(dat)[idx] <- paste0(colnames(dat)[idx], "_", suffix)
  dat
}

epithelium_raw <- read_se_events(epithelium_raw_path, "Epithelium complete JC file")
stroma_raw <- read_se_events(stroma_raw_path, "Stroma complete JC file")
epithelium_das <- read_se_events(epithelium_das_path, "Epithelium DAS list", retain_values = FALSE)
stroma_das <- read_se_events(stroma_das_path, "Stroma DAS list", retain_values = FALSE)

if (!all(epithelium_das$event_key %in% epithelium_raw$event_key)) {
  stop("Some epithelium DAS events are absent from its complete JC file.")
}
if (!all(stroma_das$event_key %in% stroma_raw$event_key)) {
  stop("Some stroma DAS events are absent from its complete JC file.")
}

union_keys <- union(stroma_das$event_key, epithelium_das$event_key)
das_union <- data.frame(event_key = union_keys, stringsAsFactors = FALSE)
das_union$stroma_DAS <- das_union$event_key %in% stroma_das$event_key
das_union$epithelium_DAS <- das_union$event_key %in% epithelium_das$event_key
category_levels <- c("Significant in stroma only", "Significant in epithelium only", "Both significant")
das_union$DAS_category <- ifelse(das_union$stroma_DAS & das_union$epithelium_DAS, "Both significant",
  ifelse(das_union$stroma_DAS, "Significant in stroma only", "Significant in epithelium only")
)

stroma_raw$present_in_JC <- TRUE
epithelium_raw$present_in_JC <- TRUE
stroma_raw <- suffix_columns(stroma_raw, "stroma")
epithelium_raw <- suffix_columns(epithelium_raw, "epithelium")

das_union <- merge(das_union, stroma_raw, by = "event_key", all.x = TRUE, sort = FALSE)
das_union <- merge(das_union, epithelium_raw, by = "event_key", all.x = TRUE, sort = FALSE)
das_union$present_in_JC_stroma[is.na(das_union$present_in_JC_stroma)] <- FALSE
das_union$present_in_JC_epithelium[is.na(das_union$present_in_JC_epithelium)] <- FALSE

# Numeric zeroes in the source are retained; absent values are never imputed.
das_union$stroma_delta_inc_for_plot <- das_union$IncLevelDifference_stroma
das_union$epithelium_delta_inc_for_plot <- das_union$IncLevelDifference_epithelium
complete_pair <- das_union$present_in_JC_stroma & das_union$present_in_JC_epithelium & is.finite(das_union$stroma_delta_inc_for_plot) &
  is.finite(das_union$epithelium_delta_inc_for_plot)
for (tissue in c("stroma", "epithelium")) {
  for (group in c("IncLevel1", "IncLevel2")) {
    n <- das_union[[paste0(group, "_n_finite_", tissue)]]
    complete_pair <- complete_pair & !is.na(n) & n >= 1
  }
}
# This is a minimum estimability criterion, not a coverage/power guarantee.
das_union$included_in_scatter <- complete_pair
das_union$exclusion_reason <- ifelse(!das_union$present_in_JC_stroma, "Exact event absent from stroma SE.MATS.JC.txt",
  ifelse(!das_union$present_in_JC_epithelium, "Exact event absent from epithelium SE.MATS.JC.txt",
    ifelse(!complete_pair, "Nonfinite effect or no finite PSI in at least one group", NA_character_)
  )
)
plot_data <- das_union[complete_pair, , drop = FALSE]
unplottable <- das_union[!complete_pair, , drop = FALSE]
if (nrow(plot_data) == 0) stop("No DAS events have paired estimates.")
stopifnot(
  all(is.finite(plot_data$stroma_delta_inc_for_plot)), all(is.finite(plot_data$epithelium_delta_inc_for_plot)),
  nrow(plot_data) + nrow(unplottable) == nrow(das_union)
)

plot_data$DAS_category <- factor(plot_data$DAS_category, levels = category_levels)
plot_data <- plot_data[order(plot_data$DAS_category), , drop = FALSE]
union_category_counts <- table(factor(das_union$DAS_category, levels = category_levels))
plotted_category_counts <- table(factor(plot_data$DAS_category, levels = category_levels))

qa <- data.frame(
  event_type = "SE", stroma_DAS_events = nrow(stroma_das), epithelium_DAS_events = nrow(epithelium_das),
  both_DAS_exact_events = unname(union_category_counts["Both significant"]), DAS_event_union = nrow(das_union),
  measured_values_in_both_JC_files = nrow(plot_data), plotted_without_imputation = nrow(plot_data),
  excluded_events = nrow(unplottable), absent_from_stroma = sum(!das_union$present_in_JC_stroma), absent_from_epithelium = sum(!das_union$present_in_JC_epithelium),
  matched_but_unquantifiable = sum(das_union$present_in_JC_stroma & das_union$present_in_JC_epithelium &
    !complete_pair), stroma_only_DAS_union = unname(union_category_counts["Significant in stroma only"]),
  epithelium_only_DAS_union = unname(union_category_counts["Significant in epithelium only"]), stroma_only_plotted = unname(plotted_category_counts["Significant in stroma only"]),
  epithelium_only_plotted = unname(plotted_category_counts["Significant in epithelium only"]), both_DAS_plotted = unname(plotted_category_counts["Both significant"]),
  stringsAsFactors = FALSE
)

if (qa$stroma_DAS_events != 698 || qa$epithelium_DAS_events != 450 || qa$both_DAS_exact_events != 13 ||
  qa$DAS_event_union != 1135) {
  stop("DAS counts do not match the existing SE overlap analysis.")
}

write.csv(plot_data, file.path(results_dir, "SE_DAS_union_delta_inc_level_scatter_data.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(das_union, file.path(results_dir, "SE_DAS_union_all_events_audit.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(unplottable, file.path(results_dir, "SE_DAS_union_unplottable_events.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(qa, file.path(results_dir, "SE_DAS_union_delta_inc_level_scatter_QA.csv"),
  row.names = FALSE,
  quote = FALSE
)

if (tables_only) {
  print(qa)
  cat("Tables only; no figures generated.\n")
  quit(save = "no", status = 0)
}

category_colors <- c(
  `Significant in stroma only` = "#D98F27", `Significant in epithelium only` = "#8F5A8F",
  `Both significant` = "#009E73"
)
category_shapes <- c(`Significant in stroma only` = 15, `Significant in epithelium only` = 17, `Both significant` = 18)
category_sizes <- c(`Significant in stroma only` = 1.9, `Significant in epithelium only` = 2.1, `Both significant` = 3)
category_alphas <- c(`Significant in stroma only` = 1, `Significant in epithelium only` = 1, `Both significant` = 1)
category_labels <- c(`Significant in stroma only` = paste0(
  "Significant in stroma only (n = ", unname(plotted_category_counts["Significant in stroma only"]),
  ")"
), `Significant in epithelium only` = paste0(
  "Significant in epithelium only (n = ", unname(plotted_category_counts["Significant in epithelium only"]),
  ")"
), `Both significant` = paste0(
  "Both significant (n = ", unname(plotted_category_counts["Both significant"]),
  ")"
))

p <- ggplot(plot_data, aes(x = stroma_delta_inc_for_plot, y = epithelium_delta_inc_for_plot)) +
  geom_hline(
    yintercept = 0,
    color = "#9EA4AA", linewidth = 0.55
  ) +
  geom_vline(xintercept = 0, color = "#9EA4AA", linewidth = 0.55) +
  geom_hline(yintercept = c(-0.1, 0.1), color = "#C7CBD0", linewidth = 0.45, linetype = "dashed") +
  geom_vline(xintercept = c(-0.1, 0.1), color = "#C7CBD0", linewidth = 0.45, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, color = "#4C5157", linewidth = 0.55, linetype = "dotted") +
  geom_point(aes(color = DAS_category, shape = DAS_category, size = DAS_category, alpha = DAS_category)) +
  scale_color_manual(values = category_colors, labels = category_labels, drop = FALSE, name = NULL) +
  scale_shape_manual(values = category_shapes, labels = category_labels, drop = FALSE, name = NULL) +
  scale_size_manual(values = category_sizes, guide = "none") +
  scale_alpha_manual(
    values = category_alphas,
    guide = "none"
  ) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE), shape = guide_legend(
    ncol = 1,
    byrow = TRUE
  )) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), expand = expansion(mult = 0)) +
  coord_equal(clip = "off") +
  labs(
    x = "Stroma IncLevelDifference", y = "Epithelium IncLevelDifference",
    caption = "Coordinate-matched SE events with estimates in both datasets; no missing-value imputation.\nColors denote separate within-dataset significance calls, not tissue specificity."
  ) +
  theme_minimal(base_family = "sans", base_size = 20) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major = element_line(
      color = "#ECEEF0",
      linewidth = 0.45
    ), axis.line = element_line(color = "#555B61", linewidth = 0.45), axis.ticks = element_line(
      color = "#555B61",
      linewidth = 0.4
    ), axis.text = element_text(color = "#33383D", size = 18), axis.title = element_text(
      color = "#20252A",
      size = 22
    ), legend.position = "top", legend.justification = "left", legend.box.just = "left", legend.title = element_text(
      face = "bold",
      size = 18
    ), legend.text = element_text(size = 17), plot.caption = element_text(size = 11, hjust = 0),
    plot.margin = margin(14, 28, 14, 14)
  )

png_path <- file.path(figures_dir, "SE_DAS_union_delta_inc_level_scatter.png")
pdf_path <- file.path(figures_dir, "SE_DAS_union_delta_inc_level_scatter.pdf")

ggsave(png_path, p, width = 9.2, height = 7.8, units = "in", dpi = 300, bg = "white")
ggsave(pdf_path, p, width = 9.2, height = 7.8, units = "in", device = cairo_pdf, bg = "white")

print(qa)
cat("Created:\n", png_path, "\n", pdf_path, "\n")
