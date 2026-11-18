#!/usr/bin/env Rscript
# Stroma, all 10 samples (including KC6): all events and raw-variance top-N PCA.
#
# Paths: --work-dir DIR > KC_WORK_DIR > current directory.

font_cache <- file.path(tempdir(), "font-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_file <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_file), "utils_paths.R"))
source(file.path(dirname(script_file), "utils_psi.R"))
out_dir <- work_path("results/pca_with_kc6")
project_dir <- ROOT
analysis_dir <- dirname(out_dir)
input_dir <- "data/stroma_with_kc6/rmats"
for (d in c("results", "figures", "cache")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, filename) fwrite(as.data.table(x), file.path(out_dir, "results", filename),
  na = "NA"
)
event_types <- c("SE", "RI", "A5SS", "A3SS", "MXE")
top_sizes <- c(500L, 1000L, 2000L, 5000L, 10000L)
set_names <- c(paste0("Top ", format(top_sizes, big.mark = ",", trim = TRUE)), "All events")
samples <- data.frame(Sample = c(paste0("KC_pos", 1:6), paste0("Control_pos", 1:4)), Label = c(paste0(
  "K",
  1:6
), paste0("C", 1:4)), Group = c(rep("Keratoconus", 6), rep("Control", 4)), PSI_column = c(rep(
  "IncLevel1",
  6
), rep("IncLevel2", 4)), Within_group_position = c(1:6, 1:4), stringsAsFactors = FALSE)

mat_list <- audit_list <- count_list <- manifest_list <- list()
for (event in event_types) {
  relative_file <- file.path(input_dir, paste0(event, ".MATS.JC.txt"))
  path <- file.path(project_dir, relative_file)
  d <- fread(path, sep = "\t", check.names = TRUE, na.strings = "NA", data.table = FALSE)
  stopifnot(!anyNA(d$ID), !anyDuplicated(d$ID))
  psi <- cbind(parse_psi(d$IncLevel1, 6L), parse_psi(d$IncLevel2, 4L))
  uid <- paste(event, d$ID, sep = "_")
  rownames(psi) <- uid
  colnames(psi) <- samples$Sample
  complete <- complete.cases(psi)
  varying <- rep(NA, nrow(psi))
  varying[complete] <- apply(psi[complete, , drop = FALSE], 1, function(v) any(v != v[1L]))
  selected <- complete & !is.na(varying) & varying
  raw_variance <- rep(NA_real_, nrow(psi))
  raw_variance[complete] <- apply(psi[complete, , drop = FALSE], 1, var)
  coordinate_cols <- setdiff(names(d)[seq_len(match("IJC_SAMPLE_1", names(d)) - 1L)], c(
    "ID", "ID.1",
    "GeneID", "geneSymbol", "chr", "strand"
  ))
  coords <- do.call(paste, c(d[, c("chr", "strand", coordinate_cols), drop = FALSE], sep = ":"))
  audit_list[[event]] <- data.table(
    EventType = event, EventUID = uid, SourceRow = seq_len(nrow(d)),
    rMATS_ID = d$ID, GeneID = d$GeneID, GeneSymbol = d$geneSymbol, CoordinateKey = paste(event, coords,
      sep = ":"
    ), FDR = d$FDR, DeltaPSI = d$IncLevelDifference, Raw_PSI_variance = raw_variance,
    Missing_PSI_samples = rowSums(is.na(psi)), Status = ifelse(!complete, "excluded_missing_PSI",
      ifelse(varying, "included", "excluded_zero_variance")
    )
  )
  mat_list[[event]] <- psi[selected, , drop = FALSE]
  count_list[[event]] <- data.table(
    EventType = event, Total_events = nrow(psi), Excluded_missing_PSI = sum(!complete),
    Complete_PSI_events = sum(complete), Excluded_zero_variance = sum(complete & !selected), PCA_eligible_events = sum(selected)
  )
  manifest_list[[event]] <- data.table(
    EventType = event, SourceFile = relative_file, Bytes = file.info(path)$size,
    MD5 = unname(tools::md5sum(path))
  )
}
mat <- do.call(rbind, mat_list)
audit <- rbindlist(audit_list)
metadata <- audit[Status == "included"]
stopifnot(identical(rownames(mat), metadata$EventUID), !anyNA(mat), !anyDuplicated(rownames(mat)))
counts <- rbindlist(count_list)
total <- counts[, lapply(.SD, sum), .SDcols = setdiff(names(counts), "EventType")]
total[, EventType := "ALL"]
write_result(rbindlist(list(counts, total), use.names = TRUE), "event_counts.csv")
write_result(rbindlist(manifest_list), "input_manifest.csv")
write_result(audit, "event_filter_audit.csv.gz")
raw_variance <- metadata$Raw_PSI_variance
ranked_index <- order(-raw_variance, metadata$EventUID, method = "radix")
ranking <- copy(metadata[ranked_index])
ranking[, Variance_rank := seq_len(.N)]
setcolorder(ranking, c("Variance_rank", "Raw_PSI_variance", setdiff(names(ranking), c(
  "Variance_rank",
  "Raw_PSI_variance"
))))
write_result(ranking, "all_event_variance_ranks.csv.gz")

# All-sample fit is computed anew; adding a sample changes complete cases and ranks.
all_fit <- prcomp(t(mat), center = TRUE, scale. = TRUE)
for (j in seq_len(ncol(all_fit$x))) {
  if (all_fit$rotation[which.max(abs(all_fit$rotation[, j])), j] < 0) {
    all_fit$x[, j] <- -all_fit$x[, j]
    all_fit$rotation[, j] <- -all_fit$rotation[, j]
  }
}
fits <- list(all = all_fit)
summary_rows <- score_rows <- variance_rows <- composition_rows <- validation_rows <- plot_data <- list()
for (k in seq_along(set_names)) {
  is_all <- k == length(set_names)
  n <- if (is_all)
    nrow(mat) else top_sizes[k]
  name <- set_names[k]
  idx <- if (is_all)
    seq_len(nrow(mat)) else ranked_index[seq_len(n)]
  x <- t(mat[idx, , drop = FALSE])
  fit <- if (is_all)
    all_fit else prcomp(x, center = TRUE, scale. = TRUE)
  if (!is_all) {
    for (j in seq_len(ncol(fit$x))) {
      ref <- all_fit$x[, j]
      r <- if (sd(ref) > 1e-10 && sd(fit$x[, j]) > 1e-10)
        cor(ref, fit$x[, j]) else NA_real_
      flip <- if (is.finite(r) && abs(r) > 1e-08)
        r < 0 else fit$rotation[which.max(abs(fit$rotation[, j])), j] < 0
      if (flip) {
        fit$x[, j] <- -fit$x[, j]
        fit$rotation[, j] <- -fit$rotation[, j]
      }
    }
    fits[[paste0("top", n)]] <- fit
    write_result(ranking[seq_len(n)], paste0("top", n, "_selected_events.csv.gz"))
    stopifnot(min(raw_variance[idx]) >= max(raw_variance[-idx]))
  }
  pct <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  scores <- cbind(data.frame(EventSet = name, N_events = n), samples, as.data.frame(fit$x))
  score_rows[[name]] <- scores
  variance_rows[[name]] <- data.table(
    EventSet = name, N_events = n, PC = paste0("PC", seq_along(pct)),
    Variance_percent = pct, Cumulative_percent = cumsum(pct)
  )
  cutoff <- min(raw_variance[idx])
  summary_rows[[name]] <- data.table(
    EventSet = name, N_events = n, Eligible_events = nrow(mat), KC_samples = 6L,
    Control_samples = 4L, Minimum_raw_PSI_variance = cutoff, Median_raw_PSI_variance = median(raw_variance[idx]),
    Events_tied_at_cutoff = sum(raw_variance == cutoff), Tied_events_selected = sum(raw_variance[idx] ==
      cutoff), PC1_percent = pct[1], PC2_percent = pct[2], PC1_PC2_percent = sum(pct[1:2])
  )
  composition_rows[[name]] <- cbind(data.table(EventSet = name, N_events = n), metadata[idx, .N, by = EventType])
  z <- scale(x, center = fit$center, scale = fit$scale)
  reconstruction_error <- max(abs(z - fit$x %*% t(fit$rotation)))
  variance_error <- abs(sum(fit$sdev^2) - ncol(z))
  stopifnot(ncol(x) == n, !anyNA(x), reconstruction_error < 1e-08, variance_error < 1e-06 * n)
  validation_rows[[name]] <- data.table(
    EventSet = name, Max_reconstruction_error = reconstruction_error,
    Total_standardized_variance_error = variance_error
  )
  plot_data[[name]] <- list(scores = scores, pct = pct, n = n, title = name)
  message(name, ": ", n, " events; PC1 ", round(pct[1], 3), "%; PC2 ", round(pct[2], 3), "%")
}
write_result(rbindlist(summary_rows), "PCA_comparison_summary.csv")
write_result(rbindlist(score_rows), "PCA_scores_all_sets.csv")
write_result(rbindlist(variance_rows), "PCA_variance_all_sets.csv")
write_result(rbindlist(composition_rows), "event_type_composition.csv")
write_result(rbindlist(validation_rows), "numerical_checks.csv")
saveRDS(list(pcas = fits, sample_info = samples, event_metadata = metadata, parameters = list(
  input = input_dir,
  sample_scope = "All samples; KC6 included; KC n=6 / Control n=4", rank_by = "raw PSI sample variance over all 10 samples, before scaling",
  top_sizes = top_sizes, tie_break = "EventUID ascending, radix", center = TRUE, scale = TRUE, significance_filter = FALSE,
  complete_cases_recomputed = TRUE
)), file.path(out_dir, "cache", "stroma_with_KC6_PCAs.rds"), compress = "gzip")

palette <- c(Control = "#222222", Keratoconus = "#E69F00")
panel <- function(a) {
  d <- a$scores
  d$Group <- factor(d$Group, levels = names(palette))
  ggplot(d, aes(PC1, PC2, color = Group, shape = Group)) +
    geom_hline(
      yintercept = 0, color = "#DDDDDD",
      linewidth = 0.3
    ) +
    geom_vline(xintercept = 0, color = "#DDDDDD", linewidth = 0.3) +
    geom_point(size = 2.9) +
    geom_text_repel(aes(label = Label),
      seed = 20261014, size = 2.9, family = "Arial", color = "#333333",
      box.padding = 0.4, point.padding = 0.2, min.segment.length = 0, segment.color = "#AAAAAA",
      segment.size = 0.25, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_color_manual(values = palette) +
    scale_shape_manual(values = c(Control = 16, Keratoconus = 17)) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(
      title = a$title,
      subtitle = paste0(format(a$n, big.mark = ","), " events"), x = sprintf("PC1 (%.2f%%)", a$pct[1]),
      y = sprintf("PC2 (%.2f%%)", a$pct[2]), color = NULL, shape = NULL
    ) +
    theme_classic(
      base_size = 11,
      base_family = "Arial"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(
        t = 5,
        b = 4
      )), plot.subtitle = element_text(size = 9, color = "#555555"), axis.text = element_text(color = "#222222"),
      axis.line = element_line(linewidth = 0.4), legend.position = "top", plot.margin = margin(
        18,
        16, 12, 12
      )
    )
}
annotate <- function(p, title, subtitle, caption) {
  p + plot_annotation(title = title, subtitle = subtitle, caption = caption, theme = theme(
    plot.title = element_text(
      family = "Arial",
      size = 19, face = "bold"
    ), plot.subtitle = element_text(family = "Arial", size = 11, color = "#555555"),
    plot.caption = element_text(family = "Arial", size = 9, color = "#555555", hjust = 0), plot.margin = margin(
      18,
      18, 18, 18
    )
  )) & theme(legend.position = "top")
}
export <- function(p, stem, width, height) {
  if (TABLES_ONLY)
    return(invisible(NULL))
  ggsave(file.path(out_dir, "figures", paste0(stem, ".pdf")), p,
    width = width, height = height, device = cairo_pdf,
    bg = "white"
  )
  ggsave(file.path(out_dir, "figures", paste0(stem, ".png")), p,
    width = width, height = height, dpi = 300,
    device = ragg::agg_png, bg = "white"
  )
}
caption <- paste0(
  "PSI complete cases and raw-variance ranks recomputed over all 10 samples; event-wise standardization; no significance filter.\n",
  "Separate PCAs and axis ranges. K/C labels = within-group rMATS input positions."
)
grid <- annotate(
  wrap_plots(lapply(plot_data, panel), ncol = 3, guides = "collect"), "Stroma: PCA including KC6",
  "Yamashita2026 | KC n = 6, Control n = 4 | All five event classes", caption
)
export(grid, "stroma_with_KC6_all_variants", 15, 11.6)

# Comparison panels reuse previous fits; no previous output is modified.
old_all_path <- file.path(analysis_dir, "pca_all/cache/stroma_PCA.rds")
old_top_path <- file.path(analysis_dir, "pca_variable/cache/stroma_top_variable_PCAs.rds")
old_all <- readRDS(old_all_path)
old_top <- readRDS(old_top_path)
write_result(data.table(ReferenceFile = c(old_all_path, old_top_path), MD5 = unname(tools::md5sum(c(
  old_all_path,
  old_top_path
)))), "comparison_reference_manifest.csv")
comparison_panels <- list()
comparison_summary <- list()
for (name in c("All events", "Top 1,000", "Top 5,000")) {
  old_fit <- if (name == "All events")
    old_all$pca else old_top$pcas[[if (name == "Top 1,000")
    "top1000" else "top5000"]]
  old_pct <- 100 * old_fit$sdev^2 / sum(old_fit$sdev^2)
  old_data <- list(
    scores = cbind(old_all$sample_info, as.data.frame(old_fit$x)), pct = old_pct, n = nrow(old_fit$rotation),
    title = paste0("KC6 excluded | ", name)
  )
  new_data <- plot_data[[name]]
  new_data$title <- paste0("KC6 included | ", name)
  comparison_panels[[paste(name, "excluded")]] <- panel(old_data)
  comparison_panels[[paste(name, "included")]] <- panel(new_data)
  comparison_summary[[name]] <- data.table(
    EventSet = name, Events_KC6_excluded = old_data$n, Events_KC6_included = new_data$n,
    PC1_excluded = old_pct[1], PC2_excluded = old_pct[2], PC1_included = new_data$pct[1], PC2_included = new_data$pct[2]
  )
}
write_result(rbindlist(comparison_summary), "KC6_inclusion_comparison.csv")
comparison <- annotate(
  wrap_plots(comparison_panels, ncol = 2, guides = "collect"), "Stroma: KC6 excluded versus included",
  "Excluded: KC 5 / Control 4 | Included: KC 6 / Control 4", paste0(
    "Separate complete-case universes and variance rankings; changes reflect both sample inclusion and event selection.\n",
    "Separate PCAs, axis ranges and arbitrary PC signs. K/C labels = within-group rMATS input positions."
  )
)
export(comparison, "stroma_KC6_inclusion_comparison", 12, 16)
writeLines(c(
  paste("Run date:", Sys.Date()), paste("Script:", script_file), paste("Script MD5:", unname(tools::md5sum(script_file))),
  capture.output(sessionInfo())
), file.path(out_dir, "results", "session_info.txt"))
message("Saved outputs to ", out_dir)
