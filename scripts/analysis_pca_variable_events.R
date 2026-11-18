#!/usr/bin/env Rscript
# Compare predefined top-N sets ranked by RAW PSI variance across all samples.  No group,
# significance, effect-size or coverage criterion selects events.  Paths: --work-dir DIR, then
# KC_WORK_DIR, then the current directory.
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
out_dir <- work_path("results/pca_variable")
project_dir <- ROOT
baseline_dir <- file.path(dirname(out_dir), "pca_all")
for (d in c("results", "figures", "cache")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}
top_sizes <- c(500L, 1000L, 2000L, 5000L, 10000L)
set_names <- c(paste0("Top ", format(top_sizes, big.mark = ",", trim = TRUE)), "All events")
event_types <- c("SE", "RI", "A5SS", "A3SS", "MXE")
write_result <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, "results", filename), na = "NA")
}


manifest <- fread(file.path(baseline_dir, "results", "input_manifest.csv"))
input_md5 <- unname(tools::md5sum(file.path(project_dir, manifest$SourceFile)))
stopifnot(identical(input_md5, manifest$MD5))
write_result(manifest, "input_manifest.csv")

summary_rows <- list()
score_rows <- list()
variance_rows <- list()
composition_rows <- list()
validation_rows <- list()
plot_data <- list()
for (tissue in c("stroma", "epithelium")) {
  message("Processing ", tissue)
  baseline_file <- file.path(baseline_dir, "cache", paste0(tissue, "_PCA.rds"))
  baseline <- readRDS(baseline_file)
  samples <- baseline$sample_info
  n_kc <- sum(samples$PSI_column == "IncLevel1")
  n_control <- sum(samples$PSI_column == "IncLevel2")
  mat_list <- list()
  for (event in event_types) {
    path <- file.path(project_dir, baseline$parameters$input, paste0(event, ".MATS.JC.txt"))
    d <- fread(path, sep = "\t", check.names = TRUE, na.strings = "NA", data.table = FALSE)
    psi <- cbind(parse_psi(d$IncLevel1, n_kc), parse_psi(d$IncLevel2, n_control))
    rownames(psi) <- paste(event, d$ID, sep = "_")
    colnames(psi) <- samples$Sample
    psi <- psi[complete.cases(psi), , drop = FALSE]
    varying <- apply(psi, 1, function(v) any(v != v[1L]))
    mat_list[[event]] <- psi[varying, , drop = FALSE]
  }
  mat <- do.call(rbind, mat_list)
  stopifnot(identical(rownames(mat), rownames(baseline$pca$rotation)))
  stopifnot(identical(colnames(mat), rownames(baseline$pca$x)))
  # Compute sample variance BEFORE scaling events to unit SD.
  raw_variance <- apply(mat, 1, var)
  stopifnot(max(abs(raw_variance - baseline$pca$scale^2)) < 1e-12)
  metadata <- copy(as.data.table(baseline$event_metadata))
  stopifnot(identical(metadata$EventUID, rownames(mat)))
  metadata[, Raw_PSI_variance := raw_variance]
  # Equal variances are resolved deterministically by event ID, never labels.
  ranked_index <- order(-raw_variance, metadata$EventUID, method = "radix")
  ranking <- metadata[ranked_index]
  ranking[, Variance_rank := seq_len(.N)]
  setcolorder(ranking, c("Variance_rank", "Raw_PSI_variance", setdiff(names(ranking), c(
    "Variance_rank",
    "Raw_PSI_variance"
  ))))
  write_result(ranking, paste0(tissue, "_all_event_variance_ranks.csv.gz"))
  selected_sets <- c(top_sizes, nrow(mat))
  fits <- list()
  plot_data[[tissue]] <- list()
  for (k in seq_along(selected_sets)) {
    n <- selected_sets[k]
    set_name <- set_names[k]
    is_all <- k == length(selected_sets)
    # Use original order for the all-event baseline; ranking only defines top-N.
    selected <- if (is_all)
      seq_len(nrow(mat)) else ranked_index[seq_len(n)]
    x <- t(mat[selected, , drop = FALSE])
    fit <- if (is_all)
      baseline$pca else prcomp(x, center = TRUE, scale. = TRUE)
    if (!is_all) {
      for (j in seq_len(ncol(fit$x))) {
        # Align arbitrary PC signs to the same PC of the all-event reference.  No rotation,
        # reordering of PCs, or group-based orientation is applied.
        ref <- baseline$pca$x[, j]
        r <- if (sd(ref) > 1e-10 && sd(fit$x[, j]) > 1e-10)
          cor(ref, fit$x[, j]) else NA_real_
        flip <- if (is.finite(r) && abs(r) > 1e-08)
          r < 0 else {
          fit$rotation[which.max(abs(fit$rotation[, j])), j] < 0
        }
        if (flip) {
          fit$x[, j] <- -fit$x[, j]
          fit$rotation[, j] <- -fit$rotation[, j]
        }
      }
      fits[[paste0("top", n)]] <- fit
      write_result(ranking[seq_len(n)], paste0(tissue, "_top", n, "_selected_events.csv.gz"))
    }
    var_pct <- 100 * fit$sdev^2 / sum(fit$sdev^2)
    cutoff <- min(raw_variance[selected])
    score <- cbind(data.frame(Tissue = tissue, EventSet = set_name, N_events = n), samples, as.data.frame(fit$x))
    score_rows[[paste(tissue, n)]] <- score
    variance_rows[[paste(tissue, n)]] <- data.frame(
      Tissue = tissue, EventSet = set_name, N_events = n,
      PC = paste0("PC", seq_along(var_pct)), Variance_percent = var_pct, Cumulative_percent = cumsum(var_pct)
    )
    summary_rows[[paste(tissue, n)]] <- data.frame(
      Tissue = tissue, EventSet = set_name, N_events = n,
      Eligible_events = nrow(mat), KC_samples = n_kc, Control_samples = n_control, Minimum_raw_PSI_variance = cutoff,
      Median_raw_PSI_variance = median(raw_variance[selected]), Percent_total_raw_variance = 100 *
        sum(raw_variance[selected]) / sum(raw_variance), Events_tied_at_cutoff = sum(raw_variance ==
        cutoff), Tied_events_selected = sum(raw_variance[selected] == cutoff), PC1_percent = var_pct[1],
      PC2_percent = var_pct[2], PC1_PC2_percent = sum(var_pct[1:2])
    )
    composition <- metadata[selected, .N, by = EventType]
    composition_rows[[paste(tissue, n)]] <- cbind(data.table(
      Tissue = tissue, EventSet = set_name,
      N_events = n
    ), composition)
    z <- scale(x, center = fit$center, scale = fit$scale)
    reconstruction_error <- max(abs(z - fit$x %*% t(fit$rotation)))
    variance_error <- abs(sum(fit$sdev^2) - ncol(z))
    # Top-N sets are nested by construction; validate boundary and exact size.
    stopifnot(ncol(x) == n, !anyNA(x), !anyDuplicated(colnames(x)))
    if (!is_all)
      stopifnot(min(raw_variance[selected]) >= max(raw_variance[-selected]))
    stopifnot(reconstruction_error < 1e-08, variance_error < 1e-06 * n)
    validation_rows[[paste(tissue, n)]] <- data.frame(
      Tissue = tissue, EventSet = set_name, N_events = n,
      Max_reconstruction_error = reconstruction_error, Total_standardized_variance_error = variance_error
    )
    plot_data[[tissue]][[set_name]] <- list(
      scores = score, variance = var_pct, n = n, tissue = tissue,
      set_name = set_name
    )
    message(
      tissue, " ", set_name, ": PC1 ", round(var_pct[1], 3), "%; PC2 ", round(var_pct[2], 3),
      "%"
    )
  }
  saveRDS(list(pcas = fits, sample_info = samples, parameters = list(
    rank_by = "raw PSI sample variance across all samples, before scaling",
    tie_break = "EventUID ascending using radix order", top_sizes = top_sizes, center = TRUE, scale = TRUE,
    significance_filter = FALSE, sample_scope = baseline$parameters$sample_scope, baseline_path = baseline_file,
    baseline_md5 = unname(tools::md5sum(baseline_file))
  )), file.path(out_dir, "cache", paste0(
    tissue,
    "_top_variable_PCAs.rds"
  )), compress = "gzip")
}
write_result(rbindlist(summary_rows), "PCA_comparison_summary.csv")
write_result(rbindlist(score_rows, fill = TRUE), "PCA_scores_all_sets.csv")
write_result(rbindlist(variance_rows), "PCA_variance_all_sets.csv")
write_result(rbindlist(composition_rows), "event_type_composition.csv")
write_result(rbindlist(validation_rows), "numerical_checks.csv")

palette <- c(Control = "#222222", Keratoconus = "#E69F00")
panel <- function(a, tissue_title = FALSE) {
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
      box.padding = 0.35, point.padding = 0.2, min.segment.length = 0, segment.color = "#AAAAAA",
      segment.size = 0.25, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_color_manual(values = palette) +
    scale_shape_manual(values = c(Control = 16, Keratoconus = 17)) +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(title = paste0(if (tissue_title)
      paste0(if (a$tissue == "stroma")
        "Stroma" else "Epithelium", " | ") else "", a$set_name), subtitle = paste0(format(a$n, big.mark = ","), " events"), x = sprintf(
      "PC1 (%.2f%%)",
      a$variance[1]
    ), y = sprintf("PC2 (%.2f%%)", a$variance[2]), color = NULL, shape = NULL) +
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
        16,
        15, 12, 12
      )
    )
}
caption <- paste0(
  "Events ranked by raw PSI variance across all samples, before event-wise standardization. No significance filter.\n",
  "Separate PCAs and axis ranges. K/C labels = rMATS input positions, not verified donor IDs."
)
annotate <- function(p, title, subtitle) {
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
grids <- list()
for (tissue in names(plot_data)) {
  p <- wrap_plots(lapply(plot_data[[tissue]], panel), ncol = 3, guides = "collect")
  p <- annotate(p, paste0(if (tissue == "stroma")
    "Stroma" else "Epithelium", ": PCA by event variability"), if (tissue == "stroma")
    "Yamashita2026 | KC n = 5, Control n = 4 | KC6 excluded" else "Sun2022 | KC n = 7, Control n = 7")
  grids[[tissue]] <- p
  export(p, paste0(tissue, "_top_variable_comparison"), 15, 11.6)
}
focus_sets <- c("Top 1,000", "Top 5,000")
focus <- wrap_plots(c(lapply(plot_data$stroma[focus_sets], panel, tissue_title = TRUE), lapply(plot_data$epithelium[focus_sets],
  panel,
  tissue_title = TRUE
)), ncol = 2, guides = "collect")
focus <- annotate(focus, "Splicing PCA: top 1,000 and top 5,000 variable events", "Stroma: KC 5 / Control 4 (KC6 excluded) | Epithelium: KC 7 / Control 7")
export(focus, "top1000_top5000_both_tissues", 12, 12)
cairo_pdf(file.path(out_dir, "figures", "all_variants_both_tissues.pdf"), width = 15, height = 11.6)
for (p in grids) print(p)
dev.off()
writeLines(c(
  paste("Run date:", Sys.Date()), paste("Script:", script_file), paste("Script MD5:", unname(tools::md5sum(script_file))),
  capture.output(sessionInfo())
), file.path(out_dir, "results", "session_info.txt"))
message("Saved outputs to ", out_dir)
