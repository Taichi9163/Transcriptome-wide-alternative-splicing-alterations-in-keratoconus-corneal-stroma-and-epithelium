#!/usr/bin/env Rscript
# PCA of all PSI-complete rMATS JC events, independently within each tissue.
# No filtering on FDR, PValue, effect size, coverage, or variable-event rank.
font_cache <- file.path(tempdir(), "font-cache")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg))
out_dir <- "results/pca_all"
project_dir <- normalizePath(".")
for (d in c("results", "figures", "cache")) {
  dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
}
event_types <- c("SE", "RI", "A5SS", "A3SS", "MXE")
configs <- list(
  stroma = list(title = "Stroma", cohort = "Yamashita2026", n_kc = 5L, n_control = 4L,
    input = "data/stroma/rmats",
    sample_scope = "KC6 excluded"),
  epithelium = list(title = "Epithelium", cohort = "Sun2022", n_kc = 7L, n_control = 7L,
    input = "data/epithelium/rmats",
    sample_scope = "All 14 samples")
)

# Names indicate within-group input position; no claim about donor/BAM identity.
# Group assignment follows the existing Heatmap_PCA scripts (IncLevel1 = KC).
parse_psi <- function(values, n_expected, context) {
  tokens <- strsplit(values, ",", fixed = TRUE)
  if (any(lengths(tokens) != n_expected)) stop("Unexpected PSI vector length: ", context)
  flat <- trimws(unlist(tokens, use.names = FALSE))
  parsed <- suppressWarnings(as.numeric(flat))
  missing_token <- is.na(flat) | flat %in% c("NA", "NaN", "nan", "")
  if (any(is.na(parsed) & !missing_token)) stop("Non-numeric PSI token: ", context)
  if (any(is.infinite(parsed))) stop("Infinite PSI: ", context)
  if (any(parsed < 0 | parsed > 1, na.rm = TRUE)) stop("PSI outside [0, 1]: ", context)
  matrix(parsed, ncol = n_expected, byrow = TRUE)
}

save_table <- function(x, filename) {
  fwrite(as.data.table(x), file.path(out_dir, "results", filename), na = "NA")
}

all_counts <- list()
all_manifest <- list()
analyses <- list()
for (tissue in names(configs)) {
  cfg <- configs[[tissue]]
  message("Processing ", tissue, " (", cfg$cohort, ")")
  sample_info <- data.frame(
    Sample = c(paste0("KC_pos", seq_len(cfg$n_kc)), paste0("Control_pos", seq_len(cfg$n_control))),
    Label = c(paste0("K", seq_len(cfg$n_kc)), paste0("C", seq_len(cfg$n_control))),
    Group = c(rep("Keratoconus", cfg$n_kc), rep("Control", cfg$n_control)),
    PSI_column = c(rep("IncLevel1", cfg$n_kc), rep("IncLevel2", cfg$n_control)),
    Within_group_position = c(seq_len(cfg$n_kc), seq_len(cfg$n_control)),
    stringsAsFactors = FALSE
  )
  mat_list <- list()
  audit_list <- list()
  for (event in event_types) {
    relative_file <- file.path(cfg$input, paste0(event, ".MATS.JC.txt"))
    input_file <- file.path(project_dir, relative_file)
    # rMATS has two columns both named ID. The first is its event identifier.
    dat <- fread(input_file, sep = "\t", header = TRUE, na.strings = "NA",
                 check.names = TRUE, data.table = FALSE)
    required <- c("ID", "GeneID", "geneSymbol", "chr", "strand", "IncLevel1", "IncLevel2", "FDR", "IncLevelDifference")
    stopifnot(all(required %in% names(dat)))
    if (anyNA(dat$ID) || anyDuplicated(dat$ID)) stop("Missing/duplicate event ID in ", relative_file)
    psi <- cbind(parse_psi(dat$IncLevel1, cfg$n_kc, paste(tissue, event, "IncLevel1")),
                 parse_psi(dat$IncLevel2, cfg$n_control, paste(tissue, event, "IncLevel2")))
    colnames(psi) <- sample_info$Sample
    uid <- paste(event, dat$ID, sep = "_")
    rownames(psi) <- uid
    complete <- complete.cases(psi) # na.omit on sample PSI, not on annotations/FDR.
    varying <- rep(NA, nrow(psi))
    varying[complete] <- apply(psi[complete, , drop = FALSE], 1, function(v) any(v != v[1L]))
    included <- complete & !is.na(varying) & varying
    status <- ifelse(!complete, "excluded_missing_PSI",
                     ifelse(varying, "included", "excluded_zero_variance"))
    # Keep exact coordinates for lookup; no gene-level collapsing or deduplication.
    coordinate_cols <- setdiff(names(dat)[seq_len(match("IJC_SAMPLE_1", names(dat)) - 1L)],
                               c("ID", "ID.1", "GeneID", "geneSymbol", "chr", "strand"))
    coord_key <- do.call(paste, c(dat[, c("chr", "strand", coordinate_cols), drop = FALSE],
                                 sep = ":"))
    audit_list[[event]] <- data.frame(
      Tissue = tissue, EventType = event, EventUID = uid, SourceRow = seq_len(nrow(dat)),
      rMATS_ID = dat$ID, GeneID = dat$GeneID, GeneSymbol = dat$geneSymbol,
      CoordinateKey = paste(event, coord_key, sep = ":"),
      FDR = dat$FDR, DeltaPSI = dat$IncLevelDifference,
      Missing_PSI_samples = rowSums(is.na(psi)), Status = status,
      stringsAsFactors = FALSE
    )
    mat_list[[event]] <- psi[included, , drop = FALSE]
    all_counts[[paste(tissue, event)]] <- data.frame(
      Tissue = tissue, Cohort = cfg$cohort, EventType = event,
      Total_events = nrow(psi), Excluded_missing_PSI = sum(!complete),
      Complete_PSI_events = sum(complete), Excluded_zero_variance = sum(complete & !included),
      PCA_events = sum(included), KC_samples = cfg$n_kc, Control_samples = cfg$n_control
    )
    all_manifest[[paste(tissue, event)]] <- data.frame(
      Tissue = tissue, EventType = event, SourceFile = relative_file,
      Bytes = file.info(input_file)$size, MD5 = unname(tools::md5sum(input_file)),
      stringsAsFactors = FALSE
    )
  }
  event_audit <- rbindlist(audit_list)
  event_matrix <- do.call(rbind, mat_list)
  stopifnot(!anyNA(event_matrix), all(is.finite(event_matrix)), !anyDuplicated(rownames(event_matrix)))
  stopifnot(identical(rownames(event_matrix), event_audit$EventUID[event_audit$Status == "included"]))
  # Center each event and scale it to unit sample standard deviation.
  # Samples are observations; events are variables. Groups never enter the fit.
  pca <- prcomp(t(event_matrix), center = TRUE, scale. = TRUE)
  # Fix arbitrary axis signs using the largest absolute loading (label-independent).
  for (j in seq_len(ncol(pca$rotation))) {
    anchor <- which.max(abs(pca$rotation[, j]))
    if (pca$rotation[anchor, j] < 0) {
      pca$rotation[, j] <- -pca$rotation[, j]
      pca$x[, j] <- -pca$x[, j]
    }
  }
  variance <- data.frame(
    PC = paste0("PC", seq_along(pca$sdev)), Eigenvalue = pca$sdev^2,
    Variance_percent = 100 * pca$sdev^2 / sum(pca$sdev^2),
    Cumulative_percent = 100 * cumsum(pca$sdev^2) / sum(pca$sdev^2)
  )
  scores <- cbind(Tissue = tissue, sample_info, as.data.frame(pca$x))
  # Meaningful numerical checks: standardized variance and full reconstruction.
  z <- scale(t(event_matrix), center = pca$center, scale = pca$scale)
  variance_error <- abs(sum(pca$sdev^2) - ncol(z))
  reconstruction_error <- max(abs(z - pca$x %*% t(pca$rotation)))
  stopifnot(variance_error < 1e-6 * ncol(z), reconstruction_error < 1e-8)
  qc <- data.frame(Tissue = tissue, Max_reconstruction_error = reconstruction_error,
                   Total_standardized_variance_error = variance_error,
                   PCA_events = ncol(z), Samples = nrow(z))
  save_table(scores, paste0(tissue, "_PCA_scores.csv"))
  save_table(variance, paste0(tissue, "_PCA_variance.csv"))
  save_table(event_audit, paste0(tissue, "_event_filter_audit.csv.gz"))
  save_table(qc, paste0(tissue, "_numerical_checks.csv"))
  loadings <- data.table(EventUID = rownames(pca$rotation), as.data.table(pca$rotation))
  save_table(loadings, paste0(tissue, "_PCA_loadings.csv.gz"))
  saveRDS(list(pca = pca, sample_info = sample_info,
               event_metadata = event_audit[Status == "included"],
               parameters = list(input = cfg$input, center = TRUE, scale = TRUE,
                 na_omit_columns = "sample PSI only", significance_filter = FALSE,
                 zero_variance_excluded = TRUE, sample_scope = cfg$sample_scope)),
          file.path(out_dir, "cache", paste0(tissue, "_PCA.rds")), compress = "gzip")
  analyses[[tissue]] <- list(cfg = cfg, scores = scores, variance = variance,
                            n_events = nrow(event_matrix))
  message(tissue, ": ", nrow(event_matrix), " events; PC1 ",
          round(variance$Variance_percent[1], 3), "%; PC2 ",
          round(variance$Variance_percent[2], 3), "%")
}

counts <- rbindlist(all_counts)
totals <- counts[, lapply(.SD, sum), by = .(Tissue, Cohort, KC_samples, Control_samples),
                 .SDcols = c("Total_events", "Excluded_missing_PSI", "Complete_PSI_events",
                             "Excluded_zero_variance", "PCA_events")]
totals[, EventType := "ALL"]
counts <- rbindlist(list(counts, totals), use.names = TRUE)
save_table(counts, "event_counts.csv")
save_table(rbindlist(all_manifest), "input_manifest.csv")

palette <- c(Control = "#222222", Keratoconus = "#E69F00")
plot_pca <- function(a, labels = TRUE) {
  plot_data <- a$scores
  plot_data$Group <- factor(plot_data$Group, levels = c("Control", "Keratoconus"))
  p <- ggplot(plot_data, aes(PC1, PC2, color = Group, shape = Group)) +
    geom_hline(yintercept = 0, color = "#DADADA", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "#DADADA", linewidth = 0.35) +
    geom_point(size = 3.5, stroke = 0.8) +
    scale_color_manual(values = palette) +
    scale_shape_manual(values = c(Control = 16, Keratoconus = 17)) +
    scale_x_continuous(expand = expansion(mult = 0.17)) +
    scale_y_continuous(expand = expansion(mult = 0.17)) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(title = paste0(a$cfg$title, " (", a$cfg$cohort, ")"),
         subtitle = sprintf("%s events  |  KC n = %d, Control n = %d%s",
           format(a$n_events, big.mark = ","), a$cfg$n_kc, a$cfg$n_control,
           if (a$cfg$title == "Stroma") "  |  KC6 excluded" else ""),
         x = sprintf("PC1 (%.2f%%)", a$variance$Variance_percent[1]),
         y = sprintf("PC2 (%.2f%%)", a$variance$Variance_percent[2]),
         color = NULL, shape = NULL) +
    theme_classic(base_size = 12, base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 10, color = "#555555", margin = margin(b = 10)),
          axis.line = element_line(linewidth = 0.45),
          axis.text = element_text(color = "#222222"),
          legend.position = "top", legend.justification = "left",
          plot.margin = margin(10, 18, 10, 12))
  if (labels) p <- p + geom_text_repel(aes(label = Label), seed = 20261014,
      size = 3.3, family = "Arial", color = "#333333", box.padding = 0.5,
      point.padding = 0.35, min.segment.length = 0, segment.color = "#AAAAAA",
      segment.size = 0.3, max.overlaps = Inf, show.legend = FALSE)
  p
}
export_plot <- function(p, name, width, height) {
  ggsave(file.path(out_dir, "figures", paste0(name, ".pdf")), p,
         width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, "figures", paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, device = ragg::agg_png, bg = "white")
}
for (tissue in names(analyses)) {
  single <- plot_pca(analyses[[tissue]], TRUE) +
    labs(caption = "K/C labels indicate within-group rMATS input positions;\nindividual donor IDs are not verified.") +
    theme(plot.caption = element_text(size = 9, color = "#555555", hjust = 0),
          plot.title = element_text(margin = margin(t = 12, b = 8)),
          plot.margin = margin(24, 18, 16, 12))
  export_plot(single, paste0(tissue, "_all_events_PCA_labeled"), 6.8, 7.0)
}
for (with_labels in c(TRUE, FALSE)) {
  combined <- (plot_pca(analyses$stroma, with_labels) +
               plot_pca(analyses$epithelium, with_labels)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Splicing PCA across all quantified event classes",
      caption = paste0("SE, RI, A5SS, A3SS and MXE; no significance filter. Complete PSI; constant events excluded; event-wise standardization.\n",
                       "Separate PCAs with independent axis ranges.",
                       if (with_labels) " Labels indicate within-group rMATS input positions, not verified donor IDs." else ""),
      theme = theme(plot.title = element_text(family = "Arial", size = 17, face = "bold"),
                    plot.caption = element_text(family = "Arial", size = 9, hjust = 0,
                                               color = "#555555"),
                    plot.margin = margin(14, 14, 14, 14))) &
    theme(legend.position = "top")
  export_plot(combined, if (with_labels) "both_tissues_all_events_PCA_labeled" else "both_tissues_all_events_PCA", 13.6, 6.6)
}

writeLines(c(paste("Run date:", Sys.Date()),
             paste("Script:", script_file),
             paste("Script MD5:", unname(tools::md5sum(script_file))),
             capture.output(sessionInfo())), file.path(out_dir, "results", "session_info.txt"))
message("Saved outputs to ", out_dir)
