#!/usr/bin/env Rscript
# Figure 1D: gene-symbol UpSet with GeneID fallback, preserving the supplied analysis.
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(patchwork)
})

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]

source(file.path(dirname(normalizePath(sub("^--file=", "", arg))), "common.R"))

rmats_dir <- file.path(ROOT, "data/stroma/rmats")

out_dir <- file.path(ROOT, "results/stroma/figure1")

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

rmats_mode <- "JC"

top_y_max <- 500

setsize_x_max <- 600

fdr_cutoff <- 0.05

delta_psi_cutoff <- 0.1

gene_column_preference <- "geneSymbol"

event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")

event_order <- c("SE", "A5SS", "A3SS", "MXE", "RI")

as_colors <- c(SE = "#FFD700", A5SS = "#32CD32", A3SS = "#0000FF", MXE = "#85CEEB", RI = "#FF6347")

files <- setNames(file.path(rmats_dir, paste0(event_types, ".MATS.", rmats_mode, ".txt")), event_types)

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
    stop("Missing rMATS input files:\n", paste(missing_files, collapse = "\n"))
}

cat("All rMATS files were found.\n")

print(files)

read_filter_rmats <- function(file, event_type, fdr_cutoff = 0.05, delta_psi_cutoff = 0.1, gene_column_preference = "geneSymbol") {
    cat("Reading:", event_type, "\n")
    df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    colnames(df) <- make.unique(colnames(df))
    required_cols <- c("FDR", "IncLevelDifference")
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop(event_type, " file is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    if (!("GeneID" %in% colnames(df)) && !("geneSymbol" %in% colnames(df))) {
        stop(event_type, " file does not contain GeneID or geneSymbol column.")
    }
    df <- df %>% mutate(FDR = as.numeric(FDR), IncLevelDifference = as.numeric(IncLevelDifference))
    sig_df <- df %>% filter(!is.na(FDR), !is.na(IncLevelDifference), FDR < fdr_cutoff, abs(IncLevelDifference) >
        delta_psi_cutoff) %>% mutate(EventType = event_type)
    if (gene_column_preference %in% colnames(sig_df)) {
        gene_col <- gene_column_preference
    }
    else if ("geneSymbol" %in% colnames(sig_df)) {
        gene_col <- "geneSymbol"
    }
    else {
        gene_col <- "GeneID"
    }
    if (gene_col == "geneSymbol" && "GeneID" %in% colnames(sig_df)) {
        sig_df <- sig_df %>% mutate(gene_for_upset = ifelse(is.na(geneSymbol) | geneSymbol == "" |
            geneSymbol == "NA", GeneID, geneSymbol))
    }
    else {
        sig_df <- sig_df %>% mutate(gene_for_upset = .data[[gene_col]])
    }
    sig_df <- sig_df %>% filter(!is.na(gene_for_upset), gene_for_upset != "", gene_for_upset !=
        "NA")
    return(sig_df)
}

sig_event_list <- map2(files, names(files), ~read_filter_rmats(file = .x, event_type = .y, fdr_cutoff = fdr_cutoff,
    delta_psi_cutoff = delta_psi_cutoff, gene_column_preference = gene_column_preference))

names(sig_event_list) <- event_types

summary_table <- tibble(EventType = event_types, Significant_events = map_int(sig_event_list, nrow),
    Significant_genes = map_int(sig_event_list, ~length(unique(.x$gene_for_upset))))

print(summary_table)

write.csv(summary_table, file = file.path(out_dir, "rMATS_significant_AS_summary.csv"), row.names = FALSE)

for (event_type in event_types) {
    write.table(sig_event_list[[event_type]], file = file.path(out_dir, paste0(event_type, "_significant_AS_events.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE)
}

gene_list <- map(sig_event_list, ~unique(.x$gene_for_upset))

for (event_type in event_types) {
    gene_df <- tibble(Gene = gene_list[[event_type]])
    write.table(gene_df, file = file.path(out_dir, paste0(event_type, "_significant_AS_genes.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

gene_membership_long <- imap_dfr(gene_list, ~tibble(Gene = .x, EventType = .y))

write.table(gene_membership_long, file = file.path(out_dir, "Gene_AS_event_membership_long.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

gene_membership_wide <- gene_membership_long %>% mutate(value = 1) %>% pivot_wider(names_from = EventType,
    values_from = value, values_fill = 0)

for (ev in event_types) {
    if (!(ev %in% colnames(gene_membership_wide))) {
        gene_membership_wide[[ev]] <- 0
    }
}

gene_membership_wide <- gene_membership_wide %>% mutate(Number_of_AS_types = SE + A5SS + A3SS +
    MXE + RI) %>% arrange(desc(Number_of_AS_types), Gene)

write.table(gene_membership_wide, file = file.path(out_dir, "Gene_AS_event_membership_wide.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)

intersection_summary <- gene_membership_wide %>% select(Gene, all_of(event_types)) %>% group_by(across(all_of(event_types))) %>%
    summarise(IntersectionSize = n(), Genes = paste(Gene, collapse = ";"), .groups = "drop") %>%
    rowwise() %>% mutate(Number_of_AS_types = sum(c_across(all_of(event_types)))) %>% ungroup() %>%
    filter(Number_of_AS_types > 0) %>% arrange(desc(IntersectionSize), desc(Number_of_AS_types)) %>%
    mutate(IntersectionID = row_number())

write.table(intersection_summary, file = file.path(out_dir, "UpSet_intersection_summary.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)

set_size_data <- gene_membership_long %>% count(EventType, name = "SetSize") %>% complete(EventType = event_types,
    fill = list(SetSize = 0))

event_order <- set_size_data %>% arrange(desc(SetSize)) %>% pull(EventType) %>% as.character()

cat("Event order from bottom to top:\n")

print(event_order)

set_size_data <- set_size_data %>% mutate(EventType = factor(EventType, levels = event_order), y = as.numeric(EventType))

matrix_data <- intersection_summary %>% select(IntersectionID, all_of(event_types)) %>% pivot_longer(cols = all_of(event_types),
    names_to = "EventType", values_to = "Present") %>% mutate(EventType = factor(EventType, levels = event_order),
    y = as.numeric(EventType), Present = as.integer(Present))

line_data <- matrix_data %>% filter(Present == 1) %>% group_by(IntersectionID) %>% summarise(ymin = min(y),
    ymax = max(y), .groups = "drop")

row_bg <- tibble(y = seq_along(event_order), fill_group = rep(c("white", "gray95"), length.out = length(event_order)))

n_intersections <- nrow(intersection_summary)

p_top <- ggplot(intersection_summary, aes(x = IntersectionID, y = IntersectionSize)) + geom_col(fill = "gray30",
    width = 0.6) + geom_text(aes(label = IntersectionSize), vjust = -0.35, size = 4) + scale_x_continuous(limits = c(0.5,
    n_intersections + 0.5), expand = c(0, 0)) + scale_y_continuous(limits = c(0, top_y_max), breaks = c(0,
    100, 300, 500), expand = expansion(mult = c(0, 0.03))) + labs(y = "Gene Intersections", x = NULL) +
    theme_classic(base_size = 14) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 13), plot.margin = margin(5,
        5, 0, 5))

p_matrix <- ggplot() + geom_rect(data = row_bg, aes(xmin = 0.5, xmax = n_intersections + 0.5, ymin = y -
    0.5, ymax = y + 0.5, fill = fill_group), inherit.aes = FALSE) + scale_fill_identity() + geom_point(data = matrix_data %>%
    filter(Present == 0), aes(x = IntersectionID, y = y), color = "gray85", size = 3.2) + geom_segment(data = line_data,
    aes(x = IntersectionID, xend = IntersectionID, y = ymin, yend = ymax), linewidth = 1.5, color = "gray20",
    lineend = "round") + geom_point(data = matrix_data %>% filter(Present == 1), aes(x = IntersectionID,
    y = y), color = "gray20", size = 3.6) + scale_x_continuous(limits = c(0.5, n_intersections +
    0.5), expand = c(0, 0)) + scale_y_continuous(breaks = seq_along(event_order), labels = event_order,
    limits = c(0.5, length(event_order) + 0.5), expand = c(0, 0)) + labs(x = NULL, y = NULL) + theme_void(base_size = 14) +
    theme(axis.text.y = element_text(size = 14, color = "black"), plot.margin = margin(0, 5, 5,
        5))

p_setsize <- ggplot(set_size_data) + geom_rect(aes(xmin = 0, xmax = SetSize, ymin = y - 0.32, ymax = y +
    0.32, fill = EventType)) + scale_fill_manual(values = as_colors, guide = "none") + scale_y_continuous(breaks = seq_along(event_order),
    labels = NULL, limits = c(0.5, length(event_order) + 0.5), expand = c(0, 0)) + scale_x_reverse(limits = c(setsize_x_max,
    0), breaks = c(500, 300, 100, 0), expand = c(0, 0)) + labs(x = "Significant AS genes", y = NULL) +
    theme_classic(base_size = 14) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 15), axis.text.x = element_text(size = 12), plot.margin = margin(0,
        5, 5, 5))

p_empty <- ggplot() + theme_void()

p_final <- (p_empty + p_top + plot_layout(widths = c(0.25, 0.75)))/(p_setsize + p_matrix + plot_layout(widths = c(0.25,
    0.75))) + plot_layout(heights = c(0.58, 0.42))

save_panel(p_final, file.path(out_dir, "Figure1D"), 10, 6)

stopifnot(sum(summary_table$Significant_events) == 1073L, nrow(gene_membership_wide) == 807L, sum(intersection_summary$IntersectionSize) ==
    807L)
