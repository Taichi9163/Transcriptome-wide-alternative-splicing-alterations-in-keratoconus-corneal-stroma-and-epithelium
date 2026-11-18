#!/usr/bin/env Rscript
# Four directional sets per event class, calculated from the full rMATS inputs.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "utils_plotting.R"))
source(file.path(script_dir, "utils_splicing_events.R"))
source(file.path(script_dir, "utils_venn_plots.R"))
if (any(!commandArgs(TRUE) %in% "--tables-only")) stop("Unknown argument")
expected_shared <- c(SE = 13L, A3SS = 1L, A5SS = 4L, MXE = 0L, RI = 5L)
summaries <- list()
for (event in names(RMATS_COORDS)) {
    data <- four_set_memberships(cross_layer_events(event))
    regions <- four_set_regions(data)
    both <- data$significant_stroma & data$significant_epithelium
    stopifnot(sum(both) == expected_shared[[event]], sum(regions$count) == nrow(data))
    out <- file.path(ROOT, if (event == "SE") "results/cross_layer/figure4A" else "results/cross_layer/suppfigure8")
    write_tsv(data, file.path(out, paste0(event, "_event_memberships.tsv")))
    write_tsv(regions, file.path(out, paste0(event, "_exclusive_regions.tsv")))
    # Retain the previous coordinate/direction table interface.
    coordinates <- as.data.frame(do.call(rbind, strsplit(data$event_key, ":", fixed = TRUE)), stringsAsFactors = FALSE)
    names(coordinates) <- RMATS_COORDS[[event]]
    coordinates$stroma <- ifelse(data$significant_stroma, ifelse(data$delta_stroma > 0, "inclusion", "skipping"), NA)
    coordinates$epithelium <- ifelse(data$significant_epithelium, ifelse(data$delta_epithelium > 0, "inclusion", "skipping"), NA)
    write_tsv(coordinates, file.path(ROOT, "results/cross_layer", paste0(event, "_union.tsv")))
    summaries[[event]] <- data.frame(Event = event, Stroma = sum(data$significant_stroma),
        Epithelium = sum(data$significant_epithelium), Shared = sum(both), Union = nrow(data),
        Shared_same_direction = sum(data$delta_stroma[both] * data$delta_epithelium[both] > 0),
        Shared_opposite_direction = sum(data$delta_stroma[both] * data$delta_epithelium[both] < 0))
    plot <- make_four_venn(data, paste("Cross-layer", event, "events"), show_genes = event == "SE")
    save_venn(plot, file.path(out, paste0(event, "_four_set_Venn")))
}
summary <- do.call(rbind, unname(summaries))
write_tsv(summary, file.path(ROOT, "results/cross_layer/directional_overlap_counts.tsv"))
print(summary, row.names = FALSE)
