#!/usr/bin/env Rscript
# KC6 sensitivity: event-level overlap, not overlap of run-local numeric IDs.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", arg)))
source(file.path(script_dir, "Utils_plotting.R"))
source(file.path(script_dir, "Utils_splicing_events.R"))
source(file.path(script_dir, "Utils_Venn_plots.R"))
if (any(!commandArgs(TRUE) %in% "--tables-only")) stop("Unknown argument")
out <- file.path(ROOT, "results/stroma/suppfigure2")
calls <- list(); class_counts <- list()
for (cohort in c("stroma", "stroma_with_kc6")) {
    hits <- lapply(names(RMATS_COORDS), function(event) {
        raw <- read_rmats(cohort, event)
        data <- raw[is_das(raw), ]
        data.frame(Event = event, Event_key = paste(event, data$event_key, sep = "|"), GeneID = data$GeneID)
    })
    for (i in seq_along(hits)) {
        class_counts[[length(class_counts) + 1L]] <- data.frame(Cohort = cohort,
            Event = names(RMATS_COORDS)[i], Events = nrow(hits[[i]]), Genes = length(unique(hits[[i]]$GeneID)))
    }
    calls[[cohort]] <- do.call(rbind, hits)
}
a <- calls$stroma$Event_key; b <- calls$stroma_with_kc6$Event_key
common <- intersect(a, b)
counts <- data.frame(KC6_excluded = length(a), KC6_included = length(b),
    Excluded_only = length(setdiff(a, b)), Shared = length(common), Included_only = length(setdiff(b, a)))
stopifnot(identical(as.integer(counts[1, ]), c(1073L, 1074L, 184L, 889L, 185L)),
          length(unique(calls$stroma_with_kc6$GeneID)) == 802L)
members <- data.frame(Event_key = union(a, b))
members$KC6_excluded <- members$Event_key %in% a
members$KC6_included <- members$Event_key %in% b
write_tsv(counts, file.path(out, "event_overlap_counts.tsv"))
write_tsv(members, file.path(out, "event_memberships.tsv"))
write_tsv(do.call(rbind, class_counts), file.path(out, "event_class_counts.tsv"))
check_venn_package()
# ind=FALSE builds the plot object without opening a graphics device.
plot <- VennDiagram::draw.pairwise.venn(length(a), length(b), length(common),
    category = c("Keratoconus 6 excluded", "Keratoconus 6 included"),
    fill = c("#FFC000", "#BFBFBF"), alpha = c(0.6, 0.6), scaled = FALSE, euler.d = FALSE,
    fontfamily = FONT, cat.fontfamily = FONT, cex = 1.5, cat.cex = 0.95, ind = FALSE)
save_venn(grid::grobTree(children = plot), file.path(out, "SuppFigure2B_KC6_event_overlap"), 7, 4.5)
print(counts, row.names = FALSE)
