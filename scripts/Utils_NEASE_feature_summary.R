# Summarize the event-preserving NEASE output for Figure 4C / Supplementary Figure 10A.
# Load Utils_plotting.R first. All counts and text labels come from event memberships.
# The lower bar counts event-class assignments, not unique mapped events.

read_nease_feature_summary <- function(cohort, panel) {
    stopifnot(cohort %in% c("stroma", "epithelium"))
    args <- setdiff(commandArgs(trailingOnly = TRUE), "--tables-only")
    event_path <- file.path(ROOT, "data", cohort, "nease", "SE_events.tsv")
    if (length(args)) {
        if (length(args) != 2L || args[1] != "--events") {
            stop("Usage: Rscript script.R [--events SE_events.tsv] [--tables-only]")
        }
        event_path <- args[2]
        if (!grepl("^(/|[A-Za-z]:[/\\\\])", event_path)) event_path <- file.path(ROOT, event_path)
    }
    events <- read_tsv(event_path)
    flags <- c("Mapped_Domains", "Mapped_Residue_evidence", "Mapped_Motifs")
    bool_columns <- c(flags, "Mapped_any_feature", "Excluded_by_NEASE_ORF_filter")
    required <- c("Event_ID", "ID", "Event_coordinates", "GeneID", "geneSymbol", "FDR",
                  "IncLevelDifference", bool_columns)
    if (!all(required %in% names(events))) stop("Incomplete NEASE event export: ", event_path)
    if (!nrow(events) || anyNA(events$Event_ID) || any(events$Event_ID == "") ||
        anyDuplicated(events$Event_ID) || anyDuplicated(events$Event_coordinates)) {
        stop("NEASE input must contain one row per distinct SE event")
    }
    for (name in bool_columns) {
        values <- tolower(as.character(events[[name]]))
        if (anyNA(values) || any(!values %in% c("true", "false"))) {
            stop("Expected TRUE/FALSE memberships in ", name)
        }
        events[[name]] <- values == "true"
    }
    if (any(!is.finite(events$FDR)) || any(!is.finite(events$IncLevelDifference)) ||
        any(events$FDR < 0 | events$FDR >= 0.05) || any(abs(events$IncLevelDifference) <= 0.1) ||
        any(events$Excluded_by_NEASE_ORF_filter)) {
        stop("Use the primary significant SE export without the optional ORF exclusion")
    }

    # Check identities and coordinates against this cohort's complete rMATS JC table.
    raw <- read.delim(file.path(ROOT, "data", cohort, "rmats", "SE.MATS.JC.txt"),
                      check.names = FALSE, stringsAsFactors = FALSE)
    selected <- raw[which(is.finite(raw$FDR) & raw$FDR < 0.05 &
                         is.finite(raw$IncLevelDifference) & abs(raw$IncLevelDifference) > 0.1), ]
    if (anyDuplicated(selected$ID) || !setequal(as.character(events$ID), as.character(selected$ID))) {
        stop("NEASE event IDs differ from the significant rMATS SE calls")
    }
    selected <- selected[match(as.character(events$ID), as.character(selected$ID)), ]
    coordinates <- c("chr", "strand", "exonStart_0base", "exonEnd", "upstreamES",
                     "upstreamEE", "downstreamES", "downstreamEE")
    keys <- do.call(paste, c(selected[, coordinates], sep = ":"))
    if (!identical(unname(keys), events$Event_coordinates) ||
        !identical(sub("\\..*$", "", selected$GeneID), events$GeneID) ||
        !isTRUE(all.equal(selected$FDR, events$FDR, tolerance = 1e-12)) ||
        !isTRUE(all.equal(selected$IncLevelDifference, events$IncLevelDifference, tolerance = 1e-12))) {
        stop("NEASE coordinates, genes or statistics differ from the rMATS input")
    }
    membership <- as.matrix(events[, flags])
    colnames(membership) <- c("Domains", "Residues", "Motifs")
    mapped <- rowSums(membership) > 0L
    if (!identical(unname(mapped), events$Mapped_any_feature)) stop("Mapped-any flag disagrees with feature classes")
    counts <- colSums(membership)
    pairwise <- c(Domains_Residues = sum(membership[, 1] & membership[, 2]),
                  Domains_Motifs = sum(membership[, 1] & membership[, 3]),
                  Residues_Motifs = sum(membership[, 2] & membership[, 3]))
    all_three <- sum(rowSums(membership) == 3L)
    total <- nrow(events)
    mapped_count <- sum(mapped)
    assignments <- sum(counts)
    stopifnot(assignments - sum(pairwise) + all_three == mapped_count)
    # Reproduction checks, separate from the data-derived counts used by the plot.
    expected <- if (cohort == "stroma") c(698, 276, 264, 50, 7, 41, 4, 3, 3) else
        c(450, 135, 129, 15, 2, 10, 1, 0, 0)
    if (!identical(as.numeric(c(total, mapped_count, counts, pairwise, all_three)), expected)) {
        stop("Feature counts differ from the manuscript for ", cohort)
    }
    out_dir <- file.path(ROOT, "results", cohort,
                         if (cohort == "stroma") "figure4" else "suppfigure10")
    write_tsv(data.frame(category = c("Mapped", "Unmapped"),
                        count = c(mapped_count, total - mapped_count), denominator = total,
                        percent = 100 * c(mapped_count, total - mapped_count) / total),
              file.path(out_dir, paste0(panel, "_coverage.tsv")))
    write_tsv(data.frame(feature = names(counts), count = unname(counts),
                        denominator = assignments, percent = 100 * unname(counts) / assignments),
              file.path(out_dir, paste0(panel, "_feature_classes.tsv")))
    write_tsv(data.frame(overlap = c(names(pairwise), "All_three"),
                        count = c(unname(pairwise), all_three)),
              file.path(out_dir, paste0(panel, "_overlaps.tsv")))
    write_tsv(cbind(events[, c("Event_ID", "GeneID", "geneSymbol", "Event_coordinates")], membership),
              file.path(out_dir, paste0(panel, "_event_memberships.tsv")))
    message(panel, ": ", mapped_count, "/", total, " mapped events; ", assignments, " class assignments.")
    list(total = total, mapped = mapped_count, counts = counts, pairwise = pairwise,
         all_three = all_three, assignments = assignments, out_dir = out_dir)
}

save_nease_feature_panel <- function(draw, out_dir, panel, width, height) {
    if (TABLES_ONLY) return(invisible(NULL))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    png(file.path(out_dir, paste0(panel, ".png")), width = width, height = height,
        res = 200, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    tryCatch(draw(), finally = dev.off())
    svg(file.path(out_dir, paste0(panel, ".svg")), width = width / 200, height = height / 200)
    tryCatch(draw(), finally = dev.off())
}
