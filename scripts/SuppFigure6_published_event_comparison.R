#!/usr/bin/env Rscript
# Supplementary Figure 6: exact-coordinate comparison with Sun et al. (2022).
# The reference coordinates were recovered from the supplied membership tables;
# see docs/suppfigure6.md. No genome-build conversion is applied.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
source(file.path(dirname(normalizePath(sub("^--file=", "", arg))), "Utils_plotting.R"))
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--tables-only")) stop("Only --tables-only is supported; set paths with KC_WORK_DIR.")

out <- file.path(ROOT, "results/epithelium/suppfigure6")
reference <- read_tsv(file.path(ROOT, "data/reference/Sun2022_published_event_coordinates.tsv"))

# All boundaries are required, including flanking exons. Gene names and event
# IDs alone do not establish equivalence between separate rMATS analyses.
coordinates <- list(
    SE = c("chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
           "downstreamES", "downstreamEE"),
    A3SS = c("chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
             "flankingES", "flankingEE"),
    A5SS = c("chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
             "flankingES", "flankingEE"),
    MXE = c("chr", "strand", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base",
            "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
    RI = c("chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE",
           "downstreamES", "downstreamEE")
)

event_keys <- function(data, event) {
    columns <- coordinates[[event]]
    if (!all(columns %in% names(data))) stop(event, ": missing coordinate columns")
    if (!nrow(data)) return(character())
    values <- lapply(columns, function(column) {
        value <- as.character(data[[column]])
        if (anyNA(value) || any(!nzchar(value))) stop(event, ": missing ", column)
        if (column == "chr") {
            value <- sub("^chr", "", value)
            if (any(!grepl("^[A-Za-z0-9_.-]+$", value))) stop("Invalid chromosome")
        } else if (column == "strand") {
            if (any(!value %in% c("+", "-"))) stop("Invalid strand")
        } else {
            number <- suppressWarnings(as.numeric(value))
            if (any(!is.finite(number) | number < 0 | number != floor(number))) {
                stop(event, ": invalid boundary in ", column)
            }
            value <- format(number, scientific = FALSE, trim = TRUE)
        }
        value
    })
    do.call(paste, c(values, sep = "|"))
}

if (!"Event" %in% names(reference) || anyNA(reference$Event) ||
    any(!reference$Event %in% names(coordinates))) stop("Invalid reference event classes")
summaries <- memberships <- audits <- list()
for (event in names(coordinates)) {
    published <- reference[reference$Event == event, , drop = FALSE]
    raw <- read_tsv(file.path(ROOT, "data/epithelium/rmats", paste0(event, ".MATS.JC.txt")))
    names(raw) <- sub("^X([12](st|nd)Exon)", "\\1", names(raw))
    if (!"FDR" %in% names(raw) || !is.numeric(raw$FDR) ||
        any(!is.na(raw$FDR) & (!is.finite(raw$FDR) | raw$FDR < 0 | raw$FDR > 1))) {
        stop(event, ": invalid or missing FDR")
    }
    raw_keys <- event_keys(raw, event)
    significant <- !is.na(raw$FDR) & raw$FDR < 0.05
    # This comparison deliberately has NO absolute delta-PSI cutoff.
    ours <- unique(raw_keys[significant])
    theirs <- event_keys(published, event)
    if (anyDuplicated(theirs)) stop(event, ": duplicate published coordinates")
    common <- intersect(ours, theirs)
    summaries[[event]] <- data.frame(Event = event, Our_total = length(ours),
        Sun_total = length(theirs), Our_only = length(setdiff(ours, theirs)),
        Overlap = length(common), Sun_only = length(setdiff(theirs, ours)))
    detail <- data.frame(Event = event, Coordinate_key = union(ours, theirs))
    detail$Membership <- ifelse(detail$Coordinate_key %in% common, "Both",
        ifelse(detail$Coordinate_key %in% ours, "Our reanalysis only", "Sun published only"))
    memberships[[event]] <- detail
    # Audit all published events, including those outside the significant set.
    audits[[event]] <- data.frame(Event = event, Coordinate_key = theirs,
        Status = ifelse(theirs %in% ours, "Significant (FDR < 0.05)",
            ifelse(theirs %in% raw_keys, "Exact match outside significant set", "No exact match")))
}
summary <- do.call(rbind, unname(summaries))
audit <- do.call(rbind, unname(audits))
stopifnot(identical(summary$Sun_total, c(33L, 3L, 6L, 11L, 14L)),
          identical(summary$Overlap, c(16L, 1L, 5L, 4L, 10L)),
          identical(summary$Our_total, c(1214L, 274L, 227L, 154L, 746L)))
for (event in names(memberships)) {
    write_tsv(memberships[[event]], file.path(out, paste0(event, "_event_membership.tsv")))
}
write_tsv(summary, file.path(out, "venn_counts.tsv"))
write_tsv(audit, file.path(out, "published_event_audit.tsv"))

# Preserve the supplied five-panel drawing specifications. Circle areas are not
# proportional to set size. No graphics device opens with --tables-only.
purple <- "#E3C4DB"
blue <- "#72C8E8"
panel <- function(i) {
    r <- summary[i, ]
    par(mar = c(1, 1, 2, 1))
    plot.new()
    plot.window(xlim = c(-1.6, 1.6), ylim = c(-1.3, 1.4), asp = 1)
    theta <- seq(0, 2 * pi, length.out = 500)
    polygon(-.48 + .9 * cos(theta), .9 * sin(theta), col = purple, border = "#000000", lwd = 1.5)
    polygon(.48 + .9 * cos(theta), .9 * sin(theta), col = adjustcolor(blue, alpha.f = .55), border = blue, lwd = 1.5)
    text(c(-.9, 0, .9), 0, c(r$Our_only, r$Overlap, r$Sun_only), cex = 1.6, font = 2)
    title(r$Event, cex.main = 1.6)
    text(-.75, -1.1, sprintf("Our reanalysis\nn = %s", r$Our_total), cex = .85)
    text(.75, -1.1, sprintf("Sun published\nn = %s", r$Sun_total), cex = .85)
}
draw_all <- function() {
    par(mfrow = c(2, 3), oma = c(2, 0, 2, 0))
    for (i in seq_len(nrow(summary))) panel(i)
    plot.new()
    legend("center", legend = c("Our reanalysis", "Sun published"), fill = c(purple, blue), border = NA, bty = "n", cex = 1.1)
    mtext("Sun2022: rMATS event overlap", outer = TRUE, side = 3, line = .4, font = 2, cex = 1.3)
    mtext("Reanalysis: FDR < 0.05; no delta-PSI cutoff. Full event coordinates matched. Circles not to scale.",
          outer = TRUE, side = 1, line = .5, cex = .8)
}
if (!TABLES_ONLY) {
    if (!requireNamespace("ragg", quietly = TRUE)) stop("Install ragg to create figure images")
    ragg::agg_png(file.path(out, "Supplementary6.png"), width = 12, height = 7.33, units = "in", res = 600)
    par(family = FONT); draw_all(); dev.off()
    pdf_device(file.path(out, "Supplementary6.pdf"), width = 12, height = 7.33)
    par(family = FONT); draw_all(); dev.off()
    for (i in seq_len(nrow(summary))) {
        pdf_device(file.path(out, paste0("Supplementary6_", summary$Event[i], ".pdf")), width = 4, height = 3.6)
        par(family = FONT); panel(i); dev.off()
    }
}
print(summary, row.names = FALSE)
