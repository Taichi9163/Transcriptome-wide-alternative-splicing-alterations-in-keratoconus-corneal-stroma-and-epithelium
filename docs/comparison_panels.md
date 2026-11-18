# Event-comparison panels

These entry scripts run directly from the bundled rMATS JC tables, without a GTF download, workbook import or new alignment. Paths are resolved under `KC_WORK_DIR`; scripts may be called by absolute path from another directory. Install `VennDiagram` in R in addition to the packages listed in the README. `ragg` is needed only when exporting images.

```sh
python3 scripts/prepare_inputs.py --work-dir work
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/SuppFigure2B_KC6_event_overlap.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/Figure4A_SuppFigure8_DAS_overlap.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/SuppFigure9AB_SE_junction_support_plots.R --tables-only
```

`--tables-only` writes numeric results and builds plot objects without opening graphics devices. Remove the flag to create PNG/PDF panels. Use a separate work directory for each new run, as result tables may be replaced within that directory. The Venn layouts are portable reconstructions, not copies of the original PowerPoint geometry. Areas are not proportional to set sizes. Four-set gene labels are listed below the diagram rather than in the original callout boxes. Supplementary Figure 8 is exported as four individual panels (A3SS, A5SS, MXE and RI, corresponding to A–D).

| Panel | Output directory under `work/results/` | Numeric checks |
| --- | --- | --- |
| Supplementary 2B | `stroma/suppfigure2/` | KC6 excluded: 1,073; included: 1,074; exclusive/shared/exclusive: 184/889/185 |
| Figure 4A | `cross_layer/figure4A/` | 698 stromal and 450 epithelial SE events; 13 shared, seven concordant and six discordant |
| Supplementary 8A–D | `cross_layer/suppfigure8/` | Shared A3SS/A5SS/MXE/RI events: 1/4/0/5 |
| Supplementary 9A,B | `cross_layer/suppfigure9/` | 15,343 qualifying SE events; 427 stromal-only, 56 epithelial-only, eight shared significant and 14,852 nonsignificant in both |

## Event identity and Venn regions

All comparisons use strict FDR < 0.05 and |IncLevelDifference| > 0.10. Identity uses chromosome, strand and all event-specific rMATS exon boundaries, including flanking exons. Numeric event IDs and gene symbols are not identity keys. An initial `chr` prefix is standardized; no coordinates are shifted and no genome-build conversion is applied. Duplicate coordinate keys within one event class and cohort are rejected. The KC6 comparison also includes event class in each key.

Four-set membership strings use this fixed bit order: stromal inclusion, stromal skipping, epithelial inclusion, epithelial skipping. For example, `1001` means increased inclusion in stroma and increased skipping in epithelium. Every panel exports all 15 exclusive region counts, including zeros, and a complete event-membership table. Within-cohort opposite directions cannot overlap, so triple/four-way regions are zero. Multiple events in one gene remain multiple events; shared gene labels are deduplicated only for display. Missing gene symbols are omitted from gene lists without dropping their events.

The Figure 4A / Supplementary Figure 8 command also retains the earlier `cross_layer/<event>_union.tsv` coordinate/direction outputs. These contain the same memberships as the previous release; row order and missing-value serialization may differ. Inputs are now read from the complete JC tables, so precomputed `results/<cohort>/das/` files are no longer required by this command.

## Supplementary Figure 9 support criterion

Start with coordinate-matched SE events regardless of significance. An event must have finite PSI and at least ten inclusion-plus-skipping junction reads in the same specimen for at least half of **each** case/control group in **both** cohorts: at least 3/5 stromal cases, 2/4 stromal controls, 4/7 epithelial cases and 4/7 epithelial controls. Both cohort effects must be finite. Zero-valued measured effects are retained.

The Venn uses the 491 qualifying events significant in either cohort; the scatter uses all 15,343 qualifying events. The eight shared significant events comprise five with the same direction and three with opposite directions. Descriptive correlations are Pearson r = 0.1011075 and Spearman rho = 0.08498896. The command exports these values; the manuscript's rounded correlation labels were separate PowerPoint annotations and are not inserted into the scatter by this script.

The scatter specification comes from `outputs/supplementary_figure8_font_v3_20261015/build/scatter.R`, whose earlier figure number became Supplementary Figure 9 in the supplied November revision. Point coordinates, color, shape, size and alpha agree for all 15,343 events. `FIGURE_FONT` replaces the hard-coded Arial dependency. The old PDF-editing implementation of panel A depended on a temporary embedded figure and fixed PDF drawing-operation indices; the public implementation calculates its counts and labels from event memberships.

The existing `SuppFigure9_SE_junction_support.R` remains the broader numeric audit, including additional sensitivity thresholds and reconciliation with the Figure 4B exports. The new `SuppFigure9AB` command is a standalone panel workflow and does not require running that audit first.

## Input provenance and validation

Five complete KC6-inclusive JC tables were added from the local `Yamashita2026/after_rMATS/allsamples/results` analysis. Lossless gzip files and their uncompressed hashes are in the input manifest. This increases the bundled input count to 26. Per-class included-cohort counts are SE 713, A3SS 78, A5SS 85, MXE 44 and RI 154, in 802 unique genes overall. The class-count output also supplies the numeric labels for Supplementary Figure 2A; its exon schematics remain a separate composition.

Validation compared all five union tables with the previous release and all 40 directional counts with the saved cross-layer matrix. Every support-qualified event key, effect, FDR, significance flag and per-group support count matched the original audit. All 15 support-filtered Venn regions and their 491 event memberships matched the saved region audit. The reconstructed Venn layouts have not been visually validated against the manuscript graphics.

Run the focused numeric/plot-object tests with:

```sh
Rscript --vanilla tests/test_comparison_panels.R
python3 -m unittest discover -s tests -v
```
