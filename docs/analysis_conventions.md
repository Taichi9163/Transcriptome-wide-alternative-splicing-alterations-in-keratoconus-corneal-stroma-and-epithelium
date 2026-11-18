# Analysis conventions and input provenance

[Repository overview](../README.md) · [Documentation](README.md) · [Figure-to-code index](figures.md)

## Data and thresholds

The bundled numeric tables come from the supplied analysis package and the recovered local analyses for the manuscript revision of 18 November 2026; see [input provenance](../inputs/README.md). Compression is lossless; `inputs/manifest.json` records the SHA-256 checksum and byte count of each uncompressed source. The complete gene table preserves its original row-name field (the header has one fewer field than the data rows); the transcript table has explicit `TranscriptID` and `GeneID` columns.

DAS requires rMATS JC FDR < 0.05 and |IncLevelDifference| > 0.10. DE genes/isoforms require DESeq2 adjusted P < 0.05 and |log2 fold change| > 1. All comparisons use keratoconus minus control. Missing effects are never replaced by zero.

DAS gene counts and expression integration use Ensembl GeneID. Figure 1D and the reconstructed Supplementary Figure 7D follow the supplied gene-symbol UpSet, with GeneID as fallback when a symbol is absent. The two identifiers give the same totals of 807 stromal and 609 epithelial genes for these inputs; they should not be assumed equivalent for other datasets.

Figure 1B uses centered, scaled event-level PSI in `prcomp`. Figure 1C clusters Euclidean distances of unscaled PSI using Ward.D2, then uses row Z scores for colors. The significant-event PCA describes a set selected using the group labels. The optional all-event and variance-ranked PCA modules avoid significance-based event selection.

For Figure 3B, select the largest absolute ΔPSI per DAS gene, breaking ties by smaller FDR and then original row order across SE, RI, A5SS, A3SS and MXE. Of 807 DAS genes, 805 have finite matched gene-expression effects. DEG-only genes are retained in exported source tables but omitted from the scatter.

## Event–transcript correspondence

`map_events.py` uses [Ensembl release 114](https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/). It converts GTF exon starts to 0-based half-open intervals and keeps rMATS coordinates in that same convention. Gene, chromosome and strand must match. A transcript in the same gene alone is insufficient.

| Event | Inclusion-compatible structure | Skipping-compatible structure |
| --- | --- | --- |
| SE | Exact target exon with both flanking junctions | Direct upstream-to-downstream junction |
| MXE | First genomic alternative exon on `+`, second on `−`, with both flanking junctions | The other alternative exon and both junctions |
| A3SS/A5SS | Long alternative boundary joined to the flanking exon | Short alternative boundary joined to the flanking exon |
| RI | A single uninterrupted exon spanning the retained-intron exon model | Exact intron-excision junction |

MXE orientation follows [rMATS's strand-dependent inclusion definition](https://github.com/Xinglab/rmats-turbo). For A3SS/A5SS, the alternative boundary and shared sequence must match; distal exon boundaries may vary. Structures compatible with both forms are excluded as ambiguous.

All compatible annotated pairs enter `event_transcript_mapping.tsv`. Figure 3D retains pairs with finite isoform log2 fold change and raw P value. An unavailable adjusted P value does not prevent plotting but cannot establish DE status. The event audit separates no compatible annotated transcript from compatible transcripts without tested expression effects. Inclusion/skipping is shown by point shape. **ΔPSI is not sign-reversed for skipping-compatible transcripts.**

Figure 3C distinguishes the 55-gene overlap at gene level from the subset of 29 genes with a structurally matched DE isoform. Ten gene callouts are selected by decreasing absolute Wald statistic among matched DE transcripts, then adjusted P value and transcript ID, taking the first transcript per gene.

## Sensitivity and scope

`leave_one_out.py` holds the primary 1,073-event DAS set fixed and recalculates mean case PSI minus mean control PSI after omitting each of the nine samples. It compares effects with the reported rMATS ΔPSI; it does not rerun rMATS or claim new FDR values. `KC_pos1` etc. denote original PSI-vector positions. The median changes can include small rounding differences because rMATS reports rounded PSI and ΔPSI values.

`stromal_sensitivity.py` separately evaluates junction support and the KC6-inclusive rMATS calls. Cross-run event matching uses exact class-specific coordinates, chromosome and strand, rather than run-local numeric IDs. See the additional-analysis guide for these inputs.

Figure generation uses portable PDF devices and an installed font (`FIGURE_FONT`, default `sans`). Font metrics and plotting-library versions may change appearance. The submitted PowerPoint composition, gel measurements and conceptual schematic are outside this computational reproduction.
