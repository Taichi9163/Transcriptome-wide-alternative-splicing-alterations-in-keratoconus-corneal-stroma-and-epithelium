# Processed input provenance

[Repository overview](../README.md) · [Reproduction guide](../docs/reproduce.md) · [Checksum manifest](manifest.json)

All bundled files are losslessly gzip-compressed. `manifest.json` records each destination, uncompressed byte count and SHA-256 checksum. `scripts/prepare_inputs.py` verifies these values and refuses to replace different existing inputs.

| Input | Source and purpose |
| --- | --- |
| `stroma/rmats/*.MATS.JC.txt.gz` | Five complete stromal rMATS JC tables from the supplied Figure 1/Figure 3 archive; primary five-case/four-control comparison |
| `stroma_with_kc6/rmats/*.MATS.JC.txt.gz` | Five complete JC tables from the local Yamashita2026 `after_rMATS/allsamples/results` analysis; six cases/four controls, 1,074 significant events in 802 genes; Supplementary Figure 2 and stromal sensitivity input |
| `stroma/expression/*_deseq2.tsv.gz` | Gene and transcript DESeq2 result tables from the same supplied archive |
| `epithelium/rmats/*.MATS.JC.txt.gz` | Five complete JC tables from the local Sun2022 `after_rMATS/results` analysis, seven samples per group; inputs to Supplementary Figures 6 and 7D |
| `<cohort>/rmaps/pVal.*.vs.bg.RNAmap.txt.gz` | Original stromal (Yamashita2026) and epithelial (Sun2022) rMAPS3 output tables, 114 motifs × eight regions per direction; candidate selection is recomputed from these P values |
| `<cohort>/rbp/tpm.tsv.gz` | Candidate-gene TPM values and Ensembl IDs from the saved outputs of the recovered RBP heatmap scripts; 19 genes × nine samples for stroma and 23 × 14 for epithelium |
| `<cohort>/nease/SE_events.tsv.gz` | Complete event-preserving NEASE exports from the local `python_revision_20261010/results` analyses; 698 stromal and 450 epithelial SE events, with the optional ORF exclusion disabled |
| `<cohort>/nease/Reactome_Python.tsv.gz` | Complete corrected 1,823-pathway exports from each `python_revision_20261010/results` analysis; primary analysis without ORF exclusion, for Fig. 4D / Supplementary Fig. 10B |

The additional `reference/Sun2022_published_event_coordinates.tsv.gz` contains 67 published event coordinates recovered from the supplied Supplementary Figure 6 membership tables. The original reference workbook was not included in the archive. [Supplementary Figure 6 provenance](../docs/suppfigure6.md) documents the extraction, source citation, attribution, matching rules and the lack of independent verification against that workbook.

The recovered source scripts were `stroma_rMAPS_hit_RNA_heatmap_BH.R` and `epithelium_Sun2022_rMAPS_hit_RNA_heatmap_BH.R` in the project's `RBP_expression_heatmap/scripts/` directory. Their saved `*_TPM.csv` values and `*_metadata.csv` gene IDs are combined into each bundled `tpm.tsv`. Sample labels and their order are retained from those exports. These matrices contain the candidate genes only; they do not replace full RSEM quantifications or support a new genome-wide expression analysis.

The original stromal heatmap read the `Gene Matrix` sheet of the marker-expression workbook. The epithelial heatmap read the Sun2022 `DataMatrix_ImportSampleList.txt` abundance columns. Gene-symbol lookups came from the study DESeq2 workbooks, including the original stromal MSI1/HNRNPCL1 lookup entries. Bundling the saved numeric matrices removes those private workbook/path dependencies without changing the plotted values.

Supplementary Figure 3 candidate selection preserves BH correction across 912 motif–region P values separately per cohort and direction, the strict cutoff of 5e-5, and the HuR → ELAVL1 and SF2-ASF → SRSF1 aliases. The heatmaps preserve log2(TPM + 1), row Z-scores, Euclidean/Ward.D2 row clustering, and fixed control-then-case column order. Row Z-scores are unchanged by using log10 instead of log2. No DESeq2 adjusted P values are inferred from TPM; the separate expression-statistics script requires full gene-model results.

Supplementary Figure 7D is reconstructed from the Figure 1D code and epithelial inputs. It uses FDR < 0.05, |ΔPSI| > 0.10, and one gene per event class (gene symbol with GeneID fallback). The total is 805 events in 609 genes, including 314 genes affected only by SE. The original epithelial UpSet source was not located. This release validates numeric data and plot construction without regenerating that figure.

## Recovered NEASE feature panels

The Figure 4C plotting source is `outputs/nease_overlap_preview_20261012/plot_overlap_extra_large_text.R`. The Supplementary Figure 10A source is `plot_epithelium_with_axis_labels.R` in the same directory. Their existing PNG outputs have identical decoded RGB pixels to the embedded images in slides 5 and 15, respectively, of the supplied 18 November 2026 figure deck. This comparison identifies the source versions; no images were regenerated.

The public scripts retain those drawing specifications and calculate their labels from the bundled event exports instead of copying fixed counts. Original NEASE domain, motif and residue-evidence annotation tables independently confirmed all memberships. The stromal counts are 276 mapped of 698 events, with 264 domain, 50 residue and seven motif assignments; the epithelial counts are 135 mapped of 450, with 129, 15 and two assignments. Inclusive domain/residue, domain/motif and residue/motif intersections are 41/4/3 in stroma and 10/1/0 in epithelium. The three-class intersections are three and zero. These data describe annotated events and do not measure functional protein changes.

The bundled NEASE tables are saved processed mappings, not a rerun of the reference annotation or pathway enrichment. `nease_analysis.py` provides the upstream calculation and exports the same event-table schema. The two feature-panel scripts check event IDs, coordinates, genes, FDR and ΔPSI against the bundled rMATS calls, require the primary analysis without ORF exclusion, and distinguish unique mapped events from event-class assignments.

The KC6-inclusive inputs are lossless copies of existing analysis outputs, not a new rMATS run. [Event-comparison panels](../docs/comparison_panels.md) documents the newly added Venn/scatter commands and their validation.

The corrected Reactome tables replace the older 1,750-row CSV inputs referenced by the recovered R sources. They are saved enrichment results, not a new NEASE run. See [Reactome selection and provenance](../docs/reactome.md). The rMAPS3 audit reuses the two already-bundled primary SE JC tables; their uncompressed bytes match the original motif-run inputs. Saved exon-coordinate checksums and scientific settings are in `config/rmaps3.json`; the external motif library and hg38 genome are not redistributed. See [rMAPS3 execution](../docs/rmaps3.md).
