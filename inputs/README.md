# Processed input provenance

All bundled files are losslessly gzip-compressed. `manifest.json` records each destination, uncompressed byte count and SHA-256 checksum. `scripts/prepare_inputs.py` verifies these values and refuses to replace different existing inputs.

| Input | Source and purpose |
| --- | --- |
| `stroma/rmats/*.MATS.JC.txt.gz` | Five complete stromal rMATS JC tables from the supplied Figure 1/Figure 3 archive; primary five-case/four-control comparison |
| `stroma/expression/*_deseq2.tsv.gz` | Gene and transcript DESeq2 result tables from the same supplied archive |
| `epithelium/rmats/*.MATS.JC.txt.gz` | Five complete JC tables from the local Sun2022 `after_rMATS/results` analysis, seven samples per group; input to Supplementary Figure 7D |
| `<cohort>/rmaps/pVal.*.vs.bg.RNAmap.txt.gz` | Original stromal (Yamashita2026) and epithelial (Sun2022) rMAPS3 output tables, 114 motifs × eight regions per direction; candidate selection is recomputed from these P values |
| `<cohort>/rbp/tpm.tsv.gz` | Candidate-gene TPM values and Ensembl IDs from the saved outputs of the recovered RBP heatmap scripts; 19 genes × nine samples for stroma and 23 × 14 for epithelium |

The recovered source scripts were `stroma_rMAPS_hit_RNA_heatmap_BH.R` and `epithelium_Sun2022_rMAPS_hit_RNA_heatmap_BH.R` in the project's `RBP_expression_heatmap/scripts/` directory. Their saved `*_TPM.csv` values and `*_metadata.csv` gene IDs are combined into each bundled `tpm.tsv`. Sample labels and their order are retained from those exports. These matrices contain the candidate genes only; they do not replace full RSEM quantifications or support a new genome-wide expression analysis.

The original stromal heatmap read the `Gene Matrix` sheet of the marker-expression workbook. The epithelial heatmap read the Sun2022 `DataMatrix_ImportSampleList.txt` abundance columns. Gene-symbol lookups came from the study DESeq2 workbooks, including the original stromal MSI1/HNRNPCL1 lookup entries. Bundling the saved numeric matrices removes those private workbook/path dependencies without changing the plotted values.

Supplementary Figure 3 candidate selection preserves BH correction across 912 motif–region P values separately per cohort and direction, the strict cutoff of 5e-5, and the HuR → ELAVL1 and SF2-ASF → SRSF1 aliases. The heatmaps preserve log2(TPM + 1), row Z-scores, Euclidean/Ward.D2 row clustering, and fixed control-then-case column order. Row Z-scores are unchanged by using log10 instead of log2. No DESeq2 adjusted P values are inferred from TPM; the separate expression-statistics script requires full gene-model results.

Supplementary Figure 7D is reconstructed from the Figure 1D code and epithelial inputs. It uses FDR < 0.05, |ΔPSI| > 0.10, and one gene per event class (gene symbol with GeneID fallback). The total is 805 events in 609 genes, including 314 genes affected only by SE. The original epithelial UpSet source was not located. This release validates numeric data and plot construction without regenerating or visually comparing figure images.
