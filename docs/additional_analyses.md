# Additional analyses

These modules complement the bundled manuscript workflow. Supplementary Figures 3A,B and 7D can run directly from the bundled processed inputs, as shown in the main README. Other modules require the paper's supplementary workbooks, additional native outputs, or RSEM quantifications.

## Dependencies and native inputs

In addition to the main workflow packages, install the dependencies for the modules you use:

```r
install.packages(c("readxl", "data.table", "ggrepel", "reshape2", "stringr", "BiocManager"))
BiocManager::install(c("DESeq2", "tximport"))
```

| Directory | Inputs |
| --- | --- |
| `data/stroma/rmats/`, `data/epithelium/rmats/`, `data/stroma_with_kc6/rmats/` | Five complete `<event>.MATS.JC.txt` files per comparison |
| `data/stroma/rmaps/`, `data/epithelium/rmaps/` | `pVal.up.vs.bg.RNAmap.txt`, `pVal.dn.vs.bg.RNAmap.txt`; include `log.motifMap.txt` when available |
| Paths supplied in expression sample sheets | RSEM `.genes.results` or `.isoforms.results` files |

The manuscript reports rMAPS3 3.1.0, hg38 SE workflow, human/mouse motifs, a 50 nt window and a 1 nt step. Preserve those settings when recreating the motif inputs.

## Supplementary workbooks

Run these commands from the repository root, or use absolute script paths from your work directory. Keep each analysis run in its own work directory. Place the paper's Supplementary Data 1, 2 and 3 workbooks in `data/` with the names shown below. The importer converts the worksheets to tab-delimited files and preserves the order of the sample-level PSI vectors. Native inputs use the directory layout below.

```sh
Rscript --vanilla scripts/Setup_import_supplementary_data.R stroma "data/Supplementary Data 1.xlsx"
Rscript --vanilla scripts/Setup_import_supplementary_data.R epithelium "data/Supplementary Data 2.xlsx"
Rscript --vanilla scripts/Setup_import_supplementary_data.R motifs "data/Supplementary Data 3.xlsx"

Rscript --vanilla scripts/Figure1A_DAS_event_summary.R stroma
Rscript --vanilla scripts/Figure1A_DAS_event_summary.R epithelium
Rscript --vanilla scripts/Figure1A_DAS_event_summary.R stroma_with_kc6
Rscript --vanilla scripts/Figure4A_SuppFigure8_DAS_overlap.R
Rscript --vanilla scripts/Figure4B_cross_layer_deltaPSI.R --tables-only
Rscript --vanilla scripts/SuppFigure9_SE_junction_support.R
python3 scripts/stromal_sensitivity.py
Rscript --vanilla scripts/Figure1EF_RBP_motif_statistics.R
```

The five rMATS event classes are SE, A5SS, A3SS, MXE and RI. Each call requires FDR < 0.05 and |IncLevelDifference| > 0.10 in the JC table. Unique gene counts use rMATS GeneID. Numeric rMATS event IDs apply within one run and event class. Cross-run matching uses chromosome, strand and all coordinates specific to that event class.

| Check from study inputs | Expected result |
| --- | --- |
| Primary stroma, all classes | 1,073 events in 807 genes |
| Epithelium, all classes | 805 events in 609 genes |
| Stroma including KC6 | 1,074 events in 802 genes, with 889 events shared with the primary calls |
| Junction support in primary stroma | 768 events have finite PSI and at least 10 inclusion + skipping junction reads in ≥3 cases and ≥2 controls |
| Stromal expression | 414 DE genes (8 with DAS) and 448 DE isoforms in 410 genes (55 with DAS) |
| SE union across tissues | 1,135 significant in either tissue, 891 estimable in both, 13 significant in both |
| rMAPS3 sensitivity | 3,648 motif–region tests and 14 pairs below the display cutoff after window correction |

The SE scatter retains measured zeroes. An event needs a coordinate match, a finite effect and at least one finite PSI in each group of each cohort. The rMAPS3 figure uses BH across the 912 regional minima within each tissue and direction. The sensitivity calculation first corrects for the 50 exonic or 250 intronic window positions, then applies BH or BY across those 912 tests.

### Gene and isoform expression

Use tab-delimited sample sheets with columns `sample`, `group`, and `path`. Each row identifies one biological sample, its group (`Control` or `Keratoconus`), and an RSEM result file. Paths are relative to the current work directory or absolute. Use the study sample annotation to identify the nine primary stromal specimens. Group membership and the KC6 exclusion must follow the study labels. For example, a row for study sample KC1 can point to `data/stroma/rsem/KC1.genes.results` if that is how the quantification file was named.

```sh
Rscript --vanilla scripts/Figure3ABCD_DESeq2_analysis.R gene data/stroma/genes.tsv results/stroma/gene
Rscript --vanilla scripts/Figure3ABCD_DESeq2_analysis.R isoform data/stroma/isoforms.tsv results/stroma/isoform
Rscript --vanilla scripts/Figure3ABCD_DESeq2_analysis.R gene data/epithelium/genes.tsv results/epithelium/gene
Rscript --vanilla scripts/Figure3ABC_expression_splicing_overlap.R
```

The model is `~ condition` with a keratoconus versus control Wald contrast. Isoforms require at least 10 counts summed across the nine stromal samples. The DESeq2 independent filtering target is 0.10 for genes and 0.05 for isoforms. Both call sets use adjusted P < 0.05 and |log2 fold change| > 1. The study analysis used DESeq2 1.52.0 and reported tximport 1.34.0. The release check used DESeq2 1.52.0 and tximport 1.40.0. All 414 gene IDs and all 448 transcript IDs matched Supplementary Data 4. Individual isoform log2 fold changes differed by up to 0.000028 from the workbook.

For the recovered candidate-RBP heatmaps, use `SuppFigure3AB_RBP_expression_heatmaps.R` with the bundled inputs. It preserves the original heatmap selection, scaling and clustering; `--tables-only` writes numeric matrices and row order without opening graphics devices. See [input provenance](../inputs/README.md) for the source matrices and their scope.

For the separate candidate-expression table with DESeq2 statistics, run `Figure1EF_RBP_motif_statistics.R` and the corresponding gene-level model first. Supply a tab-separated Ensembl release 114 annotation with columns `gene_id` and `gene_symbol`, using IDs matching the RSEM reference. The output contains DESeq2 adjusted P values and row-wise Z scores of log10(TPM + 1).

```sh
Rscript --vanilla scripts/SuppFigure3AB_RBP_expression.R stroma data/gene_symbols.tsv
Rscript --vanilla scripts/SuppFigure3AB_RBP_expression.R epithelium data/gene_symbols.tsv
```

### NEASE and Reactome

Use Python 3.12 and the package versions recorded for the NEASE analysis. Run after the two primary `Figure1A_DAS_event_summary.R` commands. The code converts the rMATS exon start to a 1-based position, preserves the event identity during feature mapping, and disambiguates repeated pathway display names by reference row during the NEASE calculation. Reactome P values are BH-adjusted within each tissue.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-nease.txt
.venv/bin/python scripts/nease_analysis.py stroma
.venv/bin/python scripts/nease_analysis.py epithelium
.venv/bin/python scripts/reactome_top10.py
```

The primary run keeps all input SE events for mapping. Pass `--orf-filter` to run the additional ORF-filtered sensitivity analysis in a separate output directory.

## Optional plots and PCA sensitivity

These commands generate figure files; the all-event PCA also writes larger caches. The three sensitivity scripts depend on the preceding scripts in the order shown.

```sh
Rscript --vanilla scripts/SuppFigure7BC_PSI_PCA_heatmap.R epithelium --tables-only
Rscript --vanilla scripts/Figure1EF_RBP_motif_heatmaps.R
Rscript --vanilla scripts/Analysis_PCA_all_events.R
Rscript --vanilla scripts/Analysis_PCA_variable_events.R
Rscript --vanilla scripts/Analysis_PCA_including_KC6.R
```

Remove `--tables-only` from `SuppFigure7BC_PSI_PCA_heatmap.R` to generate its PCA and heatmap PDFs. Use the dedicated manuscript workflow for the current stromal Figure 1B–D panels.
