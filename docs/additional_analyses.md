# Additional analyses

These modules complement the bundled manuscript workflow. They require the paper's supplementary workbooks, native rMAPS3 outputs, or RSEM quantifications.

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
Rscript --vanilla scripts/import_supplement.R stroma "data/Supplementary Data 1.xlsx"
Rscript --vanilla scripts/import_supplement.R epithelium "data/Supplementary Data 2.xlsx"
Rscript --vanilla scripts/import_supplement.R motifs "data/Supplementary Data 3.xlsx"

Rscript --vanilla scripts/call_das.R stroma
Rscript --vanilla scripts/call_das.R epithelium
Rscript --vanilla scripts/call_das.R stroma_with_kc6
Rscript --vanilla scripts/compare_das.R
Rscript --vanilla scripts/compare_se.R --tables-only
Rscript --vanilla scripts/compare_se_support.R
python3 scripts/stromal_sensitivity.py
Rscript --vanilla scripts/motif_tests.R
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
Rscript --vanilla scripts/differential_expression.R gene data/stroma/genes.tsv results/stroma/gene
Rscript --vanilla scripts/differential_expression.R isoform data/stroma/isoforms.tsv results/stroma/isoform
Rscript --vanilla scripts/differential_expression.R gene data/epithelium/genes.tsv results/epithelium/gene
Rscript --vanilla scripts/expression_overlap.R
```

The model is `~ condition` with a keratoconus versus control Wald contrast. Isoforms require at least 10 counts summed across the nine stromal samples. The DESeq2 independent filtering target is 0.10 for genes and 0.05 for isoforms. Both call sets use adjusted P < 0.05 and |log2 fold change| > 1. The study analysis used DESeq2 1.52.0 and reported tximport 1.34.0. The release check used DESeq2 1.52.0 and tximport 1.40.0. All 414 gene IDs and all 448 transcript IDs matched Supplementary Data 4. Individual isoform log2 fold changes differed by up to 0.000028 from the workbook.

For candidate RBP expression, run `motif_tests.R` and the corresponding gene-level model first. Supply a tab-separated Ensembl release 114 annotation with columns `gene_id` and `gene_symbol`, using IDs matching the RSEM reference. The output contains DESeq2 adjusted P values and row-wise Z scores of log10(TPM + 1).

```sh
Rscript --vanilla scripts/rbp_expression.R stroma data/gene_symbols.tsv
Rscript --vanilla scripts/rbp_expression.R epithelium data/gene_symbols.tsv
```

### NEASE and Reactome

Use Python 3.12 and the package versions recorded for the NEASE analysis. Run after the two primary `call_das.R` commands. The code converts the rMATS exon start to a 1-based position, preserves the event identity during feature mapping, and disambiguates repeated pathway display names by reference row during the NEASE calculation. Reactome P values are BH-adjusted within each tissue.

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
Rscript --vanilla scripts/pca_significant.R epithelium --tables-only
Rscript --vanilla scripts/motif_heatmaps.R
Rscript --vanilla scripts/pca_all.R
Rscript --vanilla scripts/pca_variable.R
Rscript --vanilla scripts/pca_with_kc6.R
```

Remove `--tables-only` from `pca_significant.R` to generate its PCA and heatmap PDFs. Use the dedicated manuscript workflow for the current stromal Figure 1B–D panels.
