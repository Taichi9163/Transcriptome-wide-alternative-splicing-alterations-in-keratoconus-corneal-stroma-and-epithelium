# Alternative splicing in keratoconus corneal stroma and epithelium

Analysis code for *Transcriptome-wide alternative splicing alterations in keratoconus corneal stroma and epithelium*. The primary stromal comparison uses five keratoconus specimens and four controls. A sensitivity analysis adds the sixth keratoconus specimen (KC6). The independent epithelial comparison uses seven specimens per group. All effect sizes use keratoconus minus control.

## Start with the supplementary data

Run the commands from the repository root in a new working copy. Place the paper's Supplementary Data 1, 2 and 3 workbooks in `data/` with the names shown below. The importer converts the worksheets to tab-delimited files and preserves the order of the sample-level PSI vectors. If you already have the native rMATS or rMAPS3 outputs, use the directory layout in the next section.

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

## Native inputs and software

The stromal dataset is identified by [GEO GSE336375](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336375). The epithelial dataset is identified by [BioProject PRJNA799648](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA799648).

| Directory | Files |
| --- | --- |
| `data/stroma/rmats/`, `data/epithelium/rmats/`, `data/stroma_with_kc6/rmats/` | Five complete `<event>.MATS.JC.txt` tables per comparison |
| `data/stroma/rmaps/`, `data/epithelium/rmaps/` | `pVal.up.vs.bg.RNAmap.txt` and `pVal.dn.vs.bg.RNAmap.txt` per tissue. Include the native `log.motifMap.txt` if available to verify scan settings. |
| User-supplied paths in sample sheets | One RSEM `*.genes.results` or `*.isoforms.results` file per sample. Use a separate sample sheet for each cohort and analysis level. |

Upstream processing used fastp 0.20.0, STAR 2.7.11a, GRCh38 with Ensembl release 114, RSEM 1.3.1 and rMATS-turbo 4.3.0. The rMATS runs used paired-end data, 150 nt read length, variable-length reads, clipped alignments, fr-unstranded libraries and the default 1 nt anchor. Supply case BAMs as group 1 and controls as group 2. rMAPS3 3.1.0 used the hg38 SE workflow, the human/mouse motif library, a 50 nt window and a 1 nt step. The repository starts from processed tables. It does not contain the original read alignment jobs, sashimi commands, gel images or the original Figure 3D event-to-transcript matching program.

R 4.6.1 was used to check the table analyses. DAS calling and the SE comparison in tables-only mode use base R. Install the additional packages for the modules you intend to run.

```r
install.packages(c("readxl", "data.table", "ggplot2", "gplots", "ggrepel", "patchwork",
                   "ragg", "RColorBrewer", "tidyr", "dplyr", "stringr", "reshape2",
                   "BiocManager"))
BiocManager::install(c("DESeq2", "tximport"))
```

## Analysis by figure

| Paper component | Script and output |
| --- | --- |
| Fig. 1A–D and epithelial DAS | `call_das.R` writes per-class calls, gene memberships and counts to `results/<cohort>/das/`. `pca_significant.R` writes a PCA and heatmap of selected events. Clustering uses Euclidean distances and Ward.D2 on PSI, with row-scaled PSI used for heatmap colors. |
| Fig. 1E,F and Supplementary Data 3 | `motif_heatmaps.R` renders the exploratory motif maps. `motif_tests.R` writes both the original and window-adjusted results under `results/motif_tests/`. |
| Fig. 2D–F and Supplementary Data 6 | `selected_exons.R` exports sample-level PSI and rMATS FDR for the three stromal SE targets. `stromal_sensitivity.py` exports the junction-support and KC6 comparisons by event and class. |
| Fig. 3A–C and Supplementary Data 4 | `differential_expression.R` fits the gene and transcript models. `expression_overlap.R` exports gene memberships and one largest-absolute-effect event per DAS gene for the gene-level scatter. Structural event–transcript pairing in Fig. 3D is outside this module. |
| Fig. 4A,B | `compare_das.R` exports exact-coordinate unions for all classes. `compare_se.R` exports SE audit tables and, when called without `--tables-only`, the scatter plot. `compare_se_support.R` checks the effect values and support thresholds. |
| Fig. 4C,D and Supplementary Data 5 | `nease_analysis.py` exports event–feature assignments, summary counts and pathway tests. `reactome_top10.py` ranks significant Reactome pathways by NEASE score and draws the top ten per tissue. |

The significant-event PCA describes the same events selected using the group labels. For PCA on complete, varying events independent of significance, run the three optional scripts in the order below. They write figure files and larger PCA caches under `results/`.

```sh
Rscript --vanilla scripts/pca_significant.R stroma
Rscript --vanilla scripts/pca_significant.R epithelium
Rscript --vanilla scripts/selected_exons.R
Rscript --vanilla scripts/motif_heatmaps.R
Rscript --vanilla scripts/pca_all.R
Rscript --vanilla scripts/pca_variable.R
Rscript --vanilla scripts/pca_with_kc6.R
```

### Gene and isoform expression

Use tab-delimited sample sheets with columns `sample`, `group`, and `path`. Each row identifies one biological sample, its group (`Control` or `Keratoconus`), and an RSEM result file. Paths are relative to the repository root or absolute. Use the study sample annotation to identify the nine primary stromal specimens. Group membership and the KC6 exclusion must follow the study labels. For example, a row for study sample KC1 can point to `data/stroma/rsem/KC1.genes.results` if that is how the quantification file was named.

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
