# Transcriptome-wide alternative splicing alterations in keratoconus corneal stroma and epithelium

**Analysis code and processed inputs accompanying the manuscript by Yuasa, Ueda and colleagues.**

[Quick start](#quick-start) · [Figure-to-code index](docs/figures.md) · [Reproduction guide](docs/reproduce.md) · [Data](#data) · [Citation](#citation)

This repository provides the RNA-seq and alternative-splicing analyses for the manuscript version dated **18 November 2026**: read processing, differential splicing, RBP motif enrichment, gene/isoform expression integration, cross-layer comparisons, and NEASE/Reactome analysis.

## Study design

| Analysis | Keratoconus | Control | Purpose |
| --- | ---: | ---: | --- |
| Corneal stroma | 5 | 4 | Primary analysis |
| Corneal epithelium | 7 | 7 | Independent cohort |
| Stroma including KC6 | 6 | 4 | Sensitivity analysis |

Effects are defined as **keratoconus minus control**. The reference is GRCh38 with Ensembl release 114. Event definitions, filtering, sample-label conventions and transcript matching are described in the [analysis conventions](docs/analysis_conventions.md) and [read-processing guide](docs/upstream.md).

## Quick start

The following workflow reproduces the numerical source tables and checks plot layers for **Figure 1B–D and Figure 3A–D**, using the bundled processed inputs. Python 3.9+ and R are required; release checks used R 4.6.1.

```sh
git clone https://github.com/Taichi9163/Transcriptome-wide-alternative-splicing-alterations-in-keratoconus-corneal-stroma-and-epithelium.git
cd Transcriptome-wide-alternative-splicing-alterations-in-keratoconus-corneal-stroma-and-epithelium
```

Install the R packages used by the manuscript panels:

```r
install.packages(c("ggplot2", "ggforce", "ragg", "gplots", "RColorBrewer",
                   "dplyr", "tidyr", "purrr", "patchwork", "VennDiagram", "scales"))
```

Prepare the inputs and run without creating figure files:

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work --download-reference
python3 scripts/workflow_run_manuscript.py --work-dir work --tables-only
```

The reference download is approximately 64 MB. To use an existing GTF, replace `--download-reference` with `--gtf /path/to/Homo_sapiens.GRCh38.114.gtf.gz`. Inputs are checksum-verified; outputs are written beneath `work/results/`.

The **[reproduction guide](docs/reproduce.md)** gives commands for the other figures, workflows that do not require the GTF, optional figure export, and software setup for native analyses. Scripts are intended to reproduce these study inputs and include manuscript-specific checks.

## Find the code for a figure

| Manuscript component | Analysis | Code index |
| --- | --- | --- |
| Figure 1 | Stromal DAS, PSI PCA/heatmap, UpSet and RBP motifs | [Figure 1](docs/figures.md#figure-1) |
| Figure 2 | Selected exon events, sashimi plots and RNA-seq PSI source tables | [Figure 2](docs/figures.md#figure-2) |
| Figure 3 | Gene/isoform expression and splicing integration | [Figure 3](docs/figures.md#figure-3) |
| Figure 4 | Cross-layer concordance, protein features and Reactome pathways | [Figure 4](docs/figures.md#figure-4) |
| Supplementary figures and tables | QC, sensitivity analyses, epithelial analyses and published-study comparison | [Supplementary material](docs/figures.md#supplementary-figures-and-tables) |

The index distinguishes analysis scripts, source-table exports and assembled or experimental panels. Figures 2D–F and 2J–L were plotted in GraphPad Prism. RT-PCR statistics, gel-intensity measurements and manual figure assembly are outside the scope of this RNA-seq code repository.

## Data

**28 compressed processed input tables** are included in [inputs/](inputs/README.md), with destinations, byte counts and SHA-256 checksums in the [manifest](inputs/manifest.json). They include complete rMATS JC tables, stromal DESeq2 results, motif-enrichment results, RBP expression matrices, NEASE mappings, corrected Reactome results and the published-event coordinate subset.

| Dataset | Accession |
| --- | --- |
| Stromal RNA-seq | [GEO GSE336375](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336375) |
| Epithelial RNA-seq | [BioProject PRJNA799648](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA799648) |

Raw FASTQ/BAM files, genome references, manuscript documents and experimental images are not bundled. Read-level analyses require the external tools and inputs described in [upstream processing](docs/upstream.md). Input provenance and the scope of each processed table are documented in [inputs/README.md](inputs/README.md).

## Repository guide

| Location | Contents |
| --- | --- |
| [scripts/](scripts/) | Numbered figure scripts, shared helpers and native-tool wrappers; [naming convention](docs/script_naming.md) |
| [config/](config/) | Sample annotations, paths and recorded analysis settings |
| [inputs/](inputs/) | Processed inputs and checksum manifest |
| [docs/](docs/README.md) | Reproduction instructions, methods and provenance |
| [tests/](tests/) | Computational checks and numeric-only validation utilities |

The [validation record](docs/validation.md) reports agreement with saved study results and distinguishes numerical checks, plot-object checks and unexecuted upstream analyses. It also records reconstruction and environment limitations; full raw-read and native motif analyses were not rerun for this repository release.

## Citation

Please cite the associated manuscript and record the Git commit used in your analysis. The complete author list and citation metadata are provided in [CITATION.cff](CITATION.cff), also available through GitHub's **Cite this repository** menu.

> Yuasa T, Ueda N, Yamashita Y, et al. *Transcriptome-wide alternative splicing alterations in keratoconus corneal stroma and epithelium.* Manuscript, 2026.

The manuscript is currently recorded as unpublished; no article DOI is assigned in this repository. For code questions or reproducibility issues, use [GitHub Issues](https://github.com/Taichi9163/Transcriptome-wide-alternative-splicing-alterations-in-keratoconus-corneal-stroma-and-epithelium/issues).

## License

Analysis code is distributed under the [MIT license](LICENSE). External tools and third-party data retain their own licenses and attribution requirements; see [input provenance](inputs/README.md).
