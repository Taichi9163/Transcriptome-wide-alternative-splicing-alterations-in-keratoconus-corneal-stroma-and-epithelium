# Transcriptome-wide alternative splicing alterations in keratoconus corneal stroma and epithelium

Analysis code accompanying the manuscript by Yuasa, Ueda and colleagues, version dated **18 November 2026**. The primary stromal comparison includes five keratoconus specimens and four controls; a sensitivity analysis adds KC6. The independent epithelial cohort includes seven specimens per group. Effects are defined as keratoconus minus control.

This repository combines the downstream analysis modules, the supplied Figure 1/Figure 3 code, and the Ubuntu read-processing commands. Paths are configurable. Original FASTQ/BAM files, manuscript documents, presentation files and gel images are not distributed here.

## Reproduce the manuscript tables and panels

Seven compressed processed input tables are included under `inputs/stroma/`: five complete rMATS JC tables and gene/transcript DESeq2 tables. Their uncompressed SHA-256 checksums are recorded in [inputs/manifest.json](inputs/manifest.json). The reference annotation is downloaded separately from Ensembl, with its checksum verified before event–transcript mapping.

Use Python 3.9 or later and R with these packages:

```r
install.packages(c("ggplot2", "ggforce", "ragg", "gplots", "RColorBrewer",
                   "dplyr", "tidyr", "purrr", "patchwork"))
```

From a clone of this repository:

```sh
# Use a separate directory for inputs expanded from the release and analysis outputs.
python3 scripts/prepare_inputs.py --work-dir work --download-reference

# Validate numeric tables and plot layers without creating images.
python3 scripts/run_manuscript.py --work-dir work --tables-only
python3 scripts/leave_one_out.py --work-dir work

# To generate Figure 1B-D and Figure 3A-D (PDF, PNG and TIFF):
python3 scripts/run_manuscript.py --work-dir work
```

To use an existing reference instead of downloading it, replace `--download-reference` with `--gtf /path/to/Homo_sapiens.GRCh38.114.gtf.gz`. Input preparation reuses identical files and refuses to overwrite different files. Set `FIGURE_FONT` to select an installed font; the default is `sans`. All generated files stay in the selected work directory. Scripts can be called by absolute path from another directory.

The manuscript workflow deliberately checks the expected study counts. It is a reproduction workflow for these inputs; adapting it to another cohort requires changing sample counts and the manuscript-specific checks.

| Output under `work/results/` | Content |
| --- | --- |
| `stroma/manuscript/` | DAS tables, PSI matrix, transcript matches and unmatched-event audit |
| `stroma/figure1/` | PCA scores/variance, clustering order and gene intersections for Fig. 1B–D |
| `stroma/figure3/` | Gene memberships, scatter source tables and Fig. 3A–D |
| `stromal_leave_one_out/` | ΔPSI after each single-sample omission; no new significance testing |
| `provenance/` | Installed R package versions and session information |

The complete transcript mapping retains 7,287 compatible event–transcript pairs. Of these, 4,047 pairs from 995 events have finite expression effects and raw P values and enter Fig. 3D. Forty-two pairs, representing 29 genes, involve DE isoforms. Both inclusion-compatible and skipping-compatible pairs retain the **original rMATS ΔPSI**; their shapes distinguish the two forms. See [analysis conventions](docs/analysis_conventions.md) for coordinate, filtering and matching details.

## Process reads on Ubuntu

The configurable workflow covers FastQC, MultiQC, fastp, STAR, RSEM, rMATS and sashimi plotting. Edit [config/upstream.json](config/upstream.json) and [config/samples_stroma.tsv](config/samples_stroma.tsv) for your files and tools.

```sh
# Preview all processing commands; no tools are run and no files are created.
python3 scripts/preprocess_reads.py all --work-dir /path/to/analysis

# Run one stage in an environment containing the required software.
python3 scripts/preprocess_reads.py trim --work-dir /path/to/analysis --execute
```

**Study sample KC6 corresponds to the original read prefix `KC_1`, whereas study sample KC3 corresponds to `KC_6`.** The supplied sample sheet records this mapping and includes all ten specimens for preprocessing. Only the nine primary specimens enter the primary rMATS comparison. The second rMATS comparison includes KC6. Read the [Ubuntu workflow guide](docs/upstream.md) before execution.

## Analyses by manuscript component

| Manuscript component | Scripts |
| --- | --- |
| Fig. 1A; DAS calls and class counts | `call_das.R` |
| Fig. 1B–D | `prepare_data.R`, `figure1BC.R`, `figure1D.R` |
| Fig. 1E,F; motif sensitivity | `motif_heatmaps.R`, `motif_tests.R` |
| Fig. 2A–C; sashimi source events | `selected_exons.R`, `preprocess_reads.py sashimi` |
| Fig. 2D–F; sample-level PSI | `selected_exons.R` |
| Fig. 3A–D | `figure3A.R`–`figure3D.R`, `map_events.py` |
| Gene/transcript model fitting from RSEM | `differential_expression.R`, `expression_overlap.R` |
| Fig. 4A,B; coordinate-matched cross-tissue comparison | `compare_das.R`, `compare_se.R`, `compare_se_support.R` |
| Fig. 4C,D; NEASE and Reactome | `nease_analysis.py`, `reactome_top10.py` |
| Supplementary sensitivity analyses | `stromal_sensitivity.py`, `leave_one_out.py`, `pca_all.R`, `pca_variable.R`, `pca_with_kc6.R` |
| Epithelial significant-event PCA/heatmap; RBP expression | `pca_significant.R`, `rbp_expression.R` |

The RT-PCR gel panels, their densitometry panels (Fig. 2G–L), and the conceptual summary (Fig. 4E) are not regenerated by this repository. The uploaded code bundle does not contain the original gel-intensity measurements. Source tables for the remaining modules can be imported from the paper's supplementary workbooks as described in [additional analyses](docs/additional_analyses.md).

## Validation and citation

The release check reproduced the primary DAS counts, expression-overlap counts, Figure 1 PCA variance and all Figure 3 event–transcript assignments. No figure images or raw-read analyses were regenerated during this check. [Validation details and limitations](docs/validation.md) distinguish numeric checks, regression tests and unexecuted upstream jobs.

```sh
python3 -m unittest discover -s tests -v
```

Use [CITATION.cff](CITATION.cff) for the manuscript authors and title, and record the Git commit used. The manuscript is currently described as unpublished; no DOI is asserted. The stromal dataset is identified by [GEO GSE336375](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336375) and the epithelial dataset by [BioProject PRJNA799648](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA799648).

Analysis code is distributed under the [MIT license](LICENSE).
