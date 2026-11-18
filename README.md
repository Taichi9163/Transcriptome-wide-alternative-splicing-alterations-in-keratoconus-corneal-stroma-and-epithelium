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

R filenames start with the corresponding `Figure` or `SuppFigure` number, followed by the analysis or panel content. Click a filename to open its code. Scripts that supply source tables are identified below; not every script draws the final assembled panel.

| Manuscript component | Scripts |
| --- | --- |
| Fig. 1A; DAS calls and class counts | [Figure1A_DAS_event_summary.R](scripts/Figure1A_DAS_event_summary.R) |
| Fig. 1B,C; PSI PCA and heatmap | [Figure1BC_PSI_PCA_heatmap.R](scripts/Figure1BC_PSI_PCA_heatmap.R) |
| Fig. 1D; gene-level UpSet | [Figure1D_DAS_gene_UpSet.R](scripts/Figure1D_DAS_gene_UpSet.R) |
| Fig. 1E,F; RBP motif heatmaps | [Figure1EF_RBP_motif_heatmaps.R](scripts/Figure1EF_RBP_motif_heatmaps.R) |
| Fig. 1E,F; motif statistics and sensitivity tables | [Figure1EF_RBP_motif_statistics.R](scripts/Figure1EF_RBP_motif_statistics.R) |
| Fig. 2A–F; source tables for sashimi and sample-level PSI | [Figure2DEF_selected_exon_PSI.R](scripts/Figure2DEF_selected_exon_PSI.R); `preprocess_reads.py sashimi` renders A–C |
| Fig. 3A; DEG/DAS gene Venn | [Figure3A_DEG_DAS_Venn.R](scripts/Figure3A_DEG_DAS_Venn.R) |
| Fig. 3B; gene expression versus ΔPSI | [Figure3B_gene_expression_vs_deltaPSI.R](scripts/Figure3B_gene_expression_vs_deltaPSI.R) |
| Fig. 3C; isoform/DAS gene Venn | [Figure3C_isoform_DAS_Venn.R](scripts/Figure3C_isoform_DAS_Venn.R) |
| Fig. 3D; matched isoform expression versus ΔPSI | [Figure3D_isoform_expression_vs_deltaPSI.R](scripts/Figure3D_isoform_expression_vs_deltaPSI.R), `map_events.py` |
| Fig. 3A–D; gene/transcript model fitting from RSEM | [Figure3ABCD_DESeq2_analysis.R](scripts/Figure3ABCD_DESeq2_analysis.R) |
| Fig. 3A–C; expression/splicing overlap source tables | [Figure3ABC_expression_splicing_overlap.R](scripts/Figure3ABC_expression_splicing_overlap.R) |
| Fig. 4A and Supplementary Fig. 8; coordinate-matched overlap tables | [Figure4A_SuppFigure8_DAS_overlap.R](scripts/Figure4A_SuppFigure8_DAS_overlap.R) |
| Fig. 4B; stromal versus epithelial ΔPSI | [Figure4B_cross_layer_deltaPSI.R](scripts/Figure4B_cross_layer_deltaPSI.R) |
| Fig. 4C,D and Supplementary Fig. 10; NEASE and Reactome | `nease_analysis.py`, `reactome_top10.py` |
| Supplementary Fig. 2 and Data 6; stromal sensitivity tables | `stromal_sensitivity.py`, `leave_one_out.py` |
| Supplementary Fig. 3A,B; candidate RBP expression tables | [SuppFigure3AB_RBP_expression.R](scripts/SuppFigure3AB_RBP_expression.R) |
| Supplementary Fig. 7B,C; epithelial PSI PCA and heatmap | [SuppFigure7BC_PSI_PCA_heatmap.R](scripts/SuppFigure7BC_PSI_PCA_heatmap.R) |
| Supplementary Fig. 9; SE junction-support sensitivity tables | [SuppFigure9_SE_junction_support.R](scripts/SuppFigure9_SE_junction_support.R) |

The DAS summary and motif scripts also support the epithelial analyses in Supplementary Fig. 7. The significant-event PCA helper accepts either cohort; the dedicated `Figure1BC` script reproduces the current stromal panels.

Shared preparation uses [Setup_check_R_packages.R](scripts/Setup_check_R_packages.R), [Setup_import_supplementary_data.R](scripts/Setup_import_supplementary_data.R), [Setup_prepare_Figure1_Figure3_data.R](scripts/Setup_prepare_Figure1_Figure3_data.R), and [Utils_plotting.R](scripts/Utils_plotting.R). Additional PCA sensitivity analyses use [Analysis_PCA_all_events.R](scripts/Analysis_PCA_all_events.R), [Analysis_PCA_variable_events.R](scripts/Analysis_PCA_variable_events.R), and [Analysis_PCA_including_KC6.R](scripts/Analysis_PCA_including_KC6.R); these optional analyses are not assigned to a numbered panel in the current manuscript.

The RT-PCR gel panels, their densitometry panels (Fig. 2G–L), and the conceptual summary (Fig. 4E) are not regenerated by this repository. The uploaded code bundle does not contain the original gel-intensity measurements. Source tables for the remaining modules can be imported from the paper's supplementary workbooks as described in [additional analyses](docs/additional_analyses.md).

## Validation and citation

The release check reproduced the primary DAS counts, expression-overlap counts, Figure 1 PCA variance and all Figure 3 event–transcript assignments. No figure images or raw-read analyses were regenerated during this check. [Validation details and limitations](docs/validation.md) distinguish numeric checks, regression tests and unexecuted upstream jobs.

```sh
python3 -m unittest discover -s tests -v
```

Use [CITATION.cff](CITATION.cff) for the manuscript authors and title, and record the Git commit used. The manuscript is currently described as unpublished; no DOI is asserted. The stromal dataset is identified by [GEO GSE336375](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336375) and the epithelial dataset by [BioProject PRJNA799648](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA799648).

Analysis code is distributed under the [MIT license](LICENSE).
