# Figure-to-code index

[Repository overview](../README.md) · [Reproduction guide](reproduce.md) · [Input provenance](../inputs/README.md) · [Validation](validation.md)

Computational sources for the **18 November 2026 manuscript** are listed below. Source-table scripts provide values used in a panel; they do not necessarily export the final assembled layout. Figure filenames follow the manuscript numbering and the [shared naming convention](script_naming.md).

[Figure 1](#figure-1) · [Figure 2](#figure-2) · [Figure 3](#figure-3) · [Figure 4](#figure-4) · [Supplementary material](#supplementary-figures-and-tables) · [Supplementary data](#supplementary-data)

## Figure 1

| Manuscript component | Scripts |
| --- | --- |
| Fig. 1A; DAS calls and class counts | [figure_01a_das_event_summary.R](../scripts/figure_01a_das_event_summary.R) |
| Fig. 1B,C; PSI PCA and heatmap | [figure_01bc_psi_pca_heatmap.R](../scripts/figure_01bc_psi_pca_heatmap.R) |
| Fig. 1D; gene-level UpSet | [figure_01d_das_gene_upset.R](../scripts/figure_01d_das_gene_upset.R) |
| Fig. 1E,F; RBP motif heatmaps | [figure_01ef_rbp_motif_heatmaps.R](../scripts/figure_01ef_rbp_motif_heatmaps.R) |
| Fig. 1E,F and Supplementary Fig. 7E,F; native rMAPS3 execution and exon-selection audit | [workflow_run_rmaps3.py](../scripts/workflow_run_rmaps3.py), [settings](../config/rmaps3.json) |
| Fig. 1E,F; motif statistics and sensitivity tables | [figure_01ef_rbp_motif_statistics.R](../scripts/figure_01ef_rbp_motif_statistics.R) |

Figure 1B–D is included in `workflow_run_manuscript.py`. Motif modules are separate; see the [rMAPS3 guide](rmaps3.md), including the native ΔPSI boundary and external-engine requirements.

## Figure 2

| Manuscript component | Scripts |
| --- | --- |
| Fig. 2A–C; sashimi plots | [workflow_preprocess_reads.py](../scripts/workflow_preprocess_reads.py), `sashimi` stage; event selection is supplied by [figure_02def_selected_exon_psi.R](../scripts/figure_02def_selected_exon_psi.R) |
| Fig. 2D–F; RNA-seq PSI graphs prepared in GraphPad Prism | [figure_02def_selected_exon_psi.R](../scripts/figure_02def_selected_exon_psi.R) exports sample-level PSI source tables; final graphing uses Prism |
| Fig. 2G–I; RT-PCR gel panels | Experimental images assembled separately |
| Fig. 2J–L; RT-PCR inclusion-fraction graphs prepared in GraphPad Prism | GraphPad Prism workflow; no custom plotting script required |

Sashimi plotting requires STAR-aligned BAM files. The selected-event script supplies the three events and their RNA-seq PSI values. RT-PCR measurements and statistics are outside this repository’s scope.

## Figure 3

| Manuscript component | Scripts |
| --- | --- |
| Fig. 3A; DEG/DAS gene Venn | [figure_03a_deg_das_venn.R](../scripts/figure_03a_deg_das_venn.R) |
| Fig. 3B; gene expression versus ΔPSI | [figure_03b_gene_expression_vs_delta_psi.R](../scripts/figure_03b_gene_expression_vs_delta_psi.R) |
| Fig. 3C; isoform/DAS gene Venn | [figure_03c_isoform_das_venn.R](../scripts/figure_03c_isoform_das_venn.R) |
| Fig. 3D; matched isoform expression versus ΔPSI | [figure_03d_isoform_expression_vs_delta_psi.R](../scripts/figure_03d_isoform_expression_vs_delta_psi.R), [figure_03d_event_transcript_mapping.py](../scripts/figure_03d_event_transcript_mapping.py) |
| Fig. 3A–D; gene/transcript model fitting from RSEM | [figure_03abcd_deseq2_analysis.R](../scripts/figure_03abcd_deseq2_analysis.R) |
| Fig. 3A–C; expression/splicing overlap source tables | [figure_03abc_expression_splicing_overlap.R](../scripts/figure_03abc_expression_splicing_overlap.R) |

`workflow_run_manuscript.py` reproduces Figure 3A–D from the bundled DESeq2 tables and the release-114 GTF. Refitting expression models requires RSEM quantifications and sample annotations; see [gene and isoform expression](additional_analyses.md#gene-and-isoform-expression).

## Figure 4

| Manuscript component | Scripts |
| --- | --- |
| Fig. 4A and Supplementary Fig. 8; coordinate-matched overlap tables and four-set Venn plots | [figure_04a_supp_figure_08_das_overlap.R](../scripts/figure_04a_supp_figure_08_das_overlap.R) |
| Fig. 4B; stromal versus epithelial ΔPSI | [figure_04b_cross_layer_delta_psi.R](../scripts/figure_04b_cross_layer_delta_psi.R) |
| Fig. 4C; stromal protein-feature coverage and overlaps | [figure_04c_nease_feature_mapping.R](../scripts/figure_04c_nease_feature_mapping.R) |
| Fig. 4D; stromal Reactome top-ten R panel | [figure_04d_reactome_enrichment.R](../scripts/figure_04d_reactome_enrichment.R) |
| Fig. 4C,D and Supplementary Fig. 10; upstream NEASE calculation and revised Python plot implementation | [analysis_nease.py](../scripts/analysis_nease.py), [analysis_reactome_top10.py](../scripts/analysis_reactome_top10.py) |
| Fig. 4E; conceptual summary | Separately assembled artwork |

Figure 4C,D R scripts use bundled corrected NEASE outputs. New NEASE calculations require the [external environment](additional_analyses.md#nease-and-reactome). The retained [analysis_reactome_top10.py](../scripts/analysis_reactome_top10.py) is the revised Python plotting implementation; see [Reactome provenance](reactome.md) for the R adaptation.

## Supplementary figures and tables

| Manuscript component | Scripts |
| --- | --- |
| Supplementary Fig. 1A; QC report generation | [workflow_preprocess_reads.py](../scripts/workflow_preprocess_reads.py): `raw-qc`, `trim`, `trimmed-qc`, `align`; see [QC coverage](qc.md) for the final panels |
| Supplementary Fig. 1B; raw, filtered and mapped read pairs | [supp_figure_01b_qc_read_pairs.py](../scripts/supp_figure_01b_qc_read_pairs.py); supplied reconstruction, with a numeric-only mode |
| Supplementary Fig. 2 and Data 6; stromal sensitivity tables | [supp_data_06_stromal_sensitivity.py](../scripts/supp_data_06_stromal_sensitivity.py), [supp_data_06_leave_one_out.py](../scripts/supp_data_06_leave_one_out.py) |
| Supplementary Fig. 2B; KC6-inclusive versus primary event Venn | [supp_figure_02b_kc6_event_overlap.R](../scripts/supp_figure_02b_kc6_event_overlap.R); also exports class counts for Supplementary Fig. 2A |
| Supplementary Fig. 3A,B; candidate RBP expression heatmaps | [supp_figure_03ab_rbp_expression_heatmaps.R](../scripts/supp_figure_03ab_rbp_expression_heatmaps.R) |
| Supplementary Fig. 3A,B; expression tables with DESeq2 adjusted P values | [supp_figure_03ab_rbp_expression.R](../scripts/supp_figure_03ab_rbp_expression.R) |
| Supplementary Figs. 4 and 5; RT-PCR gels and controls | Experimental panels assembled separately |
| Supplementary Fig. 6; exact-coordinate comparison with published Sun et al. events | [supp_figure_06_published_event_comparison.R](../scripts/supp_figure_06_published_event_comparison.R) |
| Supplementary Fig. 7A; epithelial DAS calls and class counts | [figure_01a_das_event_summary.R](../scripts/figure_01a_das_event_summary.R), `epithelium` argument |
| Supplementary Fig. 7B,C; epithelial PSI PCA and heatmap | [supp_figure_07bc_psi_pca_heatmap.R](../scripts/supp_figure_07bc_psi_pca_heatmap.R) |
| Supplementary Fig. 7D; epithelial gene-level UpSet | [supp_figure_07d_das_gene_upset.R](../scripts/supp_figure_07d_das_gene_upset.R) |
| Supplementary Fig. 7E,F; epithelial motif analysis | [workflow_run_rmaps3.py](../scripts/workflow_run_rmaps3.py), [motif statistics](../scripts/figure_01ef_rbp_motif_statistics.R), [motif heatmaps](../scripts/figure_01ef_rbp_motif_heatmaps.R) |
| Supplementary Fig. 8A–D; A3SS/A5SS/MXE/RI cross-layer overlaps | [figure_04a_supp_figure_08_das_overlap.R](../scripts/figure_04a_supp_figure_08_das_overlap.R) |
| Supplementary Fig. 9; additional SE junction-support sensitivity tables | [supp_figure_09_se_junction_support.R](../scripts/supp_figure_09_se_junction_support.R) |
| Supplementary Fig. 9A,B; support-filtered SE Venn and scatter | [supp_figure_09ab_se_junction_support_plots.R](../scripts/supp_figure_09ab_se_junction_support_plots.R) |
| Supplementary Fig. 10A; epithelial protein-feature coverage and overlaps | [supp_figure_10a_nease_feature_mapping.R](../scripts/supp_figure_10a_nease_feature_mapping.R) |
| Supplementary Fig. 10B; epithelial Reactome top-ten R panel | [supp_figure_10b_reactome_enrichment.R](../scripts/supp_figure_10b_reactome_enrichment.R) |
| Supplementary Table 1; RT-PCR primers | Experimental metadata |
| Supplementary Table 2; STAR unique-mapping and junction-support columns | [supp_table_02_star_junction_qc.py](../scripts/supp_table_02_star_junction_qc.py), [sample metadata](../config/samples_stroma_qc.tsv) |

Figure S1B is the supplied reconstruction; Figure S7D adapts the Figure 1D implementation. The Venn layouts are portable reconstructions. Numerical checks and provenance are described in [validation](validation.md), [comparison panels](comparison_panels.md) and [published-event comparison](suppfigure6.md); visual identity is not asserted. Supplementary Figure 2A uses the class-count output of `supp_figure_02b_kc6_event_overlap.R` with separately assembled exon schematics. The [QC guide](qc.md#coverage) identifies which Supplementary Table 2 columns are calculated by the supplied scripts.

## Supplementary data

| Data | Analysis source |
| --- | --- |
| 1, 2 | Complete stromal/epithelial rMATS JC outputs; [DAS summary](../scripts/figure_01a_das_event_summary.R) and [upstream commands](upstream.md) |
| 3 | [Motif multiplicity and sensitivity](../scripts/figure_01ef_rbp_motif_statistics.R) |
| 4 | [DESeq2](../scripts/figure_03abcd_deseq2_analysis.R), [expression/splicing overlap](../scripts/figure_03abc_expression_splicing_overlap.R), [event–transcript matching](../scripts/figure_03d_event_transcript_mapping.py) |
| 5 | [NEASE](../scripts/analysis_nease.py), [feature summaries](../scripts/utils_nease_feature_summary.R), [Reactome selection](../scripts/utils_reactome_enrichment.R) |
| 6 | [Stromal sensitivity](../scripts/supp_data_06_stromal_sensitivity.py), [leave-one-out effects](../scripts/supp_data_06_leave_one_out.py) |

These modules export analysis tables underlying the supplementary material; they do not assemble the final formatted Excel workbooks.

## Shared preparation and optional analyses

The primary workflow uses [workflow_run_manuscript.py](../scripts/workflow_run_manuscript.py), [setup_prepare_inputs.py](../scripts/setup_prepare_inputs.py), [setup_check_r_packages.R](../scripts/setup_check_r_packages.R), [setup_prepare_manuscript_data.R](../scripts/setup_prepare_manuscript_data.R), and [utils_plotting.R](../scripts/utils_plotting.R). [setup_import_supplementary_data.R](../scripts/setup_import_supplementary_data.R) imports supplementary workbooks for additional modules. Shared counting and plot helpers remain in `scripts/`.

Optional PCA analyses are [analysis_pca_all_events.R](../scripts/analysis_pca_all_events.R), [analysis_pca_variable_events.R](../scripts/analysis_pca_variable_events.R), and [analysis_pca_including_kc6.R](../scripts/analysis_pca_including_kc6.R). They are not assigned to numbered panels in this manuscript; their dependencies and execution order are in [additional analyses](additional_analyses.md#optional-plots-and-pca-sensitivity).
