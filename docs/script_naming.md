# Script naming convention

[Repository overview](../README.md) · [Figure-to-code index](figures.md) · [Reproduction guide](reproduce.md)

All script basenames use **lowercase snake_case**. Figure, supplementary figure, table and data numbers use **two digits**, so files sort in numerical order. Panel letters follow the figure number in lowercase (`01a`, `01bc`). Extensions identify the language: `.R` for R and `.py` for Python.

| Role | Pattern | Example |
| --- | --- | --- |
| Main figure | `figure_<number><panels>_<purpose>` | `figure_01a_das_event_summary.R` |
| Supplementary figure | `supp_figure_<number><panels>_<purpose>` | `supp_figure_10b_reactome_enrichment.R` |
| Supplementary table | `supp_table_<number>_<purpose>` | `supp_table_02_star_junction_qc.py` |
| Supplementary data | `supp_data_<number>_<purpose>` | `supp_data_06_leave_one_out.py` |
| Input preparation and dependency checks | `setup_<purpose>` | `setup_prepare_inputs.py` |
| Multi-step or external-tool execution | `workflow_<purpose>` | `workflow_run_manuscript.py` |
| Analyses shared across panels or optional analyses | `analysis_<purpose>` | `analysis_nease.py` |
| Shared helpers | `utils_<purpose>` | `utils_plotting.R` |
| Tests | `test_<purpose>` | `test_reactome_panels.R` |

Names describe the analysis purpose. Acronyms in filenames are lowercase (`psi`, `pca`, `deseq2`, `nease`, `qc`); their conventional spelling is retained in manuscript text and plot labels. A script serving two numbered figures includes both identifiers, for example `figure_04a_supp_figure_08_das_overlap.R`. Shared analyses retain one implementation rather than duplicate files for each panel.

The rename updates script imports, helper sources, workflow entry points, tests and documented commands. Scientific calculations, input tables, result-directory names and figure-output names are retained. Original filenames used solely to identify historical source material are kept in the provenance documentation.

## Migration from previous script names

Use the current paths below in new commands. Old-name wrapper scripts are not retained; previous versions remain available in Git history.

<details>
<summary>Show the complete filename mapping</summary>

| Previous filename | Current script |
| --- | --- |
| `Analysis_PCA_all_events.R` | [analysis_pca_all_events.R](../scripts/analysis_pca_all_events.R) |
| `Analysis_PCA_including_KC6.R` | [analysis_pca_including_kc6.R](../scripts/analysis_pca_including_kc6.R) |
| `Analysis_PCA_variable_events.R` | [analysis_pca_variable_events.R](../scripts/analysis_pca_variable_events.R) |
| `Figure1A_DAS_event_summary.R` | [figure_01a_das_event_summary.R](../scripts/figure_01a_das_event_summary.R) |
| `Figure1BC_PSI_PCA_heatmap.R` | [figure_01bc_psi_pca_heatmap.R](../scripts/figure_01bc_psi_pca_heatmap.R) |
| `Figure1D_DAS_gene_UpSet.R` | [figure_01d_das_gene_upset.R](../scripts/figure_01d_das_gene_upset.R) |
| `Figure1EF_RBP_motif_heatmaps.R` | [figure_01ef_rbp_motif_heatmaps.R](../scripts/figure_01ef_rbp_motif_heatmaps.R) |
| `Figure1EF_RBP_motif_statistics.R` | [figure_01ef_rbp_motif_statistics.R](../scripts/figure_01ef_rbp_motif_statistics.R) |
| `Figure2DEF_selected_exon_PSI.R` | [figure_02def_selected_exon_psi.R](../scripts/figure_02def_selected_exon_psi.R) |
| `Figure3ABCD_DESeq2_analysis.R` | [figure_03abcd_deseq2_analysis.R](../scripts/figure_03abcd_deseq2_analysis.R) |
| `Figure3ABC_expression_splicing_overlap.R` | [figure_03abc_expression_splicing_overlap.R](../scripts/figure_03abc_expression_splicing_overlap.R) |
| `Figure3A_DEG_DAS_Venn.R` | [figure_03a_deg_das_venn.R](../scripts/figure_03a_deg_das_venn.R) |
| `Figure3B_gene_expression_vs_deltaPSI.R` | [figure_03b_gene_expression_vs_delta_psi.R](../scripts/figure_03b_gene_expression_vs_delta_psi.R) |
| `Figure3C_isoform_DAS_Venn.R` | [figure_03c_isoform_das_venn.R](../scripts/figure_03c_isoform_das_venn.R) |
| `Figure3D_isoform_expression_vs_deltaPSI.R` | [figure_03d_isoform_expression_vs_delta_psi.R](../scripts/figure_03d_isoform_expression_vs_delta_psi.R) |
| `Figure4A_SuppFigure8_DAS_overlap.R` | [figure_04a_supp_figure_08_das_overlap.R](../scripts/figure_04a_supp_figure_08_das_overlap.R) |
| `Figure4B_cross_layer_deltaPSI.R` | [figure_04b_cross_layer_delta_psi.R](../scripts/figure_04b_cross_layer_delta_psi.R) |
| `Figure4C_NEASE_feature_mapping.R` | [figure_04c_nease_feature_mapping.R](../scripts/figure_04c_nease_feature_mapping.R) |
| `Figure4D_Reactome_enrichment.R` | [figure_04d_reactome_enrichment.R](../scripts/figure_04d_reactome_enrichment.R) |
| `Setup_check_R_packages.R` | [setup_check_r_packages.R](../scripts/setup_check_r_packages.R) |
| `Setup_import_supplementary_data.R` | [setup_import_supplementary_data.R](../scripts/setup_import_supplementary_data.R) |
| `Setup_prepare_Figure1_Figure3_data.R` | [setup_prepare_manuscript_data.R](../scripts/setup_prepare_manuscript_data.R) |
| `SuppFigure10A_NEASE_feature_mapping.R` | [supp_figure_10a_nease_feature_mapping.R](../scripts/supp_figure_10a_nease_feature_mapping.R) |
| `SuppFigure10B_Reactome_enrichment.R` | [supp_figure_10b_reactome_enrichment.R](../scripts/supp_figure_10b_reactome_enrichment.R) |
| `SuppFigure1B_QC_read_pairs.py` | [supp_figure_01b_qc_read_pairs.py](../scripts/supp_figure_01b_qc_read_pairs.py) |
| `SuppFigure2B_KC6_event_overlap.R` | [supp_figure_02b_kc6_event_overlap.R](../scripts/supp_figure_02b_kc6_event_overlap.R) |
| `SuppFigure3AB_RBP_expression.R` | [supp_figure_03ab_rbp_expression.R](../scripts/supp_figure_03ab_rbp_expression.R) |
| `SuppFigure3AB_RBP_expression_heatmaps.R` | [supp_figure_03ab_rbp_expression_heatmaps.R](../scripts/supp_figure_03ab_rbp_expression_heatmaps.R) |
| `SuppFigure6_published_event_comparison.R` | [supp_figure_06_published_event_comparison.R](../scripts/supp_figure_06_published_event_comparison.R) |
| `SuppFigure7BC_PSI_PCA_heatmap.R` | [supp_figure_07bc_psi_pca_heatmap.R](../scripts/supp_figure_07bc_psi_pca_heatmap.R) |
| `SuppFigure7D_DAS_gene_UpSet.R` | [supp_figure_07d_das_gene_upset.R](../scripts/supp_figure_07d_das_gene_upset.R) |
| `SuppFigure9AB_SE_junction_support_plots.R` | [supp_figure_09ab_se_junction_support_plots.R](../scripts/supp_figure_09ab_se_junction_support_plots.R) |
| `SuppFigure9_SE_junction_support.R` | [supp_figure_09_se_junction_support.R](../scripts/supp_figure_09_se_junction_support.R) |
| `SuppTable2_STAR_junction_QC.py` | [supp_table_02_star_junction_qc.py](../scripts/supp_table_02_star_junction_qc.py) |
| `Utils_DAS_gene_UpSet.R` | [utils_das_gene_upset.R](../scripts/utils_das_gene_upset.R) |
| `Utils_NEASE_feature_summary.R` | [utils_nease_feature_summary.R](../scripts/utils_nease_feature_summary.R) |
| `Utils_Reactome_enrichment.R` | [utils_reactome_enrichment.R](../scripts/utils_reactome_enrichment.R) |
| `Utils_Venn_plots.R` | [utils_venn_plots.R](../scripts/utils_venn_plots.R) |
| `Utils_plotting.R` | [utils_plotting.R](../scripts/utils_plotting.R) |
| `Utils_splicing_events.R` | [utils_splicing_events.R](../scripts/utils_splicing_events.R) |
| `leave_one_out.py` | [supp_data_06_leave_one_out.py](../scripts/supp_data_06_leave_one_out.py) |
| `map_events.py` | [figure_03d_event_transcript_mapping.py](../scripts/figure_03d_event_transcript_mapping.py) |
| `nease_analysis.py` | [analysis_nease.py](../scripts/analysis_nease.py) |
| `prepare_inputs.py` | [setup_prepare_inputs.py](../scripts/setup_prepare_inputs.py) |
| `preprocess_reads.py` | [workflow_preprocess_reads.py](../scripts/workflow_preprocess_reads.py) |
| `reactome_top10.py` | [analysis_reactome_top10.py](../scripts/analysis_reactome_top10.py) |
| `run_manuscript.py` | [workflow_run_manuscript.py](../scripts/workflow_run_manuscript.py) |
| `run_rmaps3.py` | [workflow_run_rmaps3.py](../scripts/workflow_run_rmaps3.py) |
| `stromal_sensitivity.py` | [supp_data_06_stromal_sensitivity.py](../scripts/supp_data_06_stromal_sensitivity.py) |
| `tests/no_graphics.R` | [tests/utils_no_graphics.R](../tests/utils_no_graphics.R) |

</details>
