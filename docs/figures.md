# Figure-to-code index

[Repository overview](../README.md) · [Reproduction guide](reproduce.md) · [Input provenance](../inputs/README.md) · [Validation](validation.md)

Computational sources for the **18 November 2026 manuscript** are listed below. Source-table scripts provide values used in a panel; they do not necessarily export the final assembled layout. Figure filenames follow the manuscript numbering.

[Figure 1](#figure-1) · [Figure 2](#figure-2) · [Figure 3](#figure-3) · [Figure 4](#figure-4) · [Supplementary material](#supplementary-figures-and-tables) · [Supplementary data](#supplementary-data)

## Figure 1

| Manuscript component | Scripts |
| --- | --- |
| Fig. 1A; DAS calls and class counts | [Figure1A_DAS_event_summary.R](../scripts/Figure1A_DAS_event_summary.R) |
| Fig. 1B,C; PSI PCA and heatmap | [Figure1BC_PSI_PCA_heatmap.R](../scripts/Figure1BC_PSI_PCA_heatmap.R) |
| Fig. 1D; gene-level UpSet | [Figure1D_DAS_gene_UpSet.R](../scripts/Figure1D_DAS_gene_UpSet.R) |
| Fig. 1E,F; RBP motif heatmaps | [Figure1EF_RBP_motif_heatmaps.R](../scripts/Figure1EF_RBP_motif_heatmaps.R) |
| Fig. 1E,F and Supplementary Fig. 7E,F; native rMAPS3 execution and exon-selection audit | [run_rmaps3.py](../scripts/run_rmaps3.py), [settings](../config/rmaps3.json) |
| Fig. 1E,F; motif statistics and sensitivity tables | [Figure1EF_RBP_motif_statistics.R](../scripts/Figure1EF_RBP_motif_statistics.R) |

Figure 1B–D is included in `run_manuscript.py`. Motif modules are separate; see the [rMAPS3 guide](rmaps3.md), including the native ΔPSI boundary and external-engine requirements.

## Figure 2

| Manuscript component | Scripts |
| --- | --- |
| Fig. 2A–C; sashimi plots | [preprocess_reads.py](../scripts/preprocess_reads.py), `sashimi` stage; event selection is supplied by [Figure2DEF_selected_exon_PSI.R](../scripts/Figure2DEF_selected_exon_PSI.R) |
| Fig. 2D–F; RNA-seq PSI graphs prepared in GraphPad Prism | [Figure2DEF_selected_exon_PSI.R](../scripts/Figure2DEF_selected_exon_PSI.R) exports sample-level PSI source tables; final graphing uses Prism |
| Fig. 2G–I; RT-PCR gel panels | Experimental images assembled separately |
| Fig. 2J–L; RT-PCR inclusion-fraction graphs prepared in GraphPad Prism | GraphPad Prism workflow; no custom plotting script required |

Sashimi plotting requires STAR-aligned BAM files. The selected-event script supplies the three events and their RNA-seq PSI values. RT-PCR measurements and statistics are outside this repository’s scope.

## Figure 3

| Manuscript component | Scripts |
| --- | --- |
| Fig. 3A; DEG/DAS gene Venn | [Figure3A_DEG_DAS_Venn.R](../scripts/Figure3A_DEG_DAS_Venn.R) |
| Fig. 3B; gene expression versus ΔPSI | [Figure3B_gene_expression_vs_deltaPSI.R](../scripts/Figure3B_gene_expression_vs_deltaPSI.R) |
| Fig. 3C; isoform/DAS gene Venn | [Figure3C_isoform_DAS_Venn.R](../scripts/Figure3C_isoform_DAS_Venn.R) |
| Fig. 3D; matched isoform expression versus ΔPSI | [Figure3D_isoform_expression_vs_deltaPSI.R](../scripts/Figure3D_isoform_expression_vs_deltaPSI.R), [map_events.py](../scripts/map_events.py) |
| Fig. 3A–D; gene/transcript model fitting from RSEM | [Figure3ABCD_DESeq2_analysis.R](../scripts/Figure3ABCD_DESeq2_analysis.R) |
| Fig. 3A–C; expression/splicing overlap source tables | [Figure3ABC_expression_splicing_overlap.R](../scripts/Figure3ABC_expression_splicing_overlap.R) |

`run_manuscript.py` reproduces Figure 3A–D from the bundled DESeq2 tables and the release-114 GTF. Refitting expression models requires RSEM quantifications and sample annotations; see [gene and isoform expression](additional_analyses.md#gene-and-isoform-expression).

## Figure 4

| Manuscript component | Scripts |
| --- | --- |
| Fig. 4A and Supplementary Fig. 8; coordinate-matched overlap tables and four-set Venn plots | [Figure4A_SuppFigure8_DAS_overlap.R](../scripts/Figure4A_SuppFigure8_DAS_overlap.R) |
| Fig. 4B; stromal versus epithelial ΔPSI | [Figure4B_cross_layer_deltaPSI.R](../scripts/Figure4B_cross_layer_deltaPSI.R) |
| Fig. 4C; stromal protein-feature coverage and overlaps | [Figure4C_NEASE_feature_mapping.R](../scripts/Figure4C_NEASE_feature_mapping.R) |
| Fig. 4D; stromal Reactome top-ten R panel | [Figure4D_Reactome_enrichment.R](../scripts/Figure4D_Reactome_enrichment.R) |
| Fig. 4C,D and Supplementary Fig. 10; upstream NEASE calculation and revised Python plot implementation | [nease_analysis.py](../scripts/nease_analysis.py), [reactome_top10.py](../scripts/reactome_top10.py) |
| Fig. 4E; conceptual summary | Separately assembled artwork |

Figure 4C,D R scripts use bundled corrected NEASE outputs. New NEASE calculations require the [external environment](additional_analyses.md#nease-and-reactome). The retained [reactome_top10.py](../scripts/reactome_top10.py) is the revised Python plotting implementation; see [Reactome provenance](reactome.md) for the R adaptation.

## Supplementary figures and tables

| Manuscript component | Scripts |
| --- | --- |
| Supplementary Fig. 1A; QC report generation | [preprocess_reads.py](../scripts/preprocess_reads.py): `raw-qc`, `trim`, `trimmed-qc`, `align`; see [QC coverage](qc.md) for the final panels |
| Supplementary Fig. 1B; raw, filtered and mapped read pairs | [SuppFigure1B_QC_read_pairs.py](../scripts/SuppFigure1B_QC_read_pairs.py); supplied reconstruction, with a numeric-only mode |
| Supplementary Fig. 2 and Data 6; stromal sensitivity tables | [stromal_sensitivity.py](../scripts/stromal_sensitivity.py), [leave_one_out.py](../scripts/leave_one_out.py) |
| Supplementary Fig. 2B; KC6-inclusive versus primary event Venn | [SuppFigure2B_KC6_event_overlap.R](../scripts/SuppFigure2B_KC6_event_overlap.R); also exports class counts for Supplementary Fig. 2A |
| Supplementary Fig. 3A,B; candidate RBP expression heatmaps | [SuppFigure3AB_RBP_expression_heatmaps.R](../scripts/SuppFigure3AB_RBP_expression_heatmaps.R) |
| Supplementary Fig. 3A,B; expression tables with DESeq2 adjusted P values | [SuppFigure3AB_RBP_expression.R](../scripts/SuppFigure3AB_RBP_expression.R) |
| Supplementary Figs. 4 and 5; RT-PCR gels and controls | Experimental panels assembled separately |
| Supplementary Fig. 6; exact-coordinate comparison with published Sun et al. events | [SuppFigure6_published_event_comparison.R](../scripts/SuppFigure6_published_event_comparison.R) |
| Supplementary Fig. 7A; epithelial DAS calls and class counts | [Figure1A_DAS_event_summary.R](../scripts/Figure1A_DAS_event_summary.R), `epithelium` argument |
| Supplementary Fig. 7B,C; epithelial PSI PCA and heatmap | [SuppFigure7BC_PSI_PCA_heatmap.R](../scripts/SuppFigure7BC_PSI_PCA_heatmap.R) |
| Supplementary Fig. 7D; epithelial gene-level UpSet | [SuppFigure7D_DAS_gene_UpSet.R](../scripts/SuppFigure7D_DAS_gene_UpSet.R) |
| Supplementary Fig. 7E,F; epithelial motif analysis | [run_rmaps3.py](../scripts/run_rmaps3.py), [motif statistics](../scripts/Figure1EF_RBP_motif_statistics.R), [motif heatmaps](../scripts/Figure1EF_RBP_motif_heatmaps.R) |
| Supplementary Fig. 8A–D; A3SS/A5SS/MXE/RI cross-layer overlaps | [Figure4A_SuppFigure8_DAS_overlap.R](../scripts/Figure4A_SuppFigure8_DAS_overlap.R) |
| Supplementary Fig. 9; additional SE junction-support sensitivity tables | [SuppFigure9_SE_junction_support.R](../scripts/SuppFigure9_SE_junction_support.R) |
| Supplementary Fig. 9A,B; support-filtered SE Venn and scatter | [SuppFigure9AB_SE_junction_support_plots.R](../scripts/SuppFigure9AB_SE_junction_support_plots.R) |
| Supplementary Fig. 10A; epithelial protein-feature coverage and overlaps | [SuppFigure10A_NEASE_feature_mapping.R](../scripts/SuppFigure10A_NEASE_feature_mapping.R) |
| Supplementary Fig. 10B; epithelial Reactome top-ten R panel | [SuppFigure10B_Reactome_enrichment.R](../scripts/SuppFigure10B_Reactome_enrichment.R) |
| Supplementary Table 1; RT-PCR primers | Experimental metadata |
| Supplementary Table 2; STAR unique-mapping and junction-support columns | [SuppTable2_STAR_junction_QC.py](../scripts/SuppTable2_STAR_junction_QC.py), [sample metadata](../config/samples_stroma_qc.tsv) |

Figure S1B is the supplied reconstruction; Figure S7D adapts the Figure 1D implementation. The Venn layouts are portable reconstructions. Numerical checks and provenance are described in [validation](validation.md), [comparison panels](comparison_panels.md) and [published-event comparison](suppfigure6.md); visual identity is not asserted. Supplementary Figure 2A uses the class-count output of `SuppFigure2B_KC6_event_overlap.R` with separately assembled exon schematics. The [QC guide](qc.md#coverage) identifies which Supplementary Table 2 columns are calculated by the supplied scripts.

## Supplementary data

| Data | Analysis source |
| --- | --- |
| 1, 2 | Complete stromal/epithelial rMATS JC outputs; [DAS summary](../scripts/Figure1A_DAS_event_summary.R) and [upstream commands](upstream.md) |
| 3 | [Motif multiplicity and sensitivity](../scripts/Figure1EF_RBP_motif_statistics.R) |
| 4 | [DESeq2](../scripts/Figure3ABCD_DESeq2_analysis.R), [expression/splicing overlap](../scripts/Figure3ABC_expression_splicing_overlap.R), [event–transcript matching](../scripts/map_events.py) |
| 5 | [NEASE](../scripts/nease_analysis.py), [feature summaries](../scripts/Utils_NEASE_feature_summary.R), [Reactome selection](../scripts/Utils_Reactome_enrichment.R) |
| 6 | [Stromal sensitivity](../scripts/stromal_sensitivity.py), [leave-one-out effects](../scripts/leave_one_out.py) |

These modules export analysis tables underlying the supplementary material; they do not assemble the final formatted Excel workbooks.

## Shared preparation and optional analyses

The primary workflow uses [run_manuscript.py](../scripts/run_manuscript.py), [prepare_inputs.py](../scripts/prepare_inputs.py), [Setup_check_R_packages.R](../scripts/Setup_check_R_packages.R), [Setup_prepare_Figure1_Figure3_data.R](../scripts/Setup_prepare_Figure1_Figure3_data.R), and [Utils_plotting.R](../scripts/Utils_plotting.R). [Setup_import_supplementary_data.R](../scripts/Setup_import_supplementary_data.R) imports supplementary workbooks for additional modules. Shared counting and plot helpers remain in `scripts/`.

Optional PCA analyses are [Analysis_PCA_all_events.R](../scripts/Analysis_PCA_all_events.R), [Analysis_PCA_variable_events.R](../scripts/Analysis_PCA_variable_events.R), and [Analysis_PCA_including_KC6.R](../scripts/Analysis_PCA_including_KC6.R). They are not assigned to numbered panels in this manuscript; their dependencies and execution order are in [additional analyses](additional_analyses.md#optional-plots-and-pca-sensitivity).
