# Documentation

[Repository overview](../README.md) · [Figure-to-code index](figures.md) · [Reproduction guide](reproduce.md)

## Start here

| Guide | Use it to |
| --- | --- |
| [Script naming convention](script_naming.md) | Understand filename prefixes and find renamed scripts |
| [Figure-to-code index](figures.md) | Locate scripts for a main figure, supplementary figure or source table |
| [Reproduction guide](reproduce.md) | Set up the bundled inputs and run the numbered workflows |
| [Input provenance](../inputs/README.md) | Identify the source and scope of each included dataset |
| [Analysis conventions](analysis_conventions.md) | Check thresholds, coordinate matching, expression integration and effect direction |
| [Validation record](validation.md) | Review reproduced values, software versions and limits of the checks |

## Analysis-specific guides

| Analysis | Guide |
| --- | --- |
| FASTQ processing, alignment, quantification, rMATS and sashimi plots | [Ubuntu workflow](upstream.md) |
| Supplementary Figure 1 and Table 2 | [QC definitions and commands](qc.md) |
| RBP motif enrichment and exon-selection audit | [rMAPS3](rmaps3.md) |
| Figure 4D and Supplementary Figure 10B | [Reactome R panels](reactome.md) |
| Figure 4A and Supplementary Figures 2B, 8 and 9 | [Event-comparison panels](comparison_panels.md) |
| Supplementary Figure 6 | [Published-event comparison](suppfigure6.md) |
| DESeq2, RBP expression, NEASE and additional PCA workflows | [Additional analyses](additional_analyses.md) |

Commands that use `--tables-only` export numerical results without figure files. Native tools and older modules have their own interfaces; consult the relevant guide before running them. The repository's scope is RNA-seq and alternative-splicing analysis; experimental RT-PCR statistics and manual figure assembly are separate.

See [code organization and portable paths](code_organization.md) for the shared work-directory interface and helper modules.
