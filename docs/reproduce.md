# Reproduction guide

[Overview](../README.md) · [Figure-to-code index](figures.md) · [Input provenance](../inputs/README.md) · [Validation](validation.md)

Run the commands below from the repository root after cloning.

Twenty-eight compressed processed input tables are included under `inputs/`: fifteen complete rMATS JC tables covering the primary stromal, KC6-inclusive stromal and epithelial comparisons, stromal gene/transcript DESeq2 tables, four rMAPS3 motif tables, two candidate-RBP TPM matrices, two NEASE event-mapping tables, two corrected complete Reactome pathway tables, and a recovered 67-event published-coordinate table for Supplementary Figure 6. Their sources and scope are described in [input provenance](../inputs/README.md), and their uncompressed SHA-256 checksums are recorded in [inputs/manifest.json](../inputs/manifest.json). The reference annotation is downloaded separately from Ensembl, with its checksum verified before event–transcript mapping.

Use Python 3.9 or later and R with these packages:

```r
install.packages(c("ggplot2", "ggforce", "ragg", "gplots", "RColorBrewer",
                   "dplyr", "tidyr", "purrr", "patchwork", "VennDiagram", "scales"))
```

From a clone of this repository:

```sh
# Use a separate directory for inputs expanded from the release and analysis outputs.
python3 scripts/setup_prepare_inputs.py --work-dir work --download-reference

# Validate numeric tables and plot layers without creating images.
python3 scripts/workflow_run_manuscript.py --work-dir work --tables-only
python3 scripts/supp_data_06_leave_one_out.py --work-dir work

# Optional figure export: creates Figure 1B-D and Figure 3A-D (PDF, PNG and TIFF).
python3 scripts/workflow_run_manuscript.py --work-dir work
```

To use an existing reference instead of downloading it, replace `--download-reference` with `--gtf /path/to/Homo_sapiens.GRCh38.114.gtf.gz`. Input preparation reuses identical files and refuses to overwrite different files. Set `FIGURE_FONT` to select an installed font; the default is `sans`. The commands on this page write analysis outputs in the selected work directory. The primary workflow and the numbered commands using `KC_WORK_DIR` accept absolute script paths. Older modules listed in [additional analyses](manuscript_methods.md) use the current working directory.

Use a separate work directory for each new run: analysis result tables may be replaced within that directory.

The manuscript workflow deliberately checks the expected study counts. It is a reproduction workflow for these inputs; adapting it to another cohort requires changing sample counts and the manuscript-specific checks.

| Output under `work/results/` | Content |
| --- | --- |
| `stroma/manuscript/` | DAS tables, PSI matrix, transcript matches and unmatched-event audit |
| `stroma/figure1/` | PCA scores/variance, clustering order and gene intersections for Fig. 1B–D |
| `stroma/figure3/` | Gene memberships, scatter source tables and Fig. 3A–D |
| `stromal_leave_one_out/` | ΔPSI after each single-sample omission; no new significance testing |
| `provenance/` | Installed R package versions and session information |

The complete transcript mapping retains 7,287 compatible event–transcript pairs. Of these, 4,047 pairs from 995 events have finite expression effects and raw P values and enter Fig. 3D. Forty-two pairs, representing 29 genes, involve DE isoforms. Both inclusion-compatible and skipping-compatible pairs retain the **original rMATS ΔPSI**; their shapes distinguish the two forms. See [analysis conventions](analysis_conventions.md) for coordinate, filtering and matching details.

### Supplementary Figures 3A,B and 7D

These panels use the bundled inputs and do not need the GTF or a new DESeq2 fit. From the repository root:

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_03ab_rbp_expression_heatmaps.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_07d_das_gene_upset.R --tables-only
```

Remove `--tables-only` to generate the figures. The heatmap script processes both cohorts, or accepts `stroma` or `epithelium` to select one panel. Outputs are written to `work/results/<cohort>/suppfigure3/` and `work/results/epithelium/suppfigure7/`. The recovered heatmap scripts select 19 stromal and 23 epithelial RBP genes. The epithelial UpSet adapts the supplied Figure 1D implementation and checks 805 significant events, 609 genes and 314 SE-only genes; the original Supplementary Figure 7D script was not located.

### Supplementary Figure 6

The supplied comparison script matches all event coordinates against the published Sun et al. event set, using FDR < 0.05 without a ΔPSI cutoff.

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_06_published_event_comparison.R --tables-only
```

Outputs are written to `work/results/epithelium/suppfigure6/`. The reproduced overlap is 36 of 67 published events. The reference coordinates were recovered from the supplied membership tables because the original S8 workbook was absent from the archive; they have not been independently checked against that workbook. See [comparison methods and provenance](suppfigure6.md). Remove `--tables-only` to draw the supplied five-panel comparison.

### Figure 4A and Supplementary Figures 2B, 8 and 9A,B

These panel commands use the bundled JC tables, including the newly bundled KC6-inclusive comparison:

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_02b_kc6_event_overlap.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/figure_04a_supp_figure_08_das_overlap.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_09ab_se_junction_support_plots.R --tables-only
```

Remove `--tables-only` to export PNG/PDF panels. The Venn layouts are portable reconstructions calculated from event memberships; the Supplementary Figure 9 scatter preserves the recovered point specifications with a configurable font. [Methods, output locations and validation](comparison_panels.md) describe the full-coordinate matching and junction-support rules. No external QC logs or published-reference workbook are needed for these commands.

### Figure 4C and Supplementary Figure 10A

The recovered R plotting scripts read the bundled NEASE event memberships and calculate coverage, feature-class counts and pairwise overlaps. No new NEASE annotation run is required for these plots.

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/figure_04c_nease_feature_mapping.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_10a_nease_feature_mapping.R --tables-only
```

Remove `--tables-only` to generate PNG and SVG panels. Outputs are written to `work/results/stroma/figure4/` and `work/results/epithelium/suppfigure10/`. Each command writes event memberships, coverage, feature counts and overlaps even in tables-only mode. The upper bars use all significant SE events as denominators (698 stromal; 450 epithelial). The lower bars use event–feature assignments (321; 146), because an event can map to more than one class. Pairwise overlaps include events present in all three classes.

To use an output from `analysis_nease.py`, add `--events results/<cohort>/nease/results/SE_events.tsv`; relative paths are resolved within `KC_WORK_DIR`. The script checks event identities against the cohort's rMATS input and the expected manuscript counts. The stromal plot preserves the recovered panel body; its outer axis titles and numeric tick labels were separate PowerPoint annotations. The epithelial source already contains those axis labels.

### rMAPS3 execution and exon audit (Figure 1E,F; Supplementary Figure 7E,F)

The recovered native commands are available through `workflow_run_rmaps3.py`. The audit verifies exon selection without a genome scan or any figure output:

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work
python3 scripts/workflow_run_rmaps3.py audit --work-dir work

# Preview native commands without running rMAPS3:
python3 scripts/workflow_run_rmaps3.py run --work-dir work \
  --rmaps-dir /path/to/rMAPS3 --fasta-root /path/to/genomedata \
  --python /path/to/rmaps-env/bin/python
```

The audit reproduces all saved coordinates: 367/171/12,955 stromal and 163/136/15,719 epithelial inclusion/skipping/background exons. Native execution requires `--execute` and also creates the engine's figures. The [rMAPS3 guide](rmaps3.md) documents setup, source-version limits, output paths, and the native inclusive ΔPSI boundary used in the original runs.

### Figure 4D and Supplementary Figure 10B

The Reactome R scripts use the bundled corrected pathway tables and do not require a new NEASE calculation:

```sh
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/figure_04d_reactome_enrichment.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_10b_reactome_enrichment.R --tables-only
```

Each script checks BH adjustment over all 1,823 Reactome pathways, selects the top ten by NEASE score at adjusted P ≤ 0.05, and exports pathway identities, scores, P values and distinct affected-gene counts. Remove `--tables-only` to export figures. [Reactome provenance and methods](reactome.md) explains the corrected inputs, the R adaptation and the optional `--pathways` override.

## Process reads on Ubuntu

The configurable workflow covers FastQC, MultiQC, fastp, STAR, RSEM, rMATS and sashimi plotting. Edit [config/upstream.json](../config/upstream.json) and [config/samples_stroma.tsv](../config/samples_stroma.tsv) for your files and tools.

```sh
# Preview all processing commands; no tools are run and no files are created.
python3 scripts/workflow_preprocess_reads.py all --work-dir /path/to/analysis

# Run one stage in an environment containing the required software.
python3 scripts/workflow_preprocess_reads.py trim --work-dir /path/to/analysis --execute
```

**Study sample KC6 corresponds to the original read prefix `KC_1`, whereas study sample KC3 corresponds to `KC_6`.** The supplied sample sheet records this mapping and includes all ten specimens for preprocessing. Only the nine primary specimens enter the primary rMATS comparison. The second rMATS comparison includes KC6. Read the [Ubuntu workflow guide](upstream.md) before execution.

For Supplementary Figure 1B, [supp_figure_01b_qc_read_pairs.py](../scripts/supp_figure_01b_qc_read_pairs.py) reads paired-end fastp JSON reports and STAR logs; add `--tables-only` to export counts without images. For the STAR unique-mapping and splice-junction columns of Supplementary Table 2, use [supp_table_02_star_junction_qc.py](../scripts/supp_table_02_star_junction_qc.py) with the original `Log.final.out` and `SJ.out.tab` files. The [QC guide](qc.md) describes both QC scripts, study/source filename mapping, required reports and validation limits. The supplied Figure 1B script identifies itself as a reconstruction of the plotting workflow.

## Other numeric analyses and dependencies

The [figure-to-code index](figures.md) lists every manuscript component covered by this repository. The [manuscript-methods guide](manuscript_methods.md) covers DAS summaries, the Figure 4B comparison, gene/isoform model fitting, RBP expression statistics and NEASE calculation. Every analysis entry point accepts `--work-dir`; follow the documented order and input requirements.

The bundled-table workflow uses Python 3.9+ and R. NEASE uses a separate Python 3.12 environment with [recorded requirements](../requirements-nease.txt). The native rMAPS3 engine, Ubuntu read-processing tools and their references require separate setup; they are not installed by `setup_prepare_inputs.py`.

## Check the implementation

```sh
python3 -m unittest discover -s tests -v
Rscript --vanilla tests/test_comparison_panels.R
Rscript --vanilla tests/test_reactome_panels.R
```

These checks cover computational rules and plot construction; they are not a replacement for rerunning all native tools or visually comparing rendered figures. The [validation record](validation.md) reports the study comparisons, software versions and remaining limits. Its historical input counts describe the release stage being checked; the current manifest has 28 entries.

The optional [graphics-device guard](../tests/utils_no_graphics.R) can be supplied through `R_PROFILE_USER` for direct numeric-only R checks. That startup guard is ignored by `--vanilla`; follow the comments in the guard file when using it. The ordinary examples above rely on each script's explicit `--tables-only` mode.
