# Release validation

[Repository overview](../README.md) · [Documentation](README.md) · [Figure-to-code index](figures.md)

Checked on 18 November 2026 against the supplied analysis archive and manuscript revision. The checks below distinguish numeric reproduction, plot-object construction and rendered-figure validation. Original manuscript and presentation files were not modified.

The current repository is restricted to reported manuscript results and their required support code. Earlier validation sections document the code present at the time of each check; they do not imply that exploratory modules remain in the current release.

## Reproduced results

| Check | Result |
| --- | --- |
| Primary stromal DAS | 1,073 events in 807 genes |
| SE / RI / A5SS / A3SS / MXE | 698 / 152 / 88 / 82 / 53 events |
| Figure 1 PCA | PC1 46.316%, PC2 16.153% (rounded) |
| Gene-level expression overlap | 414 DEGs; 8 genes shared with DAS |
| Figure 3B | 805 displayed DAS genes: 797 DAS-only and 8 DAS + DEG |
| Isoform-level gene sets | 448 DE isoforms in 410 genes; 55 genes shared with DAS |
| GTF coverage | All 110,618 tested isoforms present in release 114 |
| All structurally compatible event–transcript pairs | 7,287 |
| Figure 3D pairs | 4,047 from 995 events and 765 genes |
| Inclusion / skipping pairs | 1,954 / 2,093 |
| DE pairs / DE transcripts / DE genes | 42 / 31 / 29 |
| Event mapping audit | 995 mapped/tested; 13 mapped without tested isoform effects; 65 without a compatible annotated transcript |
| Figure 3C shared genes | 29 event-matched; 26 other gene-level overlaps |
| Single-sample omission | 1,060 of 1,073 retain direction across all nine omissions |

Twenty-four regenerated source tables were compared with the supplied archive, including the event–transcript mapping, plotted pairs, Figure 1 PCA scores/variance, UpSet intersections, Figure 3 membership and scatter tables. Identifiers, memberships, assignments, plotting directions and counts agreed. Floating-point values were compared with relative tolerance 1e-12 and absolute tolerance 1e-14; a few serialized FDR values differed at approximately 1e-16.

All 1,073 × 9 single-sample omission effects were also compared with the revised Supplementary Data 6. The maximum absolute difference was 1.39e-16. All case omissions retained all 1,073 directions; the four control omissions retained 1,073, 1,072, 1,070 and 1,064, respectively.

The R scripts passed syntax checks; Python scripts passed compilation checks. Regression tests cover local splice junction requirements, MXE strand orientation, R-normalized coordinate headers, retained intron structure, alternative splice boundaries, gene/chromosome/strand matching, study-label exclusion, safe paths and grouping, and execution/resume/failure behavior using mock commands. Plot layers were built without graphics devices, including the Figure 3D ΔPSI values, shapes and colors.

## Software used for these checks

| Component | Version |
| --- | --- |
| Python | 3.9.6 |
| R | 4.6.1 |
| ggplot2 | 4.0.3 |
| ggforce | 0.5.0 |
| ragg | 1.5.2 |
| gplots | 3.3.0 |
| RColorBrewer | 1.1.3 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| purrr | 1.2.2 |
| patchwork | 1.3.2 |

`setup_check_r_packages.R` records the versions installed in each reproduction environment. These are the release-check versions, not a claim that all original analyses used the same environment.

The earlier downstream release check reproduced the 414 gene and 448 isoform DE call sets using DESeq2 1.52.0 and tximport 1.40.0. The manuscript reports tximport 1.34.0; isoform fold changes differed by at most 0.000028 in that check. The current bundled-table run uses the supplied DESeq2 tables directly and does not refit those models.

## Recovered Supplementary Figure 3A,B and reconstructed Figure 7D

The additional check on 18 November 2026 ran both numbered supplementary scripts with `--tables-only` from a separate work directory, using absolute script paths and `KC_WORK_DIR`. All 18 bundled inputs passed their checksum checks. No figure images were created.

- Supplementary Figure 3A: 19 genes × nine samples; Figure 3B: 23 genes × 14 samples. Recomputed motif-selection metadata, TPM, log2(TPM + 1) and row Z-scores matched the saved outputs of the recovered scripts, with no numeric differences in the eight compared tables. Row clustering retains Euclidean distance and Ward.D2; sample columns retain the original control-then-case order.
- Supplementary Figure 7D: 805 events in 609 genes, including 314 SE-only genes. An independent Python calculation from the complete epithelial JC tables matched all 609 gene memberships and all 21 nonempty intersections. The shared plotting code built its ggplot layers without graphics devices. This is an adaptation of the supplied stromal UpSet code; it does not establish pixel-level equivalence to the original epithelial panel.
- Figure 1D regression: the previous published script and the shared implementation were run on the same stromal data. All 14 output tables were byte-for-byte identical, retaining 1,073 events, 807 genes and 487 SE-only genes.
- All R scripts passed parsing, and the existing eight automated analysis/workflow tests passed.

## Recovered Figure 4C and Supplementary Figure 10A

The R sources were recovered from `outputs/nease_overlap_preview_20261012`: `plot_overlap_extra_large_text.R` for Figure 4C and `plot_epithelium_with_axis_labels.R` for Supplementary Figure 10A. Their existing PNG outputs match the decoded RGB pixels of the embedded figures in slides 5 and 15 of the supplied 18 November 2026 deck. This identifies the actual source versions, including the feature-overlap labels.

The public scripts ran with `--tables-only` using the 20 checksum-verified bundled inputs. They calculate plot counts from the saved NEASE event memberships and check event identity, coordinates, genes, FDR and ΔPSI against the significant rMATS SE calls. Independently deduplicating the original NEASE domain, residue-evidence and motif annotation tables confirmed every exported event membership in both cohorts.

| Check | Stroma / Figure 4C | Epithelium / Supplementary Figure 10A |
| --- | --- | --- |
| Mapped / all SE events | 276 / 698 (39.5%) | 135 / 450 (30.0%) |
| Unmapped events | 422 | 315 |
| Domain / residue / motif assignments | 264 / 50 / 7 | 129 / 15 / 2 |
| Lower-bar denominator | 321 assignments | 146 assignments |
| Inclusive domain–residue / domain–motif / residue–motif overlaps | 41 / 4 / 3 | 10 / 1 / 0 |
| Three-class overlap | 3 | 0 |

The 35 stromal and 45 epithelial drawing calls, including all labels and coordinates, matched the recovered sources when recorded through non-rendering graphics stubs. No graphics device was opened. Invalid inputs for the wrong cohort, duplicate events, inconsistent mapped-any flags, and the ORF-filtered analysis were rejected. The `--events` override reproduced the same four numeric tables from another working directory. All R scripts passed parsing.

No new NEASE reference mapping, pathway-enrichment run or figure image was generated. The existing-image comparison establishes provenance; it is not a visual test of newly rendered public-script output. Figure 4C's outer axis titles and numeric tick labels remain separate annotations in the original PowerPoint composition.

## Supplied STAR junction-QC script

`STAR_junction_QC.py` was incorporated as `supp_table_02_star_junction_qc.py`. The definitions and output columns are retained; sample labels/RIN are in `config/samples_stroma_qc.tsv`, with explicit filename options for original source prefixes and the repository's per-sample alignment layout.

On synthetic STAR files for all ten samples, the adapted script matched every output value and the complete TSV bytes of the supplied script. Thirty quartiles/medians from ten vectors (including singletons and fractional interpolation) matched R `quantile(type=7)`. Tests cover zero-support exclusion, separation of unique and multimapper support, unchanged pair counts, sample-ID mapping including KC6, missing/malformed inputs and protection of existing results. At this junction-only check, all 14 automated tests passed (eight existing analysis/workflow tests and six QC tests); the subsequent Figure 1B addition brings the suite to 20.

The original `Log.final.out` and `SJ.out.tab` files were not found in the inspected local locations. No study QC values were recomputed or filled in from synthetic data, and no QC images or read alignments were regenerated. The script supplies only the unique-mapping/junction columns of Supplementary Table 2; see [QC coverage](qc.md) for Supplementary Figure 1 and other table columns.

## Supplied Supplementary Figure 1B reconstruction

`Keratoconus_QC_ReadPairs.py` was adapted as `supp_figure_01b_qc_read_pairs.py`, with explicit study/source IDs and a `--tables-only` option. The supplied source says it reconstructs the plot from the reported methods and has not been verified against the original plotting code; this release retains that distinction.

On synthetic reports for all ten samples, every TSV value and the entire TSV file matched the supplied code. The plotting function's parsed Python syntax tree is unchanged. Tests check division of fastp reads by two, unchanged STAR pair counts, exclusion of too-many-loci mappings, all sample mappings, odd/invalid counts, inconsistent input totals, duplicate log metrics and output protection. Tables-only mode ran from another directory with third-party Python packages disabled. All 20 automated tests passed, including the six additional read-pair tests.

No original fastp JSON/STAR log set was supplied or found in the inspected locations. The manuscript's actual sample-level read counts were not recomputed, and no image was rendered. Unlike the supplied script's warning, mismatched fastp filtered-pair and STAR input-pair totals now stop the calculation; the valid-input arithmetic is unchanged.

## Supplied Supplementary Figure 6 comparison

All 21 bundled inputs passed their checksums. The numbered script ran with `--tables-only` from a separate working directory. The five membership tables and the count summary matched the six supplied TSV files byte for byte. Reanalysis totals (SE/A3SS/A5SS/MXE/RI) were 1,214/274/227/154/746; published totals were 33/3/6/11/14; intersections were 16/1/5/4/10. The additional audit reproduced 36 shared significant events, 16 exact matches outside the significant set (all with FDR ≥ 0.05), and 15 without an exact match.

All 50 drawing calls matched the supplied R functions through recording stubs, without opening a graphics device. The original reference workbook was absent: the 67 published coordinate records were recovered from the supplied comparison memberships, not independently verified against the workbook. See [provenance and methods](suppfigure6.md). No liftOver, new rMATS analysis or figure rendering was performed.

## Additional event-comparison panels

The Supplementary Figure 2B, Figure 4A / Supplementary Figure 8, and Supplementary Figure 9A,B commands ran from absolute paths with `KC_WORK_DIR`, using 26 checksum-verified input files. The five added inputs are complete KC6-inclusive JC tables copied from the existing analysis, not a new alignment or rMATS run.

- Supplementary Figure 2B: 1,073 primary and 1,074 inclusive events, with 184/889/185 exclusive/shared/exclusive memberships. Inclusive-cohort class counts are 713/78/85/44/154 (SE/A3SS/A5SS/MXE/RI), matching the supplied figure labels.
- Figure 4A / Supplementary Figure 8: every coordinate/direction membership in all five union tables matched the previous released implementation, and all 40 directional counts matched the saved original matrix. Shared counts are 13/1/4/0/5; SE has seven concordant and six discordant shared events.
- Supplementary Figure 9: all 15,343 selected event keys, effects, FDRs, significance flags and four-group support counts matched the saved support audit. All 15 Venn regions and all 491 significant-event memberships matched the original region audit. Missing `NA` gene symbols are omitted only from display lists; shared-gene labels are unchanged. Pearson r = 0.1011075 and Spearman rho = 0.08498896.
- All 15,343 scatter points and their colors, shapes, sizes and alpha values matched the recovered scatter implementation. The Venn code uses VennDiagram 1.8.2 with `ind=FALSE` to construct objects. Focused R tests cover complete coordinate identity, strict cutoffs, all directional memberships, numeric Venn labels and per-group support thresholds. The existing 20 Python tests also pass.

The Venn panels are portable layout reconstructions; pixel-level equivalence is not claimed. A comparison-wrapper error generated temporary Supplementary Figure 9A/B PNG/PDF files during validation. Those files are not distributed, were not used to change the manuscript figures, and were not visually validated. The comparison was subsequently repeated with graphics-device functions disabled, producing tables only. The normal `--tables-only` entry commands also completed without images.

Independent checks against original QC logs and the published S8 workbook are outside the agreed scope. Their previously documented validation limits remain in place.

## rMAPS3 execution wrapper and Reactome R scripts

All 28 bundled inputs passed their SHA-256 checks. The new wrapper's numeric-only audit reproduced every row and the complete file hashes of all six saved native exon-coordinate files: stromal up/down/background counts 367/171/12,955 and epithelial counts 163/136/15,719. The two complete SE input files and the inspected upstream 114-motif library also match the original rMAPS3 inputs byte for byte.

The native rMAPS3 selection uses inclusive ΔPSI cutoffs at ±0.10, unlike the general DAS summaries' strict boundary. Tests preserve that distinction, missing-PSI behavior, strict background boundaries, last-qualifying-row deduplication, strand identity and exclusion of cross-group conflicts. Preview mode created no files. Mock execution checked command arguments, paths with spaces, Fisher-test environment, native-coordinate verification, output protection and failure status. All 25 Python tests passed. No native motif scan or P-value recalculation was executed; the original source/environment revision remains unverified, as described in [rMAPS3 provenance](rmaps3.md).

Both numbered Reactome R entry scripts ran with `--tables-only` from another working directory. Graphics-device functions were disabled through `tests/utils_no_graphics.R` before startup; the startup guard was confirmed active. Each complete table has 1,823 pathways, with 64 stromal and 26 epithelial pathways at adjusted P ≤ 0.05. Independent BH recomputation differed from the saved adjusted values by at most 2.23e-16.

All twenty selected pathway IDs, ranks, names, gene lists and distinct-gene counts matched the saved corrected top-ten tables exactly. NEASE scores and raw/adjusted P values, including −log10 values, matched at relative tolerance 1e-12 and absolute tolerance 1e-14 (maximum absolute serialized difference below 4.85e-14). Both plot objects retained all ten finite points, including the 14-gene epithelial pathway omitted by the old R size cap. Focused R tests also checked equal-name pathways, tied scores, the inclusive 0.05 boundary, duplicate genes, invalid IDs, incomplete BH tables and invalid scores.

No figures were generated during this addition, and no Word or PowerPoint file was modified. The R scripts are adaptations of the older plotting sources using the corrected numerical results; rendered visual equivalence was not assessed. The complete NEASE enrichment calculation was not repeated.

## Script-name standardization

All 49 analysis, workflow and helper scripts now follow the [shared naming convention](script_naming.md). A comparison with the preceding versions confirmed that script contents changed only where filenames, imports or workflow entry-point names needed updating. All 33 input and configuration files retained their SHA-256 hashes.

All 41 R files and 15 Python files, including tests, passed syntax checks. The 25 Python tests and both focused R test scripts passed. The renamed `workflow_run_manuscript.py --tables-only` command completed from a separate working directory with graphics devices disabled, exercising all nine steps for Figure 1B–D and Figure 3A–D. It retained 4,047 Figure 3D pairs from 995 events and 765 genes, including 42 DE pairs. All checked local documentation links and anchors resolved. No figure files were generated, and no Word or PowerPoint file was modified.

## Code readability and portable-path refactor

All 49 existing scripts were reviewed and organized, with three shared helpers added for R/Python paths and PSI parsing. Every analysis entry point accepts `--work-dir`; the same `KC_WORK_DIR` fallback is available in both languages. The [code organization guide](code_organization.md) defines relative-path resolution and the helper responsibilities. No personal machine path is required in scripts or configuration files.

The preceding commit (`4065e23`) was retained as the comparison baseline. Its 22 numerical commands ran in a separate work directory. The revised commands ran from outside their work directory with explicit `--work-dir` arguments; the motif-heatmap table export was checked as an additional command. Graphics devices were disabled throughout these R checks.

- All 158 unchanged result files matched byte for byte. Three provenance tables differed only in the recorded work-directory paths; their other fields and input checksums matched. No preceding result file was missing.
- The shared PSI parser matched each of the three original implementations across 9,000 tested sample vectors. No full optional PCA rerun or large PCA cache regeneration was performed.
- All four newly exported motif adjusted-P matrices matched the original selection function, including motif order and region names, with numeric tolerance 1e-12.
- Both Python Reactome top-ten selections matched the original numerical statements exactly. Its new `--tables-only` command ran with matplotlib imports blocked and wrote no images. NEASE event-annotation helper results matched the original helpers on synthetic inclusive-boundary and excluded-event cases; the full native enrichment calculation was not rerun.
- All 29 Python tests and three R test scripts passed. These include work-directory precedence, paths containing spaces, absolute and home-relative paths, malformed arguments, PSI validation and all existing analysis/QC regressions. Every Python entry point supports `--help` without optional analysis packages; importing the modules created no outputs or analysis state.
- All 44 R files and 17 Python files, including tests and helpers, passed parsing. Formatting changes were checked against parsed expressions. All bundled inputs and configuration files remained unchanged.

This work generated numerical tables only. No figure image, Word file or PowerPoint file was created or modified. Native read processing, DESeq2 refitting, full NEASE enrichment and the optional large PCA workflows retain the previously documented validation limits.

## Restriction to manuscript results

The 18 November 2026 manuscript and Supplementary Data 3 and 6 were inspected without modification. The current `scripts/` directory contains 46 files, each linked from the [figure-to-code index](figures.md) to a reported result or a required preprocessing/helper role.

The three exploratory PCA variants and their exclusive PSI helper were removed, together with the duplicate Python Reactome plotting implementation and the broader junction-support audit. The significant-event PCA implementations for Figure 1B and Supplementary Figure 7B remain. The numbered Supplementary Figure 9 command already supplies the reported support-filtered memberships and correlations. NEASE now exports only the primary feature mapping and Reactome results; the unreported KEGG and ORF/coordinate-comparison branches were removed.

The window-corrected BH and BY columns are present in Supplementary Data 3 and were retained. The coverage, KC6-inclusion and single-sample-omission calculations in Supplementary Data 6 were also retained, as were the candidate-RBP expression statistics reported in the text.

All remaining R/Python scripts passed parsing, the 29 Python tests and three R test scripts passed, and local documentation links and helper references resolved. A mocked NEASE engine compared the previous and restricted entry points: all ten retained primary TSV exports were byte-identical, primary analysis settings agreed, and the retained summary fields agreed. This checks dispatch and export behavior; it is not a new native enrichment run. Bundled inputs and configuration files were unchanged, no images were generated, and the manuscript files retained their original SHA-256 fingerprints.

## Limits of this validation

- FastQC, MultiQC, fastp, STAR, RSEM, rMATS and sashimi commands were planned and regression-tested, but no raw reads were processed on Ubuntu for this release.
- Newly rendered outputs were not visually checked; see the temporary Figure 9 rendering noted above. Existing NEASE images were compared to identify source scripts. Platform fonts and plotting versions can affect layout.
- Figures 4A,C,D and Supplementary Figures 2B, 3A,B, 6, 7D, 8, 9A,B and 10A,B have separate commands using bundled inputs; they are not part of `workflow_run_manuscript.py`. Native rMAPS3 execution requires an external engine and genome; upstream NEASE and other expression modules require their separately documented dependencies or inputs.
- Figures 2D–F and 2J–L were prepared in GraphPad Prism, as confirmed by the authors. No custom plotting scripts are required for those panels. Their Prism workflows were not rerun during code validation; the RNA-seq PSI source-table export for Figures 2D–F remains included.
- The original PowerPoint composition, RT-PCR gel-intensity measurements and conceptual schematic are outside the checked workflow.
