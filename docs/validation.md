# Release validation

Checked on 18 November 2026 against the supplied analysis archive and manuscript revision. Only numeric tables and plot specifications were regenerated; no PNG, PDF or TIFF figure files were generated.

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

`Setup_check_R_packages.R` records the versions installed in each reproduction environment. These are the release-check versions, not a claim that all original analyses used the same environment.

The earlier downstream release check reproduced the 414 gene and 448 isoform DE call sets using DESeq2 1.52.0 and tximport 1.40.0. The manuscript reports tximport 1.34.0; isoform fold changes differed by at most 0.000028 in that check. The current bundled-table run uses the supplied DESeq2 tables directly and does not refit those models.

## Recovered Supplementary Figure 3A,B and reconstructed Figure 7D

The additional check on 18 November 2026 ran both numbered supplementary scripts with `--tables-only` from a separate work directory, using absolute script paths and `KC_WORK_DIR`. All 18 bundled inputs passed their checksum checks. No figure images were created.

- Supplementary Figure 3A: 19 genes × nine samples; Figure 3B: 23 genes × 14 samples. Recomputed motif-selection metadata, TPM, log2(TPM + 1) and row Z-scores matched the saved outputs of the recovered scripts, with no numeric differences in the eight compared tables. Row clustering retains Euclidean distance and Ward.D2; sample columns retain the original control-then-case order.
- Supplementary Figure 7D: 805 events in 609 genes, including 314 SE-only genes. An independent Python calculation from the complete epithelial JC tables matched all 609 gene memberships and all 21 nonempty intersections. The shared plotting code built its ggplot layers without graphics devices. This is an adaptation of the supplied stromal UpSet code; it does not establish pixel-level equivalence to the original epithelial panel.
- Figure 1D regression: the previous published script and the shared implementation were run on the same stromal data. All 14 output tables were byte-for-byte identical, retaining 1,073 events, 807 genes and 487 SE-only genes.
- All R scripts passed parsing, and the existing eight automated analysis/workflow tests passed.

## Limits of this validation

- FastQC, MultiQC, fastp, STAR, RSEM, rMATS and sashimi commands were planned and regression-tested, but no raw reads were processed on Ubuntu for this release.
- Figure images were not regenerated or visually compared. Platform fonts and plotting versions can affect layout.
- Supplementary Figures 3A,B and 7D have separate commands using bundled inputs; they are not part of `run_manuscript.py`. Other epithelial, rMAPS3 and NEASE modules may require the separately described supplementary workbooks, logs or quantifications.
- The original PowerPoint composition, RT-PCR gel-intensity measurements and conceptual schematic are outside the checked workflow.
