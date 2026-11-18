# Supplementary Figure 6: comparison with published events

[SuppFigure6_published_event_comparison.R](../scripts/SuppFigure6_published_event_comparison.R) adapts `supplementary6/scripts/supplementary6.R` from the authors' supplied archive. It compares the epithelial reanalysis with the 67 corneal rMATS events attributed to Supplementary Table S8 of:

Sun X, Zhang H, Shan M, Dong Y, Zhang L, Chen L and Wang Y (2022). *Comprehensive Transcriptome Analysis of Patients With Keratoconus Highlights the Regulation of Immune Responses and Inflammatory Processes*. Frontiers in Genetics 13:782709. [doi:10.3389/fgene.2022.782709](https://doi.org/10.3389/fgene.2022.782709).

## Run

From the repository root, using base R and Python 3.9 or later:

```sh
python3 scripts/prepare_inputs.py --work-dir work
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/SuppFigure6_published_event_comparison.R --tables-only
```

Outputs are written to `work/results/epithelium/suppfigure6/`: `venn_counts.tsv`, five event-membership tables, and `published_event_audit.tsv`. The comparison needs no GTF download or new rMATS run. Remove `--tables-only` to generate the supplied composite PNG/PDF and five individual PDFs; PNG generation also requires `ragg`. Circle areas are not proportional to counts. Font selection uses `FIGURE_FONT` (default `sans`).

## Matching rules

- Select reanalysis events with **FDR < 0.05**, without an absolute ΔPSI cutoff. This is different from the FDR plus |ΔPSI| > 0.10 filter used for the main DAS analyses.
- Match within each event class by chromosome, strand and every rMATS event boundary, including the flanking exons. The script lists the required coordinate columns explicitly for SE, A3SS, A5SS, MXE and RI.
- Strip an initial `chr` from chromosome names. Keep the supplied coordinate convention; do not add/subtract one, infer a genome build or apply liftOver. This reproduces the manuscript's direct, as-tabulated comparison and does not establish genome-build equivalence.
- Collapse repeated reanalysis coordinate keys within each class. Reject duplicate published coordinates, missing boundaries and invalid FDR values. Gene symbols and event IDs are not matching keys.

| Event | Reanalysis, FDR < 0.05 | Published | Shared | Reanalysis only | Published only |
| --- | ---: | ---: | ---: | ---: | ---: |
| SE | 1,214 | 33 | 16 | 1,198 | 17 |
| A3SS | 274 | 3 | 1 | 273 | 2 |
| A5SS | 227 | 6 | 5 | 222 | 1 |
| MXE | 154 | 11 | 4 | 150 | 7 |
| RI | 746 | 14 | 10 | 736 | 4 |
| Total | 2,615 | 67 | 36 | 2,579 | 31 |

The additional audit searches all 118,057 reanalysis events. Of the 31 published events outside the shared significant set, 16 have an exact coordinate match with FDR ≥ 0.05 and 15 have no exact match. “No exact match” describes this coordinate comparison, not biological absence.

## Reference provenance and limits

The supplied archive includes code and saved membership tables, but omits `Sun2022_Table1.xlsx` and the shared `common.R` file referenced by that code. The original code identifies sheet `S8`, columns M:AG, with data rows 5–37 (SE), 40–42 (A3SS), 45–50 (A5SS), 53–63 (MXE) and 66–79 (RI).

For a self-contained release, `inputs/reference/Sun2022_published_event_coordinates.tsv.gz` contains the 67 reference coordinate keys recovered from rows labelled `Both` or `Sun published only` in the supplied `*_event_membership.tsv` files. Keys are split into the named rMATS coordinate columns in the exact order used by the supplied code. Unused class-specific columns are blank. No FDR values, gene annotations or genome builds are inferred for the published events. This file is a recovered coordinate subset, **not a copy independently checked against the original S8 workbook**. The full epithelial reanalysis inputs were already bundled and are independently filtered and matched on every run.

Source archive: `supplementary6.zip`, SHA-256 `125cf47f2b0eb97ff9e38a9a9dc2eb3d656bff0acc599f4f21f736d744300547`. The manifest records the recovered table's uncompressed checksum. The original Sun et al. publication is distributed under [CC BY](https://creativecommons.org/licenses/by/4.0/); credit that study when using its event coordinates. The repository's MIT license applies to the analysis code and does not replace attribution to the published data source.

The release check reproduced all six supplied comparison tables byte for byte and the 36/16/15 audit counts. All 50 drawing calls matched the supplied R functions using non-rendering recording stubs. No figure image was regenerated or visually compared during this check.

Independent verification against the original reference workbook is outside the agreed scope of this release. This does not change the documented provenance of the recovered coordinate subset.
