# Reactome R panels

[Repository overview](../README.md) · [Documentation](README.md) · [Figure-to-code index](figures.md)

The numbered entry scripts are `figure_04d_reactome_enrichment.R` (stroma) and `supp_figure_10b_reactome_enrichment.R` (epithelium). Both use `utils_reactome_enrichment.R`. Their source was recovered from each cohort's `after_NEASE/scripts/エンリッチメント.R`; local paths, `setwd()`, workspace clearing and unguarded interactive plot printing were removed.

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/figure_04d_reactome_enrichment.R --tables-only
KC_WORK_DIR="$PWD/work" Rscript --vanilla scripts/supp_figure_10b_reactome_enrichment.R --tables-only
```

`--tables-only` writes numeric tables and builds the point layers without opening a graphics device. Omit it to export PDF, PNG and TIFF panels. Packages are ggplot2 and scales, plus ragg for raster export. `FIGURE_FONT` defaults to `sans`. Script paths can be absolute and the working directory can be elsewhere.

## Input and selection

The bundled `data/<cohort>/nease/Reactome_Python.tsv` files are the complete corrected exports from the saved `python_revision_20261010/results` analyses: 1,823 pathways in each cohort, with the primary analysis's ORF exclusion disabled. They match the pathway schema produced by `analysis_nease.py`. The older R sources read 1,750-row `Reactome.csv` exports, predating the correction for repeated pathway display names; those older tables are not used here.

For a separately calculated NEASE result, add `--pathways results/<cohort>/nease/results/Reactome_Python.tsv`. Relative overrides resolve within `KC_WORK_DIR`. The code verifies required columns, unique pathway IDs, finite values, Reactome-only membership, and BH correction over the complete table. It retains the supplied adjusted P values, filters at **adjusted P ≤ 0.05**, ranks by decreasing NEASE score, breaks ties by pathway ID, and selects ten pathways. There are 64 eligible stromal and 26 eligible epithelial pathways in the bundled exports.

Point size counts distinct splicing-affected genes, after removing parenthetical interaction counts; it does not sum interactions. Point color is −log10(adjusted P). Pathway IDs define the plot positions, so equal display names cannot collapse two pathways onto one row. In particular, the old epithelial size limit of seven was removed: the corrected Vesicle-mediated transport pathway has 14 genes and must remain visible.

## Outputs and scope

Tables are `Reactome_top10_for_figure.tsv`, `Reactome_summary.tsv` and `Reactome_input_provenance.tsv`, under `results/stroma/figure4/` or `results/epithelium/suppfigure10/`. Input provenance includes a path and MD5 fingerprint; the bundled input manifest separately provides SHA-256 checksums.

The selected pathway IDs, order, names, gene lists/counts, NEASE scores, raw/adjusted P values and transformed color values match the saved corrected `Reactome_top10_for_figure.tsv` tables. Floating-point serialization is compared with relative tolerance 1e-12 and absolute tolerance 1e-14. All ten points in each cohort have finite positions/sizes.

This is a portable R adaptation using the corrected numerical results. It retains the original R scripts' bubble-plot conventions and cohort colors. These numbered R scripts are the retained implementations for Figure 4D and Supplementary Figure 10B. No pixel-level equivalence to the assembled paper panels is claimed. Neither NEASE enrichment nor figure rendering was rerun during this addition; manuscript and presentation files were not changed.
