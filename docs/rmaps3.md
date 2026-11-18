# rMAPS3 execution and exon-selection audit

[Repository overview](../README.md) · [Documentation](README.md) · [Figure-to-code index](figures.md)

`scripts/workflow_run_rmaps3.py` supplies the upstream SE motif workflow for Figure 1E,F and Supplementary Figure 7E,F. It reconstructs the commands recorded in both study `log.motifMap.txt` files on 30 June 2026, replacing workstation paths with command-line options. It does not vendor or modify the external rMAPS3 engine.

## Verify selection without generating figures

```sh
python3 scripts/setup_prepare_inputs.py --work-dir work
python3 scripts/workflow_run_rmaps3.py audit --work-dir work
```

This mode needs Python 3.9+ only. It checks the complete SE input checksum, reconstructs all exon groups, and compares every coordinate and row order with checksums of the saved native coordinate files. It writes three `*.coord.txt` files and `selection_audit.json` under `work/results/<cohort>/rmaps3_audit/`. Identical files can be reused; different files are not overwritten. It does not read genome sequences, calculate motif P values, or run any plotting code.

| Cohort | Increased inclusion (`up`) | Increased skipping (`dn`) | Background (`bg`) |
| --- | --- | --- | --- |
| Stroma | 367 | 171 | 12,955 |
| Epithelium | 163 | 136 | 15,719 |

The input is the complete JC table, including nonsignificant events needed for the background, with keratoconus as rMATS group 1. Exon identity is chromosome, strand, target-exon start and target-exon end. The last qualifying row within each group supplies its flanking coordinates; an exon assigned to more than one group is removed from every group. Coordinates remain in genomic order on both strands.

The native rMAPS3 threshold is **FDR < 0.05 and |ΔPSI| ≥ 0.10**. This inclusive ΔPSI boundary differs from the strict `> 0.10` rule used by the repository's general DAS summaries. Eight stromal and three epithelial significant input rows are exactly at ±0.10 before deduplication. Keeping the native boundary reproduces the saved rMAPS3 coordinate files exactly; do not substitute the general DAS table. Background requires FDR > 0.5, the lower group-mean PSI < 0.85, and the higher > 0.15. Missing PSI entries are omitted within each group; a group with no numeric PSI is ineligible. No fallback background is allowed.

## Prepare the external engine

Use the [official rMAPS3 repository](https://github.com/UofLBioinformatics/rMAPS3). The manuscript identifies version 3.1.0; the saved log reports its internal motifTools engine version as 2.0.1. The original Git revision and environment were not recorded. These labels refer to different layers and do not establish an exact source revision.

For a concrete source reference, the wrapper interface and SE selection were inspected at [revision b9a9dce](https://github.com/UofLBioinformatics/rMAPS3/tree/b9a9dce9ad8596fdef21c3a705fbc5f37d9fa45d). This is an inspected revision, **not a claim that it was used for the original analysis**. The original native motif scan has not been rerun for this release.

```sh
git clone https://github.com/UofLBioinformatics/rMAPS3.git /path/to/rMAPS3
git -C /path/to/rMAPS3 checkout b9a9dce9ad8596fdef21c3a705fbc5f37d9fa45d
python3.12 -m venv /path/to/rmaps-env
/path/to/rmaps-env/bin/python -m pip install -r /path/to/rMAPS3/requirements.txt
```

Follow the upstream [genome and plotting setup](https://github.com/UofLBioinformatics/rMAPS3/blob/b9a9dce9ad8596fdef21c3a705fbc5f37d9fa45d/README.md) for hg38, PyX/TeX and optional rendering dependencies. This wrapper uses the single-FASTA layout `/path/to/genomedata/hg38/hg38.fa`, with UCSC chromosome names matching the rMATS tables. It does not download a genome automatically. The external `data/knownMotifs.human.mouse.txt` must match the saved 114-motif library; its SHA-256 is in `config/rmaps3.json`. The inspected upstream library is byte-identical to both saved study inputs. External software and data retain their own terms.

## Preview or execute

```sh
# Preview only: no external tools run and no files are created.
python3 scripts/workflow_run_rmaps3.py run --work-dir work \
  --rmaps-dir /path/to/rMAPS3 --fasta-root /path/to/genomedata \
  --python /path/to/rmaps-env/bin/python
```

Add `--cohort stroma` or `--cohort epithelium` for one cohort. Add `--execute` only when intending to run the complete native analysis: **the external engine also generates figures**. Execution is deliberately separate from the numeric-only audit. The wrapper sets `RMAPS_STAT_METHOD=fisher` explicitly, matching the manuscript's one-sided Fisher test, and adds the source checkout to `PYTHONPATH`.

The logged settings are hg38, 250 intronic nucleotides, 50 exonic nucleotides, a 50-nt window, a 1-nt step, FDR 0.05, ΔPSI 0.10, and label `RBP`. Outputs go to `work/results/<cohort>/rmaps3/`; `--output-root` selects a different base. Existing nonempty output directories are refused. Each run records its command, Python version, engine/motif hashes, FASTA path/size, selection audit and completion status. It also checks the native coordinate outputs against the saved study sets. The FASTA path/size record is not a genome checksum or a dependency lockfile.

The native output includes `pVal.up.vs.bg.RNAmap.txt`, `pVal.dn.vs.bg.RNAmap.txt` and `log.motifMap.txt`. The repository's downstream R motif scripts read these names under `data/<cohort>/rmaps/`. Keep new native runs separate while assessing reproducibility; the bundled files in that directory are the original published-analysis inputs and are not replaced automatically. Downstream processing applies BH to 912 tests separately in each cohort and direction; it does not pool all four groups.
