# Ubuntu read-processing workflow

[Repository overview](../README.md) · [Documentation](README.md) · [Figure-to-code index](figures.md)

This is a path-safe reconstruction of the supplied commands, not a record of a newly executed FASTQ analysis. The manuscript reports FastQC 0.12.1, MultiQC 1.29, fastp 0.20.0, STAR 2.7.11a, RSEM 1.3.1, rMATS-turbo 4.3.0 and rmats2sashimiplot 3.0.0, with GRCh38/Ensembl release 114. Record the versions of all executables when running the workflow. The original environment was named `shortread`; its environment lock was not supplied, so an untested replacement lock is not presented as the original environment.

## Set up paths

1. Activate an Ubuntu environment containing the tools required by the selected stage.
2. Copy and edit `config/upstream.json`. Set `genome_fasta`, `annotation_gtf`, `output_dir`, `threads`, `rmats_command` and `sashimi_command`.
3. Copy and edit `config/samples_stroma.tsv`. FASTQ paths are relative to `--work-dir`, or may be absolute. Do not change the study labels to match the original file prefixes.
4. Supply the GRCh38 primary-assembly FASTA and **uncompressed** release-114 GTF. The `.gtf.gz` downloaded by `setup_prepare_inputs.py` is used directly for downstream mapping; extract a separate `.gtf` copy for STAR/RSEM/rMATS.

The default rMATS command is `rmats.py` on `PATH`. For a source installation set `rmats_command` to `["python", "/path/to/rmats-turbo/rmats.py"]`. Commands are argument lists, not shell strings. Spaces in directory names are supported. Because rMATS BAM lists are comma-separated, BAM paths must not contain commas or line breaks.

## Run stages

```sh
# Run from the repository root, after editing copies of the settings as needed.
python3 scripts/workflow_preprocess_reads.py all --config "$PWD/config/upstream.json" \
  --samples "$PWD/config/samples_stroma.tsv" --work-dir /path/to/analysis
```

This previews the commands without creating files. Add `--execute` to run them. For a first run, execute stages individually in this order:

| Stage | Processing and output beneath `output_dir` |
| --- | --- |
| `raw-qc` | FastQC for each specimen, then MultiQC separately for cases and controls |
| `trim` | Paired-end fastp with HTML/JSON reports |
| `trimmed-qc` | FastQC on trimmed reads, then group-specific MultiQC |
| `star-index` | STAR reference index |
| `align` | Coordinate-sorted and transcriptome BAMs, using distinct temporary paths per specimen |
| `rsem-index` | RSEM transcript reference |
| `quantify` | Gene and isoform estimates; primary-cohort sample sheets |
| `rmats` | Primary 5-versus-4 comparison and KC6-inclusive 6-versus-4 comparison |
| `sashimi` | Selected SE panels; requires `--events`; omitted from `all` |

Each job writes its command, input size/modification-time fingerprints, run log and success marker into its output directory. A completed job with the same inputs and command is skipped on rerun. An incomplete or conflicting output directory stops execution; preserve it for diagnosis and choose a new `output_dir` for a fresh run. Fingerprints are a practical resume check, not content hashes of the FASTQ/BAM files. Use a new output directory after changing software or reference versions.

The wrapper preserves the supplied analysis options: fastp defaults; STAR `--outSAMtype BAM SortedByCoordinate` and `--quantMode TranscriptomeSAM`; RSEM `--paired-end --estimate-rspd --output-genome-bam`; rMATS paired reads, read length 150, `--allow-clipping` and `--variable-read-length`. It leaves other settings at tool defaults, including STAR `sjdbOverhang` and rMATS library type, anchor and cstat. The downstream |ΔPSI| > 0.1 threshold is **not** substituted for rMATS's statistical cstat parameter. See the [rMATS documentation](https://github.com/Xinglab/rmats-turbo) and [STAR documentation](https://github.com/alexdobin/STAR).

## Study labels and BAM order

The source-ID mapping was checked against the study's RSEM sample annotation. In particular, study KC6 is original `KC_1`; study KC3 is original `KC_6`. The primary comparison excludes study KC6, while preprocessing and the inclusive comparison retain it.

For new rMATS runs, cases are group 1 and controls are group 2; within-group order follows the sample sheet. Each rMATS output directory includes `sample_order.tsv`, `cases.txt` and `controls.txt`. The historical bundled rMATS tables preserve their existing PSI-vector order. Their displayed numbered labels identify within-group positions; this release does not assert that the historical BAM list order has been independently recovered from those tables.

After `quantify`, `sample_sheets/genes.tsv` and `sample_sheets/isoforms.tsv` contain absolute paths for the nine primary specimens and can be passed to `figure_03abcd_deseq2_analysis.R`. These generated paths are local outputs and are not committed to GitHub.

After `align`, run `supp_table_02_star_junction_qc.py --star-dir /path/to/analysis/results/upstream/align --layout per-sample --output /path/to/analysis/results/qc/Supplementary_Table2_STAR_junction_QC.tsv` to summarize unique-mapping and junction support across all ten specimens. It reads the alignment logs, uses the QC metadata in study-label order, and leaves the alignment outputs unchanged. For Supplementary Figure 1B, `supp_figure_01b_qc_read_pairs.py` reads both the `trim` fastp JSON reports and the `align` STAR logs; use `--layout per-sample --tables-only` to export counts without images. See [Supplementary Figure 1 / Table 2 QC](qc.md) for complete commands, original flat-file layouts and the remaining report sources.

## Connect to downstream modules

Copy or link the five rMATS JC tables from `results/upstream/rmats/primary/` into `data/stroma/rmats/` **in a new work directory**. The inclusive outputs belong under `data/stroma_with_kc6/rmats/`. Keep read-derived results separate from the checksum-verified bundled inputs, because new software versions can change numerical results.

Run `figure_01a_das_event_summary.R stroma`, then `figure_02def_selected_exon_psi.R` from that work directory. The latter writes a sample-level PSI table and `results/stroma/sashimi_events.tsv` for the ZEB2, ROCK2 and TPM2 events in the manuscript. To render sashimi plots:

```sh
python3 /path/to/repository/scripts/workflow_preprocess_reads.py sashimi \
  --work-dir /path/to/analysis \
  --events /path/to/analysis/results/stroma/sashimi_events.tsv --execute
```

The generated grouping file uses cases 1–5 and controls 6–9, following the [rmats2sashimiplot group-file convention](https://github.com/Xinglab/rmats2sashimiplot). Colors and exon/intron scaling follow the supplied commands. Raw BAMs are required; they are not included in the repository.
