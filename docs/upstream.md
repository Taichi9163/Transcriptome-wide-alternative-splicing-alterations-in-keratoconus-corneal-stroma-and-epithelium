# Ubuntu read-processing workflow

This is a path-safe reconstruction of the supplied commands, not a record of a newly executed FASTQ analysis. The manuscript reports fastp 0.20.0, STAR 2.7.11a, RSEM 1.3.1 and rMATS-turbo 4.3.0, with GRCh38/Ensembl release 114. Original FastQC, MultiQC and sashimi tool versions were not supplied. Record the versions of all executables when running the workflow. The original environment was named `shortread`; its environment lock was not supplied, so an untested replacement lock is not presented as the original environment.

## Set up paths

1. Activate an Ubuntu environment containing the tools required by the selected stage.
2. Copy and edit `config/upstream.json`. Set `genome_fasta`, `annotation_gtf`, `output_dir`, `threads`, `rmats_command` and `sashimi_command`.
3. Copy and edit `config/samples_stroma.tsv`. FASTQ paths are relative to `--work-dir`, or may be absolute. Do not change the study labels to match the original file prefixes.
4. Supply the GRCh38 primary-assembly FASTA and **uncompressed** release-114 GTF. The `.gtf.gz` downloaded by `prepare_inputs.py` is used directly for downstream mapping; extract a separate `.gtf` copy for STAR/RSEM/rMATS.

The default rMATS command is `rmats.py` on `PATH`. For a source installation set `rmats_command` to `["python", "/path/to/rmats-turbo/rmats.py"]`. Commands are argument lists, not shell strings. Spaces in directory names are supported. Because rMATS BAM lists are comma-separated, BAM paths must not contain commas or line breaks.

## Run stages

```sh
python3 scripts/preprocess_reads.py all --config config/upstream.json \
  --samples config/samples_stroma.tsv --work-dir /path/to/analysis
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

After `quantify`, `sample_sheets/genes.tsv` and `sample_sheets/isoforms.tsv` contain absolute paths for the nine primary specimens and can be passed to `differential_expression.R`. These generated paths are local outputs and are not committed to GitHub.

## Connect to downstream modules

Copy or link the five rMATS JC tables from `results/upstream/rmats/primary/` into `data/stroma/rmats/` **in a new work directory**. The inclusive outputs belong under `data/stroma_with_kc6/rmats/`. Keep read-derived results separate from the checksum-verified bundled inputs, because new software versions can change numerical results.

Run `call_das.R stroma`, then `selected_exons.R` from that work directory. The latter writes a sample-level PSI table and `results/stroma/sashimi_events.tsv` for the ZEB2, ROCK2 and TPM2 events in the manuscript. To render sashimi plots:

```sh
python3 /path/to/repository/scripts/preprocess_reads.py sashimi \
  --work-dir /path/to/analysis \
  --events /path/to/analysis/results/stroma/sashimi_events.tsv --execute
```

The generated grouping file uses cases 1–5 and controls 6–9, following the [rmats2sashimiplot group-file convention](https://github.com/Xinglab/rmats2sashimiplot). Colors and exon/intron scaling follow the supplied commands. Raw BAMs are required; they are not included in the repository.
