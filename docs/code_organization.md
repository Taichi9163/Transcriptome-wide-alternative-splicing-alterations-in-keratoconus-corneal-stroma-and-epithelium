# Code organization and portable paths

[Repository overview](../README.md) · [Figure-to-code index](figures.md) · [Reproduction guide](reproduce.md)

## Choose one work directory

Every analysis entry point in `scripts/` accepts `--work-dir`. The selected directory contains `data/` and `results/`; it is independent of the location of the repository or the directory from which a command is launched.

Selection follows the same order in R and Python:

1. An explicit `--work-dir` argument.
2. The `KC_WORK_DIR` environment variable.
3. The current working directory when neither is supplied.

For example, a Linux user could keep the analysis in `/home/user/keratoconus`. This is an example, not a required location or a hard-coded default. Replace it with a directory you own; paths containing spaces should be quoted.

```sh
# Run these commands from the cloned repository.
python3 scripts/setup_prepare_inputs.py --work-dir /home/user/keratoconus --download-reference
python3 scripts/workflow_run_manuscript.py --work-dir /home/user/keratoconus --tables-only

Rscript --vanilla scripts/figure_01a_das_event_summary.R stroma --work-dir /home/user/keratoconus
Rscript --vanilla scripts/figure_04d_reactome_enrichment.R --work-dir /home/user/keratoconus --tables-only
python3 scripts/workflow_run_rmaps3.py audit --work-dir /home/user/keratoconus
```

Alternatively, set the directory once with `export KC_WORK_DIR=/home/user/keratoconus`. Explicit `--work-dir` arguments still take priority. Scripts do not change the caller's working directory.

## Input and output path rules

Relative input/output arguments, configuration paths, sample-sheet read paths and RSEM paths are resolved against the selected work directory. Absolute paths and `~` paths are also accepted. Relative paths inside sample sheets are relative to the work directory, not to the sample sheet's parent folder.

The bundled default configuration and sample annotation files are located relative to the repository itself. External reference genomes, FASTQ files, native tool checkouts and quantifications remain external inputs; the scripts do not require the authors' directory layout.

```text
/home/user/keratoconus/
├── data/
│   ├── reference/
│   ├── stroma/
│   ├── epithelium/
│   └── stroma_with_kc6/
└── results/
    ├── stroma/
    ├── epithelium/
    ├── cross_layer/
    └── upstream/
```

## Reading and modifying the code

The [filename convention](script_naming.md) distinguishes numbered figure scripts, workflow commands, setup commands and shared helpers. The [figure index](figures.md) maps manuscript panels to their implementations and inputs.

| Shared component | Responsibility |
| --- | --- |
| [utils_paths.R](../scripts/utils_paths.R), [utils_paths.py](../scripts/utils_paths.py) | Work-directory selection and path resolution |
| [utils_plotting.R](../scripts/utils_plotting.R) | Tabular output and manuscript figure export |
| [utils_psi.R](../scripts/utils_psi.R) | Validated PSI-vector parsing shared by the three optional PCA analyses |
| [utils_splicing_events.R](../scripts/utils_splicing_events.R) | Complete event-coordinate identity and directional memberships |
| [utils_nease_feature_summary.R](../scripts/utils_nease_feature_summary.R) | Protein-feature membership summaries |
| [utils_reactome_enrichment.R](../scripts/utils_reactome_enrichment.R) | Reactome selection, numerical checks and plot construction |

R uses two-space indentation and Python uses four-space indentation. Python command-line processing is inside `main()`; importing a module does not start an analysis or create result directories. NEASE-specific libraries are loaded only when its analysis is requested. The optional Python Reactome implementation separates pathway selection from plotting.

The refactor preserves the study's thresholds, sample order, event definitions, output column names and result filenames. Selection rules that differ between analyses remain distinct: in particular, the native rMAPS3 inclusive ΔPSI boundary is not replaced with the strict boundary used in general DAS summaries.

## Numerical checks without figures

Use `--tables-only` for supported plotting scripts. It also works for the motif heatmaps, the optional Python Reactome implementation and the three optional PCA scripts. Motif heatmaps export their selected adjusted-P matrices before any layout construction. Optional PCA runs still write their numerical caches, which can be large; they are separate from the main manuscript workflow.

Native rMAPS3 and upstream processing retain preview-by-default behavior; `--execute` runs the external commands. Native NEASE is a separate analysis and is not run by the manuscript table workflow. See the [validation record](validation.md) for what was executed and what was checked without rerunning the full native analysis.
