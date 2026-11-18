#!/usr/bin/env python3
"""Reproduce Figure 1B-D and Figure 3A-D, or validate only their numeric tables."""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from utils_paths import add_work_dir_argument

from setup_prepare_inputs import GTF_NAME, prepare

SCRIPTS = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_work_dir_argument(parser)
    parser.add_argument(
        "--tables-only", action="store_true", help="Build plot layers without generating images"
    )
    args = parser.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    rscript = os.environ.get("RSCRIPT", "Rscript")
    if shutil.which(rscript) is None:
        parser.error("Rscript was not found. See the installation instructions in README.md.")
    reference = work / "data/reference" / GTF_NAME
    if not reference.is_file():
        parser.error("Run setup_prepare_inputs.py with --gtf or --download-reference first.")
    prepare(work)
    environment = dict(os.environ, KC_WORK_DIR=str(work))
    steps = [
        [rscript, "--vanilla", str(SCRIPTS / "setup_check_r_packages.R")],
        [rscript, "--vanilla", str(SCRIPTS / "setup_prepare_manuscript_data.R")],
        [
            sys.executable,
            str(SCRIPTS / "figure_03d_event_transcript_mapping.py"),
            "--work-dir",
            str(work),
        ],
    ]
    panels = (
        "figure_01bc_psi_pca_heatmap",
        "figure_01d_das_gene_upset",
        "figure_03a_deg_das_venn",
        "figure_03b_gene_expression_vs_delta_psi",
        "figure_03c_isoform_das_venn",
        "figure_03d_isoform_expression_vs_delta_psi",
    )
    for name in panels:
        steps.append(
            [rscript, "--vanilla", str(SCRIPTS / f"{name}.R")]
            + (["--tables-only"] if args.tables_only else [])
        )
    for command in steps:
        print(
            "Running:",
            Path(command[2] if command[1] == "--vanilla" else command[1]).name,
            flush=True,
        )
        subprocess.run(command, cwd=work, env=environment, check=True)
    print(f"Complete: {work / 'results/stroma'}")


if __name__ == "__main__":
    main()
