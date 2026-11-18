#!/usr/bin/env python3
"""Reproduce Figure 1B-D and Figure 3A-D, or validate only their numeric tables."""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from prepare_inputs import GTF_NAME, prepare

SCRIPTS = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=Path.cwd())
    parser.add_argument("--tables-only", action="store_true", help="Build plot layers without generating images")
    args = parser.parse_args()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    rscript = os.environ.get("RSCRIPT", "Rscript")
    if shutil.which(rscript) is None:
        parser.error("Rscript was not found. See the installation instructions in README.md.")
    reference = work / "data/reference" / GTF_NAME
    if not reference.is_file():
        parser.error("Run prepare_inputs.py with --gtf or --download-reference first.")
    prepare(work)
    environment = dict(os.environ, KC_WORK_DIR=str(work))
    steps = [
        [rscript, "--vanilla", str(SCRIPTS / "Setup_check_R_packages.R")],
        [rscript, "--vanilla", str(SCRIPTS / "Setup_prepare_Figure1_Figure3_data.R")],
        [sys.executable, str(SCRIPTS / "map_events.py"), "--work-dir", str(work)],
    ]
    panels = (
        "Figure1BC_PSI_PCA_heatmap",
        "Figure1D_DAS_gene_UpSet",
        "Figure3A_DEG_DAS_Venn",
        "Figure3B_gene_expression_vs_deltaPSI",
        "Figure3C_isoform_DAS_Venn",
        "Figure3D_isoform_expression_vs_deltaPSI",
    )
    for name in panels:
        steps.append([rscript, "--vanilla", str(SCRIPTS / f"{name}.R")]
                     + (["--tables-only"] if args.tables_only else []))
    for command in steps:
        print("Running:", Path(command[2] if command[1] == "--vanilla" else command[1]).name, flush=True)
        subprocess.run(command, cwd=work, env=environment, check=True)
    print(f"Complete: {work / 'results/stroma'}")


if __name__ == "__main__":
    main()
