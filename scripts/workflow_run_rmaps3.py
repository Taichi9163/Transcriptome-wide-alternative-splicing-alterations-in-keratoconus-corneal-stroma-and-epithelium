#!/usr/bin/env python3
"""Audit SE exon selection or preview/run rMAPS3 for Fig. 1E,F / Supp. Fig. 7E,F.

The audit uses only the standard library and never opens a graphics device.
The native rMAPS3 engine runs only with `run --execute` and also creates figures.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

from utils_paths import add_work_dir_argument, resolve_path
import shlex
import subprocess
import sys

from setup_prepare_inputs import install_input

REPOSITORY = Path(__file__).resolve().parents[1]
COORD_COLUMNS = (
    "chr",
    "strand",
    "exonStart_0base",
    "exonEnd",
    "upstreamES",
    "upstreamEE",
    "downstreamES",
    "downstreamEE",
)
COORD_HEADER = (
    "chr",
    "strand",
    "exonStart",
    "exonEnd",
    "firstExonStart",
    "firstExonEnd",
    "secondExonStart",
    "secondExonEnd",
)


def sha256(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def number(value):
    """Native rMATS missing values are excluded from numeric comparisons."""
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def mean_psi(value):
    values = [number(v) for v in value.replace('"', "").split(",")]
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else None


def select_exons(rows, settings):
    """Reconstruct the logged native SE selection, preserving input order.

    Identity is the target exon (chromosome, strand, start, end). Within a group,
    the last row for an identity supplies the flanking coordinates. Exons in
    more than one group are excluded from all groups. Coordinates are kept in
    genomic order on both strands; the native engine handles strand orientation.
    """
    groups = {name: {} for name in ("up", "dn", "bg")}
    exact_boundary = 0
    for row in rows:
        coords = tuple(row[column] for column in COORD_COLUMNS)
        fdr, delta = number(row["FDR"]), number(row["IncLevelDifference"])
        if fdr is None:
            continue
        group = None
        if fdr < settings["sig_fdr"] and delta is not None:
            if abs(delta) == settings["sig_delta_psi"]:
                exact_boundary += 1
            # The logged native rMAPS selection is inclusive at +/-0.10.
            if delta >= settings["sig_delta_psi"]:
                group = "up"
            elif delta <= -settings["sig_delta_psi"]:
                group = "dn"
        elif fdr > settings["background_fdr"]:
            means = [mean_psi(row[column]) for column in ("IncLevel1", "IncLevel2")]
            if (
                all(v is not None for v in means)
                and min(means) < settings["background_min_mean_psi_upper"]
                and max(means) > settings["background_max_mean_psi_lower"]
            ):
                group = "bg"
        if group:
            groups[group][coords[:4]] = coords
    initial = {name: len(group) for name, group in groups.items()}
    conflicts = (
        set(groups["up"]) & set(groups["dn"])
        | set(groups["up"]) & set(groups["bg"])
        | set(groups["dn"]) & set(groups["bg"])
    )
    selected = {
        name: [coords for key, coords in group.items() if key not in conflicts]
        for name, group in groups.items()
    }
    # No substitute background is allowed in this study reproduction workflow.
    if not selected["bg"]:
        raise ValueError("No strict background exons; refusing a fallback background.")
    return selected, {
        "initial_counts": initial,
        "conflicting_target_exons": len(conflicts),
        "significant_rows_at_exact_delta_cutoff": exact_boundary,
    }


def coordinate_bytes(rows):
    return ("\t".join(COORD_HEADER) + "\n" + "".join("\t".join(row) + "\n" for row in rows)).encode(
        "utf-8"
    )


def audit(input_path, cohort, config):
    expected = config["cohorts"][cohort]
    actual_hash = sha256(input_path)
    if actual_hash != expected["rmats_sha256"]:
        raise ValueError(f"{cohort}: SE input differs from the recorded study input.")
    with input_path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = set(COORD_COLUMNS) | {"FDR", "IncLevelDifference", "IncLevel1", "IncLevel2"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("The input must be a complete rMATS SE JC table.")
        selected, summary = select_exons(reader, config["settings"])
    payloads = {name: coordinate_bytes(rows) for name, rows in selected.items()}
    counts = {name: len(rows) for name, rows in selected.items()}
    hashes = {name: hashlib.sha256(data).hexdigest() for name, data in payloads.items()}
    if counts != expected["selected_counts"] or hashes != expected["coordinate_sha256"]:
        raise ValueError(f"{cohort}: selected coordinates differ from the saved rMAPS3 run.")
    summary.update(
        cohort=cohort,
        input_sha256=actual_hash,
        selected_counts=counts,
        coordinate_sha256=hashes,
        matches_saved_coordinates=True,
        settings=config["settings"],
    )
    return payloads, summary


def save_bytes(path, data):
    install_input(path, hashlib.sha256(data).hexdigest(), lambda target: target.write_bytes(data))


def make_command(python, rmaps_dir, fasta_root, input_path, output_dir, config):
    settings = config["settings"]
    command = [
        str(python),
        str(rmaps_dir / "legacy/motifMapSE_MP.py"),
        "-k",
        str(rmaps_dir / "data/knownMotifs.human.mouse.txt"),
        "-m",
        "NA",
        "--fasta-root",
        str(fasta_root),
        "-g",
        "hg38",
        "-o",
        str(output_dir),
        "-r",
        str(input_path),
        "-mi",
        "NA",
        "-u",
        "NA",
        "-d",
        "NA",
        "-b",
        "NA",
        "--label",
        "RBP",
    ]
    for flag, key in (
        ("--intron", "intron"),
        ("--exon", "exon"),
        ("--window", "window"),
        ("--step", "step"),
        ("--sigFDR", "sig_fdr"),
        ("--sigDeltaPSI", "sig_delta_psi"),
    ):
        command.extend([flag, str(settings[key])])
    return command


def execute(command, rmaps_dir, fasta_root, output_dir, config, summary):
    engine = rmaps_dir / "legacy/motifMapSE_MP.py"
    if not engine.is_file():
        raise ValueError(f"Missing rMAPS3 engine: {engine}")
    motif = rmaps_dir / "data/knownMotifs.human.mouse.txt"
    if not motif.is_file() or sha256(motif) != config["motif_sha256"]:
        raise ValueError("The human/mouse motif library does not match the recorded 114 motifs.")
    # Use the documented single-FASTA layout, with UCSC chr names.
    fasta = fasta_root / "hg38/hg38.fa"
    if not fasta.is_file():
        raise ValueError(f"Missing hg38 FASTA: {fasta}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError(f"Output directory is not empty; choose a fresh directory: {output_dir}")
    env = os.environ.copy()
    env["RMAPS_STAT_METHOD"] = "fisher"
    env["PYTHONPATH"] = str(rmaps_dir) + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    # The actual interpreter, dependency versions, source and genome belong to
    # each new run. This does not assert identity to the June 2026 environment.
    version = subprocess.run([command[0], "--version"], check=True, capture_output=True, text=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance = {
        "command": command,
        "python_version": (version.stdout + version.stderr).strip(),
        "engine_sha256": sha256(engine),
        "motif_sha256": sha256(motif),
        "fasta": str(fasta),
        "fasta_bytes": fasta.stat().st_size,
        "stat_method": "fisher",
        "selection_audit": summary,
        "status": "started",
    }
    provenance_path = output_dir / "run_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    try:
        subprocess.run(command, cwd=rmaps_dir, env=env, check=True)
        for group, expected in summary["coordinate_sha256"].items():
            if sha256(output_dir / "exon" / f"{group}.coord.txt") != expected:
                raise ValueError("Native rMAPS3 exon selection differs from the saved study run.")
    except (OSError, ValueError, subprocess.CalledProcessError):
        provenance["status"] = "failed"
        raise
    else:
        provenance["status"] = "completed"
    finally:
        provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("audit", "run"))
    parser.add_argument("--cohort", choices=("stroma", "epithelium", "both"), default="both")
    add_work_dir_argument(parser)
    parser.add_argument("--rmaps-dir", type=Path, help="External rMAPS3 source checkout")
    parser.add_argument("--fasta-root", type=Path, help="Directory containing hg38/hg38.fa")
    parser.add_argument(
        "--python", default=sys.executable, help="Python interpreter with rMAPS3 dependencies"
    )
    parser.add_argument("--output-root", type=Path, help="Run output base; default: WORK/results")
    parser.add_argument(
        "--execute", action="store_true", help="Run native rMAPS3 (also generates its figures)"
    )
    args = parser.parse_args()
    if args.mode == "audit" and args.execute:
        parser.error("audit never executes rMAPS3; remove --execute")
    if args.mode == "run" and (args.rmaps_dir is None or args.fasta_root is None):
        parser.error("run requires --rmaps-dir and --fasta-root")
    config = json.loads((REPOSITORY / "config/rmaps3.json").read_text())
    work = args.work_dir.expanduser().resolve()
    cohorts = ("stroma", "epithelium") if args.cohort == "both" else (args.cohort,)
    for cohort in cohorts:
        input_path = work / "data" / cohort / "rmats/SE.MATS.JC.txt"
        if args.mode == "audit" or args.execute:
            payloads, summary = audit(input_path, cohort, config)
        if args.mode == "audit":
            target = work / "results" / cohort / "rmaps3_audit"
            for group, data in payloads.items():
                save_bytes(target / f"{group}.coord.txt", data)
            save_bytes(
                target / "selection_audit.json", (json.dumps(summary, indent=2) + "\n").encode()
            )
            print(f"{cohort}: matched saved exon coordinates {summary['selected_counts']}")
            continue
        rmaps_dir, fasta_root = resolve_path(args.rmaps_dir, work), resolve_path(
            args.fasta_root, work
        )
        output_root = resolve_path(args.output_root, work) if args.output_root else work / "results"
        output_dir = output_root / cohort / "rmaps3"
        python = str(resolve_path(args.python, work)) if "/" in args.python else args.python
        command = make_command(python, rmaps_dir, fasta_root, input_path, output_dir, config)
        print(
            f"[{cohort}] RMAPS_STAT_METHOD=fisher PYTHONPATH={shlex.quote(str(rmaps_dir))} "
            + shlex.join(command)
        )
        if args.execute:
            execute(command, rmaps_dir, fasta_root, output_dir, config, summary)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
