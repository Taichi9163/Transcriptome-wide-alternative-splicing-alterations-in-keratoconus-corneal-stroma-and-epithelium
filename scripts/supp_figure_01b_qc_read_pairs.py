#!/usr/bin/env python3
"""Supplementary Figure 1B: raw, filtered and mapped read pairs.

Adapted from the supplied Keratoconus_QC_ReadPairs.py. That source describes
itself as a reconstruction, not a verified copy of the original plotting code.
See docs/qc.md for provenance, required logs and the study/source sample mapping.

fastp total_reads counts both mates and is divided by two. STAR counts are
paired-end fragments and are not divided. Mapped = unique + multiple loci;
too-many-loci mappings are excluded. --tables-only never imports plotting tools.
"""
import argparse
import csv
import io
import json
import sys
from pathlib import Path

from supp_table_02_star_junction_qc import DEFAULT_SAMPLES, read_samples, sample_paths


def get_nonnegative_int(value, description: str):
    if isinstance(value, bool):
        raise ValueError(f"{description} is not an integer: {value!r}")
    try:
        count = int(str(value).replace(",", "").strip())
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{description} is not an integer: {value!r}") from exc
    if count < 0:
        raise ValueError(f"{description} must be nonnegative")
    return count


def read_fastp(path: Path):
    with path.open(encoding="utf-8") as f:
        obj = json.load(f)
    try:
        before = obj["summary"]["before_filtering"]["total_reads"]
        after = obj["summary"]["after_filtering"]["total_reads"]
    except (KeyError, TypeError) as exc:
        raise ValueError(f"Missing fastp summary read counts in {path}") from exc
    raw = get_nonnegative_int(before, f"fastp before_filtering in {path}")
    filtered = get_nonnegative_int(after, f"fastp after_filtering in {path}")
    if raw % 2 or filtered % 2:
        raise ValueError(
            f"Odd fastp read count in {path}: raw={raw}, filtered={filtered}. "
            "The script assumes paired-end processing and counts both mates."
        )
    if filtered > raw:
        raise ValueError(f"Filtered reads exceed raw reads in {path}")
    return raw // 2, filtered // 2


def read_star_log(path: Path):
    fields = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            if "|" in line:
                k, v = line.split("|", 1)
                key = k.strip()
                if key in fields:
                    raise ValueError(f"Duplicate STAR metric {key}: {path}")
                fields[key] = v.strip()
    keys = (
        "Number of input reads",
        "Uniquely mapped reads number",
        "Number of reads mapped to multiple loci",
    )
    try:
        input_pairs, unique_pairs, multi_pairs = [
            get_nonnegative_int(fields[k], f"STAR {k} in {path}") for k in keys
        ]
    except KeyError as exc:
        raise ValueError(f"Missing STAR Log.final.out field {exc} in {path}") from exc
    if unique_pairs + multi_pairs > input_pairs:
        raise ValueError(f"Mapped reads exceed STAR input reads in {path}")
    return input_pairs, unique_pairs, multi_pairs


def make_plot(rows, outdir: Path, dpi: int):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib and numpy: pip install matplotlib numpy") from exc

    x = np.arange(len(rows))
    width = 0.23
    fig, ax = plt.subplots(figsize=(11.5, 5.5))
    for i, (key, caption, color) in enumerate(
        (("Raw read pairs", "Raw Reads", "#080808"),
         ("Filtered read pairs", "Filtered Reads", "#818181"),
         ("Mapped read pairs", "Mapped Reads", "#bcbcbc"))
    ):
        ax.bar(x + (i - 1) * width, [r[key] / 1_000_000 for r in rows],
               width=width, label=caption, color=color)
    ax.set_ylabel("Number of read pairs (million)", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels([r["Sample"] for r in rows], rotation=43, ha="right")
    ax.set_xlim(-0.6, len(rows) - 0.4)
    ax.set_ylim(bottom=0)
    ax.yaxis.grid(True, linewidth=0.55, color="0.75")
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", frameon=False, ncol=3)
    fig.tight_layout()
    png = outdir / "Supplementary_Figure_1B_ReadPairs.png"
    pdf = outdir / "Supplementary_Figure_1B_ReadPairs.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def fastp_path(directory, sample, layout="flat", file_labels="study"):
    prefix = sample["sample_id" if file_labels == "study" else "source_id"]
    if layout == "per-sample":
        return directory / prefix / "report.json"
    return directory / f"{prefix}_report.json"


def collect_rows(fastp_dir, star_dir, samples, layout="flat", file_labels="study"):
    """Use explicit sample identities; never guess a study label from KC_ numbers."""
    rows = []
    for sample in samples:
        report = fastp_path(fastp_dir, sample, layout, file_labels)
        log, _ = sample_paths(star_dir, sample, layout, file_labels)
        raw, filtered = read_fastp(report)
        star_input, unique, multiple = read_star_log(log)
        # These quantities must describe the same paired-end input. The supplied
        # script only warned; stop here to avoid an internally inconsistent graph.
        if star_input != filtered:
            raise ValueError(
                f"{sample['sample_id']}: STAR input pairs ({star_input:,}) differ from "
                f"fastp filtered pairs ({filtered:,}); check the input reports and sample mapping."
            )
        rows.append({
            "Sample": sample["sample_name"],
            "Raw read pairs": raw,
            "Filtered read pairs": filtered,
            "Mapped read pairs": unique + multiple,
            "Uniquely mapped read pairs": unique,
            "Multimapped read pairs": multiple,
            "STAR input read pairs": star_input,
        })
    return rows


def write_table(path, rows):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    payload = buffer.getvalue().encode("utf-8")
    if path.exists():
        if path.read_bytes() != payload:
            raise FileExistsError(f"Existing output differs: {path}. Choose a new --output-dir.")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as handle:
        handle.write(payload)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastp-dir", type=Path, required=True)
    parser.add_argument("--star-dir", type=Path, required=True)
    parser.add_argument("--samples", type=Path, default=DEFAULT_SAMPLES, help="QC metadata TSV in output order")
    parser.add_argument("--layout", choices=("flat", "per-sample"), default="flat")
    parser.add_argument("--file-labels", choices=("study", "source"), default="study",
                        help="Study labels (KC1/Control1) or original prefixes (KC_3/N_1)")
    parser.add_argument("--output-dir", type=Path, default=Path("results/qc/suppfigure1"))
    parser.add_argument("--tables-only", action="store_true", help="Write numeric TSV without generating images")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if args.dpi < 50:
        parser.error("--dpi must be at least 50")
    fastp_dir, star_dir = args.fastp_dir.expanduser(), args.star_dir.expanduser()
    if not fastp_dir.is_dir() or not star_dir.is_dir():
        parser.error("Both --fastp-dir and --star-dir must exist")
    samples = read_samples(args.samples.expanduser())
    rows = collect_rows(fastp_dir, star_dir, samples, args.layout, args.file_labels)
    outdir = args.output_dir.expanduser()
    if not args.tables_only:
        for extension in ("png", "pdf"):
            path = outdir / f"Supplementary_Figure_1B_ReadPairs.{extension}"
            if path.exists():
                raise FileExistsError(f"Existing figure preserved: {path}. Choose a new --output-dir.")
    table = outdir / "Keratoconus_Raw_Filtered_Mapped_Reads.tsv"
    write_table(table, rows)
    print("Sample\tRaw (M pairs)\tFiltered (M pairs)\tMapped (M pairs)")
    for row in rows:
        print(f"{row['Sample']}\t{row['Raw read pairs']/1e6:.3f}\t"
              f"{row['Filtered read pairs']/1e6:.3f}\t{row['Mapped read pairs']/1e6:.3f}")
    print(f"Saved: {table}")
    if not args.tables_only:
        png, pdf = make_plot(rows, outdir, args.dpi)
        print(f"Saved: {png}\nSaved: {pdf}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
