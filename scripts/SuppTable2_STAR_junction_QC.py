#!/usr/bin/env python3
"""Summarize the STAR mapping/junction columns of Supplementary Table 2.

Adapted from the supplied STAR_junction_QC.py. Uses Python's standard library;
does not run STAR or draw Supplementary Figure 1. All ten stromal samples,
including the specimen excluded from the primary comparison, are retained.

Study-labelled flat inputs are <sample_id>_Log.final.out and <sample_id>_SJ.out.tab.
Use --file-labels source for the original KC_*/N_* prefixes, or --layout per-sample
for <prefix>/<prefix>_Log.final.out and <prefix>/<prefix>_SJ.out.tab.

In this paired-end study, STAR's uniquely mapped read count represents read pairs;
it is not divided by two. Junction support is column 7, retaining values >= 1.
Quartiles use R type-7 interpolation. RIN is provided metadata, not a STAR metric.
"""

import argparse
import csv
import io
import math
import re
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
DEFAULT_SAMPLES = REPOSITORY / "config/samples_stroma_qc.tsv"
FIELDS = [
    "Sample", "RIN", "Unique mapped read pairs", "Unique mapped (%)",
    "Supported splice junctions", "Median junction support", "Q1", "Q3",
    "Median [Q1-Q3]",
]


def read_samples(path):
    """Keep the metadata row order and explicit study/source ID mapping."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "source_id", "sample_name", "rin"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Sample sheet requires columns: {sorted(required)}")
        samples = list(reader)
    if not samples:
        raise ValueError("Sample sheet is empty")
    for sample in samples:
        for key in required:
            if not sample.get(key):
                raise ValueError(f"Missing {key} in sample sheet")
        for key in ("sample_id", "source_id"):
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample[key]):
                raise ValueError(f"Invalid {key}: {sample[key]}")
        sample["rin"] = float(sample["rin"])
        if not math.isfinite(sample["rin"]) or not 0 <= sample["rin"] <= 10:
            raise ValueError(f"Invalid RIN for {sample['sample_id']}")
    for key in ("sample_id", "source_id", "sample_name"):
        if len({row[key] for row in samples}) != len(samples):
            raise ValueError(f"Duplicate {key} in sample sheet")
    return samples


def percentile_type7(sorted_values, p):
    """R quantile(type=7); input must already be sorted."""
    if not sorted_values:
        raise ValueError("No splice junctions with >=1 uniquely mapping read")
    if not math.isfinite(p) or not 0 <= p <= 1:
        raise ValueError("Quantile probability must be between zero and one")
    h = (len(sorted_values) - 1) * p
    lower, upper = math.floor(h), math.ceil(h)
    return sorted_values[lower] + (h - lower) * (sorted_values[upper] - sorted_values[lower])


def read_star_log(path):
    """Extract STAR's unique mapping count and percentage without rescaling."""
    keys = {"Uniquely mapped reads number", "Uniquely mapped reads %"}
    values = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if "|" in line:
                key, value = (part.strip() for part in line.split("|", 1))
                if key in keys:
                    if key in values:
                        raise ValueError(f"Duplicate STAR metric {key}: {path}")
                    values[key] = value
    try:
        count = int(values["Uniquely mapped reads number"].replace(",", ""))
        percent = float(values["Uniquely mapped reads %"].rstrip("%"))
    except (KeyError, ValueError) as exc:
        raise ValueError(f"Cannot parse uniquely mapped reads from {path}") from exc
    if count < 0 or not math.isfinite(percent) or not 0 <= percent <= 100:
        raise ValueError(f"Invalid uniquely mapped count/percentage: {path}")
    return count, percent


def read_unique_junction_support(path):
    """Read SJ.out.tab column 7; exclude zero unique-support rows before sorting."""
    counts = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected nine SJ.out.tab columns: {path}, line {line_number}")
            try:
                support = int(fields[6])
            except ValueError as exc:
                raise ValueError(f"Invalid unique junction count: {path}, line {line_number}") from exc
            if support < 0:
                raise ValueError(f"Negative unique junction count: {path}, line {line_number}")
            if support >= 1:
                counts.append(support)
    counts.sort()
    return counts


def format_quantile(number):
    """Preserve the supplied script's display precision."""
    return f"{number:.2f}".rstrip("0").rstrip(".")


def sample_paths(star_dir, sample, layout="flat", file_labels="study"):
    prefix = sample["sample_id" if file_labels == "study" else "source_id"]
    directory = star_dir / prefix if layout == "per-sample" else star_dir
    return directory / f"{prefix}_Log.final.out", directory / f"{prefix}_SJ.out.tab"


def generate_rows(star_dir, samples, layout="flat", file_labels="study"):
    rows = []
    for sample in samples:
        log_file, sj_file = sample_paths(star_dir, sample, layout, file_labels)
        for path in (log_file, sj_file):
            if not path.is_file():
                raise FileNotFoundError(f"Required file not found: {path}")
        unique_count, unique_percent = read_star_log(log_file)
        supports = read_unique_junction_support(sj_file)
        q1, median, q3 = [format_quantile(percentile_type7(supports, p)) for p in (0.25, 0.5, 0.75)]
        rows.append({
            "Sample": sample["sample_name"], "RIN": sample["rin"],
            "Unique mapped read pairs": unique_count, "Unique mapped (%)": f"{unique_percent:.2f}",
            "Supported splice junctions": len(supports), "Median junction support": median,
            "Q1": q1, "Q3": q3, "Median [Q1-Q3]": f"{median} [{q1}-{q3}]",
        })
    return rows


def write_rows(output, rows):
    """Reuse identical output; never replace a different existing result or input."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    payload = buffer.getvalue().encode("utf-8")
    if output.exists():
        if output.read_bytes() != payload:
            raise FileExistsError(f"Existing output differs: {output}. Choose a new --output path.")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("xb") as handle:
        handle.write(payload)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--star-dir", type=Path, required=True, help="Directory containing existing STAR outputs")
    parser.add_argument("--samples", type=Path, default=DEFAULT_SAMPLES, help="TSV sample metadata in output order")
    parser.add_argument("--layout", choices=("flat", "per-sample"), default="flat")
    parser.add_argument("--file-labels", choices=("study", "source"), default="study",
                        help="Study labels (KC1/Control1) or original prefixes (KC_3/N_1)")
    parser.add_argument("--output", type=Path,
                        default=Path("results/qc/Supplementary_Table2_STAR_junction_QC.tsv"))
    args = parser.parse_args()
    try:
        star_dir = args.star_dir.expanduser()
        if not star_dir.is_dir():
            raise ValueError(f"STAR directory does not exist: {star_dir}")
        samples = read_samples(args.samples.expanduser())
        rows = generate_rows(star_dir, samples, args.layout, args.file_labels)
        write_rows(args.output.expanduser(), rows)
    except (OSError, ValueError) as exc:
        parser.exit(status=1, message=f"ERROR: {exc}\n")
    print(f"Wrote STAR junction QC for {len(rows)} samples: {args.output.expanduser()}")


if __name__ == "__main__":
    main()
