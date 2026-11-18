"""Summarize read support and the effect of adding KC6, using complete JC tables."""

import argparse
import csv
import math

from utils_paths import add_work_dir_argument
from statistics import mean, median

CLASSES = ("SE", "A5SS", "A3SS", "MXE", "RI")


def number(value):
    return math.nan if value == "NA" else float(value)


def read_events(work_dir, cohort, event):
    path = work_dir / "data" / cohort / "rmats" / f"{event}.MATS.JC.txt"
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        # rMATS places its full event-specific coordinate block between strand
        # and the second ID field. IDs are local to each rMATS run.
        fields = reader.fieldnames
        coords = fields[fields.index("chr") : fields.index("IJC_SAMPLE_1") - 1]
        assert coords[:2] == ["chr", "strand"] and len(coords) >= 8
        table = {}
        for row in reader:
            key = tuple(row[c] for c in coords)
            if key in table:
                raise ValueError(f"Duplicate event coordinates in {path}: {key}")
            table[key] = row
    return table


def significant(row):
    return number(row["FDR"]) < 0.05 and abs(number(row["IncLevelDifference"])) > 0.10


def qualifying_psi(row, group, expected):
    psi, inc, skip = [
        [number(x) for x in row[col].split(",")]
        for col in (f"IncLevel{group}", f"IJC_SAMPLE_{group}", f"SJC_SAMPLE_{group}")
    ]
    if not len(psi) == len(inc) == len(skip) == expected:
        raise ValueError("Unexpected sample count in the primary stromal JC table")
    return [p for p, i, s in zip(psi, inc, skip) if math.isfinite(p) and i + s >= 10]


def write_rows(output_dir, name, rows):
    with (output_dir / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows, label):
    primary = [r for r in rows if r["primary_DAS"]]
    qualified = [r for r in primary if math.isfinite(r["qualified_delta"])]
    paired = [r for r in primary if math.isfinite(r["kc6_delta"])]

    def change(group, field):
        return median(abs(r[field] - r["primary_delta"]) for r in group) if group else math.nan

    return dict(
        event_class=label,
        primary=len(primary),
        qualified=len(qualified),
        qualified_same_sign=sum(r["qualified_delta"] * r["primary_delta"] > 0 for r in qualified),
        qualified_median_abs_change=change(qualified, "qualified_delta"),
        with_kc6=sum(r["kc6_DAS"] for r in rows),
        shared=sum(r["kc6_DAS"] for r in primary),
        primary_paired=len(paired),
        paired_same_sign=sum(r["kc6_delta"] * r["primary_delta"] > 0 for r in paired),
        paired_median_abs_change=change(paired, "kc6_delta"),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_work_dir_argument(parser)
    args = parser.parse_args()
    output_dir = args.work_dir / "results/stromal_sensitivity"
    output_dir.mkdir(parents=True, exist_ok=True)
    all_rows, summaries = [], []
    for event in CLASSES:
        primary = read_events(args.work_dir, "stroma", event)
        inclusive = read_events(args.work_dir, "stroma_with_kc6", event)
        keys = {k for k, r in primary.items() if significant(r)}
        keys.update(k for k, r in inclusive.items() if significant(r))
        rows = []
        for key in sorted(keys):
            a, b = primary.get(key), inclusive.get(key)
            kc = qualifying_psi(a, 1, 5) if a else []
            ctrl = qualifying_psi(a, 2, 4) if a else []
            rows.append(
                dict(
                    event_class=event,
                    coordinates=":".join(key),
                    gene=(a or b)["geneSymbol"],
                    primary_DAS=bool(a and significant(a)),
                    kc6_DAS=bool(b and significant(b)),
                    primary_delta=number(a["IncLevelDifference"]) if a else math.nan,
                    kc6_delta=number(b["IncLevelDifference"]) if b else math.nan,
                    primary_FDR=number(a["FDR"]) if a else math.nan,
                    kc6_FDR=number(b["FDR"]) if b else math.nan,
                    qualified_cases=len(kc),
                    qualified_controls=len(ctrl),
                    qualified_delta=(
                        mean(kc) - mean(ctrl) if len(kc) >= 3 and len(ctrl) >= 2 else math.nan
                    ),
                )
            )
        summaries.append(summarize(rows, event))
        all_rows.extend(rows)
    write_rows(output_dir, "per_event.tsv", all_rows)
    write_rows(output_dir, "summary.tsv", summaries + [summarize(all_rows, "All")])


if __name__ == "__main__":
    main()
