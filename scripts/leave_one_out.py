#!/usr/bin/env python3
"""Recalculate stromal delta PSI after each single-sample omission.

The primary DAS set is held fixed. This is an effect-size sensitivity analysis;
rMATS P values and FDR are not recalculated. Labels indicate PSI vector position.
"""
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean, median


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    rows = []
    labels = [f"KC_pos{index}" for index in range(1, 6)] + [f"Control_pos{index}" for index in range(1, 5)]
    for event_type in ("SE", "RI", "A5SS", "A3SS", "MXE"):
        path = args.work_dir / "data/stroma/rmats" / f"{event_type}.MATS.JC.txt"
        with path.open() as stream:
            for event in csv.DictReader(stream, delimiter="\t"):
                if event["FDR"] == "NA" or event["IncLevelDifference"] == "NA":
                    continue
                primary = float(event["IncLevelDifference"])
                if not (float(event["FDR"]) < 0.05 and abs(primary) > 0.1):
                    continue
                cases = [float(value) for value in event["IncLevel1"].split(",")]
                controls = [float(value) for value in event["IncLevel2"].split(",")]
                if len(cases) != 5 or len(controls) != 4 or not all(map(math.isfinite, cases + controls)):
                    raise ValueError("Expected five cases and four controls with complete PSI")
                effects = [mean(cases[:i] + cases[i+1:]) - mean(controls) for i in range(5)]
                effects += [mean(cases) - mean(controls[:i] + controls[i+1:]) for i in range(4)]
                row = dict(EventType=event_type, ID=event["ID"], Gene=event["geneSymbol"], PrimaryDeltaPSI=primary)
                row.update(zip(labels, effects))
                row["SameDirectionAllOmissions"] = all(effect * primary > 0 for effect in effects)
                rows.append(row)
    if not rows:
        raise ValueError("No primary DAS events")
    summaries = [dict(OmittedSample=label, Events=len(rows),
                      SameDirection=sum(row[label] * row["PrimaryDeltaPSI"] > 0 for row in rows),
                      MedianAbsoluteChange=median(abs(row[label] - row["PrimaryDeltaPSI"]) for row in rows))
                 for label in labels]
    output = args.work_dir / "results/stromal_leave_one_out"
    output.mkdir(parents=True, exist_ok=True)
    for name, records in (("per_event.tsv", rows), ("summary.tsv", summaries)):
        with (output / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, list(records[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    overall = dict(events=len(rows), same_direction_all_omissions=sum(row["SameDirectionAllOmissions"] for row in rows))
    (output / "summary.json").write_text(json.dumps(overall, indent=2) + "\n")
    print(json.dumps(overall))


if __name__ == "__main__":
    main()
