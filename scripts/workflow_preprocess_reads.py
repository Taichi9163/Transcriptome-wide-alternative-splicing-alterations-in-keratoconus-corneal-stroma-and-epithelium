#!/usr/bin/env python3
"""Plan or run the supplied Ubuntu FASTQ-to-rMATS workflow.

The default is a read-only command preview. --execute runs the selected stage.
Paths in settings and sample sheets are relative to --work-dir, not the scripts.
"""
import argparse
import csv
import hashlib
import json
import re
import shlex
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
STAGES = ("raw-qc", "trim", "trimmed-qc", "star-index", "align", "rsem-index", "quantify", "rmats", "sashimi")


@dataclass
class Job:
    name: str
    output: Path
    command: list
    inputs: list
    text_files: dict = field(default_factory=dict)


def absolute(work, value):
    return (work / value).resolve()


def read_samples(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    required = {"sample_id", "group", "read1", "read2", "include_primary"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Sample sheet requires {sorted(required)}")
    names, reads = set(), set()
    for row in rows:
        name = row["sample_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or name in names:
            raise ValueError(f"Invalid or duplicate sample_id: {name}")
        if row["group"] not in ("Keratoconus", "Control"):
            raise ValueError(f"Invalid group for {name}")
        if row["include_primary"].lower() not in ("true", "false"):
            raise ValueError(f"include_primary must be true or false: {name}")
        for read in (row["read1"], row["read2"]):
            if read in reads:
                raise ValueError(f"Read file assigned more than once: {read}")
            reads.add(read)
        names.add(name)
    for group in ("Keratoconus", "Control"):
        if sum(row["group"] == group and row["include_primary"].lower() == "true" for row in rows) < 2:
            raise ValueError(f"At least two primary specimens required in {group}")
    return rows


def plan(stage, work, settings, samples, events=None):
    output = absolute(work, settings["output_dir"])
    fasta = absolute(work, settings["genome_fasta"])
    gtf = absolute(work, settings["annotation_gtf"])
    threads = str(settings["threads"])
    if int(threads) < 1:
        raise ValueError("threads must be positive")
    star_index = output / "star_index"
    rsem_index = output / "rsem_index/reference"
    jobs = []

    def add(name, destination, command, inputs, text_files=None):
        jobs.append(Job(name, destination, [str(part) for part in command], inputs, text_files or {}))

    def bam(sample, transcriptome=False):
        suffix = "Aligned.toTranscriptome.out.bam" if transcriptome else "Aligned.sortedByCoord.out.bam"
        return output / "align" / sample["sample_id"] / (sample["sample_id"] + "_" + suffix)

    if stage in ("raw-qc", "trimmed-qc"):
        for sample in samples:
            name = sample["sample_id"]
            reads = ([absolute(work, sample[key]) for key in ("read1", "read2")] if stage == "raw-qc"
                     else [output / "trim" / name / f"{name}.R{mate}.fq.gz" for mate in (1, 2)])
            destination = output / stage / sample["group"] / name
            add(f"{stage}/{name}", destination,
                ["fastqc", *reads, "-t", threads, "-o", destination], reads)
        for group in ("Keratoconus", "Control"):
            inputs = [output / stage / group / sample["sample_id"] for sample in samples if sample["group"] == group]
            destination = output / stage / f"multiqc_{group}"
            add(f"{stage}/multiqc_{group}", destination, ["multiqc", *inputs, "-o", destination], inputs)
    elif stage == "trim":
        for sample in samples:
            name = sample["sample_id"]
            reads = [absolute(work, sample[key]) for key in ("read1", "read2")]
            destination = output / "trim" / name
            # Keep fastp's supplied defaults, including its default thread count.
            add(f"trim/{name}", destination, ["fastp", "-i", reads[0], "-I", reads[1],
                "-o", destination / f"{name}.R1.fq.gz", "-O", destination / f"{name}.R2.fq.gz",
                "-h", destination / "report.html", "-j", destination / "report.json"], reads)
    elif stage == "star-index":
        add(stage, star_index, ["STAR", "--runThreadN", threads, "--runMode", "genomeGenerate",
            "--genomeDir", star_index, "--genomeFastaFiles", fasta, "--sjdbGTFfile", gtf,
            "--outFileNamePrefix", str(star_index) + "/"], [fasta, gtf])
    elif stage == "align":
        for sample in samples:
            name = sample["sample_id"]
            destination = output / "align" / name
            reads = [output / "trim" / name / f"{name}.R{mate}.fq.gz" for mate in (1, 2)]
            add(f"align/{name}", destination, ["STAR", "--runMode", "alignReads", "--genomeDir", star_index,
                "--readFilesIn", *reads, "--readFilesCommand", "gunzip", "-c",
                "--outFileNamePrefix", destination / (name + "_"), "--runThreadN", threads,
                "--outSAMtype", "BAM", "SortedByCoordinate", "--quantMode", "TranscriptomeSAM",
                "--outTmpDir", destination / "STAR_tmp"], [*reads, star_index / "Genome"])
    elif stage == "rsem-index":
        add(stage, rsem_index.parent, ["rsem-prepare-reference", "--num-threads", threads,
            "--gtf", gtf, fasta, rsem_index], [gtf, fasta])
    elif stage == "quantify":
        for sample in samples:
            name = sample["sample_id"]
            destination = output / "quantify" / name
            add(f"quantify/{name}", destination, ["rsem-calculate-expression", "--bam", "--paired-end",
                "--estimate-rspd", "--output-genome-bam", "--num-threads", threads,
                bam(sample, True), rsem_index, destination / name], [bam(sample, True), Path(str(rsem_index) + ".grp")])
    elif stage in ("rmats", "sashimi"):
        for cohort in (("primary", "with_kc6") if stage == "rmats" else ("primary",)):
            selected = [sample for sample in samples if cohort == "with_kc6" or sample["include_primary"].lower() == "true"]
            cases = [sample for sample in selected if sample["group"] == "Keratoconus"]
            controls = [sample for sample in selected if sample["group"] == "Control"]
            case_bams, control_bams = [bam(sample) for sample in cases], [bam(sample) for sample in controls]
            if any(any(character in str(path) for character in (",", "\n", "\r"))
                   for path in case_bams + control_bams):
                raise ValueError("rMATS BAM paths cannot contain commas or line breaks")
            destination = output / stage / cohort
            b1, b2 = destination / "cases.txt", destination / "controls.txt"
            text_files = {b1: ",".join(map(str, case_bams)) + "\n", b2: ",".join(map(str, control_bams)) + "\n"}
            ordered = cases + controls
            text_files[destination / "sample_order.tsv"] = "group\tindex\tsample_id\n" + "".join(
                f"{sample['group']}\t{index}\t{sample['sample_id']}\n" for index, sample in enumerate(ordered, 1))
            if stage == "rmats":
                command = [*settings["rmats_command"], "--b1", b1, "--b2", b2, "--gtf", gtf,
                    "--nthread", threads, "-t", "paired", "--readLength", "150", "--allow-clipping",
                    "--od", destination, "--tmp", destination / "tmp", "--variable-read-length"]
                inputs = [*case_bams, *control_bams, gtf]
            else:
                if events is None:
                    raise ValueError("sashimi requires --events with a selected SE rMATS table")
                group_file = destination / "groups.gf"
                text_files[group_file] = (f"Keratoconus: 1-{len(cases)}\n"
                                          f"Control: {len(cases)+1}-{len(ordered)}\n")
                command = [*settings["sashimi_command"], "--b1", b1, "--b2", b2, "--event-type", "SE",
                    "-e", events, "--l1", "Keratoconus", "--l2", "Control", "--exon_s", "1", "--intron_s", "5",
                    "-o", destination / "plots", "--group-info", group_file, "--color", "#FFA500,#000000"]
                inputs = [*case_bams, *control_bams, events]
            add(f"{stage}/{cohort}", destination, command, inputs, text_files)
    return jobs


def execute(job):
    missing = [str(path) for path in job.inputs if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing inputs: " + ", ".join(missing))
    provenance = dict(command=job.command, auxiliary_files={str(k): v for k, v in job.text_files.items()},
                      inputs={str(path): [path.stat().st_size, path.stat().st_mtime_ns] for path in job.inputs})
    signature = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    marker = job.output / "complete.json"
    if job.output.exists():
        if marker.is_file() and json.loads(marker.read_text())["signature"] == signature:
            print("Already complete:", job.name, flush=True)
            return
        raise FileExistsError(f"Output already exists or is incomplete: {job.output}. Use a new output_dir.")
    job.output.mkdir(parents=True)
    for path, contents in job.text_files.items():
        path.write_text(contents)
    (job.output / "command.json").write_text(json.dumps(provenance, indent=2) + "\n")
    with (job.output / "run.log").open("w") as log:
        subprocess.run(job.command, stdout=log, stderr=subprocess.STDOUT, check=True)
    marker.write_text(json.dumps(dict(signature=signature), indent=2) + "\n")


def export_expression_samples(work, settings, samples):
    """Write primary-cohort RSEM sheets consumed by figure_03abcd_deseq2_analysis.R."""
    output = absolute(work, settings["output_dir"])
    directory = output / "sample_sheets"
    directory.mkdir(parents=True, exist_ok=True)
    primary = [sample for sample in samples if sample["include_primary"].lower() == "true"]
    for level in ("genes", "isoforms"):
        contents = "sample\tgroup\tpath\n"
        for sample in primary:
            name = sample["sample_id"]
            result = output / "quantify" / name / f"{name}.{level}.results"
            contents += f"{name}\t{sample['group']}\t{result}\n"
        path = directory / f"{level}.tsv"
        if path.exists() and path.read_text() != contents:
            raise FileExistsError(f"Existing sample sheet differs: {path}")
        path.write_text(contents)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=(*STAGES, "all"))
    parser.add_argument("--config", type=Path, default=REPOSITORY / "config/upstream.json")
    parser.add_argument("--samples", type=Path, default=REPOSITORY / "config/samples_stroma.tsv")
    parser.add_argument("--work-dir", type=Path, default=Path.cwd())
    parser.add_argument("--events", type=Path, help="Selected SE.MATS.JC.txt for sashimi plotting")
    parser.add_argument("--execute", action="store_true", help="Run commands; default only prints the plan")
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    samples = read_samples(args.samples)
    stages = STAGES[:-1] if args.stage == "all" else (args.stage,)
    for stage in stages:
        for job in plan(stage, args.work_dir.resolve(), settings, samples,
                        args.events.resolve() if args.events else None):
            print(f"[{job.name}] {shlex.join(job.command)}", flush=True)
            if args.execute:
                execute(job)
        if stage == "quantify" and args.execute:
            export_expression_samples(args.work_dir.resolve(), settings, samples)


if __name__ == "__main__":
    main()
