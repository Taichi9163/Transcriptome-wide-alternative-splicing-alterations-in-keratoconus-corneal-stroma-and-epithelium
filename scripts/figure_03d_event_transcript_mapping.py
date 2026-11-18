#!/usr/bin/env python3
"""Match rMATS events to complete local exon/junction paths in Ensembl transcripts.

Coordinates are 0-based and half-open internally. A shared gene is insufficient:
the transcript must contain the inclusion or skipping structure of the event.
"""
import argparse
import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from setup_prepare_inputs import GTF_NAME, GTF_SHA256, digest


def read_table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_table(path, rows):
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return math.nan


def forms(event, exons):
    """Return the uniquely compatible local form, or [] for no/ambiguous match.

    MXE inclusion is the first genomic exon on + and the second on -, following
    rMATS. The other classes use the native rMATS inclusion definition.
    """
    exons = sorted(set(exons))
    junctions = {(left[1], right[0]) for left, right in zip(exons, exons[1:])}

    def coordinate(key):
        # R may prefix column names beginning with a number with X.
        return int(event[key] if key in event else event["X" + key])

    def complete_chain(target):
        return (target in exons
                and (coordinate("upstreamEE"), target[0]) in junctions
                and (target[1], coordinate("downstreamES")) in junctions)

    event_type = event["EventType"]
    if event_type == "SE":
        inclusion = complete_chain((coordinate("exonStart_0base"), coordinate("exonEnd")))
        skipping = (coordinate("upstreamEE"), coordinate("downstreamES")) in junctions
    elif event_type == "MXE":
        first = complete_chain((coordinate("1stExonStart_0base"), coordinate("1stExonEnd")))
        second = complete_chain((coordinate("2ndExonStart_0base"), coordinate("2ndExonEnd")))
        inclusion, skipping = (first, second) if event["strand"] == "+" else (second, first)
    elif event_type in ("A3SS", "A5SS"):
        long_exon = (coordinate("longExonStart_0base"), coordinate("longExonEnd"))
        short_exon = (coordinate("shortES"), coordinate("shortEE"))
        flank = (coordinate("flankingES"), coordinate("flankingEE"))

        def connected(alternative):
            # Match the alternative boundary and shared local sequence.
            # Variation at the distal end does not alter this splice event.
            if long_exon[1] <= flank[0]:
                return ((alternative[1], flank[0]) in junctions
                        and any(start <= alternative[0] and end == alternative[1]
                                for start, end in exons))
            return ((flank[1], alternative[0]) in junctions
                    and any(start == alternative[0] and end >= alternative[1]
                            for start, end in exons))

        inclusion, skipping = connected(long_exon), connected(short_exon)
    elif event_type == "RI":
        # A single exon must span the complete retained-intron model.
        inclusion = any(start <= coordinate("riExonStart_0base")
                        and end >= coordinate("riExonEnd") for start, end in exons)
        skipping = (coordinate("upstreamEE"), coordinate("downstreamES")) in junctions
    else:
        raise ValueError(f"Unsupported event class: {event_type}")
    if inclusion == skipping:
        return []
    return ["Inclusion"] if inclusion else ["Skipping"]


def read_annotation(path, genes):
    transcripts, all_ids = {}, set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] not in ("exon", "transcript"):
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")
            if transcript_id:
                all_ids.add(transcript_id)
            if gene_id not in genes or fields[2] != "exon":
                continue
            transcript = transcripts.setdefault(transcript_id, dict(
                gene=gene_id, chrom=fields[0].removeprefix("chr"), strand=fields[6],
                exons=[], version=attributes.get("transcript_version", "")))
            transcript["exons"].append((int(fields[3]) - 1, int(fields[4])))
    return transcripts, all_ids


def match_events(events, expression, transcripts):
    by_gene = defaultdict(list)
    for transcript_id, transcript in transcripts.items():
        by_gene[transcript["gene"]].append(transcript_id)
    pairs, audit = [], []
    for event in events:
        event_key = event["EventType"] + ":" + event["ID"]
        matched = []
        for transcript_id in by_gene[event["GeneID"]]:
            transcript = transcripts[transcript_id]
            if (transcript["chrom"] != event["chr"].removeprefix("chr")
                    or transcript["strand"] != event["strand"]):
                continue
            for form in forms(event, transcript["exons"]):
                result = expression.get(transcript_id, {})
                plotted = (math.isfinite(number(result.get("log2FoldChange")))
                           and math.isfinite(number(result.get("pvalue"))))
                row = dict(
                    EventKey=event_key, EventType=event["EventType"], EventID=event["ID"],
                    GeneID=event["GeneID"], GeneSymbol=event["geneSymbol"],
                    Chromosome=event["chr"], Strand=event["strand"],
                    TranscriptID=transcript_id, TranscriptVersion=transcript["version"],
                    Form=form, DeltaPSI=number(event["IncLevelDifference"]),
                    DAS_FDR=number(event["FDR"]), PSI_KC=event["IncLevel1"],
                    PSI_Control=event["IncLevel2"],
                    log2FoldChange=result.get("log2FoldChange", ""),
                    padj=result.get("padj", ""), pvalue=result.get("pvalue", ""),
                    baseMean=result.get("baseMean", ""),
                    DE_isoform=str(result.get("DE_isoform", "FALSE")).upper(), Plotted=plotted,
                    ExonChain_0based=";".join(f"{start}-{end}" for start, end in sorted(transcript["exons"])))
                matched.append(row)
        if any(row["Plotted"] for row in matched):
            status = "mapped_tested"
        elif matched:
            status = "mapped_no_tested_isoform"
        else:
            status = "no_compatible_annotated_transcript"
        audit.append(dict(
            EventKey=event_key, GeneID=event["GeneID"], GeneSymbol=event["geneSymbol"],
            EventType=event["EventType"], DeltaPSI=event["IncLevelDifference"],
            Inclusion_annotated=sum(row["Form"] == "Inclusion" for row in matched),
            Skipping_annotated=sum(row["Form"] == "Skipping" for row in matched),
            Inclusion_tested=sum(row["Form"] == "Inclusion" and row["Plotted"] for row in matched),
            Skipping_tested=sum(row["Form"] == "Skipping" and row["Plotted"] for row in matched),
            Status=status))
        pairs.extend(matched)
    return pairs, audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=Path.cwd())
    args = parser.parse_args()
    output = args.work_dir / "results/stroma/manuscript"
    gtf = args.work_dir / "data/reference" / GTF_NAME
    reference_hash = digest(gtf)
    if reference_hash != GTF_SHA256:
        raise ValueError("Reference differs from the release-114 GTF used for the manuscript")
    events = read_table(output / "source_DAS_events.tsv")
    expression_rows = read_table(args.work_dir / "data/stroma/expression/transcript_deseq2.tsv")
    expression = {row["TranscriptID"]: row for row in expression_rows}
    if len(expression) != len(expression_rows):
        raise ValueError("Duplicate transcript IDs in the expression table")
    genes = {event["GeneID"] for event in events}
    transcripts, all_ids = read_annotation(gtf, genes)
    pairs, audit = match_events(events, expression, transcripts)
    plotted = [row for row in pairs if row["Plotted"]]
    write_table(output / "event_transcript_mapping.tsv", pairs)
    write_table(output / "event_mapping_audit.tsv", audit)
    write_table(output / "isoform_scatter_displayed_pairs.tsv", plotted)
    summary = dict(
        source_events=len(events), source_genes=len(genes), gtf_transcripts=len(all_ids),
        tested_isoforms=len(expression), tested_isoforms_in_gtf=len(set(expression) & all_ids),
        mapping_pairs=len(pairs), plotted_pairs=len(plotted),
        plotted_events=len({row["EventKey"] for row in plotted}),
        plotted_genes=len({row["GeneID"] for row in plotted}),
        mapped_DE_pairs=sum(row["DE_isoform"] == "TRUE" for row in plotted),
        mapped_DE_genes=len({row["GeneID"] for row in plotted if row["DE_isoform"] == "TRUE"}),
        audit_status=dict(Counter(row["Status"] for row in audit)),
        form_counts=dict(Counter(row["Form"] for row in plotted)), reference_sha256=reference_hash)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
