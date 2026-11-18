"""Reproducible NEASE analysis with event-preserving annotation exports.

Use the Python environment specified in requirements-nease.txt. Event-level exports
apply its gene/interval overlap rules before its gene-feature deduplication.
"""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import warnings

from utils_paths import add_work_dir_argument


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tissue", choices=["stroma", "epithelium"])
    add_work_dir_argument(parser)
    return parser.parse_args(argv)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_table(table, output_dir, name):
    """Write one named NEASE result without changing its column schema."""
    table.to_csv(output_dir / "results" / f"{name}.tsv", sep="\t", index=False)


def overlap(events, ref):
    joined = events.merge(ref, left_on="GeneID", right_on="Gene stable ID", how="inner")
    keep = joined[["Start_1based", "Genomic coding start"]].max(axis=1) <= joined[
        ["exonEnd", "Genomic coding end"]
    ].min(axis=1)
    return joined.loc[keep].drop_duplicates().reset_index(drop=True)


def annotation_tables(events, reference):
    """Recover event-level domain, motif and residue-interaction memberships."""
    eligible = events.loc[~events["Excluded_by_NEASE_ORF_filter"]].copy()
    domains = overlap(eligible, reference.database_mapping["Human"].drop(columns=["Unnamed: 0"]))
    domains = domains.loc[domains["NCBI gene ID"].notna()].copy()
    motifs = overlap(eligible, reference.elm["Human"])
    residues = overlap(eligible, reference.pdb["Human"]).dropna().copy()
    return domains, motifs, residues


def keyset(df, columns):
    return set(df[columns].astype(str).itertuples(index=False, name=None))


def run_analysis(args):
    """Run the primary feature mapping and Reactome analysis reported in the paper."""
    os.environ.setdefault("MPLBACKEND", "Agg")
    import numpy as np
    import pandas as pd
    import nease
    import nease.load as reference

    input_file = args.work_dir / "results" / args.tissue / "das/SE.tsv"
    raw_file = args.work_dir / "data" / args.tissue / "rmats/SE.MATS.JC.txt"
    output_dir = args.work_dir / "results" / args.tissue / "nease"
    for directory in ("results", "provenance"):
        (output_dir / directory).mkdir(parents=True, exist_ok=True)
    events = pd.read_csv(input_file, sep="\t", dtype={"ID": str, "GeneID": str})
    assert len(events) > 0 and events.ID.is_unique
    assert (events.FDR < 0.05).all() and (events.IncLevelDifference.abs() > 0.1).all()
    raw = pd.read_csv(raw_file, sep="\t", dtype={"ID": str})
    significant = raw.loc[(raw.FDR < 0.05) & (raw.IncLevelDifference.abs() > 0.1)]
    assert set(events.ID) == set(significant.ID)
    events["Event_ID"] = "SE_" + events.ID
    events["GeneID"] = events.GeneID.str.split(".").str[0]
    events["Start_1based"] = events.exonStart_0base + 1
    coord_cols = [
        "chr",
        "strand",
        "exonStart_0base",
        "exonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
    ]
    events["Event_coordinates"] = events[coord_cols].astype(str).agg(":".join, axis=1)
    assert events.Event_coordinates.is_unique
    events["Listed_in_NEASE_ORF_exclusion"] = [
        {g, int(s), int(e)} in reference.non_coding
        for g, s, e in zip(events.GeneID, events.Start_1based, events.exonEnd)
    ]
    events["Excluded_by_NEASE_ORF_filter"] = False
    original_cols = [
        "Event_ID",
        "ID",
        "GeneID",
        "geneSymbol",
        "chr",
        "strand",
        "exonStart_0base",
        "Start_1based",
        "exonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
        "FDR",
        "IncLevelDifference",
        "Event_coordinates",
        "Listed_in_NEASE_ORF_exclusion",
        "Excluded_by_NEASE_ORF_filter",
    ]
    events = events[original_cols].copy()
    standard = events[["GeneID", "Start_1based", "exonEnd", "IncLevelDifference"]].copy()
    standard.columns = ["Gene", "Start", "End", "dPSI"]
    save_table(standard, output_dir, "NEASE_standard_input")
    domains, motifs, residues = annotation_tables(events, reference)
    # Native annotations give PDB-supported exon/partner evidence, not individual residues.
    residues = residues.rename(columns={"entrezgene": "Interaction_partner_NCBI_ID"})
    residues["Amino_acid_residue_position"] = "Not provided in NEASE 1.3.2 reference"
    residues["PDB_structure_ID"] = "Not provided in NEASE 1.3.2 reference"

    config = dict(
        organism="Human",
        input_type="Standard",
        p_value_cutoff=0.05,
        remove_non_in_frame=False,
        only_divisible_by_3=False,
        min_delta=0.1,
        only_DDIs=False,
    )
    print("Running NEASE with explicit settings:", config, flush=True)
    # NEASE 1.3.2 looks up pathways by display name across all databases. Names
    # are not unique: the resulting multirow lookup raises inside a bare except,
    # reuses the previous pathway's degree/gene set, and assigns the first ID.
    # Supply one unique temporary label per reference row; preserve all actual
    # IDs, gene sets, degrees and statistical functions. Restore display names.
    original_pathways = reference.Pathways["Human"].copy(deep=True)
    working_pathways = original_pathways.copy(deep=True)
    labels = [f"NEASE_reference_row_{i}" for i in range(len(working_pathways))]
    display_names = dict(zip(labels, original_pathways["pathway"]))
    working_pathways["pathway"] = labels
    reference.Pathways["Human"] = working_pathways
    try:
        analysis = nease.run(standard.copy(), **config)
        assert hasattr(analysis, "enrichment") and len(analysis.enrichment) > 0
        analysis.enrichment["Pathway name"] = analysis.enrichment["Pathway name"].map(display_names)
        assert analysis.enrichment["Pathway name"].notna().all()
        analysis.path = original_pathways
    finally:
        reference.Pathways["Human"] = original_pathways
    api_domains, api_motifs, api_residues = (
        analysis.get_domains(),
        analysis.get_elm(),
        analysis.get_pdb(),
    )
    assert keyset(domains, ["Gene stable ID", "Pfam ID"]) == keyset(
        api_domains, ["Gene stable ID", "Pfam ID"]
    )
    assert keyset(motifs, ["Gene stable ID", "ELMIdentifier"]) == keyset(
        api_motifs, ["Gene stable ID", "ELMIdentifier"]
    )
    assert set(residues["Gene stable ID"]) == set(api_residues["Gene stable ID"])
    for name, table in [("Domains", domains), ("Motifs", motifs), ("Residue_evidence", residues)]:
        ids = set(table.Event_ID)
        events[f"Mapped_{name}"] = events.Event_ID.isin(ids)
    mapped_cols = ["Mapped_Domains", "Mapped_Motifs", "Mapped_Residue_evidence"]
    events["Mapped_any_feature"] = events[mapped_cols].any(axis=1)
    events["Mapping_status"] = np.select(
        [events.Excluded_by_NEASE_ORF_filter, events.Mapped_any_feature],
        ["Excluded by NEASE ORF/noncoding filter", "Mapped to annotated feature"],
        default="No mapped feature in NEASE reference",
    )
    events["Feature_classes"] = events.apply(
        lambda r: "; ".join(
            label
            for col, label in zip(
                mapped_cols, ["Domain", "Linear motif", "Residue interaction evidence"]
            )
            if r[col]
        ),
        axis=1,
    )

    # Trace native gene-level PDB aggregation back to every event/partner pair.
    native_pdb_pairs = {
        (r["Gene stable ID"], p)
        for _, r in api_residues.iterrows()
        for p in str(r["Co-resolved interactions"]).split(",")
    }
    assert native_pdb_pairs == keyset(residues, ["Gene stable ID", "Interaction_partner_NCBI_ID"])

    # Keep complete native outputs and pathway tests, including nonsignificant rows.
    for name, df in [
        ("SE_events", events),
        ("Event_domains", domains),
        ("Event_motifs", motifs),
        ("Event_residue_evidence", residues),
        ("Native_domains", api_domains),
        ("Native_motifs", api_motifs),
        ("Native_residue_evidence", api_residues),
        ("Affected_edges", analysis.get_edges()),
    ]:
        save_table(df, output_dir, name)
    # The manuscript reports Reactome enrichment only.
    reactome = analysis.enrich(database=["Reactome"], cutoff=0.05)
    assert reactome["Pathway ID"].is_unique
    assert set(reactome["Pathway ID"]) == set(
        original_pathways.loc[original_pathways.source == "Reactome", "external_id"]
    )
    save_table(reactome, output_dir, "Reactome_Python")

    gene_sets = [set(t["Gene stable ID"]) for t in [domains, motifs, residues]]
    summary = {
        "input_events": len(events),
        "input_genes": events.GeneID.nunique(),
        "unique_target_exons": len(events[["GeneID", "Start_1based", "exonEnd"]].drop_duplicates()),
        "excluded_events": int(events.Excluded_by_NEASE_ORF_filter.sum()),
        "mapped_events": int(events.Mapped_any_feature.sum()),
        "unmapped_eligible_events": int(
            (~events.Mapped_any_feature & ~events.Excluded_by_NEASE_ORF_filter).sum()
        ),
        "mapped_fraction_all_events": float(events.Mapped_any_feature.mean()),
        "domain_events": domains.Event_ID.nunique(),
        "motif_events": motifs.Event_ID.nunique(),
        "residue_evidence_events": residues.Event_ID.nunique(),
        "mapped_genes": len(set.union(*gene_sets)),
        "mapped_gene_fraction": len(set.union(*gene_sets)) / events.GeneID.nunique(),
        "domain_genes": len(gene_sets[0]),
        "motif_genes": len(gene_sets[1]),
        "residue_evidence_genes": len(gene_sets[2]),
        "domain_annotation_rows": len(domains),
        "motif_annotation_rows": len(motifs),
        "residue_evidence_rows": len(residues),
        "Reactome_significant": int((reactome["adj p_value"] <= 0.05).sum()),
    }
    refdir = Path(reference.__file__).parent
    manifest = {
        "python": sys.version,
        "NEASE_version": importlib.metadata.version("nease"),
        "settings": config,
        "input_file": str(input_file),
        "input_sha256": sha(input_file),
        "pathway_name_disambiguation": "Unique temporary labels per reference row during official NEASE calculation; original display names restored. Statistical functions and reference gene sets/degrees unchanged.",
        "coordinates": "GRCh38, exonStart_0base + 1 through exonEnd, both endpoints inclusive",
        "input_selection": "FDR < 0.05 and abs(IncLevelDifference) > 0.1; rMATS JC",
        "reference_sha256": {
            str(p.relative_to(refdir)): sha(p) for p in (refdir / "data").rglob("*") if p.is_file()
        },
        "source_sha256": {p.name: sha(p) for p in refdir.glob("*.py")},
        "packages": {d.metadata["Name"]: d.version for d in importlib.metadata.distributions()},
        "summary": summary,
    }
    (output_dir / "provenance/run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False)
    )
    (output_dir / "results/summary.json").write_text(json.dumps(summary, indent=2))
    (output_dir / "provenance/requirements-lock.txt").write_text(
        "\n".join(f"{k}=={v}" for k, v in sorted(manifest["packages"].items())) + "\n"
    )
    print(json.dumps(summary, indent=2), flush=True)


def main():
    args = parse_arguments()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        run_analysis(args)


if __name__ == "__main__":
    main()
