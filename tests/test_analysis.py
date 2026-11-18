"""Scientific edge cases and command planning; no sequencing tools or images needed."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from figure_03d_event_transcript_mapping import forms, match_events
from workflow_preprocess_reads import Job, execute, plan, read_samples, REPOSITORY


class SpliceStructureTests(unittest.TestCase):
    def test_se_requires_both_flanking_junctions(self):
        event = dict(EventType="SE", upstreamEE="100", downstreamES="400", exonStart_0base="200", exonEnd="300")
        self.assertEqual(forms(event, [(0, 100), (200, 300), (400, 500)]), ["Inclusion"])
        self.assertEqual(forms(event, [(0, 100), (400, 500)]), ["Skipping"])
        self.assertEqual(forms(event, [(0, 100), (200, 310), (400, 500)]), [])
        self.assertEqual(forms(event, [(0, 100), (200, 300), (350, 375), (400, 500)]), [])

    def test_mxe_strand_and_r_column_names(self):
        event = dict(EventType="MXE", upstreamEE="100", downstreamES="400", strand="-",
                     X1stExonStart_0base="150", X1stExonEnd="200", X2ndExonStart_0base="250", X2ndExonEnd="300")
        first = [(0, 100), (150, 200), (400, 500)]
        second = [(0, 100), (250, 300), (400, 500)]
        self.assertEqual(forms(event, first), ["Skipping"])
        self.assertEqual(forms(event, second), ["Inclusion"])
        event = {key.removeprefix("X"): value for key, value in event.items()}
        event["strand"] = "+"
        self.assertEqual(forms(event, first), ["Inclusion"])
        self.assertEqual(forms(event, second), ["Skipping"])
        self.assertEqual(forms(event, [(0, 100), (150, 200), (250, 300), (400, 500)]), [])

    def test_retention_requires_full_exon_model(self):
        event = dict(EventType="RI", riExonStart_0base="0", riExonEnd="500", upstreamEE="100", downstreamES="400")
        self.assertEqual(forms(event, [(0, 500)]), ["Inclusion"])
        self.assertEqual(forms(event, [(0, 100), (400, 500)]), ["Skipping"])
        self.assertEqual(forms(event, [(50, 450)]), [])
        self.assertEqual(forms(event, [(-1, 600), (0, 100), (400, 500)]), [])

    def test_alternative_sites_accept_distal_variation(self):
        for event_type in ("A3SS", "A5SS"):
            event = dict(EventType=event_type, longExonStart_0base="0", longExonEnd="200",
                         shortES="0", shortEE="100", flankingES="400", flankingEE="500")
            self.assertEqual(forms(event, [(-50, 200), (400, 600)]), ["Inclusion"])
            self.assertEqual(forms(event, [(0, 100), (400, 500)]), ["Skipping"])
            self.assertEqual(forms(event, [(0, 200), (401, 500)]), [])
            event.update(longExonStart_0base="200", longExonEnd="500", shortES="400", shortEE="500", flankingES="0", flankingEE="100")
            self.assertEqual(forms(event, [(0, 100), (200, 600)]), ["Inclusion"])
            self.assertEqual(forms(event, [(0, 100), (400, 500)]), ["Skipping"])

    def test_matching_requires_gene_chromosome_and_strand(self):
        event = dict(EventType="SE", ID="1", GeneID="g", geneSymbol="G", chr="chr1", strand="+",
                     upstreamEE="100", downstreamES="400", exonStart_0base="200", exonEnd="300",
                     IncLevelDifference="0.2", FDR="0.001", IncLevel1="0.8", IncLevel2="0.6")
        transcript = dict(gene="g", chrom="2", strand="+", version="1", exons=[(0, 100), (400, 500)])
        pairs, audit = match_events([event], {}, {"t": transcript})
        self.assertEqual(pairs, [])
        transcript["chrom"] = "1"
        pairs, audit = match_events([event], {}, {"t": transcript})
        self.assertEqual(pairs[0]["DeltaPSI"], 0.2)  # Skipping does not reverse the plotted effect.
        self.assertEqual(audit[0]["Status"], "mapped_no_tested_isoform")
        transcript["strand"] = "-"
        self.assertEqual(match_events([event], {}, {"t": transcript})[0], [])


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.samples = read_samples(REPOSITORY / "config/samples_stroma.tsv")
        self.config = json.loads((REPOSITORY / "config/upstream.json").read_text())
        self.work = Path("/tmp/example work with spaces")

    def test_study_kc6_is_source_kc1(self):
        samples = {row["sample_id"]: row for row in self.samples}
        self.assertEqual(samples["KC6"]["source_id"], "KC_1")
        self.assertEqual(samples["KC3"]["source_id"], "KC_6")
        primary, inclusive = plan("rmats", self.work, self.config, self.samples)
        primary_list = primary.text_files[primary.output / "cases.txt"]
        self.assertIn("/KC3/KC3_", primary_list)
        self.assertNotIn("/KC6/KC6_", primary_list)
        self.assertEqual(len(primary_list.strip().split(",")), 5)
        self.assertEqual(len(inclusive.text_files[inclusive.output / "cases.txt"].strip().split(",")), 6)
        self.assertIn("--variable-read-length", primary.command)
        self.assertIn("--allow-clipping", primary.command)
        self.assertNotIn("--cstat", primary.command)

    def test_sashimi_group_order_and_safe_paths(self):
        job = plan("sashimi", self.work, self.config, self.samples, self.work / "selected.tsv")[0]
        self.assertEqual(job.text_files[job.output / "groups.gf"], "Keratoconus: 1-5\nControl: 6-9\n")
        self.assertIn(str(self.work / "selected.tsv"), job.command)
        alignments = plan("align", self.work, self.config, self.samples)
        temporary = [job.command[job.command.index("--outTmpDir") + 1] for job in alignments]
        self.assertEqual(len(set(temporary)), 10)

    def test_success_resume_failure_and_input_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.txt"
            source.write_text("original")
            job = Job("mock", root / "output", [sys.executable, "-c", "print('ok')"], [source])
            execute(job)
            self.assertTrue((job.output / "complete.json").exists())
            execute(job)
            source.write_text("changed input")
            with self.assertRaises(FileExistsError):
                execute(job)
            failed = Job("mock failure", root / "failed", [sys.executable, "-c", "raise SystemExit(2)"], [source])
            with self.assertRaises(subprocess.CalledProcessError):
                execute(failed)
            self.assertFalse((failed.output / "complete.json").exists())
            with self.assertRaises(FileExistsError):
                execute(failed)


if __name__ == "__main__":
    unittest.main()
