"""Selection boundaries, safe command execution, and complete saved exon sets."""

import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
import workflow_run_rmaps3 as rmaps

CONFIG = json.loads((REPOSITORY / "config/rmaps3.json").read_text())


def row(start, fdr=0.01, delta=0.2, psi1=".5,.5", psi2=".4,.4", strand="+", flank="10"):
    return dict(
        zip(
            rmaps.COORD_COLUMNS,
            ("chr1", strand, str(start), str(start + 10), flank, "20", "500", "510"),
        ),
        FDR=str(fdr),
        IncLevelDifference=str(delta),
        IncLevel1=psi1,
        IncLevel2=psi2,
    )


class RmapsTests(unittest.TestCase):
    def test_saved_study_coordinates(self):
        for cohort in CONFIG["cohorts"]:
            with self.subTest(cohort=cohort), tempfile.TemporaryDirectory() as tmp:
                source = REPOSITORY / "inputs" / cohort / "rmats/SE.MATS.JC.txt.gz"
                target = Path(tmp) / "SE.txt"
                target.write_bytes(gzip.decompress(source.read_bytes()))
                payloads, summary = rmaps.audit(target, cohort, CONFIG)
                self.assertTrue(summary["matches_saved_coordinates"])
                self.assertEqual(
                    summary["significant_rows_at_exact_delta_cutoff"],
                    8 if cohort == "stroma" else 3,
                )
                self.assertEqual(len(payloads), 3)
                target.write_text(target.read_text() + "\n")
                with self.assertRaisesRegex(ValueError, "differs"):
                    rmaps.audit(target, cohort, CONFIG)

    def test_cutoffs_and_missing_group_means(self):
        rows = [
            row(30, delta=0.1),
            row(60, delta=-0.1),
            row(90, fdr=0.05),
            row(120, fdr=0.5),
            row(150, fdr=0.6, psi1=".85", psi2=".9"),
            row(180, fdr=0.6, psi1=".1", psi2=".15"),
            row(210, fdr=0.6, psi1="NA,NA"),
            row(240, fdr=0.6, psi1="NA,.8", psi2=".2,NA"),
            row(270, fdr="NA"),
        ]
        groups, summary = rmaps.select_exons(rows, CONFIG["settings"])
        self.assertEqual({k: len(v) for k, v in groups.items()}, {"up": 1, "dn": 1, "bg": 1})
        self.assertEqual(groups["bg"][0][2], "240")
        self.assertEqual(summary["significant_rows_at_exact_delta_cutoff"], 2)
        with self.assertRaisesRegex(ValueError, "background"):
            rmaps.select_exons([row(30)], CONFIG["settings"])

    def test_duplicates_conflicts_and_strand(self):
        rows = [
            row(30),
            row(30, flank="5"),
            row(60),
            row(60, delta=-0.2),
            row(90),
            row(90, fdr=0.8),
            row(120, delta=-0.2),
            row(120, fdr=0.8),
            row(150, fdr=0.8),
            row(30, strand="-", delta=-0.2),
        ]
        groups, summary = rmaps.select_exons(rows, CONFIG["settings"])
        self.assertEqual(summary["conflicting_target_exons"], 3)
        self.assertEqual({k: len(v) for k, v in groups.items()}, {"up": 1, "dn": 1, "bg": 1})
        self.assertEqual(groups["up"][0][4], "5")
        self.assertEqual(groups["dn"][0][1], "-")
        self.assertEqual(groups["dn"][0][4:], ("10", "20", "500", "510"))

    def test_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory(prefix="rmaps preview ") as tmp:
            work = Path(tmp) / "absent work"
            command = [
                sys.executable,
                str(REPOSITORY / "scripts/workflow_run_rmaps3.py"),
                "run",
                "--work-dir",
                str(work),
                "--rmaps-dir",
                str(Path(tmp) / "source with spaces"),
                "--fasta-root",
                str(Path(tmp) / "fasta with spaces"),
            ]
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertIn("--sigDeltaPSI 0.1", result.stdout)
            self.assertIn("RMAPS_STAT_METHOD=fisher", result.stdout)
            self.assertFalse(work.exists())
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_execute_arguments_failures_and_existing_output(self):
        with tempfile.TemporaryDirectory(prefix="rmaps mock ") as tmp:
            root = Path(tmp)
            source, fasta_root, output = root / "source", root / "genomes", root / "results"
            (source / "legacy").mkdir(parents=True)
            (source / "data").mkdir()
            (source / "legacy/motifMapSE_MP.py").write_text("# Mock only; never executed\n")
            motif = source / "data/knownMotifs.human.mouse.txt"
            motif.write_text("mock motif\n")
            (fasta_root / "hg38").mkdir(parents=True)
            (fasta_root / "hg38/hg38.fa").write_text(">chr1\nACGT\n")
            config = dict(CONFIG, motif_sha256=rmaps.sha256(motif))
            payloads = {
                g: rmaps.coordinate_bytes([tuple(row(30)[c] for c in rmaps.COORD_COLUMNS)])
                for g in ("up", "dn", "bg")
            }
            summary = {
                "coordinate_sha256": {g: hashlib.sha256(b).hexdigest() for g, b in payloads.items()}
            }
            command = rmaps.make_command(
                sys.executable, source, fasta_root, root / "SE input.txt", output, config
            )

            def mock_run(argv, **kwargs):
                if argv[-1] == "--version":
                    return subprocess.CompletedProcess(argv, 0, "Python 3.12.0\n", "")
                self.assertEqual(argv, command)
                self.assertEqual(kwargs["cwd"], source)
                self.assertEqual(kwargs["env"]["RMAPS_STAT_METHOD"], "fisher")
                self.assertTrue(kwargs["env"]["PYTHONPATH"].startswith(str(source)))
                (output / "exon").mkdir()
                for group, data in payloads.items():
                    (output / "exon" / f"{group}.coord.txt").write_bytes(data)
                return subprocess.CompletedProcess(argv, 0)

            with patch.object(rmaps.subprocess, "run", side_effect=mock_run):
                rmaps.execute(command, source, fasta_root, output, config, summary)
                self.assertEqual(
                    json.loads((output / "run_provenance.json").read_text())["status"], "completed"
                )
                with self.assertRaisesRegex(ValueError, "not empty"):
                    rmaps.execute(command, source, fasta_root, output, config, summary)
            failed = root / "failed"
            command = rmaps.make_command(
                sys.executable, source, fasta_root, root / "SE input.txt", failed, config
            )
            version = subprocess.CompletedProcess(command, 0, "Python 3.12.0\n", "")
            with patch.object(
                rmaps.subprocess,
                "run",
                side_effect=[version, subprocess.CalledProcessError(1, command)],
            ):
                with self.assertRaises(subprocess.CalledProcessError):
                    rmaps.execute(command, source, fasta_root, failed, config, summary)
            self.assertEqual(
                json.loads((failed / "run_provenance.json").read_text())["status"], "failed"
            )


if __name__ == "__main__":
    unittest.main()
