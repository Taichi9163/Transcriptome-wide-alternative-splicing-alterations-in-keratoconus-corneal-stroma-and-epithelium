"""Synthetic paired-end QC fixtures; no study values or images are generated."""
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
from SuppFigure1B_QC_read_pairs import (
    collect_rows, fastp_path, read_fastp, read_star_log, write_table,
)
from SuppTable2_STAR_junction_QC import DEFAULT_SAMPLES, read_samples, sample_paths


def fixture(root, sample, index=0, layout="flat", file_labels="study"):
    report = fastp_path(root / "fastp", sample, layout, file_labels)
    log, _ = sample_paths(root / "star", sample, layout, file_labels)
    report.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps({"summary": {
        "before_filtering": {"total_reads": 2000 + 2 * index},
        "after_filtering": {"total_reads": 1800 + 2 * index},
    }}))
    log.write_text(f"Number of input reads | {900 + index}\n"
                   f"Uniquely mapped reads number | {700 + index}\n"
                   "Number of reads mapped to multiple loci | 100\n"
                   "Number of reads mapped to too many loci | 20\n")
    return report, log


class ReadPairQCTests(unittest.TestCase):
    def test_pair_units_and_too_many_loci_exclusion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = read_samples(DEFAULT_SAMPLES)[0]
            fixture(root, sample)
            row = collect_rows(root / "fastp", root / "star", [sample])[0]
            self.assertEqual([row[k] for k in ("Raw read pairs", "Filtered read pairs", "Mapped read pairs")],
                             [1000, 900, 800])
            self.assertEqual(row["Uniquely mapped read pairs"], 700)

    def test_source_ids_and_per_sample_paths_agree(self):
        samples = read_samples(DEFAULT_SAMPLES)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for i, sample in enumerate(samples):
                fixture(root / "source", sample, i, file_labels="source")
                fixture(root / "study", sample, i, layout="per-sample")
            source = collect_rows(root / "source/fastp", root / "source/star", samples, file_labels="source")
            study = collect_rows(root / "study/fastp", root / "study/star", samples, layout="per-sample")
            self.assertEqual(source, study)
            self.assertEqual(source[4]["Sample"], "Keratoconus 1")
            self.assertEqual(source[4]["Raw read pairs"], 1004)
            self.assertEqual(source[9]["Sample"], "Keratoconus 6")
            self.assertEqual(source[9]["Raw read pairs"], 1009)

    def test_invalid_fastp_counts_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "report.json"
            for raw, filtered in [(1001, 800), (1000, 801), (800, 1000), (-2, 0),
                                  (True, 0), (2.5, 0), ("bad", 0)]:
                path.write_text(json.dumps({"summary": {"before_filtering": {"total_reads": raw},
                                                        "after_filtering": {"total_reads": filtered}}}))
                with self.assertRaises(ValueError):
                    read_fastp(path)

    def test_mismatched_and_duplicate_star_logs_stop(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = read_samples(DEFAULT_SAMPLES)[0]
            _, log = fixture(root, sample)
            original = log.read_text()
            log.write_text(original.replace("| 900", "| 950"))
            with self.assertRaisesRegex(ValueError, "differ"):
                collect_rows(root / "fastp", root / "star", [sample])
            log.write_text(original + "Number of input reads | 900\n")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                read_star_log(log)
            log.write_text(original.replace("| 100", "| 500"))
            with self.assertRaisesRegex(ValueError, "exceed"):
                read_star_log(log)

    def test_cli_tables_only_needs_no_plotting_packages(self):
        with tempfile.TemporaryDirectory(prefix="QC read pairs ") as temporary:
            root = Path(temporary)
            samples = read_samples(DEFAULT_SAMPLES)
            for i, sample in enumerate(samples):
                fixture(root, sample, i)
            result = subprocess.run([
                sys.executable, "-S", str(REPOSITORY / "scripts/SuppFigure1B_QC_read_pairs.py"),
                "--fastp-dir", str(root / "fastp"), "--star-dir", str(root / "star"),
                "--output-dir", str(root / "out"), "--tables-only",
            ], cwd=root, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            files = list((root / "out").iterdir())
            self.assertEqual([p.suffix for p in files], [".tsv"])
            with files[0].open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 10)

    def test_existing_output_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = read_samples(DEFAULT_SAMPLES)[0]
            fixture(root, sample)
            rows = collect_rows(root / "fastp", root / "star", [sample])
            output = root / "result.tsv"
            write_table(output, rows)
            before = output.read_bytes()
            write_table(output, rows)
            self.assertEqual(before, output.read_bytes())
            output.write_text("existing results\n")
            with self.assertRaises(FileExistsError):
                write_table(output, rows)
            self.assertEqual(output.read_text(), "existing results\n")


if __name__ == "__main__":
    unittest.main()
