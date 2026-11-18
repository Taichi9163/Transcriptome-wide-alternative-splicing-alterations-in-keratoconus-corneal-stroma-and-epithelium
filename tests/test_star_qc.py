"""Synthetic STAR fixtures for QC arithmetic, sample identity and output protection."""
import csv
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
from SuppTable2_STAR_junction_QC import (
    DEFAULT_SAMPLES, generate_rows, percentile_type7, read_samples,
    read_star_log, read_unique_junction_support, sample_paths, write_rows,
)


def write_fixture(directory, sample, index, layout="flat", file_labels="study"):
    log, sj = sample_paths(directory, sample, layout, file_labels)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(f"Uniquely mapped reads number | {1000 + index:,}\nUniquely mapped reads % | 50.25%\n")
    # Column 8 deliberately has large multimapper support, including at unique=0.
    sj.write_text("\n".join(f"chr1\t{100 + i * 100}\t{150 + i * 100}\t1\t1\t0\t{v}\t999\t20"
                            for i, v in enumerate([0, 10, 1, 9, 2])) + "\n\n")


class StarQCTests(unittest.TestCase):
    def test_type7_interpolation_and_singleton(self):
        values = [1, 2, 9, 10]
        self.assertEqual([percentile_type7(values, p) for p in (0, .25, .5, .75, 1)],
                         [1, 1.75, 5.5, 9.25, 10])
        self.assertEqual(percentile_type7([7], .25), 7)
        with self.assertRaises(ValueError):
            percentile_type7([], .5)
        for p in (-.1, 1.1, math.nan):
            with self.assertRaises(ValueError):
                percentile_type7(values, p)

    def test_metadata_matches_upstream_mapping_and_includes_kc6(self):
        samples = read_samples(DEFAULT_SAMPLES)
        with (REPOSITORY / "config/samples_stroma.tsv").open() as f:
            upstream = {r["sample_id"]: r["source_id"] for r in csv.DictReader(f, delimiter="\t")}
        self.assertEqual({r["sample_id"]: r["source_id"] for r in samples}, upstream)
        self.assertEqual(len(samples), 10)
        self.assertEqual([r["sample_id"] for r in samples[:4]], [f"Control{i}" for i in range(1, 5)])
        indexed = {r["sample_id"]: r for r in samples}
        self.assertEqual((indexed["KC1"]["source_id"], indexed["KC1"]["rin"]), ("KC_3", 6.3))
        self.assertEqual((indexed["KC6"]["source_id"], indexed["KC6"]["rin"]), ("KC_1", 7.6))

    def test_cli_and_source_prefix_layouts_preserve_sample_identity(self):
        samples = read_samples(DEFAULT_SAMPLES)
        with tempfile.TemporaryDirectory(prefix="synthetic STAR QC ") as temporary:
            root = Path(temporary)
            for i, sample in enumerate(samples):
                write_fixture(root / "study", sample, i, "per-sample", "study")
                write_fixture(root / "source", sample, i, "flat", "source")
            rows = generate_rows(root / "source", samples, file_labels="source")
            self.assertEqual(rows, generate_rows(root / "study", samples, layout="per-sample"))
            for i, row in enumerate(rows):
                self.assertEqual(row["Unique mapped read pairs"], 1000 + i)  # No division by two.
                self.assertEqual(row["Supported splice junctions"], 4)
                self.assertEqual(row["Median [Q1-Q3]"], "5.5 [1.75-9.25]")
            output = root / "output/qc.tsv"
            subprocess.run([sys.executable, str(REPOSITORY / "scripts/SuppTable2_STAR_junction_QC.py"),
                            "--star-dir", str(root / "source"), "--file-labels", "source",
                            "--output", str(output)], cwd=root, check=True, capture_output=True)
            with output.open() as f:
                actual = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(len(actual), 10)
            self.assertEqual(actual[4]["Sample"], "Keratoconus 1")
            self.assertEqual(actual[4]["Unique mapped read pairs"], "1004")

    def test_malformed_logs_and_junctions_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.txt"
            for text in ["", "Uniquely mapped reads number | -1\nUniquely mapped reads % | 10%\n",
                         "Uniquely mapped reads number | 1\nUniquely mapped reads % | nan%\n",
                         "Uniquely mapped reads number | 1\nUniquely mapped reads % | 101%\n",
                         "Uniquely mapped reads number | 1\nUniquely mapped reads number | 2\n"]:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_star_log(path)
            for text in ["chr1\t1\t2\t1\t1\t0\t5\n", "chr1\t1\t2\t1\t1\t0\t-1\t0\t20\n",
                         "chr1\t1\t2\t1\t1\t0\tbad\t0\t20\n"]:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_unique_junction_support(path)

    def test_missing_input_and_zero_supported_junctions_fail(self):
        sample = read_samples(DEFAULT_SAMPLES)[0]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(FileNotFoundError):
                generate_rows(root, [sample])
            write_fixture(root, sample, 0)
            _, sj = sample_paths(root, sample)
            sj.write_text("chr1\t100\t150\t1\t1\t0\t0\t500\t20\n")
            with self.assertRaisesRegex(ValueError, "No splice junctions"):
                generate_rows(root, [sample])

    def test_existing_outputs_are_preserved(self):
        sample = read_samples(DEFAULT_SAMPLES)[0]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            write_fixture(root, sample, 0)
            rows = generate_rows(root, [sample])
            output = root / "qc.tsv"
            write_rows(output, rows)
            before = output.read_bytes()
            write_rows(output, rows)
            self.assertEqual(output.read_bytes(), before)
            output.write_text("preserve existing result\n")
            with self.assertRaises(FileExistsError):
                write_rows(output, rows)
            self.assertEqual(output.read_text(), "preserve existing result\n")


if __name__ == "__main__":
    unittest.main()
