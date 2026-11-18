"""Check portability and import safety without running native analyses."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
from utils_paths import add_work_dir_argument, resolve_path


class PortablePathsTests(unittest.TestCase):
    def test_explicit_work_directory_overrides_environment(self):
        with tempfile.TemporaryDirectory(prefix="kc paths ") as temporary:
            root = Path(temporary)
            with patch.dict(os.environ, {"KC_WORK_DIR": str(root / "environment")}):
                parser = argparse.ArgumentParser()
                add_work_dir_argument(parser)
                self.assertEqual(parser.parse_args([]).work_dir, (root / "environment").resolve())
                self.assertEqual(
                    parser.parse_args(["--work-dir", str(root / "explicit run")]).work_dir,
                    (root / "explicit run").resolve(),
                )

    def test_relative_absolute_and_home_paths(self):
        with tempfile.TemporaryDirectory(prefix="kc paths ") as temporary:
            root = Path(temporary).resolve()
            self.assertEqual(resolve_path("data/input.tsv", root), root / "data/input.tsv")
            self.assertEqual(
                resolve_path(root / "external.tsv", root / "other"), root / "external.tsv"
            )
            self.assertEqual(resolve_path("~/input.tsv", root), Path.home() / "input.tsv")
            with patch.dict(os.environ, {"KC_WORK_DIR": ""}):
                parser = argparse.ArgumentParser()
                add_work_dir_argument(parser)
                self.assertEqual(parser.parse_args([]).work_dir, Path.cwd().resolve())

    def test_all_python_entry_points_offer_help_without_optional_packages(self):
        with tempfile.TemporaryDirectory(prefix="kc help ") as temporary:
            for script in sorted((REPOSITORY / "scripts").glob("*.py")):
                if script.name.startswith("utils_"):
                    continue
                result = subprocess.run(
                    [sys.executable, "-B", "-S", str(script), "--help"],
                    cwd=temporary,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, (script.name, result.stderr))
                self.assertIn("--work-dir", result.stdout)
            self.assertEqual(list(Path(temporary).iterdir()), [])

    def test_imports_do_not_run_analyses_or_load_optional_packages(self):
        modules = [path.stem for path in sorted((REPOSITORY / "scripts").glob("*.py"))]
        code = (
            "import importlib, sys; "
            f"sys.path.insert(0, {str(REPOSITORY / 'scripts')!r}); "
            f"[importlib.import_module(name) for name in {modules!r}]; "
            "assert not any(name.split('.')[0] in "
            "{'pandas', 'numpy', 'matplotlib', 'nease'} for name in sys.modules)"
        )
        with tempfile.TemporaryDirectory(prefix="kc import ") as temporary:
            result = subprocess.run(
                [sys.executable, "-B", "-S", "-c", code],
                cwd=temporary,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(list(Path(temporary).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
