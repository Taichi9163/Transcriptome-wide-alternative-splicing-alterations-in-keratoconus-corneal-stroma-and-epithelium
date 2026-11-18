"""Portable work-directory selection shared by all Python command-line scripts.

Priority: --work-dir, then KC_WORK_DIR, then the current working directory.
Input and output paths are relative to that directory unless already absolute.
"""

import os
from pathlib import Path


def resolve_path(value, root=None):
    """Expand a user path and resolve relative values beneath the chosen root."""
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (Path.cwd() if root is None else Path(root)) / path
    return path.resolve()


def add_work_dir_argument(parser):
    """Add the same work-directory option to each public entry point."""
    parser.add_argument(
        "--work-dir",
        type=resolve_path,
        default=resolve_path(os.environ.get("KC_WORK_DIR") or Path.cwd()),
        help="Directory containing data/ and results/ (default: KC_WORK_DIR or current directory)",
    )
