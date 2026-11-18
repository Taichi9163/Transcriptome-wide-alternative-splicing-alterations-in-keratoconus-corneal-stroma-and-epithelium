#!/usr/bin/env python3
"""Expand the published numeric inputs without overwriting different local data."""
import argparse
import gzip
import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path

from utils_paths import add_work_dir_argument, resolve_path

REPOSITORY = Path(__file__).resolve().parents[1]
GTF_NAME = "Homo_sapiens.GRCh38.114.gtf.gz"
GTF_URL = "https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/" + GTF_NAME
GTF_SHA256 = "75ec5b7a1bbb8ce566622dd9abf974ad37d390130325dae6d197d1a590854c50"


def digest(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def install_input(destination, expected_hash, write):
    """Write atomically; identical existing files are reusable."""
    if destination.exists():
        if digest(destination) != expected_hash:
            raise ValueError(f"Existing file differs: {destination}. Use a new work directory.")
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(dir=destination.parent)
    os.close(handle)
    temporary = Path(temporary)
    try:
        write(temporary)
        if digest(temporary) != expected_hash:
            raise ValueError(f"Checksum mismatch for {destination.name}")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def prepare(work_dir, gtf=None, download_reference=False):
    manifest = json.loads((REPOSITORY / "inputs/manifest.json").read_text())
    for entry in manifest:

        def unpack(temporary, source=REPOSITORY / entry["file"]):
            with gzip.open(source, "rb") as incoming, temporary.open("wb") as outgoing:
                shutil.copyfileobj(incoming, outgoing)

        install_input(work_dir / entry["destination"], entry["sha256"], unpack)
    destination = work_dir / "data/reference" / GTF_NAME
    if gtf:
        install_input(destination, GTF_SHA256, lambda temporary: shutil.copyfile(gtf, temporary))
    elif download_reference:

        def download(temporary):
            with urllib.request.urlopen(GTF_URL, timeout=120) as response, temporary.open(
                "wb"
            ) as stream:
                shutil.copyfileobj(response, stream)

        install_input(destination, GTF_SHA256, download)
    print(f"Verified and prepared {len(manifest)} analysis inputs in {work_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_work_dir_argument(parser)
    reference = parser.add_mutually_exclusive_group()
    reference.add_argument("--gtf", type=Path, help="Existing release-114 .gtf.gz from Ensembl")
    reference.add_argument(
        "--download-reference", action="store_true", help="Download about 64 MB from Ensembl"
    )
    args = parser.parse_args()
    prepare(
        args.work_dir,
        resolve_path(args.gtf, args.work_dir) if args.gtf else None,
        args.download_reference,
    )


if __name__ == "__main__":
    main()
