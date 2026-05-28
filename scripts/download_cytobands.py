#!/usr/bin/env python3
"""Download UCSC cytoBand files for GRCh37 (hg19) and GRCh38 (hg38).

Files are saved as ``cytoBand_GRCh37.txt.gz`` and ``cytoBand_GRCh38.txt.gz``
in the destination directory (default: ``src/iscn_parser/data/``).

Usage
-----
    python scripts/download_cytobands.py
    python scripts/download_cytobands.py --dest path/to/dir
    python scripts/download_cytobands.py --builds GRCh38
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

BUILDS: dict[str, str] = {
    "GRCh37": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",
    "GRCh38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
}

# Script lives in <repo>/scripts/; data dir is relative to repo root.
_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_DEST = _REPO_ROOT / "src" / "iscn_parser" / "data"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def download(build: str, dest_dir: Path, overwrite: bool = False) -> Path:
    url = BUILDS[build]
    out = dest_dir / f"cytoBand_{build}.txt.gz"
    if out.exists() and not overwrite:
        print(f"  {out.name} already exists, skipping (use --overwrite to re-download).")
        return out
    print(f"  Downloading {build} from {url} …", end=" ", flush=True)
    urllib.request.urlretrieve(url, out)
    size_kb = out.stat().st_size // 1024
    print(f"done ({size_kb} KB, sha256={_sha256(out)[:12]}…)")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download UCSC cytoBand files for use with iscn-parser --cytoband."
    )
    parser.add_argument(
        "--dest",
        default=str(_DEFAULT_DEST),
        help=f"Destination directory (default: {_DEFAULT_DEST})",
    )
    parser.add_argument(
        "--builds",
        nargs="+",
        choices=list(BUILDS),
        default=list(BUILDS),
        metavar="BUILD",
        help="Genome build(s) to download: GRCh37, GRCh38 (default: both)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download even if the file already exists.",
    )
    args = parser.parse_args()

    dest = Path(args.dest)
    dest.mkdir(parents=True, exist_ok=True)

    print(f"Saving cytoBand files to: {dest}\n")
    for build in args.builds:
        try:
            path = download(build, dest, overwrite=args.overwrite)
            print(f"  → {path}")
        except Exception as exc:  # noqa: BLE001
            print(f"  ERROR downloading {build}: {exc}", file=sys.stderr)
            sys.exit(1)

    print("\nDone. Pass a file to iscn-parser with:")
    example_build = args.builds[0]
    example_path = dest / f"cytoBand_{example_build}.txt.gz"
    print(f"  iscn-parser --cytoband {example_path} ...")


if __name__ == "__main__":
    main()
