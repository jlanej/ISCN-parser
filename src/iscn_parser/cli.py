"""Command-line interface for ISCN-parser."""

from __future__ import annotations

import argparse
import sys
from typing import Optional

from . import __version__
from .converters import make_cnv_map, to_bed, to_plink_cnv, to_plink_fam, write_plink
from .models import GenomeBuild, Severity
from .parser import parse_iscn, parse_iscn_file
from .reference import load_cytoband_file


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="iscn-parser",
        description=(
            "Parse ISCN (International System for Human Cytogenomic Nomenclature) "
            "results and convert them to PLINK .cnv or BED format."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("-i", "--input", help="Input file of ISCN records (one per line).")
    src.add_argument(
        "-s", "--string", help="A single ISCN string to parse instead of a file."
    )

    parser.add_argument(
        "-f",
        "--format",
        choices=["cnv", "bed"],
        default="cnv",
        help="Output format: PLINK 'cnv' (default) or 'bed'.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output path. For --format cnv this is a prefix and writes "
            "<prefix>.cnv/.fam/.cnv.map; for bed it is a file. Defaults to stdout."
        ),
    )
    parser.add_argument(
        "-b",
        "--build",
        help="Default genome build (GRCh37/hg19, GRCh38/hg38, NCBI36/hg18).",
    )
    parser.add_argument(
        "--sample-id",
        default="SAMPLE",
        help="Sample id used when parsing a single --string (default: SAMPLE).",
    )
    parser.add_argument(
        "--cytoband",
        help="Optional UCSC cytoBand.txt(.gz) file for band-level karyotype resolution.",
    )
    parser.add_argument(
        "--sample-column",
        type=int,
        default=0,
        help="0-based column holding the sample id in --input (default: 0).",
    )
    parser.add_argument(
        "--iscn-column",
        type=int,
        default=1,
        help="0-based column holding the ISCN string in --input (default: 1).",
    )
    parser.add_argument(
        "--delimiter", default="\t", help="Field delimiter for --input (default: tab)."
    )
    parser.add_argument(
        "--has-header", action="store_true", help="Treat the first input line as a header."
    )
    parser.add_argument(
        "--no-map",
        action="store_true",
        help="Do not write the .cnv.map file when writing PLINK output to a prefix.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero if any parse errors are encountered.",
    )
    parser.add_argument(
        "--quiet", action="store_true", help="Suppress diagnostics on stderr."
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    build = None
    if args.build:
        build = GenomeBuild.from_string(args.build)
        if build is None:
            print(f"error: unknown genome build {args.build!r}", file=sys.stderr)
            return 2

    cytoband = None
    if args.cytoband:
        try:
            cytoband = load_cytoband_file(args.cytoband)
        except OSError as exc:
            print(f"error: cannot read cytoband file: {exc}", file=sys.stderr)
            return 2

    if args.string is not None:
        result = parse_iscn(
            args.string, build=build, sample_id=args.sample_id, cytoband=cytoband
        )
    else:
        try:
            result = parse_iscn_file(
                args.input,
                build=build,
                cytoband=cytoband,
                sample_column=args.sample_column,
                iscn_column=args.iscn_column,
                delimiter=args.delimiter,
                has_header=args.has_header,
            )
        except OSError as exc:
            print(f"error: cannot read input file: {exc}", file=sys.stderr)
            return 2

    if not args.quiet:
        for message in result.messages:
            if message.severity in (Severity.WARNING, Severity.ERROR):
                print(str(message), file=sys.stderr)

    if args.format == "cnv":
        if args.output:
            written = write_plink(result.variants, args.output, make_map=not args.no_map)
            if not args.quiet:
                print("wrote " + ", ".join(written), file=sys.stderr)
        else:
            sys.stdout.write(to_plink_cnv(result.variants))
    else:  # bed
        bed = to_bed(result.variants)
        if args.output:
            with open(args.output, "w", encoding="utf-8") as handle:
                handle.write(bed)
            if not args.quiet:
                print(f"wrote {args.output}", file=sys.stderr)
        else:
            sys.stdout.write(bed)

    if not args.quiet:
        print(
            f"parsed {len(result.variants)} variant(s); "
            f"{len(result.warnings)} warning(s), {len(result.errors)} error(s)",
            file=sys.stderr,
        )

    if args.strict and result.errors:
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
