"""Convert resolved copy-number variants to PLINK ``.cnv`` and BED formats."""

from __future__ import annotations

import os
from collections.abc import Iterable

from .models import CopyNumberVariant

# PLINK .cnv header (PLINK 1.07 / 1.9 CNV segment file).
_CNV_HEADER = ["FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES"]


def _plink_chrom(chrom: str) -> str:
    """Map a chromosome label to the numeric code PLINK expects."""
    mapping = {"X": "23", "Y": "24", "XY": "25", "MT": "26", "M": "26"}
    return mapping.get(chrom.upper(), chrom)


def to_plink_cnv(variants: Iterable[CopyNumberVariant]) -> str:
    """Render variants as a PLINK ``.cnv`` segment file (string).

    ``TYPE`` holds the integer copy number, the convention used by PLINK
    (``0`` homozygous deletion, ``1`` deletion, ``3``/``4`` gains; ``2`` is
    diploid). LOH/copy-neutral regions are emitted with their copy number
    (typically ``2``).
    """
    lines = ["\t".join(_CNV_HEADER)]
    for v in variants:
        lines.append(
            "\t".join(
                [
                    v.sample_id,
                    v.sample_id,
                    _plink_chrom(v.chrom),
                    str(v.start),
                    str(v.end),
                    str(v.copy_number),
                    f"{v.score:g}",
                    str(v.sites),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def to_plink_fam(variants: Iterable[CopyNumberVariant]) -> str:
    """Render the companion ``.fam`` file (one row per unique sample).

    Columns are the standard PLINK six: FID, IID, PAT, MAT, SEX, PHENOTYPE.
    SEX and PHENOTYPE are emitted as unknown (``0``/``-9``).
    """
    seen: list[str] = []
    rows: list[str] = []
    for v in variants:
        if v.sample_id in seen:
            continue
        seen.append(v.sample_id)
        rows.append("\t".join([v.sample_id, v.sample_id, "0", "0", "0", "-9"]))
    return "\n".join(rows) + ("\n" if rows else "")


def make_cnv_map(variants: Iterable[CopyNumberVariant]) -> str:
    """Build a PLINK ``.cnv.map`` file enumerating event boundaries.

    The map lists every distinct (chromosome, position) breakpoint as a marker,
    mirroring the output of ``plink --cnv-make-map``. Columns: CHR, SNP, BP.
    """
    positions: set[tuple[str, int]] = set()
    for v in variants:
        positions.add((v.chrom, v.start))
        positions.add((v.chrom, v.end))

    def sort_key(item: tuple[str, int]) -> tuple[int, int]:
        chrom = _plink_chrom(item[0])
        try:
            chrom_num = int(chrom)
        except ValueError:
            chrom_num = 99
        return (chrom_num, item[1])

    rows: list[str] = []
    for chrom, pos in sorted(positions, key=sort_key):
        rows.append("\t".join([_plink_chrom(chrom), f"{chrom}-{pos}", str(pos)]))
    return "\n".join(rows) + ("\n" if rows else "")


def write_plink(
    variants: Iterable[CopyNumberVariant], prefix: str, make_map: bool = True
) -> list[str]:
    """Write ``prefix.cnv``, ``prefix.fam`` and optionally ``prefix.cnv.map``.

    Returns the list of written file paths.
    """
    variants = list(variants)
    written: list[str] = []
    out_dir = os.path.dirname(prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    cnv_path = f"{prefix}.cnv"
    with open(cnv_path, "w", encoding="utf-8") as handle:
        handle.write(to_plink_cnv(variants))
    written.append(cnv_path)

    fam_path = f"{prefix}.fam"
    with open(fam_path, "w", encoding="utf-8") as handle:
        handle.write(to_plink_fam(variants))
    written.append(fam_path)

    if make_map:
        map_path = f"{prefix}.cnv.map"
        with open(map_path, "w", encoding="utf-8") as handle:
            handle.write(make_cnv_map(variants))
        written.append(map_path)

    return written


# Colour coding used by genome browsers: red for losses, blue for gains.
_DEL_RGB = "255,0,0"
_DUP_RGB = "0,0,255"
_LOH_RGB = "0,170,0"


def to_bed(variants: Iterable[CopyNumberVariant], chr_prefix: bool = True) -> str:
    """Render variants as a BED9 file (string).

    BED is 0-based, half-open, so the inclusive 1-based start is decremented.
    The ``name`` column encodes the event type and copy number; ``itemRgb``
    colours losses red, gains blue, and copy-neutral LOH green.
    """
    rows: list[str] = []
    for v in variants:
        chrom = f"chr{v.chrom}" if chr_prefix and not str(v.chrom).startswith("chr") else v.chrom
        bed_start = max(0, v.start - 1)
        name_bits = [v.type_label, f"CN{v.copy_number}"]
        if v.cytoband:
            name_bits.append(v.cytoband)
        if v.mosaic_fraction is not None:
            name_bits.append(f"mosaic{v.mosaic_fraction:g}")
        name = "|".join(name_bits)
        if v.is_loh:
            rgb = _LOH_RGB
        elif v.copy_number < 2:
            rgb = _DEL_RGB
        else:
            rgb = _DUP_RGB
        score = min(1000, max(0, int(v.score)))
        rows.append(
            "\t".join(
                [
                    chrom,
                    str(bed_start),
                    str(v.end),
                    name,
                    str(score),
                    ".",
                    str(bed_start),
                    str(v.end),
                    rgb,
                ]
            )
        )
    return "\n".join(rows) + ("\n" if rows else "")
