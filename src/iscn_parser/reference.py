"""Reference genome metadata: chromosome lengths, centromeres and cytobands.

Bundled TSV files provide chromosome lengths and (approximate) centromere
boundaries for GRCh37 and GRCh38, which are sufficient to resolve whole
chromosome and chromosome-arm events from conventional karyotypes. For
band-level resolution a UCSC ``cytoBand.txt`` file can be supplied at runtime.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from functools import cache
from importlib import resources

from .models import GenomeBuild


@dataclass(frozen=True)
class ChromInfo:
    """Length and centromere boundaries for a single chromosome."""

    name: str
    length: int
    centromere_start: int
    centromere_end: int


def normalize_chrom(chrom: str) -> str:
    """Normalise a chromosome label (strip ``chr`` prefix, upper-case sex)."""
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.upper() in {"X", "Y", "MT", "M"}:
        return "MT" if c.upper() in {"MT", "M"} else c.upper()
    return c


@cache
def _load_chrom_table(build: GenomeBuild) -> dict[str, ChromInfo]:
    filename = f"chromosomes_{build.value}.tsv"
    try:
        text = resources.files("iscn_parser.data").joinpath(filename).read_text()
    except (FileNotFoundError, ModuleNotFoundError):  # pragma: no cover - defensive
        return {}
    table: dict[str, ChromInfo] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        name = normalize_chrom(parts[0])
        table[name] = ChromInfo(
            name=name,
            length=int(parts[1]),
            centromere_start=int(parts[2]),
            centromere_end=int(parts[3]),
        )
    return table


def chrom_info(chrom: str, build: GenomeBuild) -> ChromInfo | None:
    """Return :class:`ChromInfo` for ``chrom`` in ``build`` (or ``None``)."""
    return _load_chrom_table(build).get(normalize_chrom(chrom))


def arm_region(chrom: str, arm: str, build: GenomeBuild) -> tuple[int, int] | None:
    """Return inclusive (start, end) coordinates of a chromosome arm.

    ``arm`` is ``"p"`` or ``"q"``. The p arm spans from base 1 to the start of
    the centromere; the q arm spans from the end of the centromere to the
    chromosome end.
    """
    info = chrom_info(chrom, build)
    if info is None:
        return None
    if arm == "p":
        if info.centromere_start <= 1:
            return None
        return (1, info.centromere_start)
    if arm == "q":
        if info.centromere_end <= 0:
            return None
        return (info.centromere_end + 1, info.length)
    return None


# --- Optional UCSC cytoBand support ----------------------------------------

_BAND_RE = re.compile(r"^(p|q)(\d+(?:\.\d+)?)$")


@dataclass
class CytobandTable:
    """A loaded UCSC ``cytoBand`` table supporting band -> coordinate lookups."""

    # chrom -> list of (band_name, start_1based, end)
    bands: dict[str, list[tuple[str, int, int]]]

    def band_region(self, chrom: str, band: str) -> tuple[int, int] | None:
        """Resolve a single band (e.g. ``p36.33``) to inclusive coordinates."""
        chrom = normalize_chrom(chrom)
        entries = self.bands.get(chrom)
        if not entries:
            return None
        band = band.strip()
        matches = [(s, e) for (name, s, e) in entries if name == band]
        if not matches:
            # Allow prefix matching (e.g. band "q21" matching "q21.1").
            matches = [(s, e) for (name, s, e) in entries if name.startswith(band)]
        if not matches:
            return None
        return (min(s for s, _ in matches), max(e for _, e in matches))

    def range_region(self, chrom: str, band1: str, band2: str) -> tuple[int, int] | None:
        """Resolve a band range (e.g. ``p36.33`` .. ``p36.31``) to coordinates."""
        r1 = self.band_region(chrom, band1)
        r2 = self.band_region(chrom, band2)
        if r1 is None or r2 is None:
            return None
        return (min(r1[0], r2[0]), max(r1[1], r2[1]))


def load_cytoband_file(path: str) -> CytobandTable:
    """Load a UCSC ``cytoBand.txt`` (optionally gzipped) file.

    The expected columns are ``chrom start end name gieStain`` where ``start``
    is 0-based half-open (UCSC convention).
    """
    opener = gzip.open if path.endswith(".gz") else open
    bands: dict[str, list[tuple[str, int, int]]] = {}
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = normalize_chrom(parts[0])
            start = int(parts[1]) + 1  # convert 0-based half-open to 1-based inclusive
            end = int(parts[2])
            name = parts[3].strip()
            bands.setdefault(chrom, []).append((name, start, end))
    return CytobandTable(bands=bands)
