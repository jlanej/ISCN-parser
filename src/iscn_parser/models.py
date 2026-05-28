"""Data models shared across the ISCN parser and converters."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class GenomeBuild(str, Enum):
    """Supported human reference genome builds, with common aliases."""

    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    NCBI36 = "NCBI36"

    @classmethod
    def from_string(cls, value: str | None) -> GenomeBuild | None:
        """Resolve a build from a free-text token (e.g. ``hg19``), or ``None``."""
        if value is None:
            return None
        token = value.strip().lower().lstrip("[(").rstrip("])")
        aliases = {
            "grch37": cls.GRCh37,
            "hg19": cls.GRCh37,
            "b37": cls.GRCh37,
            "grch38": cls.GRCh38,
            "hg38": cls.GRCh38,
            "b38": cls.GRCh38,
            "ncbi36": cls.NCBI36,
            "hg18": cls.NCBI36,
            "b36": cls.NCBI36,
        }
        return aliases.get(token)


class Severity(str, Enum):
    """Severity levels for parser diagnostics."""

    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


@dataclass
class ParseMessage:
    """A single diagnostic emitted while parsing an ISCN string."""

    severity: Severity
    message: str
    line: int | None = None
    text: str | None = None

    def __str__(self) -> str:  # pragma: no cover - trivial formatting
        loc = f" (line {self.line})" if self.line is not None else ""
        snippet = f": {self.text!r}" if self.text else ""
        return f"[{self.severity.value}]{loc} {self.message}{snippet}"


# PLINK CNV convention: TYPE is an integer copy number. 2 == diploid (normal).
DIPLOID_COPY_NUMBER = 2


@dataclass
class CopyNumberVariant:
    """A single copy-number event resolved to genomic coordinates.

    Coordinates are 1-based, inclusive (the convention used by ISCN and by the
    PLINK ``.cnv`` format). BED output converts these to 0-based half-open.
    """

    chrom: str
    start: int
    end: int
    copy_number: int
    sample_id: str = "SAMPLE"
    build: GenomeBuild | None = None
    cytoband: str | None = None
    score: float = 0.0
    sites: int = 0
    mosaic_fraction: float | None = None
    is_loh: bool = False
    source: str | None = None
    notation: str = "array"

    @property
    def length(self) -> int:
        """Length of the event in base pairs (inclusive coordinates)."""
        return self.end - self.start + 1

    @property
    def type_label(self) -> str:
        """Human-readable DEL/DUP/LOH/NORMAL label."""
        if self.is_loh:
            return "LOH"
        if self.copy_number < DIPLOID_COPY_NUMBER:
            return "DEL"
        if self.copy_number > DIPLOID_COPY_NUMBER:
            return "DUP"
        return "NORMAL"

    def __post_init__(self) -> None:
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        if self.start < 1:
            self.start = 1


@dataclass
class ParseResult:
    """The outcome of parsing one or more ISCN strings."""

    variants: list[CopyNumberVariant] = field(default_factory=list)
    messages: list[ParseMessage] = field(default_factory=list)

    @property
    def errors(self) -> list[ParseMessage]:
        return [m for m in self.messages if m.severity is Severity.ERROR]

    @property
    def warnings(self) -> list[ParseMessage]:
        return [m for m in self.messages if m.severity is Severity.WARNING]

    @property
    def ok(self) -> bool:
        """True when no ERROR-level messages were produced."""
        return not self.errors

    def extend(self, other: ParseResult) -> None:
        self.variants.extend(other.variants)
        self.messages.extend(other.messages)
