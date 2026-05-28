"""Parser for ISCN (International System for Human Cytogenomic Nomenclature).

The parser focuses on the two practically convertible classes of ISCN results:

* **Microarray ("arr") nomenclature** which carries explicit base-pair
  coordinates and copy numbers, e.g.
  ``arr[GRCh37] 1p36.33p36.31(849466_6015936)x1``. This is converted directly
  and losslessly.
* **Conventional karyotypes**, e.g. ``47,XX,+21`` or ``46,XY,del(5)(p15.2p15.33)``.
  Whole-chromosome and chromosome-arm events are resolved using bundled
  reference data; band-level events are resolved when a cytoband table is
  supplied.

The parser is deliberately tolerant: malformed components produce diagnostics
(:class:`~iscn_parser.models.ParseMessage`) rather than exceptions, so that a
single bad record never aborts a whole file.
"""

from __future__ import annotations

import re

from .models import (
    DIPLOID_COPY_NUMBER,
    CopyNumberVariant,
    GenomeBuild,
    ParseMessage,
    ParseResult,
    Severity,
)
from .reference import CytobandTable, arm_region, chrom_info, normalize_chrom

# --- low level helpers ------------------------------------------------------


def _strip_int(value: str) -> int:
    """Parse an integer that may contain thousands separators or whitespace."""
    return int(re.sub(r"[,\s]", "", value))


def split_top_level(text: str, sep: str = ",") -> list[str]:
    """Split ``text`` on ``sep`` only at bracket depth zero.

    This keeps coordinate groups such as ``(849,466_6,015,936)`` intact even
    when they contain comma thousands-separators.
    """
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for ch in text:
        if ch in "([":
            depth += 1
            current.append(ch)
        elif ch in ")]":
            depth = max(0, depth - 1)
            current.append(ch)
        elif ch == sep and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(ch)
    if current:
        parts.append("".join(current))
    return [p.strip() for p in parts if p.strip()]


_BUILD_RE = re.compile(r"[\[(]\s*(GRCh3[78]|hg(?:18|19|38)|NCBI36|b3[678])\s*[\])]", re.I)

# Array aberration with explicit coordinates, e.g. 1p36.33p36.31(849466_6015936)x1
_ARR_COORD_RE = re.compile(
    r"""
    (?P<chrom>[0-9]{1,2}|X|Y)
    (?P<bands>(?:[pq][0-9.]+){0,2})?
    \(\s*(?P<start>[\d,\s]+)\s*[-_]\s*(?P<end>[\d,\s]+)\s*\)
    \s*x\s*(?P<cn>\d+|\?)
    (?:\s*~\s*(?P<cn2>\d+))?
    (?:\s*\[(?P<frac>[0-9.]+)\])?
    (?P<rest>.*)$
    """,
    re.X | re.I,
)

# Conventional aneuploidy, e.g. +21 / -7 / +X
_ANEUPLOIDY_RE = re.compile(r"^(?P<sign>[+-])(?P<chrom>[0-9]{1,2}|X|Y)$")

# del/dup with bands, e.g. del(5)(p15.2p15.33) or dup(1)(q21q32) or del(5)(p15.2)
_DELDUP_RE = re.compile(
    r"""^(?P<kind>del|dup)\(
        (?P<chrom>[0-9]{1,2}|X|Y)\)\(
        (?P<band1>[pq][0-9.]+)(?P<band2>[pq][0-9.]+)?\)$""",
    re.X | re.I,
)

_SEX_RE = re.compile(r"^[XY]+$")


def _resolve_cn(cn: str, cn2: str | None) -> tuple[int | None, float | None]:
    """Resolve a copy-number token (and optional mosaic upper bound)."""
    if cn == "?":
        return None, None
    base = int(cn)
    if cn2 is not None:
        # Mosaic range "x1~2"; use the value furthest from diploid as the call.
        upper = int(cn2)
        far_upper = abs(upper - DIPLOID_COPY_NUMBER) > abs(base - DIPLOID_COPY_NUMBER)
        choice = upper if far_upper else base
        return choice, None
    return base, None


# --- array notation ---------------------------------------------------------


def _parse_array_component(
    component: str,
    build: GenomeBuild | None,
    sample_id: str,
    line: int | None,
) -> tuple[list[CopyNumberVariant], list[ParseMessage]]:
    variants: list[CopyNumberVariant] = []
    messages: list[ParseMessage] = []

    match = _ARR_COORD_RE.search(component)
    if not match:
        # Components without coordinates are typically normal baselines such as
        # "(1-22)x2" or "(X)x1"; report them as informational, not errors.
        if re.search(r"x\s*\d", component):
            messages.append(
                ParseMessage(
                    Severity.INFO,
                    "Skipping component without base-pair coordinates",
                    line,
                    component,
                )
            )
        else:
            messages.append(
                ParseMessage(Severity.WARNING, "Unrecognised array component", line, component)
            )
        return variants, messages

    chrom = normalize_chrom(match.group("chrom"))
    try:
        start = _strip_int(match.group("start"))
        end = _strip_int(match.group("end"))
    except ValueError:
        messages.append(
            ParseMessage(Severity.ERROR, "Invalid coordinates", line, component)
        )
        return variants, messages

    cn, _ = _resolve_cn(match.group("cn"), match.group("cn2"))
    rest = (match.group("rest") or "").strip()
    is_loh = bool(re.search(r"\bhmz\b|\bloh\b", rest, re.I))

    if cn is None:
        if is_loh:
            cn = DIPLOID_COPY_NUMBER
        else:
            messages.append(
                ParseMessage(Severity.WARNING, "Unknown copy number (x?)", line, component)
            )
            return variants, messages

    if cn == DIPLOID_COPY_NUMBER and not is_loh:
        messages.append(
            ParseMessage(
                Severity.INFO, "Skipping copy-neutral (x2) region", line, component
            )
        )
        return variants, messages

    frac = match.group("frac")
    mosaic = float(frac) if frac else None

    bands = match.group("bands") or None
    cytoband = f"{chrom}{bands}" if bands else None

    info = chrom_info(chrom, build) if build else None
    if info and (end > info.length or start < 1):
        messages.append(
            ParseMessage(
                Severity.WARNING,
                f"Coordinates outside chromosome {chrom} bounds for {build.value}",
                line,
                component,
            )
        )

    variants.append(
        CopyNumberVariant(
            chrom=chrom,
            start=start,
            end=end,
            copy_number=cn,
            sample_id=sample_id,
            build=build,
            cytoband=cytoband,
            sites=0,
            mosaic_fraction=mosaic,
            is_loh=is_loh,
            source=component,
            notation="array",
        )
    )
    return variants, messages


# --- conventional karyotype -------------------------------------------------


def _parse_sex_complement(
    complement: str, sample_id: str, build: GenomeBuild | None, line: int | None
) -> tuple[list[CopyNumberVariant], list[ParseMessage]]:
    variants: list[CopyNumberVariant] = []
    messages: list[ParseMessage] = []
    n_x = complement.count("X")
    n_y = complement.count("Y")
    # Choose the baseline complement nearest to the observed one.
    base_x, base_y = (1, 1) if n_y >= 1 else (2, 0)
    for chrom, observed, baseline in (("X", n_x, base_x), ("Y", n_y, base_y)):
        if observed == baseline:
            continue
        info = chrom_info(chrom, build) if build else None
        if info is None:
            messages.append(
                ParseMessage(
                    Severity.WARNING,
                    f"Cannot resolve chromosome {chrom} length without a known build",
                    line,
                    complement,
                )
            )
            continue
        variants.append(
            CopyNumberVariant(
                chrom=chrom,
                start=1,
                end=info.length,
                copy_number=observed,
                sample_id=sample_id,
                build=build,
                cytoband=chrom,
                source=complement,
                notation="karyotype",
            )
        )
    return variants, messages


def _parse_karyotype_component(
    component: str,
    build: GenomeBuild | None,
    sample_id: str,
    cytoband: CytobandTable | None,
    line: int | None,
) -> tuple[list[CopyNumberVariant], list[ParseMessage]]:
    variants: list[CopyNumberVariant] = []
    messages: list[ParseMessage] = []

    aneu = _ANEUPLOIDY_RE.match(component)
    if aneu:
        chrom = normalize_chrom(aneu.group("chrom"))
        info = chrom_info(chrom, build) if build else None
        if info is None:
            messages.append(
                ParseMessage(
                    Severity.WARNING,
                    f"Cannot resolve length of chromosome {chrom} (build required)",
                    line,
                    component,
                )
            )
            return variants, messages
        cn = DIPLOID_COPY_NUMBER + (1 if aneu.group("sign") == "+" else -1)
        variants.append(
            CopyNumberVariant(
                chrom=chrom,
                start=1,
                end=info.length,
                copy_number=cn,
                sample_id=sample_id,
                build=build,
                cytoband=chrom,
                source=component,
                notation="karyotype",
            )
        )
        return variants, messages

    deldup = _DELDUP_RE.match(component)
    if deldup:
        chrom = normalize_chrom(deldup.group("chrom"))
        band1 = deldup.group("band1")
        band2 = deldup.group("band2")
        cn = 1 if deldup.group("kind").lower() == "del" else 3
        region, exact = _resolve_band_region(chrom, band1, band2, build, cytoband)
        if region is None:
            messages.append(
                ParseMessage(
                    Severity.WARNING,
                    "Cannot resolve band coordinates (supply a cytoBand file for "
                    "band-level resolution)",
                    line,
                    component,
                )
            )
            return variants, messages
        if not exact:
            messages.append(
                ParseMessage(
                    Severity.WARNING,
                    "Imprecise arm-level coordinates used (supply a cytoBand file "
                    "for exact band-level resolution)",
                    line,
                    component,
                )
            )
        start, end = region
        variants.append(
            CopyNumberVariant(
                chrom=chrom,
                start=start,
                end=end,
                copy_number=cn,
                sample_id=sample_id,
                build=build,
                cytoband=f"{chrom}{band1}{band2 or ''}",
                source=component,
                notation="karyotype",
            )
        )
        return variants, messages

    # Structural rearrangements that do not change copy number, or unsupported
    # constructs: report and move on.
    if re.match(r"^(t|inv|ins)\(", component, re.I):
        messages.append(
            ParseMessage(
                Severity.INFO,
                "Balanced/structural rearrangement carries no copy-number change",
                line,
                component,
            )
        )
    else:
        messages.append(
            ParseMessage(
                Severity.WARNING, "Unsupported karyotype component", line, component
            )
        )
    return variants, messages


def _resolve_band_region(
    chrom: str,
    band1: str,
    band2: str | None,
    build: GenomeBuild | None,
    cytoband: CytobandTable | None,
) -> tuple[tuple[int, int] | None, bool]:
    """Resolve band(s) to coordinates.

    Returns ``(region, exact)`` where ``exact`` is True when band-level
    coordinates were available and False when an imprecise arm-level fallback
    was used. ``region`` is ``None`` when nothing could be resolved.
    """
    if cytoband is not None:
        if band2:
            region = cytoband.range_region(chrom, band1, band2)
        else:
            region = cytoband.band_region(chrom, band1)
        if region is not None:
            return region, True
    # Fall back to arm-level resolution when no cytoband table is available.
    arm = band1[0]
    if build is not None:
        return arm_region(chrom, arm, build), False
    return None, False


# --- public API -------------------------------------------------------------


def parse_iscn(
    text: str,
    build: GenomeBuild | None = None,
    sample_id: str = "SAMPLE",
    cytoband: CytobandTable | None = None,
    line: int | None = None,
) -> ParseResult:
    """Parse a single ISCN string into copy-number variants.

    Parameters
    ----------
    text:
        The ISCN nomenclature string.
    build:
        Default genome build. If the string embeds a build (e.g. ``[GRCh38]``)
        the embedded value takes precedence.
    sample_id:
        Sample identifier assigned to the resulting variants.
    cytoband:
        Optional cytoband table for band-level karyotype resolution.
    line:
        Optional 1-based line number used in diagnostics.
    """
    result = ParseResult()
    raw = text.strip()
    if not raw or raw.startswith("#"):
        return result

    embedded = _BUILD_RE.search(raw)
    if embedded:
        resolved = GenomeBuild.from_string(embedded.group(1))
        if resolved is not None:
            build = resolved

    is_array = bool(re.match(r"^\s*arr\b", raw, re.I)) or "arr" in raw[:6].lower()

    if is_array:
        body = _BUILD_RE.sub("", raw)
        body = re.sub(r"^\s*arr\b", "", body, flags=re.I).strip()
        components = split_top_level(body)
        if not components:
            result.messages.append(
                ParseMessage(Severity.WARNING, "Empty array result", line, raw)
            )
        for component in components:
            variants, messages = _parse_array_component(component, build, sample_id, line)
            result.variants.extend(variants)
            result.messages.extend(messages)
        return result

    # Conventional karyotype: first field is the modal chromosome number,
    # second field is the sex complement, remaining fields are aberrations.
    components = split_top_level(raw)
    if len(components) < 2:
        result.messages.append(
            ParseMessage(Severity.ERROR, "Unrecognised ISCN string", line, raw)
        )
        return result

    # components[0] is the chromosome count (and possibly ranges like 46~48); skip.
    sex = components[1]
    if _SEX_RE.match(sex):
        variants, messages = _parse_sex_complement(sex, sample_id, build, line)
        result.variants.extend(variants)
        result.messages.extend(messages)
    else:
        result.messages.append(
            ParseMessage(
                Severity.WARNING, "Unrecognised sex chromosome complement", line, sex
            )
        )

    for component in components[2:]:
        variants, messages = _parse_karyotype_component(
            component, build, sample_id, cytoband, line
        )
        result.variants.extend(variants)
        result.messages.extend(messages)

    return result


def parse_iscn_file(
    path: str,
    build: GenomeBuild | None = None,
    cytoband: CytobandTable | None = None,
    sample_column: int = 0,
    iscn_column: int = 1,
    delimiter: str = "\t",
    has_header: bool = False,
) -> ParseResult:
    """Parse a file of ISCN records.

    Each non-empty, non-comment line is either a bare ISCN string or a
    delimited record carrying a sample identifier and an ISCN string. When a
    line cannot be split into the requested columns, the whole line is treated
    as an ISCN string with an auto-generated sample id.
    """
    result = ParseResult()
    with open(path, encoding="utf-8") as handle:
        lines = handle.readlines()

    start_index = 1 if has_header else 0
    for offset, raw_line in enumerate(lines[start_index:], start=start_index + 1):
        line_text = raw_line.rstrip("\n")
        stripped = line_text.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = line_text.split(delimiter)
        if len(fields) > max(sample_column, iscn_column):
            sample_id = fields[sample_column].strip() or f"SAMPLE_{offset}"
            iscn_text = fields[iscn_column].strip()
        else:
            sample_id = f"SAMPLE_{offset}"
            iscn_text = stripped
        result.extend(
            parse_iscn(iscn_text, build=build, sample_id=sample_id, cytoband=cytoband, line=offset)
        )
    return result
