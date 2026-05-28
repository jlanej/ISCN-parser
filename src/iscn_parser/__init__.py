"""ISCN-parser: parse ISCN cytogenomic nomenclature and convert to PLINK CNV / BED.

The public API re-exports the most commonly used helpers so that callers can do::

    from iscn_parser import parse_iscn, to_plink_cnv, to_bed

See :mod:`iscn_parser.parser` and :mod:`iscn_parser.converters` for details.
"""

from .converters import make_cnv_map, to_bed, to_plink_cnv, to_plink_fam, write_plink
from .models import CopyNumberVariant, GenomeBuild, ParseMessage, ParseResult
from .parser import parse_iscn, parse_iscn_file

__all__ = [
    "CopyNumberVariant",
    "GenomeBuild",
    "ParseMessage",
    "ParseResult",
    "parse_iscn",
    "parse_iscn_file",
    "to_plink_cnv",
    "to_plink_fam",
    "to_bed",
    "make_cnv_map",
    "write_plink",
]

__version__ = "0.1.0"
