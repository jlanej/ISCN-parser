"""Tests for the ISCN parsing logic."""

import os

import pytest

from iscn_parser import GenomeBuild, parse_iscn, parse_iscn_file
from iscn_parser.models import Severity
from iscn_parser.parser import split_top_level
from iscn_parser.reference import load_cytoband_file

DATA = os.path.join(os.path.dirname(__file__), "data")


def first(result):
    assert result.variants, f"expected at least one variant; messages={result.messages}"
    return result.variants[0]


# --- tokenisation -----------------------------------------------------------


def test_split_top_level_keeps_coordinate_commas():
    parts = split_top_level("1p36(849,466_6,015,936)x1,7q11(1_2)x3")
    assert parts == ["1p36(849,466_6,015,936)x1", "7q11(1_2)x3"]


# --- array notation ---------------------------------------------------------


def test_basic_deletion():
    v = first(parse_iscn("arr[GRCh37] 1p36.33p36.31(849466_6015936)x1"))
    assert (v.chrom, v.start, v.end, v.copy_number) == ("1", 849466, 6015936, 1)
    assert v.build is GenomeBuild.GRCh37
    assert v.type_label == "DEL"
    assert v.cytoband == "1p36.33p36.31"


def test_duplication_with_thousands_separators_and_build():
    v = first(parse_iscn("arr[GRCh38] 7q11.23(72,726,578_74,139,390)x3"))
    assert (v.chrom, v.start, v.end, v.copy_number) == ("7", 72726578, 74139390, 3)
    assert v.build is GenomeBuild.GRCh38
    assert v.type_label == "DUP"


def test_dash_separator_and_hg19_alias():
    v = first(parse_iscn("arr[hg19] 16p11.2(29581742-30176508)x1"))
    assert v.build is GenomeBuild.GRCh37
    assert (v.start, v.end) == (29581742, 30176508)


def test_homozygous_deletion():
    v = first(parse_iscn("arr[GRCh37] Xp22.31(6455151_8135644)x0"))
    assert v.chrom == "X"
    assert v.copy_number == 0
    assert v.type_label == "DEL"


def test_mosaic_range_picks_value_furthest_from_diploid():
    v = first(parse_iscn("arr[GRCh37] 8p23.1(8100000_11900000)x1~2"))
    assert v.copy_number == 1


def test_mosaic_fraction_recorded():
    v = first(parse_iscn("arr[GRCh37] 13q14.2(48800000_49100000)x1[0.45]"))
    assert v.mosaic_fraction == pytest.approx(0.45)
    assert v.copy_number == 1


def test_loh_region_is_copy_neutral():
    v = first(parse_iscn("arr[GRCh37] 6q16.1q22.31(95000000_120000000)x2 hmz"))
    assert v.is_loh is True
    assert v.copy_number == 2
    assert v.type_label == "LOH"


def test_multiple_aberrations_single_string():
    result = parse_iscn(
        "arr[GRCh37] 1p36.33p36.31(849466_6015936)x1,7q11.23(72726578_74139390)x3"
    )
    assert len(result.variants) == 2
    assert {v.copy_number for v in result.variants} == {1, 3}


def test_normal_array_produces_no_variants():
    result = parse_iscn("arr(1-22,X)x2")
    assert result.variants == []
    assert result.ok


def test_copy_neutral_x2_skipped_with_info():
    result = parse_iscn("arr[GRCh37] 1p36.33(849466_6015936)x2")
    assert result.variants == []
    assert any(m.severity is Severity.INFO for m in result.messages)


def test_embedded_build_overrides_default():
    v = first(parse_iscn("arr[GRCh38] 7q11.23(1_2)x3", build=GenomeBuild.GRCh37))
    assert v.build is GenomeBuild.GRCh38


# --- conventional karyotype -------------------------------------------------


def test_trisomy_21_whole_chromosome():
    v = first(parse_iscn("47,XX,+21", build=GenomeBuild.GRCh37))
    assert v.chrom == "21"
    assert v.copy_number == 3
    assert v.start == 1
    assert v.end == 249250621 or v.end == 48129895  # chr21 length GRCh37
    assert v.end == 48129895


def test_monosomy_7():
    v = first(parse_iscn("45,XY,-7", build=GenomeBuild.GRCh37))
    assert v.chrom == "7" and v.copy_number == 1


def test_turner_syndrome_loses_x():
    v = first(parse_iscn("45,X", build=GenomeBuild.GRCh37))
    assert v.chrom == "X" and v.copy_number == 1


def test_klinefelter_gains_x():
    v = first(parse_iscn("47,XXY", build=GenomeBuild.GRCh37))
    assert v.chrom == "X" and v.copy_number == 2


def test_xyy_gains_y():
    v = first(parse_iscn("47,XYY", build=GenomeBuild.GRCh37))
    assert v.chrom == "Y" and v.copy_number == 2


def test_normal_male_no_variants():
    result = parse_iscn("46,XY", build=GenomeBuild.GRCh37)
    assert result.variants == []


def test_deldup_arm_fallback_warns_without_cytoband():
    result = parse_iscn("46,XY,del(5)(p15.2p15.33)", build=GenomeBuild.GRCh37)
    assert result.variants and result.variants[0].copy_number == 1
    assert any("Imprecise" in m.message for m in result.warnings)


def test_deldup_with_cytoband_is_precise():
    table = load_cytoband_file(os.path.join(DATA, "cytoBand_test.txt"))
    result = parse_iscn(
        "46,XY,del(5)(p15.2p15.33)", build=GenomeBuild.GRCh37, cytoband=table
    )
    v = first(result)
    assert v.chrom == "5" and v.copy_number == 1
    # p15.33 starts at 0 (->1) and p15.2 ends at 18,400,000.
    assert v.start == 1
    assert v.end == 18400000
    assert not any("Imprecise" in m.message for m in result.warnings)


def test_dup_with_cytoband():
    table = load_cytoband_file(os.path.join(DATA, "cytoBand_test.txt"))
    v = first(parse_iscn("46,XX,dup(1)(q21q32)", build=GenomeBuild.GRCh37, cytoband=table))
    assert v.copy_number == 3
    assert v.start == 142600001
    assert v.end == 211500000


def test_translocation_has_no_copy_change():
    result = parse_iscn("46,XY,t(9;22)(q34;q11.2)", build=GenomeBuild.GRCh37)
    assert result.variants == []
    assert any(m.severity is Severity.INFO for m in result.messages)


# --- error handling ---------------------------------------------------------


def test_missing_coordinates_is_not_fatal():
    result = parse_iscn("arr[GRCh37] 1p36.33p36.31(849466_6015936)")
    assert result.variants == []
    # No exception, just no usable variant.


def test_invalid_coordinates_report_error_or_skip():
    result = parse_iscn("arr[GRCh37] 1p36.33(abc_def)x1")
    assert result.variants == []


def test_unknown_string_reports_error():
    result = parse_iscn("this is not iscn at all")
    assert not result.ok
    assert result.errors


def test_out_of_bounds_coordinates_warn():
    result = parse_iscn("arr[GRCh37] 21q22.3(1_999999999)x3")
    assert result.variants  # still emitted
    assert result.warnings


def test_unknown_copy_number_warns():
    result = parse_iscn("arr[GRCh37] 1p36.33(849466_6015936)x?")
    assert result.variants == []
    assert result.warnings


def test_empty_and_comment_lines_ignored():
    assert parse_iscn("").variants == []
    assert parse_iscn("# a comment").variants == []


# --- file parsing -----------------------------------------------------------


def test_parse_array_file():
    result = parse_iscn_file(os.path.join(DATA, "array_samples.tsv"))
    assert len(result.variants) >= 14
    sample_ids = {v.sample_id for v in result.variants}
    assert "ARR001" in sample_ids


def test_parse_karyotype_file():
    result = parse_iscn_file(
        os.path.join(DATA, "karyotype_samples.tsv"), build=GenomeBuild.GRCh37
    )
    assert any(v.chrom == "21" and v.copy_number == 3 for v in result.variants)


def test_parse_malformed_file_never_raises():
    result = parse_iscn_file(
        os.path.join(DATA, "malformed_samples.tsv"), build=GenomeBuild.GRCh37
    )
    # Some lines yield diagnostics; the call must complete without raising.
    assert isinstance(result.variants, list)
    assert result.messages
