"""Comprehensive tests for the bundled UCSC cytoBand files and cytoband resolution.

Covers:
- Loading the bundled ``cytoBand_GRCh37.txt.gz`` and ``cytoBand_GRCh38.txt.gz`` files.
- :class:`~iscn_parser.reference.CytobandTable` band/range lookups (exact, prefix,
  multi-band, missing, edge cases).
- Integration of the bundled cytoband tables with :func:`~iscn_parser.parser.parse_iscn`
  for conventional karyotype ``del``/``dup`` events on both GRCh37 and GRCh38.
- Comparison of band-level precision vs. arm-level fallback.
- CLI ``--cytoband`` flag wired to the bundled files.
- Malformed/edge-case inputs to the loader.
"""

from __future__ import annotations

import gzip
import os

import pytest

from iscn_parser import GenomeBuild, parse_iscn
from iscn_parser.cli import main
from iscn_parser.reference import CytobandTable, load_cytoband_file

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(__file__)
_REPO_ROOT = os.path.dirname(_HERE)
_DATA = os.path.join(_REPO_ROOT, "src", "iscn_parser", "data")

CYTOBAND_37 = os.path.join(_DATA, "cytoBand_GRCh37.txt.gz")
CYTOBAND_38 = os.path.join(_DATA, "cytoBand_GRCh38.txt.gz")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cb37() -> CytobandTable:
    return load_cytoband_file(CYTOBAND_37)


@pytest.fixture(scope="module")
def cb38() -> CytobandTable:
    return load_cytoband_file(CYTOBAND_38)


# ---------------------------------------------------------------------------
# 1. Loading bundled files
# ---------------------------------------------------------------------------


def test_bundled_grch37_file_exists():
    assert os.path.isfile(CYTOBAND_37), f"Missing bundled file: {CYTOBAND_37}"


def test_bundled_grch38_file_exists():
    assert os.path.isfile(CYTOBAND_38), f"Missing bundled file: {CYTOBAND_38}"


def test_load_grch37_returns_cytoband_table():
    table = load_cytoband_file(CYTOBAND_37)
    assert isinstance(table, CytobandTable)


def test_load_grch38_returns_cytoband_table():
    table = load_cytoband_file(CYTOBAND_38)
    assert isinstance(table, CytobandTable)


def test_grch37_band_count(cb37):
    # The bundled GRCh37 file has 862 rows (one per band entry).
    total = sum(len(v) for v in cb37.bands.values())
    assert total == 862


def test_grch38_band_count(cb38):
    # The bundled GRCh38 file is larger; 1549 entries.
    total = sum(len(v) for v in cb38.bands.values())
    assert total == 1549


def test_grch37_chromosomes_present(cb37):
    expected = {str(i) for i in range(1, 23)} | {"X", "Y"}
    assert expected.issubset(cb37.bands.keys())


def test_grch38_chromosomes_present(cb38):
    expected = {str(i) for i in range(1, 23)} | {"X", "Y"}
    assert expected.issubset(cb38.bands.keys())


def test_grch37_no_chr_prefix_in_keys(cb37):
    """Keys must be normalized (no 'chr' prefix)."""
    for key in cb37.bands:
        assert not key.lower().startswith("chr"), f"Unexpected key: {key!r}"


def test_grch38_no_chr_prefix_in_keys(cb38):
    for key in cb38.bands:
        assert not key.lower().startswith("chr"), f"Unexpected key: {key!r}"


def test_load_plain_text_cytoband(tmp_path):
    """load_cytoband_file handles uncompressed plain-text files."""
    plain = tmp_path / "cytoBand.txt"
    plain.write_text(
        "chr1\t0\t2300000\tp36.33\tgneg\n"
        "chr1\t2300000\t5400000\tp36.32\tgpos25\n"
    )
    table = load_cytoband_file(str(plain))
    assert "1" in table.bands
    assert len(table.bands["1"]) == 2


def test_load_skips_comment_and_blank_lines(tmp_path):
    content = (
        "# comment\n"
        "\n"
        "chr1\t0\t2300000\tp36.33\tgneg\n"
    )
    plain = tmp_path / "cb.txt"
    plain.write_text(content)
    table = load_cytoband_file(str(plain))
    assert len(table.bands["1"]) == 1


def test_load_skips_malformed_lines(tmp_path):
    """Lines with fewer than 4 columns are silently skipped."""
    content = "chr1\t0\tp36.33\n"  # only 3 columns
    plain = tmp_path / "bad.txt"
    plain.write_text(content)
    table = load_cytoband_file(str(plain))
    assert table.bands == {}


def test_load_gzipped_round_trip(tmp_path):
    """Confirm that a hand-crafted .gz file loads correctly."""
    raw = b"chr5\t0\t4500000\tp15.33\tgneg\nchr5\t4500000\t6300000\tp15.32\tgpos25\n"
    gz_path = tmp_path / "cb.txt.gz"
    with gzip.open(gz_path, "wb") as fh:
        fh.write(raw)
    table = load_cytoband_file(str(gz_path))
    assert "5" in table.bands
    assert len(table.bands["5"]) == 2


# ---------------------------------------------------------------------------
# 2. CytobandTable.band_region()
# ---------------------------------------------------------------------------


class TestBandRegion:
    """Tests for :meth:`CytobandTable.band_region`."""

    # --- GRCh37 exact matches ---

    def test_chr1_p36_33_grch37(self, cb37):
        # File: chr1 0 2300000 p36.33 → 1-based: start=1, end=2300000
        r = cb37.band_region("1", "p36.33")
        assert r == (1, 2300000)

    def test_chr1_p36_31_grch37(self, cb37):
        # File: chr1 5400000 7200000 p36.31 → start=5400001, end=7200000
        r = cb37.band_region("1", "p36.31")
        assert r == (5400001, 7200000)

    def test_chr5_p15_33_grch37(self, cb37):
        # File: chr5 0 4500000 p15.33 → start=1, end=4500000
        r = cb37.band_region("5", "p15.33")
        assert r == (1, 4500000)

    def test_chr5_p15_2_grch37(self, cb37):
        # File: chr5 9800000 15000000 p15.2 → start=9800001, end=15000000
        r = cb37.band_region("5", "p15.2")
        assert r == (9800001, 15000000)

    def test_chr7_q11_23_grch37(self, cb37):
        # File: chr7 72200000 77500000 q11.23 → start=72200001, end=77500000
        r = cb37.band_region("7", "q11.23")
        assert r == (72200001, 77500000)

    def test_chr1_q21_1_grch37(self, cb37):
        # File: chr1 142600000 147000000 q21.1 → start=142600001, end=147000000
        r = cb37.band_region("1", "q21.1")
        assert r == (142600001, 147000000)

    def test_chr8_p23_1_grch37(self, cb37):
        # File: chr8 6200000 12700000 p23.1 → start=6200001, end=12700000
        r = cb37.band_region("8", "p23.1")
        assert r == (6200001, 12700000)

    def test_chrx_p22_31_grch37(self, cb37):
        # File: chrX 6000000 9500000 p22.31 → start=6000001, end=9500000
        r = cb37.band_region("X", "p22.31")
        assert r == (6000001, 9500000)

    # --- chromosome name normalisation ---

    def test_chr_prefix_stripped(self, cb37):
        assert cb37.band_region("chr1", "p36.33") == cb37.band_region("1", "p36.33")

    def test_lowercase_x(self, cb37):
        assert cb37.band_region("x", "p22.31") == cb37.band_region("X", "p22.31")

    # --- prefix matching ---

    def test_prefix_match_q21_spans_sub_bands(self, cb37):
        # "q21" should match q21.1, q21.2, q21.3 → union span
        q21_1 = cb37.band_region("1", "q21.1")
        q21_3 = cb37.band_region("1", "q21.3")
        r = cb37.band_region("1", "q21")
        assert r is not None
        assert r[0] == q21_1[0]
        assert r[1] == q21_3[1]

    def test_prefix_match_p15_spans_sub_bands(self, cb37):
        r = cb37.band_region("5", "p15")
        assert r is not None
        # Should span from p15.33 (start=1) to p15.1 (end=18400000)
        assert r[0] == 1
        assert r[1] == 18400000

    def test_prefix_match_q32_spans_all_sub_bands(self, cb37):
        q32_1 = cb37.band_region("1", "q32.1")
        q32_3 = cb37.band_region("1", "q32.3")
        r = cb37.band_region("1", "q32")
        assert r[0] == q32_1[0]
        assert r[1] == q32_3[1]

    # --- missing band / chrom ---

    def test_unknown_band_returns_none(self, cb37):
        assert cb37.band_region("1", "p99.99") is None

    def test_unknown_chrom_returns_none(self, cb37):
        assert cb37.band_region("99", "p11.1") is None

    def test_empty_band_prefix_matches_all(self, cb37):
        # An empty string is a prefix of every band name, so all bands match.
        # The result spans the whole chromosome rather than returning None.
        r = cb37.band_region("1", "")
        assert r is not None
        assert r[0] == 1  # first band on chr1 starts at 1

    # --- GRCh38 exact matches ---

    def test_chr1_p36_33_grch38(self, cb38):
        # File: chr1 0 2300000 p36.33 → same start, same end in GRCh38
        r = cb38.band_region("1", "p36.33")
        assert r == (1, 2300000)

    def test_chr1_p36_31_grch38(self, cb38):
        # File: chr1 5300000 7100000 p36.31 → start=5300001, end=7100000 (differs from GRCh37)
        r = cb38.band_region("1", "p36.31")
        assert r == (5300001, 7100000)

    def test_chr5_p15_33_grch38(self, cb38):
        # File: chr5 0 4400000 p15.33 → start=1, end=4400000
        r = cb38.band_region("5", "p15.33")
        assert r == (1, 4400000)

    def test_chr5_p15_2_grch38(self, cb38):
        # File: chr5 9900000 15000000 p15.2 → start=9900001, end=15000000
        r = cb38.band_region("5", "p15.2")
        assert r == (9900001, 15000000)

    def test_chrx_p22_31_grch38(self, cb38):
        # File: chrX 6100000 9600000 p22.31 → start=6100001, end=9600000
        r = cb38.band_region("X", "p22.31")
        assert r == (6100001, 9600000)


# ---------------------------------------------------------------------------
# 3. CytobandTable.range_region()
# ---------------------------------------------------------------------------


class TestRangeRegion:
    """Tests for :meth:`CytobandTable.range_region`."""

    def test_chr5_p15_33_to_p15_2_grch37(self, cb37):
        # p15.33 (1, 4500000) ↔ p15.2 (9800001, 15000000)
        r = cb37.range_region("5", "p15.33", "p15.2")
        assert r == (1, 15000000)

    def test_chr5_p15_33_to_p15_1_grch37(self, cb37):
        # Spans all of p15
        r = cb37.range_region("5", "p15.33", "p15.1")
        assert r == (1, 18400000)

    def test_chr1_p36_33_to_p36_31_grch37(self, cb37):
        # p36.33 (1, 2300000) ↔ p36.31 (5400001, 7200000)
        r = cb37.range_region("1", "p36.33", "p36.31")
        assert r == (1, 7200000)

    def test_chr1_q21_to_q32_grch37(self, cb37):
        # Prefix-matched multi-band range across q21* and q32*
        q21_start = cb37.band_region("1", "q21")[0]
        q32_end = cb37.band_region("1", "q32")[1]
        r = cb37.range_region("1", "q21", "q32")
        assert r == (q21_start, q32_end)

    def test_range_inverted_order_is_unioned(self, cb37):
        # Providing bands in reverse cytogenetic order still returns the span.
        r_fwd = cb37.range_region("5", "p15.33", "p15.2")
        r_rev = cb37.range_region("5", "p15.2", "p15.33")
        assert r_fwd == r_rev

    def test_single_band_same_as_band_region(self, cb37):
        r_range = cb37.range_region("7", "q11.23", "q11.23")
        r_band = cb37.band_region("7", "q11.23")
        assert r_range == r_band

    def test_missing_band1_returns_none(self, cb37):
        assert cb37.range_region("1", "p99.99", "p36.31") is None

    def test_missing_band2_returns_none(self, cb37):
        assert cb37.range_region("1", "p36.33", "p99.99") is None

    def test_unknown_chrom_returns_none(self, cb37):
        assert cb37.range_region("99", "p11.1", "q11.1") is None

    def test_chr_prefix_accepted(self, cb37):
        r = cb37.range_region("chr5", "p15.33", "p15.2")
        assert r == (1, 15000000)

    def test_grch38_range(self, cb38):
        # GRCh38 chr5 p15.33(0-4400000)→(1,4400000) .. p15.2(9900000-15000000)→(9900001,15000000)
        r = cb38.range_region("5", "p15.33", "p15.2")
        assert r == (1, 15000000)

    def test_grch38_chr1_p36_range(self, cb38):
        # chr1 p36.33 (1, 2300000) .. p36.31 (5300001, 7100000)
        r = cb38.range_region("1", "p36.33", "p36.31")
        assert r == (1, 7100000)


# ---------------------------------------------------------------------------
# 4. Integration: parse_iscn() + bundled GRCh37 cytoband
# ---------------------------------------------------------------------------


class TestParseIscnWithGRCh37Cytoband:
    """Integration tests using the bundled GRCh37 cytoBand table."""

    def test_del_chr5_p15_33_p15_2_precise(self, cb37):
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        assert result.variants, result.messages
        v = result.variants[0]
        assert v.chrom == "5"
        assert v.copy_number == 1
        assert v.start == 1
        assert v.end == 15000000

    def test_del_chr5_no_imprecise_warning(self, cb37):
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        assert not any("Imprecise" in m.message for m in result.warnings)

    def test_dup_chr7_q11_23_single_band(self, cb37):
        result = parse_iscn(
            "46,XX,dup(7)(q11.23)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        assert result.variants
        v = result.variants[0]
        assert v.chrom == "7"
        assert v.copy_number == 3
        assert v.start == 72200001
        assert v.end == 77500000

    def test_del_chr1_p36_33_single_band(self, cb37):
        result = parse_iscn(
            "46,XY,del(1)(p36.33)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        v = result.variants[0]
        assert v.chrom == "1"
        assert v.copy_number == 1
        assert v.start == 1
        assert v.end == 2300000

    def test_dup_chr1_q21_q32_prefix_range(self, cb37):
        result = parse_iscn(
            "46,XX,dup(1)(q21q32)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        v = result.variants[0]
        assert v.chrom == "1"
        assert v.copy_number == 3
        expected_start = cb37.band_region("1", "q21")[0]
        expected_end = cb37.band_region("1", "q32")[1]
        assert v.start == expected_start
        assert v.end == expected_end

    def test_del_chr13_q14_2_q14_3(self, cb37):
        result = parse_iscn(
            "46,XY,del(13)(q14.2q14.3)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        v = result.variants[0]
        assert v.chrom == "13"
        assert v.copy_number == 1
        # q14.2 (47300001, 50900000) .. q14.3 (50900001, 55300000)
        assert v.start == 47300001
        assert v.end == 55300000

    def test_del_chr8_p23_1(self, cb37):
        result = parse_iscn(
            "46,XY,del(8)(p23.1)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        v = result.variants[0]
        assert v.chrom == "8"
        assert v.start == 6200001
        assert v.end == 12700000

    def test_dup_chrx_p22_31(self, cb37):
        result = parse_iscn(
            "46,XX,dup(X)(p22.31)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        v = result.variants[0]
        assert v.chrom == "X"
        assert v.copy_number == 3
        assert v.start == 6000001
        assert v.end == 9500000

    def test_cytoband_attribute_set_correctly(self, cb37):
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        assert result.variants[0].cytoband == "5p15.33p15.2"

    def test_multiple_events_with_cytoband(self, cb37):
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2),dup(7)(q11.23)",
            build=GenomeBuild.GRCh37,
            cytoband=cb37,
        )
        assert len(result.variants) == 2
        del_v = next(v for v in result.variants if v.chrom == "5")
        dup_v = next(v for v in result.variants if v.chrom == "7")
        assert del_v.copy_number == 1
        assert dup_v.copy_number == 3


# ---------------------------------------------------------------------------
# 5. Integration: parse_iscn() + bundled GRCh38 cytoband
# ---------------------------------------------------------------------------


class TestParseIscnWithGRCh38Cytoband:
    """Integration tests using the bundled GRCh38 cytoBand table."""

    def test_del_chr5_p15_33_p15_2_grch38(self, cb38):
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh38, cytoband=cb38
        )
        v = result.variants[0]
        # GRCh38: p15.33(0-4400000)→(1,4400000), p15.2(9900000-15000000)→(9900001,15000000)
        assert v.start == 1
        assert v.end == 15000000

    def test_del_chr1_p36_33_p36_31_grch38(self, cb38):
        result = parse_iscn(
            "46,XY,del(1)(p36.33p36.31)", build=GenomeBuild.GRCh38, cytoband=cb38
        )
        v = result.variants[0]
        # GRCh38: p36.33(0,2300000)→(1,2300000), p36.31(5300000,7100000)→(5300001,7100000)
        assert v.start == 1
        assert v.end == 7100000

    def test_dup_chrx_p22_31_grch38(self, cb38):
        result = parse_iscn(
            "46,XX,dup(X)(p22.31)", build=GenomeBuild.GRCh38, cytoband=cb38
        )
        v = result.variants[0]
        # GRCh38: chrX p22.31 (6100000-9600000) → (6100001, 9600000)
        assert v.start == 6100001
        assert v.end == 9600000

    def test_no_imprecise_warning_grch38(self, cb38):
        result = parse_iscn(
            "46,XX,dup(X)(p22.31)", build=GenomeBuild.GRCh38, cytoband=cb38
        )
        assert not any("Imprecise" in m.message for m in result.warnings)


# ---------------------------------------------------------------------------
# 6. Cytoband vs. arm-level fallback comparison
# ---------------------------------------------------------------------------


class TestCytobandVsArmFallback:
    """Confirm that cytoband resolution is more precise than the arm fallback."""

    def test_band_level_more_precise_than_arm(self, cb37):
        # Without cytoband: del(5)(p15.33p15.2) falls back to whole p-arm.
        arm_result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37
        )
        band_result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        arm_v = arm_result.variants[0]
        band_v = band_result.variants[0]
        # Band-level end should be less than the whole p-arm end.
        assert band_v.end < arm_v.end
        # Band-level start is still 1 (p arm starts at telomere).
        assert band_v.start == 1

    def test_arm_fallback_warns_imprecise(self):
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37
        )
        assert any("Imprecise" in m.message for m in result.warnings)

    def test_band_level_suppresses_imprecise_warning(self, cb37):
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        assert not any("Imprecise" in m.message for m in result.warnings)

    def test_grch37_and_grch38_give_different_coordinates(self, cb37, cb38):
        r37 = parse_iscn(
            "46,XY,del(1)(p36.33p36.31)", build=GenomeBuild.GRCh37, cytoband=cb37
        )
        r38 = parse_iscn(
            "46,XY,del(1)(p36.33p36.31)", build=GenomeBuild.GRCh38, cytoband=cb38
        )
        v37 = r37.variants[0]
        v38 = r38.variants[0]
        # The end differs between builds (GRCh37: 7200000, GRCh38: 7100000).
        assert v37.end != v38.end


# ---------------------------------------------------------------------------
# 7. CLI --cytoband with bundled files
# ---------------------------------------------------------------------------


class TestCLIWithBundledCytobands:
    """CLI integration tests using the bundled cytoBand files."""

    def test_cli_grch37_cytoband_precise_del(self, capsys):
        rc = main(
            [
                "-s",
                "46,XY,del(5)(p15.33p15.2)",
                "-b",
                "GRCh37",
                "--cytoband",
                CYTOBAND_37,
                "--quiet",
            ]
        )
        out = capsys.readouterr().out
        assert rc == 0
        # Precise end coordinate from the bundled file.
        assert "15000000" in out

    def test_cli_grch38_cytoband_dup_x(self, capsys):
        rc = main(
            [
                "-s",
                "46,XX,dup(X)(p22.31)",
                "-b",
                "GRCh38",
                "--cytoband",
                CYTOBAND_38,
                "--quiet",
            ]
        )
        out = capsys.readouterr().out
        assert rc == 0
        # GRCh38 p22.31 start → 6100001
        assert "6100001" in out

    def test_cli_cytoband_suppresses_imprecise_warning(self, capsys):
        main(
            [
                "-s",
                "46,XY,del(5)(p15.33p15.2)",
                "-b",
                "GRCh37",
                "--cytoband",
                CYTOBAND_37,
            ]
        )
        err = capsys.readouterr().err
        assert "Imprecise" not in err

    def test_cli_without_cytoband_emits_imprecise_warning(self, capsys):
        main(["-s", "46,XY,del(5)(p15.33p15.2)", "-b", "GRCh37"])
        err = capsys.readouterr().err
        assert "Imprecise" in err

    def test_cli_cytoband_bed_output(self, capsys):
        rc = main(
            [
                "-s",
                "46,XY,del(5)(p15.33p15.2)",
                "-b",
                "GRCh37",
                "--cytoband",
                CYTOBAND_37,
                "-f",
                "bed",
                "--quiet",
            ]
        )
        out = capsys.readouterr().out
        assert rc == 0
        assert out.startswith("chr5\t")
        # BED is 0-based: start 1 → 0
        assert "\t0\t15000000\t" in out

    def test_cli_cytoband_file_to_plink(self, tmp_path, capsys):
        prefix = os.path.join(str(tmp_path), "out")
        rc = main(
            [
                "-s",
                "46,XY,del(5)(p15.33p15.2)",
                "-b",
                "GRCh37",
                "--cytoband",
                CYTOBAND_37,
                "-o",
                prefix,
                "--quiet",
            ]
        )
        assert rc == 0
        assert os.path.exists(prefix + ".cnv")

    def test_cli_invalid_cytoband_path_returns_nonzero(self, capsys):
        rc = main(
            [
                "-s",
                "46,XY,del(5)(p15.33p15.2)",
                "-b",
                "GRCh37",
                "--cytoband",
                "/nonexistent/cytoBand.txt",
                "--quiet",
            ]
        )
        assert rc != 0


# ---------------------------------------------------------------------------
# 8. Edge cases and boundary conditions
# ---------------------------------------------------------------------------


class TestEdgeCases:
    """Edge cases for the cytoband loader and lookup."""

    def test_band_at_chromosome_start_has_start_1(self, cb37):
        """A band that starts at position 0 in the file must map to 1-based 1."""
        r = cb37.band_region("1", "p36.33")
        assert r[0] == 1

    def test_band_end_matches_file_value_exactly(self, cb37):
        """The end coordinate is the raw value from the file (no +1)."""
        r = cb37.band_region("1", "p36.33")
        # File: chr1 0 2300000 → end should be exactly 2300000.
        assert r[1] == 2300000

    def test_range_where_both_bands_are_identical(self, cb37):
        r1 = cb37.band_region("5", "p15.33")
        r2 = cb37.range_region("5", "p15.33", "p15.33")
        assert r1 == r2

    def test_cytoband_table_empty_chromosome_list(self):
        table = CytobandTable(bands={})
        assert table.band_region("1", "p36.33") is None
        assert table.range_region("1", "p36.33", "p36.31") is None

    def test_loading_empty_file(self, tmp_path):
        empty = tmp_path / "empty.txt"
        empty.write_text("")
        table = load_cytoband_file(str(empty))
        assert table.bands == {}

    def test_loading_comments_only_file(self, tmp_path):
        f = tmp_path / "comments.txt"
        f.write_text("# header\n# another comment\n\n")
        table = load_cytoband_file(str(f))
        assert table.bands == {}

    def test_loading_file_without_stain_column(self, tmp_path):
        """Files with exactly 4 columns (no gieStain) should load correctly."""
        f = tmp_path / "no_stain.txt"
        f.write_text("chr1\t0\t2300000\tp36.33\n")
        table = load_cytoband_file(str(f))
        assert table.band_region("1", "p36.33") == (1, 2300000)

    def test_band_region_whitespace_stripped(self, cb37):
        """Leading/trailing whitespace in band name should be handled."""
        r = cb37.band_region("1", "  p36.33  ")
        assert r == (1, 2300000)

    def test_grch37_has_more_bands_than_test_fixture(self, cb37):
        """Bundled file should cover far more bands than the minimal test fixture."""
        from iscn_parser.reference import load_cytoband_file as lcf

        test_fixture = os.path.join(_HERE, "data", "cytoBand_test.txt")
        small = lcf(test_fixture)
        small_total = sum(len(v) for v in small.bands.values())
        full_total = sum(len(v) for v in cb37.bands.values())
        assert full_total > small_total * 10

    def test_parse_iscn_with_none_cytoband_uses_arm_fallback(self):
        """Passing cytoband=None (default) must not raise."""
        result = parse_iscn(
            "46,XY,del(5)(p15.33p15.2)", build=GenomeBuild.GRCh37, cytoband=None
        )
        # Falls back to arm-level; a variant should still be emitted.
        assert result.variants

    def test_bundled_grch37_chr5_p15_33_is_first_band(self, cb37):
        """The very first band on chr5 should start at position 1."""
        r = cb37.band_region("5", "p15.33")
        assert r[0] == 1

    def test_bundled_grch37_chrm_absent(self, cb37):
        """chrMT/chrM is not present in the standard GRCh37 cytoBand table."""
        assert cb37.bands.get("MT") is None

    def test_bundled_grch38_chrm_present(self, cb38):
        """GRCh38 cytoBand includes chrM."""
        assert cb38.bands.get("MT") is not None
