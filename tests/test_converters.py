"""Tests for the output converters."""

import os

from iscn_parser import (
    GenomeBuild,
    make_cnv_map,
    parse_iscn,
    to_bed,
    to_plink_cnv,
    to_plink_fam,
    write_plink,
)


def variants(text="arr[GRCh37] 1p36.33p36.31(849466_6015936)x1,7q11.23(72726578_74139390)x3"):
    return parse_iscn(text).variants


def test_plink_cnv_header_and_rows():
    text = to_plink_cnv(variants())
    lines = text.strip().splitlines()
    assert lines[0].split("\t") == ["FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES"]
    assert len(lines) == 3
    # TYPE column holds the integer copy number.
    assert lines[1].split("\t")[5] == "1"
    assert lines[2].split("\t")[5] == "3"


def test_plink_cnv_sex_chromosome_recoded():
    v = parse_iscn("arr[GRCh37] Xp22.31(6455151_8135644)x0").variants
    row = to_plink_cnv(v).strip().splitlines()[1].split("\t")
    assert row[2] == "23"  # X -> 23


def test_plink_fam_one_row_per_sample():
    v = parse_iscn("S1\t").variants  # empty
    assert to_plink_fam(v) == ""
    multi = (
        parse_iscn("arr[GRCh37] 1p36(1_2)x1", sample_id="A").variants
        + parse_iscn("arr[GRCh37] 2q11(3_4)x3", sample_id="A").variants
        + parse_iscn("arr[GRCh37] 3p11(5_6)x1", sample_id="B").variants
    )
    rows = to_plink_fam(multi).strip().splitlines()
    assert len(rows) == 2
    assert rows[0].split("\t")[:2] == ["A", "A"]


def test_cnv_map_breakpoints_sorted():
    text = make_cnv_map(variants())
    rows = [r.split("\t") for r in text.strip().splitlines()]
    # Four breakpoints (two per event), chromosome-sorted.
    assert len(rows) == 4
    chroms = [r[0] for r in rows]
    assert chroms == sorted(chroms, key=int)


def test_bed_is_zero_based_half_open():
    text = to_bed(variants())
    rows = [r.split("\t") for r in text.strip().splitlines()]
    # 1-based inclusive start 849466 -> 0-based 849465; end unchanged.
    assert rows[0][1] == "849465"
    assert rows[0][2] == "6015936"
    assert rows[0][0] == "chr1"
    assert "DEL" in rows[0][3]
    assert rows[0][8] == "255,0,0"  # deletion = red


def test_bed_gain_is_blue():
    text = to_bed(parse_iscn("arr[GRCh37] 7q11.23(72726578_74139390)x3").variants)
    assert text.strip().split("\t")[8] == "0,0,255"


def test_bed_without_chr_prefix():
    text = to_bed(variants(), chr_prefix=False)
    assert text.startswith("1\t")


def test_write_plink_creates_files(tmp_path):
    prefix = os.path.join(str(tmp_path), "out", "result")
    written = write_plink(variants(), prefix)
    assert any(p.endswith(".cnv") for p in written)
    assert any(p.endswith(".fam") for p in written)
    assert any(p.endswith(".cnv.map") for p in written)
    for path in written:
        assert os.path.exists(path)


def test_write_plink_no_map(tmp_path):
    prefix = os.path.join(str(tmp_path), "result")
    written = write_plink(variants(), prefix, make_map=False)
    assert not any(p.endswith(".cnv.map") for p in written)


def test_empty_inputs_produce_header_only_cnv_and_empty_bed():
    assert to_plink_cnv([]).strip() == "FID\tIID\tCHR\tBP1\tBP2\tTYPE\tSCORE\tSITES"
    assert to_bed([]) == ""
    assert make_cnv_map([]) == ""


def test_loh_bed_colour_green():
    v = parse_iscn(
        "arr[GRCh37] 6q16.1q22.31(95000000_120000000)x2 hmz", build=GenomeBuild.GRCh37
    ).variants
    assert to_bed(v).strip().split("\t")[8] == "0,170,0"
