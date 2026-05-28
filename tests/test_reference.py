"""Tests for reference data and build resolution."""

import pytest

from iscn_parser import GenomeBuild
from iscn_parser.reference import arm_region, chrom_info, normalize_chrom


@pytest.mark.parametrize(
    "value,expected",
    [
        ("GRCh37", GenomeBuild.GRCh37),
        ("hg19", GenomeBuild.GRCh37),
        ("GRCh38", GenomeBuild.GRCh38),
        ("hg38", GenomeBuild.GRCh38),
        ("[GRCh37]", GenomeBuild.GRCh37),
        ("hg18", GenomeBuild.NCBI36),
        ("nonsense", None),
        (None, None),
    ],
)
def test_build_from_string(value, expected):
    assert GenomeBuild.from_string(value) is expected


@pytest.mark.parametrize(
    "value,expected",
    [("chr1", "1"), ("1", "1"), ("chrX", "X"), ("x", "X"), ("M", "MT"), ("MT", "MT")],
)
def test_normalize_chrom(value, expected):
    assert normalize_chrom(value) == expected


def test_chrom_lengths_known():
    assert chrom_info("21", GenomeBuild.GRCh37).length == 48129895
    assert chrom_info("1", GenomeBuild.GRCh38).length == 248956422
    assert chrom_info("X", GenomeBuild.GRCh37).length == 155270560


def test_unknown_chrom_returns_none():
    assert chrom_info("99", GenomeBuild.GRCh37) is None


def test_arm_region_p_and_q():
    p = arm_region("1", "p", GenomeBuild.GRCh37)
    q = arm_region("1", "q", GenomeBuild.GRCh37)
    assert p[0] == 1
    assert q[1] == chrom_info("1", GenomeBuild.GRCh37).length
    assert p[1] < q[0]
