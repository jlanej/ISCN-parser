"""Tests for the command-line interface."""

import os

import pytest

from iscn_parser.cli import main

DATA = os.path.join(os.path.dirname(__file__), "data")


def test_cli_string_to_cnv_stdout(capsys):
    rc = main(["-s", "arr[GRCh37] 1p36.33p36.31(849466_6015936)x1", "-f", "cnv", "--quiet"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "FID\tIID\tCHR" in out
    assert "849466" in out


def test_cli_string_to_bed_stdout(capsys):
    rc = main(["-s", "arr[GRCh38] 7q11.23(72726578_74139390)x3", "-f", "bed", "--quiet"])
    out = capsys.readouterr().out
    assert rc == 0
    assert out.startswith("chr7\t")


def test_cli_file_to_plink_prefix(tmp_path, capsys):
    prefix = os.path.join(str(tmp_path), "out")
    rc = main(["-i", os.path.join(DATA, "array_samples.tsv"), "-o", prefix, "--quiet"])
    assert rc == 0
    assert os.path.exists(prefix + ".cnv")
    assert os.path.exists(prefix + ".fam")
    assert os.path.exists(prefix + ".cnv.map")


def test_cli_unknown_build_returns_2(capsys):
    rc = main(["-s", "47,XX,+21", "-b", "notabuild", "--quiet"])
    assert rc == 2


def test_cli_strict_mode_fails_on_error(capsys):
    rc = main(["-s", "this is not iscn", "--strict", "--quiet"])
    assert rc == 1


def test_cli_non_strict_tolerates_error(capsys):
    rc = main(["-s", "this is not iscn", "--quiet"])
    assert rc == 0


def test_cli_diagnostics_on_stderr(capsys):
    main(["-s", "46,XY,del(5)(p15.2p15.33)", "-b", "GRCh37"])
    err = capsys.readouterr().err
    assert "WARNING" in err
    assert "parsed" in err


def test_cli_cytoband_option(tmp_path, capsys):
    rc = main(
        [
            "-s",
            "46,XY,del(5)(p15.2p15.33)",
            "-b",
            "GRCh37",
            "--cytoband",
            os.path.join(DATA, "cytoBand_test.txt"),
            "--quiet",
        ]
    )
    out = capsys.readouterr().out
    assert rc == 0
    # Precise end coordinate from the cytoband table.
    assert "18400000" in out


def test_cli_requires_input_source():
    with pytest.raises(SystemExit):
        main([])


def test_cli_version(capsys):
    with pytest.raises(SystemExit):
        main(["--version"])
    assert capsys.readouterr().out.strip()
