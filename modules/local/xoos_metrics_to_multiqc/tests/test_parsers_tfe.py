"""Unit tests for parse_tfe_results."""

import textwrap
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers import parse_tfe_results


TFE_RESULTS_TSV = textwrap.dedent("""\
    #Version=1.0.0
    #CommandLine=tumor_fraction_estimator --input test.bam
    metric_name\tvalue
    tumor_fraction\t0.0012
    mean_vaf\t0.0015
    error_rate\t0.00023
    total_alt\t5
    total_depth\t50000
    total_variants\t100
    variants_detected\t3
    passing_variants\t80
    passing_variants_detected\t2
    some_other_key\t999
""")

TFE_RESULTS_MINIMAL = textwrap.dedent("""\
    #Version=1.0.0
    metric_name\tvalue
    tumor_fraction\t0.005
    total_variants\t50
""")

EMPTY_TFE_TSV = textwrap.dedent("""\
    #Version=1.0.0
    #CommandLine=tumor_fraction_estimator --input test.bam
    metric_name\tvalue
""")


@pytest.fixture
def tfe_file(tmp_path: Path) -> Path:
    p = tmp_path / "tumor_fraction.results.tsv"
    p.write_text(TFE_RESULTS_TSV)
    return p


@pytest.fixture
def tfe_file_minimal(tmp_path: Path) -> Path:
    p = tmp_path / "tumor_fraction.results.tsv"
    p.write_text(TFE_RESULTS_MINIMAL)
    return p


@pytest.fixture
def empty_tfe_file(tmp_path: Path) -> Path:
    p = tmp_path / "tumor_fraction.results.tsv"
    p.write_text(EMPTY_TFE_TSV)
    return p


class TestParseTfeResults:
    def test_extracts_selected_metrics(self, tfe_file: Path):
        result = parse_tfe_results(tfe_file)

        assert result["tumor fraction"] == 0.0012
        assert result["mean vaf"] == 0.0015
        assert result["error rate"] == 0.00023
        assert result["total alt"] == 5
        assert result["total depth"] == 50000
        assert result["total variants"] == 100
        assert result["variants detected"] == 3
        assert result["passing variants"] == 80
        assert result["passing variants detected"] == 2

    def test_excludes_non_selected_keys(self, tfe_file: Path):
        result = parse_tfe_results(tfe_file)
        assert "some other key" not in result

    def test_replaces_underscores_with_spaces(self, tfe_file: Path):
        result = parse_tfe_results(tfe_file)
        assert "tumor_fraction" not in result
        assert "tumor fraction" in result

    def test_skips_comment_and_header_lines(self, tfe_file: Path):
        result = parse_tfe_results(tfe_file)
        assert len(result) == 9

    def test_empty_file_returns_empty_dict(self, empty_tfe_file: Path):
        result = parse_tfe_results(empty_tfe_file)
        assert result == {}

    def test_minimal_file(self, tfe_file_minimal: Path):
        result = parse_tfe_results(tfe_file_minimal)
        assert result["tumor fraction"] == 0.005
        assert result["total variants"] == 50
        assert len(result) == 2
