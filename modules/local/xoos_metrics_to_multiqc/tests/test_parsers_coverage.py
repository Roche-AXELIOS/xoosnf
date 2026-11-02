"""Unit tests for coverage metric parsers."""

import textwrap
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers import parse_coverage_metrics


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SIMPLEX_COVERAGE_STATS = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    metric_name\tany_coverage\tpost_filter_coverage
    median\t30\t25
    mean\t31.50\t26.20
    stddev\t10.20\t8.50
    min\t5\t3
    max\t60\t50
    percentile_10\t12\t10
    percentile_25\t20\t15
    percentile_50\t30\t25
    percentile_75\t40\t35
    percentile_90\t52\t45
    ratio_90_to_10_percentile\t4.333333\t4.500000
    ratio_median_to_10_percentile\t2.500000\t2.500000
""")

ALL_NA_COVERAGE_STATS = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    metric_name\tany_coverage\tpost_filter_coverage
    median\tNA\tNA
    mean\tNA\tNA
    stddev\tNA\tNA
    min\tNA\tNA
    max\tNA\tNA
    percentile_25\tNA\tNA
    percentile_75\tNA\tNA
    ratio_90_to_10_percentile\tNA\tNA
    ratio_median_to_10_percentile\tNA\tNA
""")

PARTIAL_NA_COVERAGE_STATS = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    metric_name\tany_coverage\tpost_filter_coverage
    median\t30\tNA
    mean\t31.50\tNA
    stddev\t10.20\tNA
    min\t5\tNA
    max\t60\tNA
    percentile_25\t20\tNA
    percentile_75\t40\tNA
    ratio_90_to_10_percentile\t4.333333\tNA
    ratio_median_to_10_percentile\t2.500000\tNA
""")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestParseCoverageMetrics:
    """Tests for parse_coverage_metrics."""

    def test_simplex_extracts_boxplot_stats(self, tmp_path):
        f = tmp_path / "coverage_stats.tsv"
        f.write_text(SIMPLEX_COVERAGE_STATS)
        result = parse_coverage_metrics(f)

        assert result["coverage"] == {
            "min": 5, "q1": 20, "median": 30, "q3": 40, "max": 60, "mean": 31.5,
        }
        assert result["post_filter_coverage"] == {
            "min": 3, "q1": 15, "median": 25, "q3": 35, "max": 50, "mean": 26.2,
        }

    def test_all_na_produces_empty_dicts(self, tmp_path):
        """When all values are NA (empty histogram), dicts should be empty."""
        f = tmp_path / "coverage_stats.tsv"
        f.write_text(ALL_NA_COVERAGE_STATS)
        result = parse_coverage_metrics(f)

        assert result["coverage"] == {}
        assert result["post_filter_coverage"] == {}

    def test_partial_na_only_populates_numeric_columns(self, tmp_path):
        """When post_filter_coverage is all NA but any_coverage has data."""
        f = tmp_path / "coverage_stats.tsv"
        f.write_text(PARTIAL_NA_COVERAGE_STATS)
        result = parse_coverage_metrics(f)

        assert result["coverage"] == {
            "min": 5, "q1": 20, "median": 30, "q3": 40, "max": 60, "mean": 31.5,
        }
        assert result["post_filter_coverage"] == {}

    def test_ignores_non_summary_rows(self, tmp_path):
        """Rows like stddev, ratio_*, percentile_10 should be ignored."""
        f = tmp_path / "coverage_stats.tsv"
        f.write_text(SIMPLEX_COVERAGE_STATS)
        result = parse_coverage_metrics(f)

        # Only the mapped keys should be present
        expected_keys = {"min", "q1", "median", "q3", "max", "mean"}
        assert set(result["coverage"].keys()) == expected_keys
        assert set(result["post_filter_coverage"].keys()) == expected_keys

    def test_all_values_are_numeric(self, tmp_path):
        """No string values should appear in the output."""
        f = tmp_path / "coverage_stats.tsv"
        f.write_text(SIMPLEX_COVERAGE_STATS)
        result = parse_coverage_metrics(f)

        for category in ("coverage", "post_filter_coverage"):
            for key, value in result[category].items():
                assert isinstance(value, (int, float)), (
                    f"{category}.{key} = {value!r} is not numeric"
                )
