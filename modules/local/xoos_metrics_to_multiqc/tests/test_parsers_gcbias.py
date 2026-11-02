"""Unit tests for alignment_metrics GC-bias parsers."""

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers import parse_gc_bias_metrics, parse_gc_bias_summary


def test_parse_gc_bias_metrics_returns_numeric_points_sorted_by_gc_percentage(tmp_path: Path):
    gc_bias = tmp_path / "gc_bias.tsv"
    gc_bias.write_text(
        """##RocheCommandLine=<ID=alignment_metrics>
gc_percentage\twindow_count\tread_count\tpost_filter_read_count\tnormalized_coverage\tnormalized_post_filter_coverage\terror_bar\tpost_filter_error_bar
70\t100\t10\t8\t1.2\t1.1\t0.1\t0.1
20\t100\t10\t8\t1.0\t0.9\t0.1\t0.1
21\t100\t10\t8\t1.0\tNA\t0.1\t0.1
"""
    )

    assert parse_gc_bias_metrics(gc_bias) == [[20, 1.0], [21, 1.0], [70, 1.2]]


def test_parse_gc_bias_summary_returns_unfiltered_dropout_values(tmp_path: Path):
    gc_bias_summary = tmp_path / "gc_bias_summary.tsv"
    gc_bias_summary.write_text(
        """##RocheCommandLine=<ID=alignment_metrics>
metric_name     value
at_dropout      1.390350
gc_dropout      0.266587
post_filter_at_dropout  0.762381
post_filter_gc_dropout  0.008090
"""
    )

    assert parse_gc_bias_summary(gc_bias_summary) == {
        "AT Dropout": 1.39035,
        "GC Dropout": 0.266587,
    }
