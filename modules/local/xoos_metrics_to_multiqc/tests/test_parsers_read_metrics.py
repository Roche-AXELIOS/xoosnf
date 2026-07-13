"""Unit tests for read metric parsing (Any vs Post-filter variants)."""

import sys
import textwrap
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers import parse_read_metrics

# read_metrics_summary.tsv encodes the filter variant in the metric name:
# bare names are the "Any" reads/bases; post_filter_*/*_passing_filter names are
# the post-filter subset. Column layout is metric_name/count/percentage.
READ_METRICS = textwrap.dedent("""\
    metric_name\tcount\tpercentage
    total_reads\t10000\t100.0
    duplicate_reads\t1500\t15.0
    unmapped_reads\t200\t2.0
    post_filter_read_count\t8000\t80.0
    total_bases\t1000000\t100.0
    aligned_bases\t950000\t95.0
    soft_clipped_bases\t50000\t5.0
    total_bases_passing_filter\t800000\t80.0
    aligned_bases_passing_filter\t790000\t79.0
    soft_clipped_bases_passing_filter\t10000\t1.0
""")


@pytest.fixture
def read_metrics_file(tmp_path: Path) -> Path:
    p = tmp_path / "read_metrics_summary.tsv"
    p.write_text(READ_METRICS)
    return p


def test_returns_any_and_post_filter_variants(read_metrics_file: Path):
    result = parse_read_metrics(read_metrics_file)
    assert set(result) == {"Any", "Post-filter"}


def test_any_variant_uses_bare_metric_names(read_metrics_file: Path):
    any_data = parse_read_metrics(read_metrics_file)["Any"]
    assert any_data["read_counts"]["Total Reads"] == 10000
    assert any_data["read_counts"]["Duplicate Reads"] == 1500
    assert any_data["read_counts"]["Unmapped Reads"] == 200
    assert any_data["base_counts"]["Total Bases"] == 1000000
    assert any_data["base_counts"]["Aligned Bases"] == 950000
    assert any_data["base_counts"]["Soft-clipped Bases"] == 50000


def test_post_filter_variant_uses_passing_filter_names(read_metrics_file: Path):
    pf_data = parse_read_metrics(read_metrics_file)["Post-filter"]
    # Post-filter read count maps to the same display label as "Total Reads".
    assert pf_data["read_counts"]["Total Reads"] == 8000
    assert pf_data["base_counts"]["Total Bases"] == 800000
    assert pf_data["base_counts"]["Aligned Bases"] == 790000
    assert pf_data["base_counts"]["Soft-clipped Bases"] == 10000


def test_percentages_surfaced(read_metrics_file: Path):
    any_data = parse_read_metrics(read_metrics_file)["Any"]
    assert any_data["read_counts"]["Total Reads (%)"] == 100.0
    assert any_data["base_counts"]["Aligned Bases (%)"] == 95.0


def test_empty_variants_are_dropped(tmp_path: Path):
    # Only "Any" data present -> "Post-filter" variant omitted entirely.
    p = tmp_path / "read_metrics_summary.tsv"
    p.write_text("metric_name\tcount\tpercentage\ntotal_reads\t10\t100.0\n")
    result = parse_read_metrics(p)
    assert set(result) == {"Any"}
