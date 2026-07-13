"""Unit tests for read_collapser parsers."""

import textwrap
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers import (
    parse_read_collapser_summary,
    parse_read_collapser_cluster_size_distributions,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SUMMARY_STATS_TSV = textwrap.dedent("""\
    ##RocheCommandLine=read_collapser markdup --input test.bam
    metric_name\tvalue\tpercentage\tdenominator
    input_reads\t1000\tNA\tNA
    discarded_missing_umi_reads\t0\t0.00\tinput_reads
    discarded_low_mapq_reads\t5\t0.50\tinput_reads
    discarded_by_flags_reads\t10\t1.00\tinput_reads
    discarded_total_reads\t15\t1.50\tinput_reads
    unmapped_reads\t20\t2.00\tinput_reads
    clustering_input_reads\t965\t96.50\tinput_reads
    clustering_reads\t965\tNA\tNA
    duplicate_reads\t100\t10.36\tclustering_reads
    total_clusters\t865\tNA\tNA
    singleton_clusters\t800\t92.49\ttotal_clusters
    forward_strand_clusters\t430\t49.71\ttotal_clusters
    reverse_strand_clusters\t435\t50.29\ttotal_clusters
    mixed_strand_clusters\t0\t0.00\ttotal_clusters
""")

SUMMARY_STATS_CONSENSUS_TSV = textwrap.dedent("""\
    ##RocheCommandLine=read_collapser consensus --input test.bam
    metric_name\tvalue\tpercentage\tdenominator
    input_reads\t2000\tNA\tNA
    discarded_total_reads\t30\t1.50\tinput_reads
    clustering_input_reads\t1970\t98.50\tinput_reads
    duplicate_reads\t200\t10.15\tclustering_reads
    total_clusters\t1770\tNA\tNA
    singleton_clusters\t1500\t84.75\ttotal_clusters
    consensus_clusters\t270\t15.25\ttotal_clusters
    consensus_median_cluster_size\t3\tNA\tNA
    total_consensus_reads\t810\tNA\tNA
""")

CLUSTER_SIZE_DISTRIBUTIONS_TSV = textwrap.dedent("""\
    ##RocheCommandLine=read_collapser markdup --input test.bam
    cluster_sizes\ttotal_clusters\tforward_strand_clusters\treverse_strand_clusters
    2\t50\t25\t25
    3\t30\t15\t15
    4\t10\t5\t5
    5\t5\t3\t2
""")

CLUSTER_SIZE_DISTRIBUTIONS_BINNED_TSV = textwrap.dedent("""\
    ##RocheCommandLine=read_collapser markdup --input test.bam
    cluster_sizes\ttotal_clusters\tforward_strand_clusters
    2\t50\t25
    3\t30\t15
    50+\t7\t4
""")

EMPTY_TSV = textwrap.dedent("""\
    ##RocheCommandLine=read_collapser markdup --input test.bam
    metric_name\tvalue\tpercentage\tdenominator
""")


@pytest.fixture
def summary_stats_file(tmp_path: Path) -> Path:
    p = tmp_path / "summary_stats.tsv"
    p.write_text(SUMMARY_STATS_TSV)
    return p


@pytest.fixture
def summary_stats_consensus_file(tmp_path: Path) -> Path:
    p = tmp_path / "summary_stats.tsv"
    p.write_text(SUMMARY_STATS_CONSENSUS_TSV)
    return p


@pytest.fixture
def cluster_size_file(tmp_path: Path) -> Path:
    p = tmp_path / "cluster_size_distributions.tsv"
    p.write_text(CLUSTER_SIZE_DISTRIBUTIONS_TSV)
    return p


@pytest.fixture
def empty_summary_file(tmp_path: Path) -> Path:
    p = tmp_path / "summary_stats.tsv"
    p.write_text(EMPTY_TSV)
    return p


@pytest.fixture
def cluster_size_binned_file(tmp_path: Path) -> Path:
    p = tmp_path / "cluster_size_distributions.tsv"
    p.write_text(CLUSTER_SIZE_DISTRIBUTIONS_BINNED_TSV)
    return p


@pytest.fixture
def empty_cluster_file(tmp_path: Path) -> Path:
    p = tmp_path / "cluster_size_distributions.tsv"
    p.write_text(textwrap.dedent("""\
        ##RocheCommandLine=read_collapser markdup --input test.bam
        cluster_sizes\ttotal_clusters
    """))
    return p


# ---------------------------------------------------------------------------
# parse_read_collapser_summary
# ---------------------------------------------------------------------------

class TestParseReadCollapserSummary:
    def test_extracts_selected_metrics(self, summary_stats_file: Path):
        result = parse_read_collapser_summary(summary_stats_file)

        assert result["Input Reads"] == 1000
        assert result["Discarded Reads"] == 15
        assert result["Clustering Input Reads"] == 965
        assert result["Duplicate Reads"] == 100
        assert result["Total Clusters"] == 865
        assert result["Singleton Clusters"] == 800

    def test_includes_percentages(self, summary_stats_file: Path):
        result = parse_read_collapser_summary(summary_stats_file)

        assert result["Discarded Reads (%)"] == 1.50
        assert result["Clustering Input Reads (%)"] == 96.50
        assert result["Duplicate Reads (%)"] == 10.36
        assert result["Singleton Clusters (%)"] == 92.49

    def test_omits_na_percentages(self, summary_stats_file: Path):
        result = parse_read_collapser_summary(summary_stats_file)

        # input_reads has percentage=NA, so no "(%)" entry
        assert "Input Reads (%)" not in result
        assert "Total Clusters (%)" not in result

    def test_excludes_non_selected_metrics(self, summary_stats_file: Path):
        result = parse_read_collapser_summary(summary_stats_file)

        # These metrics exist in the TSV but are not in the selection
        assert "Unmapped Reads" not in result
        assert "Forward Strand Clusters" not in result
        assert "Clustering Reads" not in result

    def test_consensus_only_metrics_present(self, summary_stats_consensus_file: Path):
        result = parse_read_collapser_summary(summary_stats_consensus_file)

        assert result["Consensus Clusters"] == 270
        assert result["Consensus Median Cluster Size"] == 3
        assert result["Total Consensus Reads"] == 810

    def test_consensus_metrics_absent_in_markdup(self, summary_stats_file: Path):
        result = parse_read_collapser_summary(summary_stats_file)

        assert "Consensus Clusters" not in result
        assert "Consensus Median Cluster Size" not in result
        assert "Total Consensus Reads" not in result

    def test_empty_file_returns_empty_dict(self, empty_summary_file: Path):
        result = parse_read_collapser_summary(empty_summary_file)
        assert result == {}


# ---------------------------------------------------------------------------
# parse_read_collapser_cluster_size_distributions
# ---------------------------------------------------------------------------

class TestParseReadCollapserClusterSizeDistributions:
    def test_parses_all_distributions(self, cluster_size_file: Path):
        result = parse_read_collapser_cluster_size_distributions(cluster_size_file)

        assert "total_clusters" in result
        assert "forward_strand_clusters" in result
        assert "reverse_strand_clusters" in result

    def test_xy_data_structure(self, cluster_size_file: Path):
        result = parse_read_collapser_cluster_size_distributions(cluster_size_file)

        total = result["total_clusters"]
        assert total == [[2, 50], [3, 30], [4, 10], [5, 5]]

    def test_preserves_per_distribution_values(self, cluster_size_file: Path):
        result = parse_read_collapser_cluster_size_distributions(cluster_size_file)

        fwd = result["forward_strand_clusters"]
        assert fwd == [[2, 25], [3, 15], [4, 5], [5, 3]]

        rev = result["reverse_strand_clusters"]
        assert rev == [[2, 25], [3, 15], [4, 5], [5, 2]]

    def test_cluster_sizes_column_excluded(self, cluster_size_file: Path):
        result = parse_read_collapser_cluster_size_distributions(cluster_size_file)
        assert "cluster_sizes" not in result

    def test_binned_values_converted_to_numeric(self, cluster_size_binned_file: Path):
        result = parse_read_collapser_cluster_size_distributions(cluster_size_binned_file)

        total = result["total_clusters"]
        assert total == [[2, 50], [3, 30], [50, 7]]

        fwd = result["forward_strand_clusters"]
        assert fwd == [[2, 25], [3, 15], [50, 4]]

    def test_empty_file_returns_empty_dict(self, empty_cluster_file: Path):
        result = parse_read_collapser_cluster_size_distributions(empty_cluster_file)
        assert result == {}
