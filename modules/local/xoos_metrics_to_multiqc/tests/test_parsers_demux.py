"""Unit tests for demux metric parsers."""

import textwrap
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers import parse_demux_run_stats, parse_demux_sample_stats


# ---------------------------------------------------------------------------
# run_stats.tsv
# ---------------------------------------------------------------------------

RUN_STATS = textwrap.dedent("""\
    #Version=1.0.0
    metrics_name\tcount\tpercentage
    total_reads\t10000\t100.0
    passing_reads\t8000\t80.0
    failed_assigned_reads\t1000\t10.0
    unassigned_reads\t1000\t10.0
    full_duplex_reads\t6000\t60.0
    partial_duplex_reads\t2000\t20.0
    no_hairpin_reads\t500\t5.0
    too_long_reads\t300\t3.0
    too_short_reads\t200\t2.0
""")

# Simplex run_stats.tsv uses a different set of field names. The disposition
# (assigned + failed + unassigned) still sums to total reads.
SIMPLEX_RUN_STATS = textwrap.dedent("""\
    #Version=1.0.0
    metrics_name\tcount\tpercentage
    total_reads\t10000\t100.0
    preassignment_passing_reads\t9000\t90.0
    failed_reads\t1000\t10.0
    too_short_reads\t1000\t10.0
    assigned_reads\t7000\t70.0
    failed_assigned_reads\t800\t8.0
    too_short_trimmed_reads\t500\t5.0
    sid_discordant_discarded_reads\t300\t3.0
    unassigned_reads\t2000\t20.0
    full_reads\t5000\t50.0
    partial_reads\t2000\t20.0
""")

# ---------------------------------------------------------------------------
# sample_stats.tsv
# ---------------------------------------------------------------------------

SAMPLE_STATS = textwrap.dedent("""\
    #Version=1.0.0
    metric\tSampleA\tSampleB
    index_sequence\tACGT\tTGCA
    assigned_reads\t5000\t4000
    passing_reads\t4000\t3500
    full_duplex_reads\t3000\t2500
    partial_duplex_reads\t1000\t1000
    failed_assigned_reads\t1000\t500
    too_many_errors_reads\t400\t200
    too_short_trimmed_reads\t300\t150
    too_long_consensus_reads\t200\t100
    no_endadapter_reads\t100\t50
    total_bases\t500000\t400000
    concordant_duplex_bases\t300000\t250000
    discordant_duplex_bases\t50000\t40000
    non_duplex_bases\t150000\t110000
""")


@pytest.fixture
def run_stats(tmp_path: Path) -> Path:
    p = tmp_path / "run_stats.tsv"
    p.write_text(RUN_STATS)
    return p


@pytest.fixture
def simplex_run_stats(tmp_path: Path) -> Path:
    p = tmp_path / "run_stats.tsv"
    p.write_text(SIMPLEX_RUN_STATS)
    return p


@pytest.fixture
def sample_stats(tmp_path: Path) -> Path:
    p = tmp_path / "sample_stats.tsv"
    p.write_text(SAMPLE_STATS)
    return p


# Simplex sample_stats.tsv: counts use full_reads/partial_reads (no _duplex
# suffix) and assigned_reads; no base composition rows.
SIMPLEX_SAMPLE_STATS = textwrap.dedent("""\
    #Version=1.0.0
    metric\tSampleX\tSampleY\tUnassigned
    index_sequence\tACGT\tTGCA\tNA
    assigned_reads\t5000\t4000\t2000
    full_reads\t3000\t2500\t0
    partial_reads\t2000\t1500\t0
    both_sid_detected_reads\t4800\t3800\t0
    sid_discordant_reads\t100\t80\t0
    index_hopping_reads\t50\t40\t0
    perfect_index_reads\t4500\t3600\t0
""")


@pytest.fixture
def simplex_sample_stats(tmp_path: Path) -> Path:
    p = tmp_path / "sample_stats.tsv"
    p.write_text(SIMPLEX_SAMPLE_STATS)
    return p


class TestParseDemuxRunStats:
    def test_duplex_run_stats(self, run_stats: Path):
        result = parse_demux_run_stats(run_stats)
        assert result["_mode"] == "duplex"
        assert result["Total Reads"] == 10000
        # Canonical disposition (sums to total)
        assert result["Passing Reads"] == 8000
        assert result["Failed Reads"] == 1000
        assert result["Unassigned Reads"] == 1000
        assert result["Passing Reads"] + result["Failed Reads"] + result["Unassigned Reads"] == result["Total Reads"]
        # Duplex unassigned breakdown
        assert result["No Hairpin"] == 500
        assert result["Too Long"] == 300
        assert result["Too Short"] == 200
        # Passing composition
        assert result["Full Reads"] == 6000
        assert result["Partial Reads"] == 2000

    def test_simplex_run_stats(self, simplex_run_stats: Path):
        result = parse_demux_run_stats(simplex_run_stats)
        assert result["_mode"] == "simplex"
        assert result["Total Reads"] == 10000
        # In simplex, assigned reads are the passing reads.
        assert result["Passing Reads"] == 7000
        assert result["Failed Reads"] == 1000
        assert result["Unassigned Reads"] == 2000
        # Simplex breakdown categories
        assert result["Too Short After Trimming"] == 500
        assert result["SID Discordant"] == 300
        # No duplex-only keys leak into simplex output
        assert "No Hairpin" not in result

    def test_empty_file_yields_empty_result(self, tmp_path: Path):
        p = tmp_path / "run_stats.tsv"
        p.write_text("#Version=1.0.0\nmetrics_name\tcount\tpercentage\n")
        assert parse_demux_run_stats(p) == {}


class TestParseDemuxSampleStats:
    def test_duplex_sample_stats(self, sample_stats: Path):
        result = parse_demux_sample_stats(sample_stats)
        assert "SampleA" in result
        assert "SampleB" in result
        assert result["SampleA"]["Assigned Reads"] == 5000
        assert result["SampleA"]["Full Reads"] == 3000
        assert result["SampleA"]["Partial Reads"] == 1000
        assert result["SampleB"]["Assigned Reads"] == 4000

    def test_duplex_base_composition_metrics(self, sample_stats: Path):
        result = parse_demux_sample_stats(sample_stats)
        assert result["SampleA"]["Concordant Bases"] == 300000
        assert result["SampleA"]["Discordant Bases"] == 50000
        assert result["SampleA"]["Non-duplex Bases"] == 150000
        assert result["SampleB"]["Concordant Bases"] == 250000
        assert result["SampleB"]["Discordant Bases"] == 40000

    def test_simplex_sample_stats(self, simplex_sample_stats: Path):
        """Simplex field names (full_reads/partial_reads) map to canonical names."""
        result = parse_demux_sample_stats(simplex_sample_stats)
        assert result["SampleX"]["Assigned Reads"] == 5000
        assert result["SampleX"]["Full Reads"] == 3000
        assert result["SampleX"]["Partial Reads"] == 2000
        assert result["SampleX"]["SID Discordant"] == 100
        # index_hopping_reads is an alias for sid_discordant_reads and must not
        # be surfaced as a separate metric (would double-count).
        assert "Index Hopping" not in result["SampleX"]
        assert result["SampleY"]["Assigned Reads"] == 4000
        # Simplex has no base composition rows
        assert "Concordant Bases" not in result["SampleX"]
