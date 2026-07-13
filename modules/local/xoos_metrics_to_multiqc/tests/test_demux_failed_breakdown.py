"""Tests for the demux per-sample Failed Reads Breakdown bar.

"No End Adapter" (no_endadapter_reads) is a passing-read attribute in the demux
source (PassingCounts: "reads that didn't have an endadapter at all"), not a
failure reason (FailedAssignedCounts). It must not appear in the Failed Reads
Breakdown bar, otherwise passing reads are miscounted as failures.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from xoos_to_multiqc import process_inputs

# Duplex sample_stats with both genuine failure reasons and the passing-read
# no_endadapter attribute. The three genuine reasons partition
# failed_assigned_reads (400 + 300 + 300 = 1000).
SAMPLE_STATS = """\
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
too_long_consensus_reads\t300\t150
no_endadapter_reads\t100\t50
"""


def _run(tmp_path: Path) -> dict:
    ss = tmp_path / "ss.tsv"
    ss.write_text(SAMPLE_STATS)
    out = tmp_path / "out"
    ns = argparse.Namespace(
        _sample_inputs={},
        output_dir=out,
        region=None,
        stage=None,
        demux_run_stats=None,
        demux_sample_stats=ss,
    )
    process_inputs(ns)
    return json.loads((out / "xoos_demux_sample_metrics_mqc.json").read_text())


def _failed_breakdown(section: dict) -> dict:
    """Return the Failed Reads Breakdown dataset (per-sample categories)."""
    labels = section["pconfig"]["data_labels"]
    datasets = section["data"]
    for label, dataset in zip(labels, datasets):
        name = label["name"] if isinstance(label, dict) else label
        if "Failed Reads Breakdown" in name:
            return dataset
    raise AssertionError("Failed Reads Breakdown dataset not found")


def test_no_end_adapter_excluded_from_failed_breakdown(tmp_path: Path):
    section = _run(tmp_path)
    failed = _failed_breakdown(section)

    for sample, categories in failed.items():
        assert "No End Adapter" not in categories, sample


def test_failed_breakdown_keeps_genuine_failure_reasons(tmp_path: Path):
    section = _run(tmp_path)
    failed = _failed_breakdown(section)

    assert failed["SampleA"] == {
        "Too Many Errors": 400,
        "Too Short After Trimming": 300,
        "Too Long": 300,
    }


def test_failed_breakdown_sums_to_failed_assigned(tmp_path: Path):
    # The three genuine reasons partition failed_assigned_reads; adding
    # no_endadapter would inflate the total beyond it.
    section = _run(tmp_path)
    failed = _failed_breakdown(section)

    assert sum(failed["SampleA"].values()) == 1000  # failed_assigned_reads
    assert sum(failed["SampleB"].values()) == 500
