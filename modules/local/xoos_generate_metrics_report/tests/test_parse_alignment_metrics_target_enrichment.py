"""Tests for target-enrichment alignment-metrics report parsing."""

import importlib.util
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parent.parent / "generate-metrics-report.py"
spec = importlib.util.spec_from_file_location("generate_metrics_report", SCRIPT_PATH)
generate_metrics_report = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = generate_metrics_report
assert spec.loader is not None
spec.loader.exec_module(generate_metrics_report)


UNIFORMITY_SUMMARY = """metric_name\tvalue
fold_80_base_penalty\t1.8
pct_bases_at_0.5x_to_2x_mean_coverage\t91.2
mean_coverage\t125.4
median_coverage\t120.0
target_region_count\t42
target_regions_with_no_coverage\t1
"""

READ_METRICS_TE = """metric_name\tcount\tdenominator\tpercentage
total_mapped_alignments_in_bam\t1000\t\t
on_target_alignments\t900\t1000\t90.0
on_target_reads_passing_filter\t850\t1000\t85.0
"""


def _values(metrics):
    return {metric.name: metric.value for metric in metrics}


def test_parse_target_enrichment_alignment_metrics(tmp_path: Path):
    coverage_dir = tmp_path / "coverage"
    read_metrics_dir = tmp_path / "read_metrics"
    coverage_dir.mkdir()
    read_metrics_dir.mkdir()

    uniformity_summary = coverage_dir / "coverage_uniformity_summary.tsv"
    read_metrics_te = read_metrics_dir / "read_metrics_te.tsv"
    uniformity_summary.write_text(UNIFORMITY_SUMMARY)
    read_metrics_te.write_text(READ_METRICS_TE)

    values = _values(
        generate_metrics_report.parse_alignment_metrics(
            "HG001", "alignment_metrics/dedup/bwa/targets", uniformity_summary
        )
        + generate_metrics_report.parse_alignment_metrics(
            "HG001", "alignment_metrics/dedup/bwa/targets", read_metrics_te
        )
    )

    assert values == {
        "alignment_metrics/dedup/bwa/targets/coverage/fold_80_base_penalty/value": "1.8",
        "alignment_metrics/dedup/bwa/targets/coverage/pct_bases_at_0.5x_to_2x_mean_coverage/value": "91.2",
        "alignment_metrics/dedup/bwa/targets/coverage/mean_coverage/value": "125.4",
        "alignment_metrics/dedup/bwa/targets/coverage/median_coverage/value": "120.0",
        "alignment_metrics/dedup/bwa/targets/coverage/target_region_count/value": "42",
        "alignment_metrics/dedup/bwa/targets/coverage/target_regions_with_no_coverage/value": "1",
        "alignment_metrics/dedup/bwa/targets/read_metrics/total_mapped_alignments_in_bam": "1000",
        "alignment_metrics/dedup/bwa/targets/read_metrics/on_target_alignments": "900",
        "alignment_metrics/dedup/bwa/targets/read_metrics/on_target_alignments%": "90.0",
        "alignment_metrics/dedup/bwa/targets/read_metrics/on_target_reads_passing_filter": "850",
        "alignment_metrics/dedup/bwa/targets/read_metrics/on_target_reads_passing_filter%": "85.0",
    }
