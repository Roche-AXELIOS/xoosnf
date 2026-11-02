"""Tests for alignment_metrics GC bias report parsing."""

import importlib.util
import json
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parent.parent / "generate-metrics-report.py"
spec = importlib.util.spec_from_file_location("generate_metrics_report", SCRIPT_PATH)
generate_metrics_report = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = generate_metrics_report
assert spec.loader is not None
spec.loader.exec_module(generate_metrics_report)


SUMMARY = """##RocheCommandLine=<ID=alignment_metrics>
metric_name\tvalue
at_dropout\t0.9
gc_dropout\t1.7
post_filter_at_dropout\t0.6
post_filter_gc_dropout\t1.4
post_filter_first_quintile_normalized_coverage\t0.8
post_filter_second_quintile_normalized_coverage\t1.0
post_filter_third_quintile_normalized_coverage\t0.9
post_filter_fourth_quintile_normalized_coverage\t0.8
post_filter_fifth_quintile_normalized_coverage\t0.7
"""

DETAIL = """##RocheCommandLine=<ID=alignment_metrics>
gc_percentage\twindow_count\tread_count\tpost_filter_read_count\tnormalized_coverage\tnormalized_post_filter_coverage\terror_bar\tpost_filter_error_bar
20\t100\t10\t8\t1.0\t0.9\t0.1\t0.1
21\t100\t10\t8\t1.0\tNA\t0.1\t0.1
70\t100\t10\t8\t1.0\t1.1\t0.1\t0.1
"""


def _values(metrics):
    return {metric.name: metric.value for metric in metrics}


def test_parse_alignment_metrics_gc_bias_uses_post_filter_metrics(tmp_path: Path):
    summary = tmp_path / "gc_bias_summary.tsv"
    summary.write_text(SUMMARY)

    values = _values(
        generate_metrics_report.parse_alignment_metrics_gc_bias_summary(
            "HG001", "gc_bias/dedup", summary
        )
    )

    assert values == {
        "gc_bias/dedup/at_dropout": "0.6",
        "gc_bias/dedup/gc_dropout": "1.4",
        "gc_bias/dedup/gc_nc_0_19": "0.8",
        "gc_bias/dedup/gc_nc_20_39": "1.0",
        "gc_bias/dedup/gc_nc_40_59": "0.9",
        "gc_bias/dedup/gc_nc_60_79": "0.8",
        "gc_bias/dedup/gc_nc_80_100": "0.7",
    }


def test_parse_alignment_metrics_gc_bias_detail_uses_post_filter_coverage(tmp_path: Path):
    detail = tmp_path / "gc_bias.tsv"
    detail.write_text(DETAIL)

    values = _values(
        generate_metrics_report.parse_alignment_metrics_gc_bias_detail(
            "HG001", "gc_bias/dedup", detail
        )
    )

    assert values["gc_bias/dedup/20/normalized_coverage"] == "0.9"
    assert values["gc_bias/dedup/70/normalized_coverage"] == "1.1"
    assert "gc_bias/dedup/21/normalized_coverage" not in values
    assert values["gc_bias/dedup/bins_20-70_normalized_coverage_0.8-1.2%"] == "3.92"


def test_alignment_metrics_gc_bias_cli_writes_report(tmp_path: Path):
    summary = tmp_path / "gc_bias_summary.tsv"
    detail = tmp_path / "gc_bias.tsv"
    output = tmp_path / "gc_bias.json"
    summary.write_text(SUMMARY)
    detail.write_text(DETAIL)

    parser = generate_metrics_report.create_parser()
    args = parser.parse_args(
        [
            "parse-alignment-metrics-gc-bias",
            "--sample-name",
            "HG001",
            "--source",
            "gc_bias/dedup",
            "--summary",
            str(summary),
            "--detail",
            str(detail),
            "--output",
            str(output),
        ]
    )
    args.func(args)

    assert {metric["name"] for metric in json.loads(output.read_text())} >= {
        "gc_bias/dedup/at_dropout",
        "gc_bias/dedup/20/normalized_coverage",
    }
