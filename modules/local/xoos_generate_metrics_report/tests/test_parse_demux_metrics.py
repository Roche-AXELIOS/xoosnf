"""Tests for generate-metrics-report demux parsing."""

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


FIXTURE_DIR = (
    Path(__file__).resolve().parent.parent
    / "parse_demux_metrics"
    / "tests"
    / "fixtures"
)


def _metric_values(metrics):
    return {metric.name: metric.value for metric in metrics}


def test_parse_run_stats_accepts_canonical_header_with_denominator():
    run_stats = FIXTURE_DIR / "canonical_run_stats.tsv"

    values = _metric_values(generate_metrics_report.parse_run_stats(run_stats))

    assert values == {
        "demux/total_reads": "100",
        "demux/total_reads%": "100.00",
        "demux/too_short_reads": "0",
        "demux/too_short_reads%": "NA",
    }


def test_parse_demux_metrics_cli_accepts_canonical_run_stats(tmp_path: Path):
    run_stats = FIXTURE_DIR / "canonical_run_stats.tsv"
    sample_stats = FIXTURE_DIR / "sample_stats.tsv"
    output = tmp_path / "demux_metrics.json"

    parser = generate_metrics_report.create_parser()
    args = parser.parse_args(
        [
            "parse-demux-metrics",
            "--run-stats",
            str(run_stats),
            "--sample-stats",
            str(sample_stats),
            "--output",
            str(output),
        ]
    )
    args.func(args)

    metrics = json.loads(output.read_text())
    assert {metric["name"] for metric in metrics} >= {
        "demux/total_reads",
        "demux/total_reads%",
        "demux/too_short_reads",
        "demux/too_short_reads%",
    }
    assert output.exists()
