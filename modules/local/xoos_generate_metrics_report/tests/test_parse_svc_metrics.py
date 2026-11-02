"""Tests for generate-metrics-report small variant caller parsing."""

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


FIXTURE = (
    Path(__file__).resolve().parent.parent
    / "parse_svc_metrics"
    / "tests"
    / "fixtures"
    / "filter_variants_metrics.tsv"
)


def test_parse_svc_metrics_reads_filter_variants_tsv():
    metrics = generate_metrics_report.parse_svc_metrics(
        "HG002", "small_variant_tumor_normal", FIXTURE
    )

    assert [(metric.name, metric.value) for metric in metrics] == [
        ("small_variant_tumor_normal/qc_status", "pass"),
        ("small_variant_tumor_normal/qc_status%", "NA"),
        ("small_variant_tumor_normal/num_pass_variants", "10"),
        ("small_variant_tumor_normal/num_pass_variants%", "50.00"),
        ("small_variant_tumor_normal/num_filter_ML_Score", "2"),
        ("small_variant_tumor_normal/num_filter_ML_Score%", "20.00"),
    ]
    assert all(metric.sample_name == "HG002" for metric in metrics)
    assert all(
        metric.sort_order == generate_metrics_report.SORT_ORDER_SVC
        for metric in metrics
    )


def test_parse_svc_metrics_cli_writes_json(tmp_path: Path):
    output = tmp_path / "svc_metrics.json"
    parser = generate_metrics_report.create_parser()
    args = parser.parse_args(
        [
            "parse-svc-metrics",
            "--sample-name",
            "HG002",
            "--source",
            "small_variant_tumor_normal",
            "--input",
            str(FIXTURE),
            "--output",
            str(output),
        ]
    )
    args.func(args)

    data = json.loads(output.read_text())
    assert data[0] == {
        "name": "small_variant_tumor_normal/qc_status",
        "value": "pass",
        "sort_order": generate_metrics_report.SORT_ORDER_SVC,
        "sample_name": "HG002",
    }
