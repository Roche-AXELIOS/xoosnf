"""Unit tests for combine_reports merge logic."""

import json
from pathlib import Path

import pytest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from combine_reports import combine_reports, _merge_data, _merge_data_labels


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_report(report_id, data_labels, data):
    """Build a minimal MultiQC JSON report dict."""
    return {
        "id": report_id,
        "section_name": "Test Section",
        "plot_type": "boxplot",
        "pconfig": {
            "id": f"{report_id}_table",
            "title": "Test",
            "data_labels": data_labels,
        },
        "data": data,
    }


def _write_report(directory, filename, report):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / filename).write_text(json.dumps(report))


def _read_combined(output_dir, filename):
    return json.loads((output_dir / filename).read_text())


# ---------------------------------------------------------------------------
# _merge_data_labels
# ---------------------------------------------------------------------------

class TestMergeDataLabels:
    def test_incoming_longer_appends_extra(self):
        result = _merge_data_labels(["A", "B"], ["A", "B", "C"])
        assert result == ["A", "B", "C"]

    def test_existing_longer_keeps_existing(self):
        result = _merge_data_labels(["A", "B", "C"], ["A", "B"])
        assert result == ["A", "B", "C"]

    def test_same_length_keeps_existing(self):
        result = _merge_data_labels(["A", "B"], ["X", "Y"])
        assert result == ["A", "B"]

    def test_empty_existing(self):
        result = _merge_data_labels([], ["A", "B"])
        assert result == ["A", "B"]

    def test_empty_incoming(self):
        result = _merge_data_labels(["A", "B"], [])
        assert result == ["A", "B"]

    def test_both_empty(self):
        result = _merge_data_labels([], [])
        assert result == []


# ---------------------------------------------------------------------------
# _merge_data (list-style)
# ---------------------------------------------------------------------------

class TestMergeDataList:
    def test_same_length_merges_elementwise(self):
        existing = [{"S1": 1}, {"S1": 2}]
        incoming = [{"S2": 3}, {"S2": 4}]
        result, _ = _merge_data(existing, incoming)
        assert result == [{"S1": 1, "S2": 3}, {"S1": 2, "S2": 4}]

    def test_incoming_longer_extends(self):
        existing = [{"S1": 1}]
        incoming = [{"S2": 3}, {"S2": 4}]
        result, _ = _merge_data(existing, incoming)
        assert len(result) == 2
        assert result[0] == {"S1": 1, "S2": 3}
        assert result[1] == {"S2": 4}

    def test_existing_longer_preserves(self):
        existing = [{"S1": 1}, {"S1": 2}, {"S1": 3}]
        incoming = [{"S2": 4}]
        result, _ = _merge_data(existing, incoming)
        assert len(result) == 3
        assert result[0] == {"S1": 1, "S2": 4}
        assert result[1] == {"S1": 2}
        assert result[2] == {"S1": 3}


class TestMergeDataLabelAware:
    """Label-aware merge for region/stage alignment tabs."""

    def test_same_label_merges_samples(self):
        existing = [{"S1": {"reads": 1}}]
        incoming = [{"S2": {"reads": 2}}]
        result, labels = _merge_data(
            existing, incoming, ["Dedup - Autosomes"], ["Dedup - Autosomes"],
        )
        assert labels == ["Dedup - Autosomes"]
        assert result == [{"S1": {"reads": 1}, "S2": {"reads": 2}}]

    def test_different_labels_union_in_order(self):
        existing = [{"S1": {"reads": 1}}]
        incoming = [{"S1": {"reads": 9}}]
        result, labels = _merge_data(
            existing, incoming, ["Dedup - Autosomes"], ["Dedup - High Confidence"],
        )
        assert labels == ["Dedup - Autosomes", "Dedup - High Confidence"]
        assert result == [{"S1": {"reads": 1}}, {"S1": {"reads": 9}}]

    def test_mixed_shared_and_new_labels(self):
        existing = [{"S1": {"reads": 1}}, {"S1": {"reads": 2}}]
        incoming = [{"S2": {"reads": 3}}]
        result, labels = _merge_data(
            existing,
            incoming,
            ["Dedup - Autosomes", "Dedup - High Confidence"],
            ["Dedup - Autosomes"],
        )
        assert labels == ["Dedup - Autosomes", "Dedup - High Confidence"]
        assert result[0] == {"S1": {"reads": 1}, "S2": {"reads": 3}}
        assert result[1] == {"S1": {"reads": 2}}


# ---------------------------------------------------------------------------
# combine_reports: data_labels sync
# ---------------------------------------------------------------------------

class TestCombineReportsDataLabels:
    """Verify data_labels stays in sync with data when tab counts differ."""

    def test_simplex_then_duplex(self, tmp_path):
        """Simplex (2 tabs) merged before duplex (3 tabs)."""
        simplex = _make_report(
            "cov", ["Coverage", "Post-filter"],
            [{"S1": {"median": 30}}, {"S1": {"median": 25}}],
        )
        duplex = _make_report(
            "cov", ["Coverage", "Post-filter", "Duplex"],
            [{"S2": {"median": 35}}, {"S2": {"median": 28}}, {"S2": {"median": 10}}],
        )

        d1, d2 = tmp_path / "a", tmp_path / "b"
        _write_report(d1, "cov_mqc.json", simplex)
        _write_report(d2, "cov_mqc.json", duplex)

        out = tmp_path / "out"
        combine_reports(out, input_dirs=[d1, d2])
        result = _read_combined(out, "cov_mqc.json")

        assert result["pconfig"]["data_labels"] == ["Coverage", "Post-filter", "Duplex"]
        assert len(result["data"]) == 3
        assert "S1" in result["data"][0] and "S2" in result["data"][0]
        assert "S1" in result["data"][1] and "S2" in result["data"][1]
        assert "S2" in result["data"][2]

    def test_duplex_then_simplex(self, tmp_path):
        """Duplex (3 tabs) merged before simplex (2 tabs)."""
        duplex = _make_report(
            "cov", ["Coverage", "Post-filter", "Duplex"],
            [{"S1": {"median": 35}}, {"S1": {"median": 28}}, {"S1": {"median": 10}}],
        )
        simplex = _make_report(
            "cov", ["Coverage", "Post-filter"],
            [{"S2": {"median": 30}}, {"S2": {"median": 25}}],
        )

        d1, d2 = tmp_path / "a", tmp_path / "b"
        _write_report(d1, "cov_mqc.json", duplex)
        _write_report(d2, "cov_mqc.json", simplex)

        out = tmp_path / "out"
        combine_reports(out, input_dirs=[d1, d2])
        result = _read_combined(out, "cov_mqc.json")

        assert result["pconfig"]["data_labels"] == ["Coverage", "Post-filter", "Duplex"]
        assert len(result["data"]) == 3

    def test_same_tab_count(self, tmp_path):
        """Both reports have the same number of tabs."""
        r1 = _make_report("cov", ["A", "B"], [{"S1": 1}, {"S1": 2}])
        r2 = _make_report("cov", ["A", "B"], [{"S2": 3}, {"S2": 4}])

        d1, d2 = tmp_path / "a", tmp_path / "b"
        _write_report(d1, "cov_mqc.json", r1)
        _write_report(d2, "cov_mqc.json", r2)

        out = tmp_path / "out"
        combine_reports(out, input_dirs=[d1, d2])
        result = _read_combined(out, "cov_mqc.json")

        assert result["pconfig"]["data_labels"] == ["A", "B"]
        assert len(result["data"]) == 2

    def test_dict_style_data_unaffected(self, tmp_path):
        """Dict-style data (single-tab tables) should merge normally."""
        r1 = {"id": "acc", "pconfig": {"id": "acc_table"}, "data": {"S1": {"f1": 0.95}}}
        r2 = {"id": "acc", "pconfig": {"id": "acc_table"}, "data": {"S2": {"f1": 0.90}}}

        d1, d2 = tmp_path / "a", tmp_path / "b"
        _write_report(d1, "acc_mqc.json", r1)
        _write_report(d2, "acc_mqc.json", r2)

        out = tmp_path / "out"
        combine_reports(out, input_dirs=[d1, d2])
        result = _read_combined(out, "acc_mqc.json")

        assert result["data"] == {"S1": {"f1": 0.95}, "S2": {"f1": 0.90}}

    def test_three_reports_incremental_tabs(self, tmp_path):
        """Three reports with 1, 2, and 3 tabs respectively."""
        r1 = _make_report("cov", ["A"], [{"S1": 1}])
        r2 = _make_report("cov", ["A", "B"], [{"S2": 2}, {"S2": 3}])
        r3 = _make_report("cov", ["A", "B", "C"], [{"S3": 4}, {"S3": 5}, {"S3": 6}])

        d1, d2, d3 = tmp_path / "a", tmp_path / "b", tmp_path / "c"
        _write_report(d1, "cov_mqc.json", r1)
        _write_report(d2, "cov_mqc.json", r2)
        _write_report(d3, "cov_mqc.json", r3)

        out = tmp_path / "out"
        combine_reports(out, input_dirs=[d1, d2, d3])
        result = _read_combined(out, "cov_mqc.json")

        assert result["pconfig"]["data_labels"] == ["A", "B", "C"]
        assert len(result["data"]) == 3
        assert result["data"][0] == {"S1": 1, "S2": 2, "S3": 4}
