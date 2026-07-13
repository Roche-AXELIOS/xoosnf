"""Tests for region/stage/variant flattening in table sections.

MultiQC custom-content tables cannot render the ``data_labels`` dataset
switcher (only plot types can). Accuracy and read-metrics are table sections,
so region/stage and the Any/Post-filter filter variant must be folded into the
table row key instead of emitted as selectable tabs. These tests pin that
behaviour so a table section is never silently reduced to its first dataset.
"""

import argparse
import json
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import xoos_to_multiqc
from xoos_to_multiqc import _row_key, _tab_label, process_inputs

ACCURACY_TSV = textwrap.dedent("""\
    metric_name\tvalue
    substitution_phred\t20.25
    indel_phred\t23.10
    overall_phred\t18.40
""")

READ_METRICS_TSV = textwrap.dedent("""\
    metric_name\tcount\tpercentage
    total_reads\t10000\t100.0
    post_filter_read_count\t8000\t80.0
    total_bases\t1000000\t100.0
    total_bases_passing_filter\t800000\t80.0
""")


# ---------------------------------------------------------------------------
# _row_key
# ---------------------------------------------------------------------------

class TestRowKey:
    def test_bare_sample_when_no_tab_or_variant(self):
        assert _row_key("HG001", None) == "HG001"

    def test_folds_tab_label(self):
        assert _row_key("HG001", "Dedup - Autosomes") == "HG001 | Dedup - Autosomes"

    def test_folds_variant_only(self):
        assert _row_key("HG001", None, "Post-filter") == "HG001 | Post-filter"

    def test_folds_tab_and_variant(self):
        assert (
            _row_key("HG001", "Dedup - Autosomes", "Post-filter")
            == "HG001 | Dedup - Autosomes - Post-filter"
        )


# ---------------------------------------------------------------------------
# _tab_label (region/stage labelling feeds the row key)
# ---------------------------------------------------------------------------

class TestTabLabel:
    def test_none_when_unset(self):
        assert _tab_label(None, None) is None

    def test_stage_and_region(self):
        assert _tab_label("autosomes", "dedup") == "Dedup - Autosomes"

    def test_unknown_values_titlecased(self):
        assert _tab_label("my_region", "my_stage") == "My Stage - My Region"


# ---------------------------------------------------------------------------
# process_inputs: table sections are flat dicts, never tabbed lists
# ---------------------------------------------------------------------------

def _args(sample_inputs, output_dir, region=None, stage=None):
    ns = argparse.Namespace()
    ns._sample_inputs = sample_inputs
    ns.output_dir = output_dir
    ns.region = region
    ns.stage = stage
    ns.demux_run_stats = None
    ns.demux_sample_stats = None
    return ns


def _write(tmp_path: Path, name: str, content: str) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


def _load(output_dir: Path, section_id: str) -> dict:
    return json.loads((output_dir / f"{section_id}_mqc.json").read_text())


def test_accuracy_is_flat_table_with_region_in_row(tmp_path: Path):
    out = tmp_path / "out"
    acc = _write(tmp_path, "acc.tsv", ACCURACY_TSV)
    args = _args(
        {"HG001": {"accuracy_metrics": [acc]}},
        out,
        region="autosomes",
        stage="dedup",
    )
    process_inputs(args)

    report = _load(out, "xoos_accuracy_metrics")
    # Table data must be a dict (flat), not a list of datasets, and must not
    # carry a data_labels selector (MultiQC tables ignore it).
    assert isinstance(report["data"], dict)
    assert "data_labels" not in report["pconfig"]
    assert list(report["data"].keys()) == ["HG001 | Dedup - Autosomes"]


def test_read_metrics_folds_region_and_variant_into_rows(tmp_path: Path):
    out = tmp_path / "out"
    rm = _write(tmp_path, "rm.tsv", READ_METRICS_TSV)
    args = _args(
        {"HG001": {"read_metrics": [rm]}},
        out,
        region="high_confidence",
        stage="dedup",
    )
    process_inputs(args)

    counts = _load(out, "xoos_read_metrics_read_counts")
    assert isinstance(counts["data"], dict)
    assert "data_labels" not in counts["pconfig"]
    assert set(counts["data"]) == {
        "HG001 | Dedup - High Confidence - Any",
        "HG001 | Dedup - High Confidence - Post-filter",
    }


def test_multiple_regions_combine_into_distinct_rows(tmp_path: Path):
    # Two regions of the same sample must produce two rows, not overwrite each
    # other (the bug that left the report showing only the first region).
    out = tmp_path / "out"
    acc = _write(tmp_path, "acc.tsv", ACCURACY_TSV)

    for region in ("autosomes", "high_confidence"):
        args = _args(
            {"HG001": {"accuracy_metrics": [acc]}},
            out,
            region=region,
            stage="dedup",
        )
        process_inputs(args)
        # Simulate per-region tasks: rename so the second run does not clobber
        # the first before the combine step would merge them.
        (out / "xoos_accuracy_metrics_mqc.json").rename(
            out / f"xoos_accuracy_metrics_{region}.json"
        )

    rows = set()
    for region in ("autosomes", "high_confidence"):
        data = json.loads((out / f"xoos_accuracy_metrics_{region}.json").read_text())["data"]
        rows.update(data.keys())
    assert rows == {
        "HG001 | Dedup - Autosomes",
        "HG001 | Dedup - High Confidence",
    }


def test_legacy_layout_without_region_uses_bare_sample(tmp_path: Path):
    out = tmp_path / "out"
    acc = _write(tmp_path, "acc.tsv", ACCURACY_TSV)
    args = _args({"HG001": {"accuracy_metrics": [acc]}}, out)
    process_inputs(args)

    report = _load(out, "xoos_accuracy_metrics")
    assert list(report["data"].keys()) == ["HG001"]
