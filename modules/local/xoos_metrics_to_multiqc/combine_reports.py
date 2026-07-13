#!/usr/bin/env python3
"""Combine multiple MultiQC JSON reports into one.

Takes directories and/or individual ``*_mqc.json`` files produced by
``xoos_to_multiqc.py`` and merges them into a single output directory.

Merge rules:
  - Files with the **same ``id``** are merged: sample-level data from all
    sources is combined into a single file. Both dict-style data
    (``{sample: metrics}``) and list-style data (``[{sample: metrics}, ...]``,
    used for multi-tab tables and boxplots) are supported.
  - Files with **different ``id``** values are written as-is.
  - If the same sample name appears in multiple inputs for the same ``id``,
    the later input overwrites the earlier one.

Usage:

    # Combine entire directories
    python3 combine_reports.py \\
        --input-dirs runA_output/ runB_output/ \\
        --output-dir combined_output/

    # Combine individual files
    python3 combine_reports.py \\
        --input-files runA/xoos_accuracy_metrics_mqc.json \\
                      runB/xoos_accuracy_metrics_mqc.json \\
        --output-dir combined_output/

    # Mix directories and files
    python3 combine_reports.py \\
        --input-dirs runA_output/ \\
        --input-files runB/xoos_accuracy_metrics_mqc.json \\
        --output-dir combined_output/

    # Then run MultiQC on the combined output
    multiqc combined_output/
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__version__ = "latest"


def _merge_data(
    existing: dict | list,
    incoming: dict | list,
    existing_labels: list[Any] | None = None,
    incoming_labels: list[Any] | None = None,
) -> tuple[dict | list, list[Any] | None]:
    """Merge incoming sample data into existing data.

    Handles both formats:
      - dict: {sample_name: {metric: value}}
      - list: [{sample_name: {metric: value}}, ...] (one per data_labels tab)

    For list-style data, tabs are merged **by label name** rather than by
    position. This is required for region/stage alignment tabs, where each
    conversion run emits a single, distinctly-named tab and the same sample
    appears across many runs. Returns ``(merged_data, merged_labels)`` where
    ``merged_labels`` is ``None`` for dict-style data.
    """
    if isinstance(existing, dict) and isinstance(incoming, dict):
        return {**existing, **incoming}, None

    if isinstance(existing, list) and isinstance(incoming, list):
        # Label-aware merge when both sides carry data_labels of matching
        # length. Align tabs by label, union the label set (existing order
        # first), and merge the per-sample dicts for shared labels.
        if (
            existing_labels is not None
            and incoming_labels is not None
            and len(existing_labels) == len(existing)
            and len(incoming_labels) == len(incoming)
        ):
            merged_map: dict[Any, dict] = {}
            order: list[Any] = []
            for label, datum in zip(existing_labels, existing):
                merged_map[label] = dict(datum) if isinstance(datum, dict) else datum
                order.append(label)
            for label, datum in zip(incoming_labels, incoming):
                if label in merged_map and isinstance(merged_map[label], dict) and isinstance(datum, dict):
                    merged_map[label].update(datum)
                else:
                    merged_map[label] = dict(datum) if isinstance(datum, dict) else datum
                    if label not in order:
                        order.append(label)
            return [merged_map[label] for label in order], order

        # Fallback: positional element-wise merge (no/uneven labels).
        max_len = max(len(existing), len(incoming))
        merged = []
        for i in range(max_len):
            e = existing[i] if i < len(existing) else {}
            inc = incoming[i] if i < len(incoming) else {}
            merged.append({**e, **inc})
        merged_labels = _merge_data_labels(existing_labels or [], incoming_labels or [])
        return merged, (merged_labels or None)

    # Mismatched types — incoming wins with a warning
    logger.warning(
        "Data type mismatch during merge (existing=%s, incoming=%s); using incoming",
        type(existing).__name__,
        type(incoming).__name__,
    )
    return incoming, incoming_labels


def _merge_data_labels(
    existing_labels: list[Any],
    incoming_labels: list[Any],
) -> list[Any]:
    """Merge positional data_labels lists, keeping the longer as superset.

    Used only for the positional fallback path. Label-aware merging in
    ``_merge_data`` returns its own ordered label list.
    """
    if len(incoming_labels) <= len(existing_labels):
        return existing_labels
    return existing_labels + incoming_labels[len(existing_labels):]


def _ingest_file(
    json_file: Path,
    merged: dict[str, dict[str, Any]],
) -> None:
    """Read a single *_mqc.json file and merge it into the accumulated state."""
    with json_file.open() as f:
        report = json.load(f)

    report_id = report.get("id")
    if report_id is None:
        logger.warning("Skipping file without id: %s", json_file)
        return

    if report_id in merged:
        existing = merged[report_id]["report"]
        existing_labels = existing.get("pconfig", {}).get("data_labels", []) or None
        incoming_labels = report.get("pconfig", {}).get("data_labels", []) or None

        merged_data, merged_labels = _merge_data(
            existing["data"], report["data"], existing_labels, incoming_labels,
        )
        existing["data"] = merged_data

        if merged_labels is not None:
            existing.setdefault("pconfig", {})["data_labels"] = merged_labels
        elif incoming_labels and not existing_labels:
            existing.setdefault("pconfig", {})["data_labels"] = incoming_labels

        logger.info("Merged %s into existing report id=%s", json_file, report_id)
    else:
        merged[report_id] = {
            "report": report,
            "filename": json_file.name,
        }
        logger.info("Added report id=%s from %s", report_id, json_file)


def combine_reports(
    output_dir: Path,
    input_dirs: list[Path] | None = None,
    input_files: list[Path] | None = None,
) -> None:
    """Read *_mqc.json files from directories and/or individual paths and merge by id."""
    output_dir.mkdir(parents=True, exist_ok=True)

    merged: dict[str, dict[str, Any]] = {}

    # Process directories
    for input_dir in (input_dirs or []):
        if not input_dir.is_dir():
            logger.warning("Skipping non-directory: %s", input_dir)
            continue
        for json_file in sorted(input_dir.glob("*_mqc.json")):
            _ingest_file(json_file, merged)

    # Process individual files
    for input_file in (input_files or []):
        if not input_file.exists():
            logger.warning("File not found: %s", input_file)
            continue
        _ingest_file(input_file, merged)

    # Write merged reports
    files_written = 0
    for report_id, entry in merged.items():
        out_file = output_dir / entry["filename"]
        with out_file.open("w") as f:
            json.dump(entry["report"], f, indent=2)
        logger.info("Wrote %s", out_file)
        files_written += 1

    if files_written == 0:
        logger.warning("No *_mqc.json files found in inputs.")
    else:
        logger.info("Wrote %d combined report(s) to %s", files_written, output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combine multiple MultiQC JSON reports into one.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )
    parser.add_argument(
        "--input-dirs",
        type=Path,
        nargs="+",
        default=[],
        metavar="DIR",
        help="Directories containing *_mqc.json files to combine.",
    )
    parser.add_argument(
        "--input-files",
        type=Path,
        nargs="+",
        default=[],
        metavar="FILE",
        help="Individual *_mqc.json files to combine.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        metavar="DIR",
        help="Output directory for combined reports.",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )

    args = parser.parse_args()

    if not args.input_dirs and not args.input_files:
        parser.error("At least one of --input-dirs or --input-files is required.")

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    combine_reports(args.output_dir, args.input_dirs, args.input_files)


if __name__ == "__main__":
    main()
