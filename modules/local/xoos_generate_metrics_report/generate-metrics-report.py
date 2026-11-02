#!/usr/bin/env python3

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass
import gzip
import io
import json
import logging
from pathlib import Path
import re
from typing import Any, Optional

__version__ = "latest"


def setup_logging(log_level: str) -> None:
    """Setup logging configuration."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


RUN_COLUMN_NAME = "Entire Run"

SORT_ORDER_RUN_STATS = 0
SORT_ORDER_SAMPLE_STATS = 1
SORT_ORDER_ALIGNMENT_METRICS = 2
SORT_ORDER_READ_COLLAPSER = 3
SORT_ORDER_VCF_EVAL = 4
SORT_ORDER_GC_BIAS = 5
SORT_ORDER_HS_METRICS = 6
SORT_ORDER_CONTAMINATION = 7
SORT_ORDER_TUMOR_FRACTION = 8
SORT_ORDER_STR = 9
SORT_ORDER_CNC = 10
SORT_ORDER_SOMALIER = 11
SORT_ORDER_SVC = 12


@dataclass
class Metric:
    name: str
    value: str
    sort_order: int
    sample_name: str = RUN_COLUMN_NAME


def create_dict_reader(
    file_handle: io.TextIOBase, delimiter: str = "\t"
) -> csv.DictReader[str]:
    return csv.DictReader(
        (line for line in file_handle if line.strip() and not line.startswith("#")),
        delimiter=delimiter,
    )


def parse_run_stats(input_file: Path) -> list[Metric]:
    logging.info(f"Parsing run metrics from: {input_file}")

    run_metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            metric_name = record["metric_name"]
            value = record["value"]
            percentage = record.get("percentage", "NA")
            run_metrics.append(
                Metric(
                    name=f"demux/{metric_name}",
                    value=value,
                    sort_order=SORT_ORDER_RUN_STATS,
                )
            )
            run_metrics.append(
                Metric(
                    name=f"demux/{metric_name}%",
                    value=percentage,
                    sort_order=SORT_ORDER_RUN_STATS,
                )
            )

    return run_metrics


def parse_sample_stats(input_file: Path) -> list[Metric]:
    logging.info(f"Parsing sample metrics from: {input_file}")

    sample_metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            metric_name = record["metric"]
            for k, v in record.items():
                if k == "metric":
                    continue
                sample_metrics.append(
                    Metric(
                        name=f"demux/{metric_name}",
                        sample_name=k,
                        value=v,
                        sort_order=SORT_ORDER_SAMPLE_STATS,
                    )
                )

    return sample_metrics


def read_collapser_metrics(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    logging.info(f"Parsing read_collapser metrics from: {input_file}")

    metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            metric_name = record.get("metric_name")
            for k, v in record.items():
                if k in ("metric_name", "denominator"):
                    continue

                if k == "value":
                    suffix = ""
                elif k == "percentage":
                    suffix = "%"
                else:
                    suffix = f"/{k}"

                metrics.append(
                    Metric(
                        name=f"{source}/{metric_name}{suffix}",
                        value=v,
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_READ_COLLAPSER,
                    )
                )

    return metrics


def parse_svc_metrics(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    """Parse small variant caller filter_variants QC metrics TSV."""
    logging.info(f"Parsing small variant caller metrics from: {input_file}")

    metrics: list[Metric] = []
    with input_file.open("r", newline="") as f:
        reader = create_dict_reader(f)
        for record in reader:
            metric_name = record.get("metric_name")
            if not metric_name:
                continue

            for column, suffix in (("value", ""), ("percentage", "%")):
                value = record.get(column)
                if value is not None:
                    metrics.append(
                        Metric(
                            name=f"{source}/{metric_name}{suffix}",
                            value=value,
                            sample_name=sample_name,
                            sort_order=SORT_ORDER_SVC,
                        )
                    )

    return metrics


def parse_alignment_metrics(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    logging.info(f"Parsing alignment metrics from: {input_file}")

    category = input_file.parent.name

    alignment_metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            metric_name = record.get("metric_name")
            if metric_name is None:
                metric_name = record.get("type")
            for k, v in record.items():
                if k in ("metric_name", "type", "denominator"):
                    continue

                if k == "count":
                    suffix = ""
                elif k == "percentage":
                    suffix = "%"
                else:
                    suffix = f"/{k}"

                alignment_metrics.append(
                    Metric(
                        name=f"{source}/{category}/{metric_name}{suffix}",
                        value=v,
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_ALIGNMENT_METRICS,
                    )
                )

    return alignment_metrics


def parse_alignment_metrics_hp_errors(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    logging.info(f"Parsing alignment metrics hp errors from: {input_file}")

    def determine_alignment_metrics_category(path: Path) -> str:
        for part in reversed(path.parts):
            if part in ("coverage", "accuracy", "read_metrics"):
                return part
        return "unknown"

    category = determine_alignment_metrics_category(input_file)
    alignment_metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            hp_base = record.get("hp_base")
            hp_length = record.get("hp_length")
            for k, v in record.items():
                if k in ("hp_base", "hp_length"):
                    continue
                alignment_metrics.append(
                    Metric(
                        name=f"{source}/{category}/hp_errors/{hp_base}/{hp_length}/{k}",
                        value=v,
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_ALIGNMENT_METRICS,
                    )
                )
    return alignment_metrics


def parse_vcf_eval_metrics(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    logging.info(f"Parsing vcf eval metrics from: {input_file}")

    if not input_file.exists():
        logging.warning(f"VCF eval metrics file does not exist: {input_file}")
        return []

    vcf_eval_metrics: list[Metric] = []
    with gzip.open(input_file, "rt") as f:
        reader = csv.reader(
            (line for line in f if not line.startswith("#")), delimiter="\t"
        )

        records = list(reader)
        metric_definitions = [
            (1, "true_positives_baseline"),
            (2, "false_positives"),
            (3, "true_positives_call"),
            (4, "false_negatives"),
            (5, "precision"),
            (6, "sensitivity"),
            (7, "f_measure"),
        ]

        if not records:
            logging.warning(
                f"VCF eval metrics file is empty or contains no score rows: {input_file}"
            )
            metrics = [
                Metric(
                    name=f"{source}/{name}",
                    value="NA",
                    sample_name=sample_name,
                    sort_order=SORT_ORDER_VCF_EVAL,
                )
                for _, name in metric_definitions
            ]
        else:
            record = records[-1]
            metrics = [
                Metric(
                    name=f"{source}/{name}",
                    value=record[idx],
                    sample_name=sample_name,
                    sort_order=SORT_ORDER_VCF_EVAL,
                )
                for idx, name in metric_definitions
            ]
        vcf_eval_metrics.extend(metrics)

    return vcf_eval_metrics


def parse_gc_bias_summary(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    logging.info(f"Parsing GC bias summary from: {input_file}")

    gc_bias_metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            if record.get("READS_USED") == "UNIQUE":
                metric_fields = [
                    "AT_DROPOUT",
                    "GC_DROPOUT",
                    "GC_NC_0_19",
                    "GC_NC_20_39",
                    "GC_NC_40_59",
                    "GC_NC_60_79",
                    "GC_NC_80_100",
                ]
                for field in metric_fields:
                    metric_name = field.lower()
                    gc_bias_metrics.append(
                        Metric(
                            name=f"{source}/{metric_name}",
                            value=record[field],
                            sample_name=sample_name,
                            sort_order=SORT_ORDER_GC_BIAS,
                        )
                    )
                break

    return gc_bias_metrics


def parse_gc_bias_detail(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    logging.info(f"Parsing GC bias detail from: {input_file}")

    bin_range = (20, 70)
    normalized_coverage_range = (0.8, 1.2)
    bins_in_range = 0

    gc_bias_metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            if record.get("READS_USED") == "UNIQUE":
                gc = record["GC"]
                normalized_coverage = record["NORMALIZED_COVERAGE"]
                gc_bias_metrics.append(
                    Metric(
                        name=f"{source}/{gc}/normalized_coverage",
                        value=normalized_coverage,
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_GC_BIAS,
                    )
                )

                if bin_range[0] <= int(gc) <= bin_range[1]:
                    if (
                        normalized_coverage_range[0]
                        <= float(normalized_coverage)
                        <= normalized_coverage_range[1]
                    ):
                        bins_in_range += 1
    total_bins = bin_range[1] - bin_range[0] + 1
    gc_bias_metrics.append(
        Metric(
            name=f"{source}/bins_20-70_normalized_coverage_0.8-1.2%",
            value=f"{(bins_in_range / total_bins) * 100:.2f}",
            sample_name=sample_name,
            sort_order=SORT_ORDER_GC_BIAS,
        )
    )

    return gc_bias_metrics


def parse_alignment_metrics_gc_bias_summary(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    """Parse the post-filter summary values from alignment_metrics gcbias."""
    logging.info(f"Parsing alignment_metrics GC bias summary from: {input_file}")

    metric_name_map = {
        "post_filter_at_dropout": "at_dropout",
        "post_filter_gc_dropout": "gc_dropout",
        "post_filter_first_quintile_normalized_coverage": "gc_nc_0_19",
        "post_filter_second_quintile_normalized_coverage": "gc_nc_20_39",
        "post_filter_third_quintile_normalized_coverage": "gc_nc_40_59",
        "post_filter_fourth_quintile_normalized_coverage": "gc_nc_60_79",
        "post_filter_fifth_quintile_normalized_coverage": "gc_nc_80_100",
    }
    gc_bias_metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            metric_name = record.get("metric_name")
            report_metric_name = metric_name_map.get(metric_name)
            if report_metric_name is not None:
                gc_bias_metrics.append(
                    Metric(
                        name=f"{source}/{report_metric_name}",
                        value=record["value"],
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_GC_BIAS,
                    )
                )

    return gc_bias_metrics


def parse_alignment_metrics_gc_bias_detail(
    sample_name: str, source: str, input_file: Path
) -> list[Metric]:
    """Parse post-filter GC coverage points from alignment_metrics gcbias."""
    logging.info(f"Parsing alignment_metrics GC bias detail from: {input_file}")

    bin_range = (20, 70)
    normalized_coverage_range = (0.8, 1.2)
    bins_in_range = 0
    gc_bias_metrics: list[Metric] = []

    with input_file.open("r") as f:
        reader = create_dict_reader(f)
        for record in reader:
            gc = record.get("gc_percentage")
            normalized_coverage = record.get("normalized_post_filter_coverage")
            if gc is None or normalized_coverage is None:
                continue
            try:
                gc_value = int(gc)
                normalized_coverage_value = float(normalized_coverage)
            except ValueError:
                continue

            gc_bias_metrics.append(
                Metric(
                    name=f"{source}/{gc}/normalized_coverage",
                    value=normalized_coverage,
                    sample_name=sample_name,
                    sort_order=SORT_ORDER_GC_BIAS,
                )
            )
            if (
                bin_range[0] <= gc_value <= bin_range[1]
                and normalized_coverage_range[0]
                <= normalized_coverage_value
                <= normalized_coverage_range[1]
            ):
                bins_in_range += 1

    total_bins = bin_range[1] - bin_range[0] + 1
    gc_bias_metrics.append(
        Metric(
            name=f"{source}/bins_20-70_normalized_coverage_0.8-1.2%",
            value=f"{(bins_in_range / total_bins) * 100:.2f}",
            sample_name=sample_name,
            sort_order=SORT_ORDER_GC_BIAS,
        )
    )

    return gc_bias_metrics


def parse_hs_metrics(sample_name: str, source: str, input_file: Path) -> list[Metric]:
    logging.info(f"Parsing HsMetrics from: {input_file}")

    lines: list[str] = []
    in_hs_metrics_section = False
    with input_file.open("r") as f:
        for line in f:
            if line.startswith("## METRICS CLASS\tpicard.analysis.directed.HsMetrics"):
                in_hs_metrics_section = True
                continue
            if in_hs_metrics_section and not line.strip():
                in_hs_metrics_section = False
                break
            if in_hs_metrics_section:
                lines.append(line)

    hs_metrics: list[Metric] = []
    reader = csv.DictReader(lines, delimiter="\t")
    for record in reader:
        bait_set = record.get("BAIT_SET", "unknown")
        for k, v in record.items():
            if k == "BAIT_SET":
                continue
            hs_metrics.append(
                Metric(
                    name=f"{source}/{bait_set}/{k}",
                    value=v,
                    sample_name=sample_name,
                    sort_order=SORT_ORDER_HS_METRICS,
                )
            )

    return hs_metrics


def parse_contamination_metrics(
    sample_name: str,
    input_file: Path,
) -> list[Metric]:
    """Parse TFE contamination detection results TSV.

    The file is a 2-column TSV (metric_name, value) with a header row.
    Lines starting with ## are metadata comments.
    """
    metrics: list[Metric] = []
    with input_file.open("r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("#") or not stripped:
                continue
            parts = stripped.split("\t")
            # Skip the header row
            if parts[0] == "metric_name":
                continue
            if len(parts) >= 2:
                metrics.append(
                    Metric(
                        name=f"contamination_detection/{parts[0]}",
                        value=parts[1],
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_CONTAMINATION,
                    )
                )
    return metrics


def parse_tumor_fraction_metrics(
    sample_name: str,
    input_file: Path,
) -> list[Metric]:
    """Parse TFE tumor fraction estimation results TSV.

    The file is a 2-column TSV (metric_name, value) with a header row.
    Lines starting with ## are metadata comments.
    """
    metrics: list[Metric] = []
    with input_file.open("r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("#") or not stripped:
                continue
            parts = stripped.split("\t")
            # Skip the header row
            if parts[0] == "metric_name":
                continue
            if len(parts) >= 2:
                metrics.append(
                    Metric(
                        name=f"tumor_fraction/{parts[0]}",
                        value=parts[1],
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_TUMOR_FRACTION,
                    )
                )
    return metrics


def parse_str_metrics(
    sample_name: str,
    input_file: Path,
) -> list[Metric]:
    """Parse STR caller metrics TSV."""
    metrics: list[Metric] = []
    with input_file.open("r") as f:
        for line in f:
            if line.startswith("##") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                metrics.append(
                    Metric(
                        name=f"str_caller/{parts[0]}",
                        value=parts[1],
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_STR,
                    )
                )
    return metrics


def parse_cnc_metrics(
    sample_name: str,
    input_file: Path,
) -> list[Metric]:
    """Parse copy number caller metrics TSV."""
    metrics: list[Metric] = []
    with input_file.open("r") as f:
        for line in f:
            if line.startswith("##") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                name, value = parts[0], parts[1]
                if name == "Metric":
                    continue
                normalized = name.replace(" ", "_").lower()
                metrics.append(
                    Metric(
                        name=f"copy_number_caller/{normalized}",
                        value=value,
                        sample_name=sample_name,
                        sort_order=SORT_ORDER_CNC,
                    )
                )
    return metrics


def _read_somalier_tsv(file_handle: io.TextIOBase) -> csv.DictReader[str]:
    """Read a somalier TSV file whose header starts with '#'.

    Strips the leading '#' from the header so DictReader field names are clean
    (e.g. 'sample_a' instead of '#sample_a').
    """
    lines: list[str] = []
    for line in file_handle:
        stripped = line.strip()
        if not stripped:
            continue
        lines.append(stripped)

    if not lines:
        return csv.DictReader(iter([]), delimiter="\t")

    # Strip leading '#' from header row
    if lines[0].startswith("#"):
        lines[0] = lines[0][1:]

    return csv.DictReader(lines, delimiter="\t")


def parse_somalier_pair_metrics(
    key_name: str,
    input_file: Path,
) -> list[Metric]:
    """Parse somalier relate pairs.tsv.

    The file is a TSV with a header row (starting with #sample_a).
    Each data row represents a pairwise comparison between two samples.

    Every column is reported as-is under the key_name column (the T/N group),
    keeping the original header names including the _a/_b suffixes. The
    sample_a and sample_b columns are emitted as metrics too, so the pair
    membership is readable straight off the report.
    """
    logging.info(f"Parsing somalier pair metrics from: {input_file}")

    metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = _read_somalier_tsv(f)
        for record in reader:
            for k, v in record.items():
                metrics.append(
                    Metric(
                        name=f"somalier/{k}",
                        value=v,
                        sample_name=key_name,
                        sort_order=SORT_ORDER_SOMALIER,
                    )
                )
    return metrics


def parse_somalier_sample_metrics(
    key_name: str,
    input_file: Path,
    sample_ids: Optional[list[str]] = None,
) -> list[Metric]:
    """Parse somalier relate samples.tsv.

    The file is a TSV with a header row (starting with #family_id).
    Each row contains per-sample QC stats (depth, het/hom counts, etc.).
    Metrics are assigned to the per-sample column derived from sample_id.
    If sample_id is missing in a row, the provided key_name is used as
    a fallback output column.

    A samples.tsv produced by `somalier relate` on a T/N pair contains rows
    for both samples, which may belong to different sequencing runs. When
    sample_ids is supplied, only rows for those samples are emitted, so a
    per-run invocation reports only the samples belonging to that run.
    """
    logging.info(f"Parsing somalier sample metrics from: {input_file}")

    run_sample_ids = set(sample_ids) if sample_ids else None

    metrics: list[Metric] = []
    with input_file.open("r") as f:
        reader = _read_somalier_tsv(f)
        for record in reader:
            sample_id = record.get("sample_id") or record.get("family_id") or key_name
            if run_sample_ids is not None and sample_id not in run_sample_ids:
                continue
            for k, v in record.items():
                if k in ("family_id", "sample_id", "paternal_id", "maternal_id"):
                    continue
                metrics.append(
                    Metric(
                        name=f"somalier/{k}",
                        value=v,
                        sample_name=sample_id,
                        sort_order=SORT_ORDER_SOMALIER,
                    )
                )
    return metrics


def merge_metrics_by_name(metrics: list[Metric]) -> list[dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = defaultdict(dict)
    for metric in metrics:
        merged[metric.name].update({metric.sample_name: metric.value})
    records: list[dict[str, Any]] = []
    for name, values in merged.items():
        record = {"metric_name": name}
        record.update(values)
        records.append(record)
    return records


def write_json_output(data: Any, output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w") as f:
        json.dump(data, f, indent=2, default=lambda k: k.__dict__)

    logging.info(f"Successfully wrote JSON output to: {output_file}")


def cmd_parse_demux_metrics(args: argparse.Namespace) -> None:
    data: list[Metric] = []
    data.extend(parse_run_stats(args.run_stats))
    data.extend(parse_sample_stats(args.sample_stats))
    write_json_output(data, args.output)

    logging.info("Successfully processed demux metrics")


def cmd_parse_read_collapser_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        if not input_file.exists():
            logging.warning(f"Read collapser metrics file does not exist: {input_file}")
            continue
        data: list[Metric] = read_collapser_metrics(
            args.sample_name,
            args.source,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed read collapser metrics")


def cmd_parse_svc_metrics(args: argparse.Namespace) -> None:
    """Parse small variant caller filter_variants metrics and write JSON."""
    output = parse_svc_metrics(args.sample_name, args.source, args.input)
    write_json_output(output, args.output)

    logging.info("Successfully processed small variant caller metrics")


def cmd_parse_alignment_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        if not input_file.exists():
            logging.warning(f"Alignment metrics file does not exist: {input_file}")
            continue
        data: list[Metric] = parse_alignment_metrics(
            args.sample_name,
            args.source,
            input_file,
        )
        output.extend(data)
    for input_file in args.input_hp_errors:
        if not input_file.exists():
            logging.warning(f"HP errors file does not exist: {input_file}")
            continue
        data: list[Metric] = parse_alignment_metrics_hp_errors(
            args.sample_name,
            args.source,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed alignment metrics")


def cmd_parse_vcf_eval_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        data: list[Metric] = parse_vcf_eval_metrics(
            args.sample_name,
            args.source,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed vcf eval metrics")


def cmd_parse_gc_bias_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []

    # Parse summary file
    if args.summary.exists():
        data: list[Metric] = parse_gc_bias_summary(
            args.sample_name,
            args.source,
            args.summary,
        )
        output.extend(data)
    else:
        logging.warning(f"GC bias summary file does not exist: {args.summary}")

    # Parse detail file
    if args.detail.exists():
        data: list[Metric] = parse_gc_bias_detail(
            args.sample_name,
            args.source,
            args.detail,
        )
        output.extend(data)
    else:
        logging.warning(f"GC bias detail file does not exist: {args.detail}")

    write_json_output(output, args.output)

    logging.info("Successfully processed GC bias metrics")


def cmd_parse_alignment_metrics_gc_bias(args: argparse.Namespace) -> None:
    output: list[Metric] = []

    if args.summary.exists():
        output.extend(
            parse_alignment_metrics_gc_bias_summary(
                args.sample_name,
                args.source,
                args.summary,
            )
        )
    else:
        logging.warning(f"GC bias summary file does not exist: {args.summary}")

    if args.detail.exists():
        output.extend(
            parse_alignment_metrics_gc_bias_detail(
                args.sample_name,
                args.source,
                args.detail,
            )
        )
    else:
        logging.warning(f"GC bias detail file does not exist: {args.detail}")

    write_json_output(output, args.output)

    logging.info("Successfully processed alignment_metrics GC bias metrics")


def cmd_parse_hs_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        if not input_file.exists():
            logging.warning(f"HsMetrics file does not exist: {input_file}")
            continue
        data: list[Metric] = parse_hs_metrics(
            args.sample_name,
            args.source,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed HsMetrics")


def cmd_parse_contamination_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        if not input_file.exists():
            logging.warning(f"Contamination results file does not exist: {input_file}")
            continue
        data: list[Metric] = parse_contamination_metrics(
            args.sample_name,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed contamination detection metrics")


def cmd_parse_tumor_fraction_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        if not input_file.exists():
            logging.warning(
                f"Tumor fraction results file does not exist: {input_file}"
            )
            continue
        data: list[Metric] = parse_tumor_fraction_metrics(
            args.sample_name,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed tumor fraction estimation metrics")


def cmd_parse_str_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        if not input_file.exists():
            logging.warning(f"STR metrics file does not exist: {input_file}")
            continue
        data: list[Metric] = parse_str_metrics(
            args.sample_name,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed STR caller metrics")


def cmd_parse_cnc_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    for input_file in args.input:
        if not input_file.exists():
            logging.warning(f"CNC metrics file does not exist: {input_file}")
            continue
        data: list[Metric] = parse_cnc_metrics(
            args.sample_name,
            input_file,
        )
        output.extend(data)
    write_json_output(output, args.output)

    logging.info("Successfully processed copy number caller metrics")


def cmd_parse_somalier_metrics(args: argparse.Namespace) -> None:
    output: list[Metric] = []
    if args.pairs_tsv and args.pairs_tsv.exists():
        output.extend(parse_somalier_pair_metrics(args.key_name, args.pairs_tsv))
    elif args.pairs_tsv:
        logging.warning(f"Somalier pairs.tsv file does not exist: {args.pairs_tsv}")

    for samples_tsv in (args.samples_tsv or []):
        if samples_tsv.exists():
            output.extend(
                parse_somalier_sample_metrics(
                    args.key_name, samples_tsv, args.sample_ids
                )
            )
        else:
            logging.warning(f"Somalier samples.tsv file does not exist: {samples_tsv}")

    write_json_output(output, args.output)

    logging.info("Successfully processed somalier metrics")


def determine_fieldnames(metrics: list[dict[str, Any]]) -> list[str]:
    fieldnames: set[str] = set()
    for metric in metrics:
        fieldnames.update(
            k for k in metric.keys() if k not in ("metric_name", RUN_COLUMN_NAME)
        )
    return ["metric_name", RUN_COLUMN_NAME] + sorted(fieldnames)


def cmd_generate_report(args: argparse.Namespace) -> None:
    all_metrics: list[Metric] = []
    if args.input is not None:
        for input_file in args.input:
            with input_file.open("r") as f:
                data = json.load(f)
                metrics = [Metric(**item) for item in data]
                all_metrics.extend(metrics)

    all_metrics.sort(key=lambda m: (m.sort_order, m.name, m.sample_name))

    # Compile regex patterns for filtering
    patterns = [re.compile(pattern) for pattern in (args.patterns or [".*"])]

    # Filter metrics based on patterns
    filtered_metrics = [
        metric
        for metric in all_metrics
        if any(pattern.match(metric.name) for pattern in patterns)
    ]

    logging.info(
        f"Filtered {len(filtered_metrics)} metrics from {len(all_metrics)} total metrics"
    )

    merged_records = merge_metrics_by_name(filtered_metrics)
    with args.output.open("w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=determine_fieldnames(merged_records),
            delimiter="\t",
            restval="NA",
        )
        writer.writeheader()
        writer.writerows(merged_records)

    logging.info("Successfully generated consolidated metrics report")


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate metrics reports from various input formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )

    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level (default: INFO)",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        help="Available commands",
        required=True,
    )

    parse_demux_parser = subparsers.add_parser(
        "parse-demux-metrics",
        help="Parse demux metrics TSV file and generate JSON report",
    )
    parse_demux_parser.add_argument(
        "--run-stats",
        type=Path,
        help="Input TSV file containing demux run metrics",
        required=True,
    )
    parse_demux_parser.add_argument(
        "--sample-stats",
        type=Path,
        help="Input TSV file containing demux sample metrics",
        required=True,
    )
    parse_demux_parser.add_argument(
        "--output", type=Path, help="Output JSON file path", required=True
    )
    parse_demux_parser.set_defaults(func=cmd_parse_demux_metrics)

    parse_read_collapser_parser = subparsers.add_parser(
        "parse-read-collapser-metrics",
        help="Parse read collapser metrics TSV file and generate JSON report",
    )
    parse_read_collapser_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the read collapser metrics",
        required=True,
    )
    parse_read_collapser_parser.add_argument(
        "--source",
        type=str,
        help="Source of the read collapser metrics",
        required=True,
    )
    parse_read_collapser_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input TSV file containing read collapser metrics",
        required=True,
    )
    parse_read_collapser_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_read_collapser_parser.set_defaults(func=cmd_parse_read_collapser_metrics)

    parse_svc_parser = subparsers.add_parser(
        "parse-svc-metrics",
        help="Parse small variant caller filter_variants metrics TSV and generate JSON report",
    )
    parse_svc_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the small variant caller metrics",
        required=True,
    )
    parse_svc_parser.add_argument(
        "--source",
        type=str,
        help="Source of the small variant caller metrics",
        required=True,
    )
    parse_svc_parser.add_argument(
        "--input",
        type=Path,
        help="Input small variant caller metrics TSV file",
        required=True,
    )
    parse_svc_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_svc_parser.set_defaults(func=cmd_parse_svc_metrics)

    parse_alignment_metrics_parser = subparsers.add_parser(
        "parse-alignment-metrics",
        help="Parse alignment metrics TSV file and generate JSON report",
    )
    parse_alignment_metrics_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the alignment metrics",
        required=True,
    )
    parse_alignment_metrics_parser.add_argument(
        "--source",
        type=str,
        help="Source of the alignment metrics",
        required=True,
    )
    parse_alignment_metrics_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input TSV file containing alignment metrics",
        required=True,
    )
    parse_alignment_metrics_parser.add_argument(
        "--input-hp-errors",
        type=Path,
        nargs="+",
        action="extend",
        help="Input TSV file containing alignment hp error metrics",
        required=True,
    )
    parse_alignment_metrics_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_alignment_metrics_parser.set_defaults(func=cmd_parse_alignment_metrics)

    parse_vcf_eval_metrics_parser = subparsers.add_parser(
        "parse-vcf-eval-metrics",
        help="Parse VCF eval metrics TSV file and generate JSON report",
    )
    parse_vcf_eval_metrics_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the VCF eval metrics",
        required=True,
    )
    parse_vcf_eval_metrics_parser.add_argument(
        "--source",
        type=str,
        help="Source of the VCF eval metrics",
        required=True,
    )
    parse_vcf_eval_metrics_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input TSV file containing VCF eval metrics",
        required=True,
    )
    parse_vcf_eval_metrics_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_vcf_eval_metrics_parser.set_defaults(func=cmd_parse_vcf_eval_metrics)

    parse_gc_bias_parser = subparsers.add_parser(
        "parse-gc-bias",
        help="Parse GC bias metrics files and generate JSON report",
    )
    parse_gc_bias_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the GC bias metrics",
        required=True,
    )
    parse_gc_bias_parser.add_argument(
        "--source",
        type=str,
        help="Source of the GC bias metrics",
        required=True,
    )
    parse_gc_bias_parser.add_argument(
        "--summary",
        type=Path,
        help="Input TSV file containing GC bias summary metrics",
        required=True,
    )
    parse_gc_bias_parser.add_argument(
        "--detail",
        type=Path,
        help="Input TSV file containing GC bias detail metrics",
        required=True,
    )
    parse_gc_bias_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_gc_bias_parser.set_defaults(func=cmd_parse_gc_bias_metrics)

    parse_alignment_metrics_gc_bias_parser = subparsers.add_parser(
        "parse-alignment-metrics-gc-bias",
        help="Parse alignment_metrics GC bias files and generate JSON report",
    )
    parse_alignment_metrics_gc_bias_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the GC bias metrics",
        required=True,
    )
    parse_alignment_metrics_gc_bias_parser.add_argument(
        "--source",
        type=str,
        help="Source of the GC bias metrics",
        required=True,
    )
    parse_alignment_metrics_gc_bias_parser.add_argument(
        "--summary",
        type=Path,
        help="Input TSV file containing GC bias summary metrics",
        required=True,
    )
    parse_alignment_metrics_gc_bias_parser.add_argument(
        "--detail",
        type=Path,
        help="Input TSV file containing GC bias detail metrics",
        required=True,
    )
    parse_alignment_metrics_gc_bias_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_alignment_metrics_gc_bias_parser.set_defaults(
        func=cmd_parse_alignment_metrics_gc_bias
    )

    parse_hs_metrics_parser = subparsers.add_parser(
        "parse-hs-metrics",
        help="Parse Picard CollectHsMetrics TSV file and generate JSON report",
    )
    parse_hs_metrics_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the HsMetrics",
        required=True,
    )
    parse_hs_metrics_parser.add_argument(
        "--source",
        type=str,
        help="Source of the HsMetrics",
        required=True,
    )
    parse_hs_metrics_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input TSV file containing HsMetrics",
        required=True,
    )
    parse_hs_metrics_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_hs_metrics_parser.set_defaults(func=cmd_parse_hs_metrics)

    parse_contamination_parser = subparsers.add_parser(
        "parse-contamination-metrics",
        help="Parse contamination detection (TFE) results CSV and generate JSON report",
    )
    parse_contamination_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the contamination metrics",
        required=True,
    )
    parse_contamination_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input contamination results TSV file (contamination_estimation_results.tsv)",
        required=True,
    )
    parse_contamination_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_contamination_parser.set_defaults(func=cmd_parse_contamination_metrics)

    parse_tumor_fraction_parser = subparsers.add_parser(
        "parse-tumor-fraction-metrics",
        help="Parse tumor fraction estimation (TFE) results TSV and generate JSON report",
    )
    parse_tumor_fraction_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the tumor fraction metrics",
        required=True,
    )
    parse_tumor_fraction_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input tumor fraction results TSV file (tumor_fraction_estimation_results.tsv)",
        required=True,
    )
    parse_tumor_fraction_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_tumor_fraction_parser.set_defaults(func=cmd_parse_tumor_fraction_metrics)

    parse_str_parser = subparsers.add_parser(
        "parse-str-metrics",
        help="Parse STR caller metrics TSV and generate JSON report",
    )
    parse_str_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the STR metrics",
        required=True,
    )
    parse_str_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input STR caller metrics TSV file",
        required=True,
    )
    parse_str_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_str_parser.set_defaults(func=cmd_parse_str_metrics)

    parse_cnc_parser = subparsers.add_parser(
        "parse-cnc-metrics",
        help="Parse copy number caller metrics TSV and generate JSON report",
    )
    parse_cnc_parser.add_argument(
        "--sample-name",
        type=str,
        help="Sample name associated with the CNC metrics",
        required=True,
    )
    parse_cnc_parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        action="extend",
        help="Input copy number caller metrics TSV file",
        required=True,
    )
    parse_cnc_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_cnc_parser.set_defaults(func=cmd_parse_cnc_metrics)

    parse_somalier_parser = subparsers.add_parser(
        "parse-somalier-metrics",
        help="Parse somalier relate pairs.tsv and samples.tsv and generate JSON report",
    )
    parse_somalier_parser.add_argument(
        "--key-name",
        type=str,
        help=(
            "Key the metrics are reported under: the T/N group id for pair "
            "metrics, or the run name for sample metrics."
        ),
        required=True,
    )
    parse_somalier_parser.add_argument(
        "--pairs-tsv",
        type=Path,
        help="Somalier relate pairs.tsv file",
        default=None,
    )
    parse_somalier_parser.add_argument(
        "--samples-tsv",
        type=Path,
        nargs='*',
        help="Somalier relate samples.tsv file(s)",
        default=None,
    )
    parse_somalier_parser.add_argument(
        "--sample-ids",
        type=str,
        nargs='*',
        help=(
            "Restrict samples.tsv parsing to these sample ids. Used to report "
            "only the samples belonging to a given sequencing run when a T/N "
            "pair spans multiple runs."
        ),
        default=None,
    )
    parse_somalier_parser.add_argument(
        "--output",
        type=Path,
        help="Output JSON file path",
        required=True,
    )
    parse_somalier_parser.set_defaults(func=cmd_parse_somalier_metrics)

    generate_report_parser = subparsers.add_parser(
        "generate-report",
        help="Generate consolidated metrics report from JSON inputs",
    )
    generate_report_parser.add_argument(
        "--input",
        type=Path,
        nargs="*",
        action="extend",
        help="Input JSON files for the report",
    )
    generate_report_parser.add_argument(
        "--pattern",
        type=str,
        action="append",
        dest="patterns",
        default=None,
        help="Regex pattern to filter metrics (can be specified multiple times, default: ['.*'] includes all)",
    )
    generate_report_parser.add_argument(
        "--output",
        type=Path,
        help="Output report file path",
        required=True,
    )
    generate_report_parser.set_defaults(func=cmd_generate_report)

    return parser


def main() -> None:
    parser = create_parser()
    args = parser.parse_args()

    setup_logging(args.log_level)

    logging.info(f"Starting generate-metrics-report with command: {args.command}")

    args.func(args)


if __name__ == "__main__":
    main()
