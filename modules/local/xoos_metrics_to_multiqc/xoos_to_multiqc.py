#!/usr/bin/env python3
"""Convert XOOS metric files into MultiQC-compatible JSON files.

Reads metric files produced by XOOS modules (demux, alignment_metrics,
read_collapser) and writes self-contained JSON files that MultiQC can
parse via its custom content feature.

Output files use the ``_mqc.json`` suffix so MultiQC discovers them
automatically without additional configuration.

Usage examples:

    # Multiple samples with all alignment metric types
    python3 xoos_to_multiqc.py \\
        --demux-run-stats demux/run_stats.tsv \\
        --demux-sample-stats demux/sample_stats.tsv \\
        --sample-name SampleA \\
        --accuracy-metrics sampleA/accuracy/base_level_accuracy_summary.tsv \\
        --hp-errors sampleA/accuracy/hp_errors.tsv \\
        --coverage-metrics sampleA/coverage/coverage_stats.tsv \\
        --read-metrics sampleA/read_metrics/read_metrics_summary.tsv \\
        --read-length-summary sampleA/read_metrics/read_length_summary.tsv \\
        --read-collapser-summary sampleA/read_collapser/summary_stats.tsv \\
        --read-collapser-cluster-sizes sampleA/read_collapser/cluster_size_distributions.tsv \\
        --sample-name SampleB \\
        --accuracy-metrics sampleB/accuracy/base_level_accuracy_summary.tsv \\
        --hp-errors sampleB/accuracy/hp_errors.tsv \\
        --coverage-metrics sampleB/coverage/coverage_stats.tsv \\
        --read-metrics sampleB/read_metrics/read_metrics_summary.tsv \\
        --read-length-summary sampleB/read_metrics/read_length_summary.tsv \\
        --read-collapser-summary sampleB/read_collapser/summary_stats.tsv \\
        --read-collapser-cluster-sizes sampleB/read_collapser/cluster_size_distributions.tsv \\
        --output-dir multiqc_input/

    # Then run MultiQC
    multiqc multiqc_input/
"""

import argparse
import copy
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

from parsers import (
    DEMUX_DISPOSITION_KEYS,
    DEMUX_DUPLEX_BREAKDOWN_KEYS,
    DEMUX_SIMPLEX_BREAKDOWN_KEYS,
    parse_accuracy_metrics,
    parse_cnc_metrics,
    parse_coverage_metrics,
    parse_demux_run_stats,
    parse_demux_sample_stats,
    parse_hp_accuracy,
    parse_tfe_results,
    parse_read_collapser_cluster_size_distributions,
    parse_read_collapser_summary,
    parse_tfe_results,
    parse_read_length_summary,
    parse_read_metrics,
    parse_str_metrics,
    parse_vcfeval_summary,
    parse_vcfeval_roc,
    parse_vcfeval_summary_from_roc
)

__version__ = "latest"

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# MultiQC section definitions
# ---------------------------------------------------------------------------

SECTION_DEFS: dict[str, dict[str, Any]] = {
    "accuracy_metrics": {
        "id": "xoos_accuracy_metrics",
        "section_name": "XOOS Accuracy Metrics",
        "description": (
            "Phred-scale accuracy metrics from XOOS alignment metrics. Each row "
            "is one sample / pipeline stage / genomic region combination "
            "(e.g. 'HG001 | Dedup - Autosomes')."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_accuracy_metrics_table",
            "title": "XOOS Accuracy Metrics",
        },
    },
    "coverage_metrics": {
        "id": "xoos_coverage_metrics",
        "section_name": "XOOS Coverage Metrics",
        "description": (
            "Coverage depth metrics from XOOS alignment metrics. Use the dataset "
            "selector to switch between pipeline stage (dedup / pre-dedup / "
            "subsample), genomic region (autosomes / high-confidence / targets) "
            "and coverage type."
        ),
        "plot_type": "boxplot",
        "pconfig": {
            "id": "xoos_coverage_metrics_table",
            "title": "XOOS Coverage Metrics",
        },
    },
    "read_metrics_read_counts": {
        "id": "xoos_read_metrics_read_counts",
        "section_name": "XOOS Read Metrics - Read Counts",
        "description": (
            "Read-level statistics from XOOS alignment metrics. Each row is one "
            "sample / pipeline stage / genomic region / filter variant "
            "combination (e.g. 'HG001 | Dedup - Autosomes - Post-filter'). "
            "The 'Any' variant counts all reads; 'Post-filter' counts only "
            "reads passing filters."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_read_metrics_read_counts_table",
            "title": "XOOS Read Metrics - Read Counts",
        },
    },
    "read_metrics_base_counts": {
        "id": "xoos_read_metrics_base_counts",
        "section_name": "XOOS Read Metrics - Base Counts",
        "description": (
            "Base-level statistics from XOOS alignment metrics. Each row is one "
            "sample / pipeline stage / genomic region / filter variant "
            "combination (e.g. 'HG001 | Dedup - Autosomes - Post-filter')."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_read_metrics_base_counts_table",
            "title": "XOOS Read Metrics - Base Counts",
        },
    },
    "demux_run_metrics": {
        "id": "xoos_demux_run_metrics",
        "section_name": "XOOS Demux - Run-level Metrics",
        "description": (
            "Run-level read breakdowns from XOOS demux. "
            "'Read Disposition' partitions all reads into passing, failed, and "
            "unassigned (sums to 100% of total reads). 'Unassigned Breakdown' "
            "details why reads were not assigned (duplex: no-hairpin / too-long / "
            "too-short; simplex: too-short-after-trimming / SID-discordant)."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "xoos_demux_run_metrics",
            "title": "XOOS Demux - Run-level Metrics",
            "stacking": "normal",
            "tt_decimals": 0,
        },
    },
    "demux_sample_counts": {
        "id": "xoos_demux_sample_counts",
        "section_name": "XOOS Demux - Per-sample Counts",
        "description": "Per-sample read and base counts from XOOS demux.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_demux_sample_counts_table",
            "title": "XOOS Demux - Per-sample Counts",
        },
        "headers": {
            "Assigned Reads": {"format": "{:,.0f}"},
            "Passing Reads": {"format": "{:,.0f}"},
            "Full Reads": {"format": "{:,.0f}"},
            "Partial Reads": {"format": "{:,.0f}"},
            "Consensus Bases": {"format": "{:,.0f}"},
        },
    },
    "demux_sample_metrics": {
        "id": "xoos_demux_sample_metrics",
        "section_name": "XOOS Demux - Per-sample Metrics",
        "description": "Per-sample read and base breakdowns from XOOS demux. Includes passing read breakdown (full vs partial), failed assigned read breakdown (failure reasons), and base composition (concordant vs discordant vs non-duplex).",
        "plot_type": "bargraph",
        "pconfig": {
            "id": "xoos_demux_sample_metrics",
            "title": "XOOS Demux - Per-sample Metrics",
            "stacking": "normal",
            "tt_decimals": 0
        },
    },
    "vcfeval_summary": {
        "id": "xoos_vcfeval_summary",
        "section_name": "Small Variant VCF Evaluation Summary",
        "description": "VCF evaluation metrics from vcfeval.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_vcfeval_summary_table",
            "title": "Small Variant VCF Evaluation Summary",
        },
        "headers": {
            "Best Threshold": {
                "format": "{:,.2f}"
            },
            "True Positives (Baseline)": {
                "format": "{:,.0f}"
            },
            "True Positives (Call)": {
                "format": "{:,.0f}"
            },
            "False Positives": {
                "format": "{:,.0f}"
            },
            "False Negatives": {
                "format": "{:,.0f}"
            },
            "Precision": {
                "format": "{:,.4f}"
            },
            "Sensitivity": {
                "format": "{:,.4f}"
            },
            "F1 Score": {
                "format": "{:,.4f}"
            }
        }
    },
    "vcfeval_summary_snp": {
        "id": "xoos_vcfeval_summary_snp",
        "section_name": "Small Variant VCF Evaluation Summary - SNPs",
        "description": "VCF evaluation metrics for SNPs from vcfeval.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_vcfeval_summary_snp_table",
            "title": "Small Variant VCF Evaluation Summary - SNPs",
        },
        "headers": {
            "Best Threshold": {
                "format": "{:,.2f}"
            },
            "True Positives (Baseline)": {
                "format": "{:,.0f}"
            },
            "True Positives (Call)": {
                "format": "{:,.0f}"
            },
            "False Positives": {
                "format": "{:,.0f}"
            },
            "False Negatives": {
                "format": "{:,.0f}"
            },
            "Precision": {
                "format": "{:,.4f}"
            },
            "Sensitivity": {
                "format": "{:,.4f}"
            },
            "F1 Score": {
                "format": "{:,.4f}"
            }
        }
    },
    "vcfeval_summary_indel": {
        "id": "xoos_vcfeval_summary_indel",
        "section_name": "Small Variant VCF Evaluation Summary - Indels",
        "description": "VCF evaluation metrics for indels from vcfeval.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_vcfeval_summary_indel_table",
            "title": "Small Variant VCF Evaluation Summary - Indels",
        },
        "headers": {
            "Best Threshold": {
                "format": "{:,.2f}"
            },
            "True Positives (Baseline)": {
                "format": "{:,.0f}"
            },
            "True Positives (Call)": {
                "format": "{:,.0f}"
            },
            "False Positives": {
                "format": "{:,.0f}"
            },
            "False Negatives": {
                "format": "{:,.0f}"
            },
            "Precision": {
                "format": "{:,.4f}"
            },
            "Sensitivity": {
                "format": "{:,.4f}"
            },
            "F1 Score": {
                "format": "{:,.4f}"
            }
        }
    },
    "vcfeval_roc": {
        "id": "xoos_vcfeval_roc",
        "section_name": "Small Variant VCF Evaluation - ROC Curve",
        "description": "ROC curve data for SNPs from vcfeval.",
        "plot_type": "linegraph",
        "pconfig": {
            "id": "xoos_vcfeval_roc",
            "title": "Small Variant VCF Evaluation - ROC Curve",
            "xlab": "False Positives",
            "ylab": "True Positives"
        },
    },
    "tfe_contamination_metrics": {
        "id": "xoos_tfe_contamination_metrics",
        "section_name": "XOOS TFE / Contamination Metrics",
        "description": "Summary metrics from the tumor fraction estimator.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_tfe_contamination_metrics_table",
            "title": "XOOS TFE / Contamination Metrics",
        },
        "headers": {
            "contamination fraction": {"format": "{:,.6f}", "min": 0},
            "mean vaf": {"format": "{:,.6f}", "min": 0},
            "error rate": {"format": "{:,.6f}", "min": 0},
            "total alt": {"format": "{:,.0f}"},
            "total depth": {"format": "{:,.0f}"},
            "total variants": {"format": "{:,.0f}"},
            "variants detected": {"format": "{:,.0f}"},
            "passing variants": {"format": "{:,.0f}"},
            "passing variants detected": {"format": "{:,.0f}"},
        },
    },
    "tfe_metrics": {
        "id": "xoos_tfe_metrics",
        "section_name": "XOOS Tumor Fraction Estimation Metrics",
        "description": "Summary metrics from the tumor fraction estimator.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_tfe_metrics_table",
            "title": "XOOS Tumor Fraction Estimation Metrics",
        },
        "headers": {
            "tumor fraction": {"format": "{:,.6f}", "min": 0},
            "mean vaf": {"format": "{:,.6f}", "min": 0},
            "error rate": {"format": "{:,.6f}", "min": 0},
            "total alt": {"format": "{:,.0f}"},
            "total depth": {"format": "{:,.0f}"},
            "total variants": {"format": "{:,.0f}"},
            "variants detected": {"format": "{:,.0f}"},
            "passing variants": {"format": "{:,.0f}"},
            "passing variants detected": {"format": "{:,.0f}"},
        },
    },
    "str_metrics": {
        "id": "xoos_str_metrics",
        "section_name": "XOOS STR Caller Metrics",
        "description": "Summary metrics from the STR caller.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_str_metrics_table",
            "title": "XOOS STR Caller Metrics",
        },
        "headers": {
            "NumSites": {"format": "{:,.0f}", "description": "Number of STR sites analyzed"},
            "LowDepth": {"format": "{:,.0f}", "description": "Sites with low depth"},
            "LowConfidence": {"format": "{:,.0f}", "description": "Sites with low confidence"},
        },
    },
    "cnc_metrics": {
        "id": "xoos_cnc_metrics",
        "section_name": "XOOS Copy Number Caller Metrics",
        "description": "Summary metrics from the copy number caller.",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_cnc_metrics_table",
            "title": "XOOS Copy Number Caller Metrics",
        },
        "headers": {
            "Number of targets": {"format": "{:,.0f}", "description": "Total number of target intervals"},
            "Sex": {"description": "Inferred sex"},
            "Number of segments": {"format": "{:,.0f}", "description": "Total segments called"},
            "Number of gain segments": {"format": "{:,.0f}", "description": "Segments with copy number gain"},
            "Number of loss segments": {"format": "{:,.0f}", "description": "Segments with copy number loss"},
            "Number of passing gain segments": {"format": "{:,.0f}", "description": "Passing gain segments"},
            "Number of passing loss segments": {"format": "{:,.0f}", "description": "Passing loss segments"},
            "Median coverage": {"format": "{:,.2f}", "description": "Median coverage across targets"},
            "Percentage of low coverage windows": {"format": "{:,.2f}", "description": "Percentage of windows with low coverage"},
        },
    },
    "read_collapser_metrics": {
        "id": "xoos_read_collapser_metrics",
        "section_name": "XOOS Read Collapser Metrics",
        "description": "Summary metrics from the XOOS read collapser (markdup or consensus mode).",
        "plot_type": "table",
        "pconfig": {
            "id": "xoos_read_collapser_metrics_table",
            "title": "XOOS Read Collapser Metrics",
        },
        "headers": {
            "Input Reads": {"format": "{:,.0f}"},
            "Discarded Reads": {"format": "{:,.0f}"},
            "Discarded Reads (%)": {"format": "{:,.2f}", "suffix": "%"},
            "Clustering Input Reads": {"format": "{:,.0f}"},
            "Clustering Input Reads (%)": {"format": "{:,.2f}", "suffix": "%"},
            "Duplicate Reads": {"format": "{:,.0f}"},
            "Duplicate Reads (%)": {"format": "{:,.2f}", "suffix": "%"},
            "Total Clusters": {"format": "{:,.0f}"},
            "Singleton Clusters": {"format": "{:,.0f}"},
            "Singleton Clusters (%)": {"format": "{:,.2f}", "suffix": "%"},
            "Consensus Clusters": {"format": "{:,.0f}"},
            "Consensus Median Cluster Size": {"format": "{:,.0f}"},
            "Total Consensus Reads": {"format": "{:,.0f}"},
        },
    },
    "read_collapser_cluster_sizes": {
        "id": "xoos_read_collapser_cluster_sizes",
        "section_name": "XOOS Read Collapser - Cluster Size Distribution",
        "description": "Cluster size distributions from the XOOS read collapser.",
        "plot_type": "linegraph",
        "pconfig": {
            "id": "xoos_read_collapser_cluster_sizes",
            "title": "XOOS Read Collapser - Cluster Size Distribution",
            "xlab": "Cluster Size",
            "ylab": "Count",
            "xmin": 1,
            "tt_decimals": 0,
        },
    },
}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

class SampleAction(argparse.Action):
    """Groups file arguments by the most recently specified --sample-name."""

    def __call__(self, parser, namespace, values, option_string=None):
        if not hasattr(namespace, "_sample_inputs") or namespace._sample_inputs is None:
            namespace._sample_inputs = defaultdict(lambda: defaultdict(list))

        current_sample = getattr(namespace, "_current_sample", None)
        if current_sample is None:
            parser.error(f"{option_string} requires a preceding --sample-name")

        key = option_string.lstrip("-").replace("-", "_")
        if isinstance(values, list):
            namespace._sample_inputs[current_sample][key].extend(values)
        else:
            namespace._sample_inputs[current_sample][key].append(values)


class SampleNameAction(argparse.Action):
    """Sets the current sample name context for subsequent file arguments."""

    def __call__(self, parser, namespace, values, option_string=None):
        namespace._current_sample = values
        if not hasattr(namespace, "_sample_inputs") or namespace._sample_inputs is None:
            namespace._sample_inputs = defaultdict(lambda: defaultdict(list))
        namespace._sample_inputs[values]  # noqa: B018


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert XOOS metric files to MultiQC-compatible JSON.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )

    parser.add_argument(
        "--sample-name",
        action=SampleNameAction,
        help="Sample name for the following metric files. Can be repeated.",
    )

    parser.add_argument(
        "--region",
        type=str,
        default=None,
        help=(
            "Genomic region these alignment metrics were computed over "
            "(e.g. 'autosomes', 'high_confidence', 'targets'). Used to label "
            "the read/coverage tab so the report can show all regions together."
        ),
    )
    parser.add_argument(
        "--stage",
        type=str,
        default=None,
        help=(
            "Pipeline stage these alignment metrics correspond to "
            "(e.g. 'pre_dedup', 'dedup', 'dedup_subsample'). Combined with "
            "--region to label the read/coverage tab."
        ),
    )

    # Per-sample metric file arguments
    parser.add_argument(
        "--accuracy-metrics",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="base_level_accuracy_summary.tsv from alignment_metrics accuracy.",
    )
    parser.add_argument(
        "--hp-errors",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="hp_errors.tsv from alignment_metrics accuracy.",
    )
    parser.add_argument(
        "--coverage-metrics",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="coverage_stats.tsv from alignment_metrics coverage.",
    )
    parser.add_argument(
        "--read-metrics",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="read_metrics_summary.tsv from alignment_metrics read.",
    )
    parser.add_argument(
        "--read-length-summary",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="read_length_summary.tsv from alignment_metrics read.",
    )
    parser.add_argument(
        "--vcfeval-snp-roc",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="snp_roc.tsv.gz from vcfeval.",
    )
    parser.add_argument(
        "--vcfeval-indel-roc",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="non_snp_roc.tsv.gz from vcfeval.",
    )
    parser.add_argument(
        "--vcfeval-summary",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="summary.tsv from vcfeval.",
    )
    parser.add_argument(
        "--tfe-contamination-results",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="TFE results CSV from tumor fraction estimator (contamination mode).",
    )
    parser.add_argument(
        "--tfe-results",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="Results TSV from tumor_fraction_estimator (e.g. tumor_fraction.results.tsv or contamination.results.tsv).",
    )
    parser.add_argument(
        "--str-metrics",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="metrics.tsv from str_caller.",
    )
    parser.add_argument(
        "--cnc-metrics",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="*_metrics.tsv from copy_number_caller.",
    )
    parser.add_argument(
        "--read-collapser-summary",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="summary_stats.tsv from read_collapser (markdup or consensus).",
    )
    parser.add_argument(
        "--read-collapser-cluster-sizes",
        type=Path,
        action=SampleAction,
        metavar="FILE",
        help="cluster_size_distributions.tsv from read_collapser.",
    )

    parser.add_argument(
        "--demux-run-stats",
        type=Path,
        metavar="FILE",
        help="Demux run_stats.tsv (run-level).",
    )
    parser.add_argument(
        "--demux-sample-stats",
        type=Path,
        metavar="FILE",
        help="Demux sample_stats.tsv (contains all samples).",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory for MultiQC JSON files (default: cwd).",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
    )

    return parser


# ---------------------------------------------------------------------------
# Region / stage tab labelling
# ---------------------------------------------------------------------------

# Human-readable names for known stages and regions, and the order in which
# stage tabs should appear. "dedup" (the post-filter / deduplicated output) is
# listed first so it becomes the default selected tab in MultiQC.
_STAGE_ORDER = ["dedup", "dedup_subsample", "pre_dedup"]
_STAGE_LABELS = {
    "dedup": "Dedup",
    "dedup_subsample": "Dedup Subsample",
    "pre_dedup": "Pre-dedup",
}
_REGION_LABELS = {
    "autosomes": "Autosomes",
    "high_confidence": "High Confidence",
    "targets": "Targets",
}


def _tab_label(region: str | None, stage: str | None) -> str | None:
    """Build a human-readable tab label from region and stage.

    Returns None when neither is provided (legacy single-table behaviour).
    """
    if not region and not stage:
        return None
    stage_label = _STAGE_LABELS.get(stage, stage.replace("_", " ").title()) if stage else None
    region_label = _REGION_LABELS.get(region, region.replace("_", " ").title()) if region else None
    parts = [p for p in (stage_label, region_label) if p]
    return " - ".join(parts)


def _row_key(sample: str, tab_label: str | None, variant: str | None = None) -> str:
    """Build a table row key that folds region/stage and filter variant in.

    MultiQC custom-content tables do not support the ``data_labels`` dataset
    switcher (only plot types do), so region/stage/variant cannot be exposed as
    selectable tabs for table sections. Instead they are folded into the row key
    (e.g. ``HG001 | Dedup - Autosomes - Post-filter``) so a single table shows
    every region/stage/variant combination as its own row.

    With no tab label and no variant this returns the bare sample name, keeping
    the legacy single-row-per-sample layout.
    """
    suffix_parts = [p for p in (tab_label, variant) if p]
    if not suffix_parts:
        return sample
    return f"{sample} | {' - '.join(suffix_parts)}"


# ---------------------------------------------------------------------------
# Bar plot data builders
# ---------------------------------------------------------------------------

def _build_bargraph(
    samples: dict[str, dict[str, Any]],
    keys: list[str],
) -> dict[str, dict[str, float]]:
    """Extract the given keys from each sample into a bar plot data dict."""
    data: dict[str, dict[str, float]] = {}
    for sample, m in samples.items():
        entry: dict[str, float] = {}
        for key in keys:
            val = m.get(key)
            if isinstance(val, (int, float)):
                entry[key] = val
        if entry:
            data[sample] = entry
    return data

# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def _parse_file(path: Path, parser_fn, target: dict[str, Any]) -> None:
    """Parse a file if it exists and merge results into target dict."""
    if path.exists():
        target.update(parser_fn(path))
    else:
        logger.warning("File not found: %s", path)


def process_inputs(args: argparse.Namespace) -> None:
    """Parse all input files and write MultiQC JSON output."""
    sample_inputs: dict[str, dict[str, list]] = getattr(args, "_sample_inputs", {}) or {}
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Tab label for region/stage-aware alignment metric sections. None falls
    # back to the legacy single-table layout.
    tab_label = _tab_label(getattr(args, "region", None), getattr(args, "stage", None))

    # Sections with only one report
    # {section_key: {sample_name: {metric: value}}}
    single_report_sections: dict[str, dict[str, dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

    # Sections with multiple reports
    # {section_key: {report_name: {sample_name: {metric: value}}}}
    multi_report_sections: dict[str, dict[str, dict[str, dict[str, Any]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    # --- Per-sample alignment metrics ---
    for sample_name, inputs in sample_inputs.items():
        # Accuracy is region/stage-aware, but it is a table section and MultiQC
        # tables cannot render a dataset selector. Fold region/stage into the
        # row key so every region/stage appears as its own row in one table
        # instead of silently dropping all but the first dataset.
        accuracy_row = _row_key(sample_name, tab_label)
        accuracy_target = single_report_sections["accuracy_metrics"][accuracy_row]

        for f in inputs.get("accuracy_metrics", []):
            _parse_file(f, parse_accuracy_metrics, accuracy_target)

        for f in inputs.get("hp_errors", []):
            _parse_file(f, parse_hp_accuracy, accuracy_target)

        for f in inputs.get("coverage_metrics", []):
            coverage_metrics = parse_coverage_metrics(f)
            # When region/stage is provided, prefix each coverage sub-type with
            # the region+stage tab so MultiQC shows one selectable dataset per
            # region/stage/sub-type. Otherwise keep the legacy fixed labels.
            for sub_key, sub_label in (
                ("coverage", "Coverage"),
                ("post_filter_coverage", "Post-filter Coverage"),
                ("concordant_duplex_coverage", "Concordant Duplex Coverage"),
            ):
                if coverage_metrics[sub_key]:
                    label = f"{tab_label} - {sub_label}" if tab_label else sub_label
                    multi_report_sections["coverage_metrics"][label][sample_name] = coverage_metrics[sub_key]

        for f in inputs.get("read_metrics", []):
            # Read metrics expose an Any vs Post-filter filter variant and are
            # region/stage-aware. These are table sections and MultiQC tables
            # cannot render a dataset selector, so region/stage and the filter
            # variant are folded into the row key (e.g.
            # "HG001 | Dedup - Autosomes - Post-filter"). Each combination then
            # appears as its own row instead of being dropped.
            read_metrics = parse_read_metrics(f)
            for variant, data in read_metrics.items():
                row = _row_key(sample_name, tab_label, variant)
                single_report_sections["read_metrics_read_counts"][row].update(data["read_counts"])
                single_report_sections["read_metrics_base_counts"][row].update(data["base_counts"])

        for f in inputs.get("read_length_summary", []):
            # Read length has no filter variant; attach it to the "Any" row so it
            # sits alongside the unfiltered read counts for the same region/stage.
            length_row = _row_key(sample_name, tab_label, "Any")
            _parse_file(f, parse_read_length_summary, single_report_sections["read_metrics_read_counts"][length_row])

        for f in inputs.get("vcfeval_summary", []):
            _parse_file(f, parse_vcfeval_summary, single_report_sections["vcfeval_summary"][sample_name])

        for f in inputs.get("vcfeval_snp_roc", []):
            multi_report_sections["vcfeval_roc"]["SNP"][sample_name] = parse_vcfeval_roc(f)
            _parse_file(f, parse_vcfeval_summary_from_roc, single_report_sections["vcfeval_summary_snp"][sample_name])

        for f in inputs.get("vcfeval_indel_roc", []):
            multi_report_sections["vcfeval_roc"]["Indel"][sample_name] = parse_vcfeval_roc(f)
            _parse_file(f, parse_vcfeval_summary_from_roc, single_report_sections["vcfeval_summary_indel"][sample_name])

        for f in inputs.get("tfe_contamination_results", []):
            _parse_file(f, parse_tfe_results, single_report_sections["tfe_contamination_metrics"][sample_name])

        for f in inputs.get("tfe_results", []):
            _parse_file(f, parse_tfe_results, single_report_sections["tfe_metrics"][sample_name])

        for f in inputs.get("str_metrics", []):
            _parse_file(f, parse_str_metrics, single_report_sections["str_metrics"][sample_name])

        for f in inputs.get("cnc_metrics", []):
            _parse_file(f, parse_cnc_metrics, single_report_sections["cnc_metrics"][sample_name])

        for f in inputs.get("read_collapser_summary", []):
            _parse_file(f, parse_read_collapser_summary, single_report_sections["read_collapser_metrics"][sample_name])

        for f in inputs.get("read_collapser_cluster_sizes", []):
            if not f.exists():
                logger.warning("File not found: %s", f)
                continue
            distributions = parse_read_collapser_cluster_size_distributions(f)
            for dist_name, xy_data in distributions.items():
                label = dist_name.replace("_", " ").title()
                multi_report_sections["read_collapser_cluster_sizes"][label][sample_name] = xy_data

    # --- Demux run-level metrics ---
    run_stats: dict[str, Any] = {}
    if args.demux_run_stats and args.demux_run_stats.exists():
        run_stats = parse_demux_run_stats(args.demux_run_stats)
        if run_stats:
            mode = run_stats.get("_mode", "duplex")

            # Read disposition: passing + failed + unassigned = 100% of total.
            # Canonical keys are identical for simplex and duplex.
            disposition_bar = _build_bargraph(
                {"Entire Run": run_stats},
                DEMUX_DISPOSITION_KEYS,
            )
            if disposition_bar:
                multi_report_sections["demux_run_metrics"]["Read Disposition"] = disposition_bar

            # Second breakdown of reads that did not pass assignment. The
            # categories differ by run type but both are shown as
            # "Unassigned Breakdown" for a consistent report layout.
            breakdown_keys = (
                DEMUX_DUPLEX_BREAKDOWN_KEYS
                if mode == "duplex"
                else DEMUX_SIMPLEX_BREAKDOWN_KEYS
            )
            unassigned_bar = _build_bargraph(
                {"Entire Run": run_stats},
                breakdown_keys,
            )
            if unassigned_bar:
                multi_report_sections["demux_run_metrics"]["Unassigned Breakdown"] = unassigned_bar

    # --- Demux per-sample metrics ---
    per_sample = {}
    if args.demux_sample_stats and args.demux_sample_stats.exists():
        per_sample = parse_demux_sample_stats(args.demux_sample_stats)

    # --- Build bar plot sections from parsed data ---
    if per_sample:
        # Per-sample counts table: the core assignment counts (and bases when
        # present for duplex). Drop the index_sequence string so the table is
        # purely numeric.
        counts_keys = [
            "Assigned Reads",
            "Passing Reads",
            "Full Reads",
            "Partial Reads",
            "Consensus Bases",
        ]
        counts_table = {
            sample: {k: m[k] for k in counts_keys if isinstance(m.get(k), (int, float))}
            for sample, m in per_sample.items()
        }
        counts_table = {s: v for s, v in counts_table.items() if v}
        if counts_table:
            single_report_sections["demux_sample_counts"] = counts_table

        # Passing read composition: full + partial
        passing_reads_bar = _build_bargraph(
            per_sample,
            ["Full Reads", "Partial Reads"],
        )
        if passing_reads_bar:
            multi_report_sections["demux_sample_metrics"]["Passing Reads Breakdown"] = passing_reads_bar

        # Failed / assignment-quality breakdown. Duplex emits failure-reason
        # categories; simplex emits assignment-quality categories. Include the
        # union so whichever set is present renders.
        failed_reads_bar = _build_bargraph(
            per_sample,
            [
                # duplex failed-assigned reasons. These mirror the demux
                # FailedAssignedCounts struct (too_many_errors,
                # trimmed_read_too_short, consensus_too_long) and partition
                # failed_assigned_reads. "No End Adapter" is intentionally
                # excluded: in the demux source it is a PassingCounts attribute
                # ("reads that didn't have an endadapter at all"), not a failure
                # reason, so it must not be shown as a failed category.
                "Too Many Errors",
                "Too Short After Trimming",
                "Too Long",
                # simplex assignment-quality categories
                # (Index Hopping is the same metric as SID Discordant; only the
                # latter is shown to avoid double-counting.)
                "SID Discordant",
            ],
        )
        if failed_reads_bar:
            multi_report_sections["demux_sample_metrics"]["Failed Reads Breakdown"] = failed_reads_bar

        # Base composition: concordant + discordant + non-duplex
        bases_bar = _build_bargraph(
            per_sample,
            ["Concordant Bases", "Discordant Bases", "Non-duplex Bases"],
        )
        if bases_bar:
            multi_report_sections["demux_sample_metrics"]["Base Composition"] = bases_bar

    # --- Write MultiQC JSON files ---
    files_written = 0
    for section_key, data in single_report_sections.items():
        if not data:
            continue

        section_def = SECTION_DEFS.get(section_key)
        if section_def is None:
            logger.warning("No section definition for key: %s", section_key)
            continue

        output = {**section_def, "data": dict(data)}
        out_file = output_dir / f"{section_def['id']}_mqc.json"
        with out_file.open("w") as f:
            json.dump(output, f, indent=2)
            logger.info("Wrote %s", out_file)
            files_written += 1

    for section_key, reports in multi_report_sections.items():
        if not reports:
            continue

        section_def = SECTION_DEFS.get(section_key)
        if section_def is None:
            logger.warning("No section definition for key: %s", section_key)
            continue

        section_def = copy.deepcopy(section_def)
        if "pconfig" in section_def:
            section_def["pconfig"]["data_labels"] = list(reports.keys())
        else:
            section_def["pconfig"] = {"data_labels": list(reports.keys())}

        output = {**section_def, "data": list(reports.values())}

        out_file = output_dir / f"{section_def['id']}_mqc.json"
        with out_file.open("w") as f:
            json.dump(output, f, indent=2)
            logger.info("Wrote %s", out_file)
            files_written += 1

    if files_written == 0:
        logger.warning("No metric files were processed. Check your inputs.")
    else:
        logger.info("Wrote %d MultiQC JSON file(s) to %s", files_written, output_dir)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    if not hasattr(args, "_sample_inputs"):
        args._sample_inputs = {}

    has_input = bool(args._sample_inputs) or args.demux_run_stats or args.demux_sample_stats
    if not has_input:
        parser.error("No input files specified. Use --sample-name with metric file "
                     "arguments, or --demux-run-stats / --demux-sample-stats.")

    process_inputs(args)


if __name__ == "__main__":
    main()
