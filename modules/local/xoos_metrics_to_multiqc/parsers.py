"""Parsers for XOOS metric files.

Each parser reads a specific XOOS metric file format and returns a flat
dict of {metric_name: value} for a single sample.
"""

import csv
import io
import logging
from pathlib import Path
from typing import Any
import re
import gzip

logger = logging.getLogger(__name__)


def _dict_reader(fh: io.TextIOBase, delimiter: str = "\t") -> csv.DictReader:
    return csv.DictReader(
        (line for line in fh if line.strip() and not line.startswith("#")),
        delimiter=delimiter,
    )


def _num(value: str) -> int | float | str:
    """Convert a string to int or float if possible, otherwise return as-is."""
    if value in ("NA", ""):
        return value
    try:
        if "." in value or "e" in value.lower():
            return float(value)
        return int(value)
    except (ValueError, TypeError):
        return value


# ---------------------------------------------------------------------------
# Accuracy metrics (base_level_accuracy_summary.tsv)
# ---------------------------------------------------------------------------

def parse_accuracy_metrics(input_file: Path) -> dict[str, Any]:
    """Extract Phred-scale accuracy values from base_level_accuracy_summary.tsv.

    Returns metrics keyed as:
      Overall Phred-scale Accuracy
      Indel Phred-scale Accuracy
      Substitution Phred-scale Accuracy
    For substitutions (which has two rows with different denominators), we take
    the row whose denominator contains 'indel_bases'.
    """
    logger.info("Parsing accuracy metrics: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        for row in _dict_reader(f):
            error_type = row.get("type", "")
            phred = _num(row.get("phred", "NA"))

            if error_type == "all_errors_with_indel_bases":
                metrics["Overall Phred-scale Accuracy"] = phred
            elif error_type == "indel_bases":
                metrics["Indel Phred-scale Accuracy"] = phred
            elif error_type == "substitutions":
                denom = row.get("denominator", "")
                # Take the row with indel_bases denominator for consistency
                if "indel_bases" in denom:
                    metrics["Substitution Phred-scale Accuracy"] = phred
    return metrics


# ---------------------------------------------------------------------------
# HP errors (hp_errors.tsv)
# ---------------------------------------------------------------------------

def parse_hp_accuracy(input_file: Path) -> dict[str, Any]:
    """Extract HP Phred-scale accuracy for hp_base='any' at lengths 8 and 10."""
    logger.info("Parsing HP accuracy: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        for row in _dict_reader(f):
            hp_base = row.get("hp_base", "")
            hp_length = row.get("hp_length", "")
            phred = _num(row.get("phred", "NA"))

            if hp_base == "any" and hp_length == "8":
                metrics["8-base HP Phred-scale Accuracy"] = phred
            elif hp_base == "any" and hp_length == "10":
                metrics["10-base HP Phred-scale Accuracy"] = phred
    return metrics


# ---------------------------------------------------------------------------
# Coverage metrics (coverage_stats.tsv)
# ---------------------------------------------------------------------------

_DISTRIBUTION_SUMMARY_METRICS_MAP = {
    "min": "min",
    "percentile_25": "q1",
    "median": "median",
    "percentile_75": "q3",
    "max": "max",
    "mean": "mean",
}


def _set_numeric_coverage(dest: dict[str, Any], key: str, raw: str) -> None:
    """Store a coverage value only if it is numeric.

    coverage_stats.tsv uses ``NA`` for metrics that cannot be computed
    (e.g. empty histograms after excluding zero-coverage bins).  MultiQC
    boxplots require all values to be numeric, so non-numeric entries are
    silently dropped.
    """
    value = _num(raw)
    if isinstance(value, (int, float)):
        dest[key] = value


def parse_coverage_metrics(input_file: Path) -> dict[str, dict[str, Any]]:
    """Extract boxplot summary statistics from coverage_stats.tsv.

    Reads the five-number summary (min, q1, median, q3, max) plus mean
    for each coverage type.  Non-numeric values (``NA``) are dropped so
    the output is safe to pass directly to MultiQC's boxplot renderer.

    Returns a dict with keys ``coverage``, ``post_filter_coverage``, and
    ``concordant_duplex_coverage``.  All three keys are always present;
    each value is a dict of numeric summary stats, or an empty dict when
    the column is missing from the input or all its values were NA.
    """
    logger.info("Parsing coverage metrics: %s", input_file)
    metrics: dict[str, dict[str, Any]] = {"coverage": {}, "concordant_duplex_coverage": {}, "post_filter_coverage": {}}
    with input_file.open() as f:
        for row in _dict_reader(f):
            name = row.get("metric_name", "")
            has_duplex = "concordant_duplex_coverage" in row

            if name in _DISTRIBUTION_SUMMARY_METRICS_MAP:
                mapped = _DISTRIBUTION_SUMMARY_METRICS_MAP[name]
                _set_numeric_coverage(metrics["coverage"], mapped, row["any_coverage"])
                _set_numeric_coverage(metrics["post_filter_coverage"], mapped, row["post_filter_coverage"])
                if has_duplex:
                    _set_numeric_coverage(metrics["concordant_duplex_coverage"], mapped, row["concordant_duplex_coverage"])
    return metrics

# ---------------------------------------------------------------------------
# Read metrics (read_metrics_summary.tsv, read_length_summary.tsv)
# ---------------------------------------------------------------------------

# read_metrics_summary.tsv encodes the filter variant in the metric name rather
# than in a column: "any" reads/bases use the bare names while the post-filter
# subset uses post_filter_*/*_passing_filter names. Map each variant separately
# so the report can offer an Any vs Post-filter selector. Coverage uses literal
# columns instead and is handled in parse_coverage_metrics; accuracy has no such
# split (it is always computed on post-filter bases).
READ_METRIC_FILTER_VARIANTS = ("Any", "Post-filter")

_READ_METRIC_NAME_MAP = {
    "Any": {
        "total_reads": "Total Reads",
        "duplicate_reads": "Duplicate Reads",
        "unmapped_reads": "Unmapped Reads",
    },
    "Post-filter": {
        "post_filter_read_count": "Total Reads",
    },
}

_BASE_METRICS_NAME_MAP = {
    "Any": {
        "total_bases": "Total Bases",
        "aligned_bases": "Aligned Bases",
        "soft_clipped_bases": "Soft-clipped Bases",
    },
    "Post-filter": {
        "total_bases_passing_filter": "Total Bases",
        "aligned_bases_passing_filter": "Aligned Bases",
        "soft_clipped_bases_passing_filter": "Soft-clipped Bases",
    },
}

def parse_read_metrics(input_file: Path) -> dict[str, dict[str, dict[str, Any]]]:
    """Extract selected read metrics from read_metrics_summary.tsv.

    Returns a nested mapping keyed by filter variant ("Any" / "Post-filter"),
    each containing "read_counts" and "base_counts" sub-dicts. Only variants
    that actually have data are returned, so reports never show an empty
    selector entry.
    """
    logger.info("Parsing read metrics: %s", input_file)
    metrics: dict[str, dict[str, dict[str, Any]]] = {
        variant: {"read_counts": {}, "base_counts": {}}
        for variant in READ_METRIC_FILTER_VARIANTS
    }
    with input_file.open() as f:
        for row in _dict_reader(f):
            name = row.get("metric_name", "")
            count = _num(row.get("count", "NA"))
            pct = _num(row.get("percentage", "NA"))
            for variant in READ_METRIC_FILTER_VARIANTS:
                if name in _READ_METRIC_NAME_MAP[variant]:
                    display = _READ_METRIC_NAME_MAP[variant][name]
                    metrics[variant]["read_counts"][display] = count
                    if pct != "NA":
                        metrics[variant]["read_counts"][f"{display} (%)"] = pct
                elif name in _BASE_METRICS_NAME_MAP[variant]:
                    display = _BASE_METRICS_NAME_MAP[variant][name]
                    metrics[variant]["base_counts"][display] = count
                    if pct != "NA":
                        metrics[variant]["base_counts"][f"{display} (%)"] = pct
    # Drop variants with no data so callers can skip empty selector entries.
    return {
        variant: data
        for variant, data in metrics.items()
        if data["read_counts"] or data["base_counts"]
    }


def parse_read_length_summary(input_file: Path) -> dict[str, Any]:
    """Extract mean read length from read_length_summary.tsv."""
    logger.info("Parsing read length summary: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        for row in _dict_reader(f):
            name = row.get("metric_name", "")
            if name == "mean":
                # Use the first value column
                for col, val in row.items():
                    if col != "metric_name":
                        metrics["Mean Read Length"] = _num(val)
                        break
    return metrics


# ---------------------------------------------------------------------------
# Demux metrics (run_stats.tsv, sample_stats.tsv)
# ---------------------------------------------------------------------------

# Demux writes different field names for simplex and duplex runs. The parsers
# below normalize both into a single canonical set of display names so the
# MultiQC report renders identically regardless of run type.
#
# Run-level read disposition (each partition sums to 100% of Total Reads):
#   duplex : Total = passing_reads + failed_assigned_reads + unassigned_reads
#   simplex: Total = assigned_reads + failed_reads        + unassigned_reads
# Second breakdown (the categories explaining reads that did not pass assignment):
#   duplex : unassigned_reads = no_hairpin + too_long + too_short
#   simplex: (subset of unassigned) = too_short_trimmed + sid_discordant_discarded

# Canonical disposition keys shared by both modes (used to build the 100% bar).
DEMUX_DISPOSITION_KEYS = ["Passing Reads", "Failed Reads", "Unassigned Reads"]

# Per-mode breakdown keys (the second run-level bar). Both are titled
# "Unassigned Breakdown" in the report; the categories differ by run type.
DEMUX_DUPLEX_BREAKDOWN_KEYS = ["No Hairpin", "Too Long", "Too Short"]
DEMUX_SIMPLEX_BREAKDOWN_KEYS = ["Too Short After Trimming", "SID Discordant"]


def _set_if_present(target: dict[str, Any], key: str, value: Any) -> None:
    """Assign value to target[key] only when value is not None."""
    if value is not None:
        target[key] = value


def parse_demux_run_stats(input_file: Path) -> dict[str, Any]:
    """Extract run-level demux metrics from run_stats.tsv.

    Normalizes simplex and duplex field names into canonical display names.
    The returned dict carries a ``_mode`` key (``"duplex"`` or ``"simplex"``)
    so callers can select the appropriate second breakdown. The disposition
    metrics (``DEMUX_DISPOSITION_KEYS``) sum to ``Total Reads`` in both modes.
    """
    logger.info("Parsing demux run stats: %s", input_file)
    raw: dict[str, Any] = {}
    with input_file.open() as f:
        for row in _dict_reader(f):
            name = row.get("metrics_name", "")
            if name:
                raw[name] = _num(row.get("count", 0))

    if not raw:
        return {}

    # Duplex emits hairpin/duplex-specific fields; simplex does not.
    is_duplex = any(
        k in raw for k in ("no_hairpin_reads", "full_duplex_reads", "passing_reads")
    )

    metrics: dict[str, Any] = {}
    _set_if_present(metrics, "Total Reads", raw.get("total_reads"))

    if is_duplex:
        metrics["_mode"] = "duplex"
        _set_if_present(metrics, "Passing Reads", raw.get("passing_reads"))
        _set_if_present(metrics, "Failed Reads", raw.get("failed_assigned_reads"))
        _set_if_present(metrics, "Unassigned Reads", raw.get("unassigned_reads"))
        _set_if_present(metrics, "Full Reads", raw.get("full_duplex_reads"))
        _set_if_present(metrics, "Partial Reads", raw.get("partial_duplex_reads"))
        # Unassigned breakdown
        _set_if_present(metrics, "No Hairpin", raw.get("no_hairpin_reads"))
        _set_if_present(metrics, "Too Long", raw.get("too_long_reads"))
        _set_if_present(metrics, "Too Short", raw.get("too_short_reads"))
    else:
        metrics["_mode"] = "simplex"
        # In simplex, assigned reads are the passing reads.
        _set_if_present(metrics, "Passing Reads", raw.get("assigned_reads"))
        _set_if_present(metrics, "Failed Reads", raw.get("failed_reads"))
        _set_if_present(metrics, "Unassigned Reads", raw.get("unassigned_reads"))
        _set_if_present(metrics, "Full Reads", raw.get("full_reads"))
        _set_if_present(metrics, "Partial Reads", raw.get("partial_reads"))
        # Breakdown of reads that passed pre-filter but were not assigned.
        _set_if_present(
            metrics, "Too Short After Trimming", raw.get("too_short_trimmed_reads")
        )
        _set_if_present(
            metrics, "SID Discordant", raw.get("sid_discordant_discarded_reads")
        )

    return metrics


# Per-sample metric names, unioned across simplex and duplex sample_stats.tsv.
# Whichever variant is present in the file is mapped to the canonical name.
_DEMUX_SAMPLE_METRIC_MAP = {
    # counts
    "index_sequence": "Index Sequence",
    "assigned_reads": "Assigned Reads",
    "passing_reads": "Passing Reads",
    "full_duplex_reads": "Full Reads",
    "full_reads": "Full Reads",
    "partial_duplex_reads": "Partial Reads",
    "partial_reads": "Partial Reads",
    # duplex failed-assigned breakdown (partitions failed_assigned_reads;
    # mirrors the demux FailedAssignedCounts struct)
    "too_many_errors_reads": "Too Many Errors",
    "too_short_trimmed_reads": "Too Short After Trimming",
    "too_long_consensus_reads": "Too Long",
    # passing-read attribute, NOT a failure reason: in the demux source
    # no_endadapter lives in PassingCounts ("reads that didn't have an
    # endadapter at all"). Kept available for the counts table but excluded
    # from the Failed Reads Breakdown bar.
    "no_endadapter_reads": "No End Adapter",
    # simplex assignment-quality breakdown
    # (index_hopping_reads is an alias for sid_discordant_reads in the demux
    # output; mapping both would double-count, so only SID Discordant is kept.)
    "sid_discordant_reads": "SID Discordant",
    "both_sid_detected_reads": "Both SID Detected",
    "perfect_index_reads": "Perfect Index",
    # bases (duplex only)
    "total_bases": "Consensus Bases",
    "concordant_duplex_bases": "Concordant Bases",
    "discordant_duplex_bases": "Discordant Bases",
    "non_duplex_bases": "Non-duplex Bases",
}


def parse_demux_sample_stats(input_file: Path) -> dict[str, dict[str, Any]]:
    """Parse demux sample_stats.tsv into per-sample canonical metrics.

    Handles both simplex and duplex field names. Returns
    {sample_name: {metric_display_name: value}}.
    """
    logger.info("Parsing demux sample stats: %s", input_file)

    # First pass: read all rows into {display_name: {sample: value}}
    raw: dict[str, dict[str, Any]] = {}
    sample_names: list[str] = []
    with input_file.open() as f:
        reader = _dict_reader(f)
        for row in reader:
            metric = row.get("metric", "")
            if not sample_names:
                sample_names = [k for k in row.keys() if k != "metric"]
            display = _DEMUX_SAMPLE_METRIC_MAP.get(metric)
            # Don't let a later duplicate (e.g. full_reads vs full_duplex_reads
            # both mapping to "Full Reads") overwrite a populated value.
            if display and display not in raw:
                raw[display] = {s: _num(row.get(s, 0)) for s in sample_names}

    per_sample: dict[str, dict[str, Any]] = {s: {} for s in sample_names}
    for display, sample_values in raw.items():
        for s, val in sample_values.items():
            per_sample[s][display] = val

    return per_sample

# ---------------------------------------------------------------------------
# TFE results (tumor_fraction.results.tsv / contamination.results.tsv)
# ---------------------------------------------------------------------------

_TFE_SUMMARY_KEYS = [
    "tumor_fraction",
    "contamination_fraction",
    "mean_vaf",
    "error_rate",
    "total_alt",
    "total_depth",
    "total_variants",
    "variants_detected",
    "passing_variants",
    "passing_variants_detected",
]


def parse_tfe_results(input_file: Path) -> dict[str, Any]:
    """Parse TFE results TSV and extract key summary metrics.

    The file is a 2-column TSV with header ``metric_name\\tvalue``.
    Lines starting with # are metadata comments.

    Returns a flat dict with display-friendly keys (underscores replaced
    with spaces).
    """
    logger.info("Parsing TFE results: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            key = parts[0]
            # Skip header row
            if key == "metric_name":
                continue
            if key in _TFE_SUMMARY_KEYS:
                metrics[key.replace("_", " ")] = _num(parts[1])
    return metrics


def parse_vcfeval_summary(input_file: Path) -> dict[str, Any]:
    """Parse vcfeval summary.tsv and extract key metrics."""
    logger.info("Parsing vcfeval summary: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        rows_converted = [re.sub(r"[ \t]+", "\t", row.strip()) for row in f if not row.startswith("-")]
        for row in _dict_reader(rows_converted):
            if row.get("Threshold", "NA") not in ("NA", "None"):
                threshold = _num(row["Threshold"])
                tp_baseline = _num(row.get("True-pos-baseline", "NA"))
                tp_call = _num(row.get("True-pos-call", "NA"))
                fp = _num(row.get("False-pos", "NA"))
                fn = _num(row.get("False-neg", "NA"))
                precision = _num(row.get("Precision", "NA"))
                sensitivity = _num(row.get("Sensitivity", "NA"))
                f1 = _num(row.get("F-measure", "NA"))
                metrics.update({
                    "Best Threshold": threshold,
                    "True Positives (Baseline)": tp_baseline,
                    "True Positives (Call)": tp_call,
                    "False Positives": fp,
                    "False Negatives": fn,
                    "Precision": precision,
                    "Sensitivity": sensitivity,
                    "F1 Score": f1
                })

    return metrics

def parse_vcfeval_roc(input_file: Path) -> list[list[Any]]:
    """Parse vcfeval ROC tsv and generate data for ROC curve."""
    logger.info("Parsing vcfeval ROC data: %s", input_file)
    roc_data = []
    with gzip.open(input_file, 'rt') as f:
        rows_converted = [re.sub(r"[ \t]+", "\t", row.strip()).replace('#', '') for row in f if not row.startswith("#") or (row.startswith("#") and "false_positives" in row)]
        for row in _dict_reader(rows_converted):
            roc_data.append([_num(row["false_positives"]), _num(row["true_positives_call"])])
    return roc_data

def parse_vcfeval_summary_from_roc(input_file: Path) -> dict[str, Any]:
    """Parse vcfeval ROC tsv and extract the best-threshold row by max F1."""
    logger.info("Parsing vcfeval summary from ROC: %s", input_file)
    current_best_f1 = 0
    current_best = {}
    with gzip.open(input_file, 'rt') as f:
        rows_converted = [re.sub(r"[ \t]+", "\t", row.strip()).replace('#', '') for row in f if not row.startswith("#") or (row.startswith("#") and "false_positives" in row)]
        for row in _dict_reader(rows_converted):
            if row.get("score", "NA") not in ("NA", "None") and row.get("f_measure", "NA") not in ("NA", "None"):
                score = _num(row["score"])
                f_measure = _num(row["f_measure"])
                if f_measure > current_best_f1:
                    current_best_f1 = f_measure
                    current_best = {
                        "Best Threshold": score,
                        "True Positives (Baseline)": _num(row.get("true_positives_baseline", "NA")),
                        "True Positives (Call)": _num(row.get("true_positives_call", "NA")),
                        "False Positives": _num(row.get("false_positives", "NA")),
                        "False Negatives": _num(row.get("false_negatives", "NA")),
                        "Precision": _num(row.get("precision", "NA")),
                        "Sensitivity": _num(row.get("sensitivity", "NA")),
                        "F1 Score": f_measure
                    }
    return current_best


# ---------------------------------------------------------------------------
# STR caller metrics (metrics.tsv)
# ---------------------------------------------------------------------------

def parse_str_metrics(input_file: Path) -> dict[str, Any]:
    """Parse STR caller metrics.tsv.

    The file is a two-column TSV (name, value) with optional comment/blank
    lines.  Returns a flat dict with metric names as keys.
    """
    logger.info("Parsing STR caller metrics: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        for line in f:
            if line.startswith("##") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                metrics[parts[0]] = _num(parts[1])
    return metrics


# ---------------------------------------------------------------------------
# Copy number caller metrics (*_metrics.tsv)
# ---------------------------------------------------------------------------

def parse_cnc_metrics(input_file: Path) -> dict[str, Any]:
    """Parse copy number caller metrics TSV.

    The file is a two-column TSV (Metric, Value) with a header row.
    Returns a flat dict with metric names as keys.
    """
    logger.info("Parsing CNC metrics: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        for line in f:
            if line.startswith("##") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2 and parts[0] != "Metric":
                metrics[parts[0]] = _num(parts[1])
    return metrics


# ---------------------------------------------------------------------------
# Read collapser metrics (summary_stats.tsv)
# ---------------------------------------------------------------------------

# Metrics selected for the MultiQC summary table. Maps raw metric_name
# to a display-friendly label. Consensus-only metrics (consensus_clusters,
# consensus_median_cluster_size, total_consensus_reads) will be absent
# for markdup runs — the parser simply omits them from the output dict,
# resulting in fewer columns in the MultiQC table for those samples.
_READ_COLLAPSER_SUMMARY_METRICS = {
    "input_reads": "Input Reads",
    "discarded_total_reads": "Discarded Reads",
    "clustering_input_reads": "Clustering Input Reads",
    "duplicate_reads": "Duplicate Reads",
    "total_clusters": "Total Clusters",
    "singleton_clusters": "Singleton Clusters",
    "consensus_clusters": "Consensus Clusters",
    "consensus_median_cluster_size": "Consensus Median Cluster Size",
    "total_consensus_reads": "Total Consensus Reads",
}


def parse_read_collapser_summary(input_file: Path) -> dict[str, Any]:
    """Extract key metrics from read_collapser summary_stats.tsv.

    Returns a flat dict with display-friendly keys for the MultiQC table.
    Percentage values are included as separate "<metric> (%)" entries.
    """
    logger.info("Parsing read collapser summary: %s", input_file)
    metrics: dict[str, Any] = {}
    with input_file.open() as f:
        for row in _dict_reader(f):
            name = row.get("metric_name", "")
            display = _READ_COLLAPSER_SUMMARY_METRICS.get(name)
            if display is None:
                continue
            value = _num(row.get("value", "NA"))
            pct = _num(row.get("percentage", "NA"))
            metrics[display] = value
            if pct != "NA":
                metrics[f"{display} (%)"] = pct
    return metrics


# ---------------------------------------------------------------------------
# Read collapser cluster size distributions (cluster_size_distributions.tsv)
# ---------------------------------------------------------------------------

def _parse_cluster_size(raw: str) -> int | float | None:
    """Convert a cluster_sizes value to a number.

    Handles binned values like ``"50+"`` by stripping the trailing ``+``
    so they become numeric (e.g. 50).  Returns ``None`` for values that
    cannot be converted.
    """
    val = _num(raw)
    if isinstance(val, (int, float)):
        return val
    # Handle "50+", "100+" style bins
    stripped = raw.rstrip("+")
    if stripped != raw:
        val = _num(stripped)
        if isinstance(val, (int, float)):
            return val
    return None


def parse_read_collapser_cluster_size_distributions(
    input_file: Path,
) -> dict[str, list[list[int | float]]]:
    """Parse cluster_size_distributions.tsv for MultiQC line graphs.

    Returns ``{distribution_name: [[cluster_size, count], ...]}`` where each
    series is a list of ``[x, y]`` pairs sorted by cluster size.

    MultiQC linegraph data is emitted to JSON, which does not preserve dict
    key order. Per the MultiQC custom-content docs, line data must be a list
    of ``[x, y]`` pairs so the x-axis is plotted in numeric order instead of
    being re-sorted lexicographically (e.g. 1, 10, 11, 2). See
    https://docs.seqera.io/multiqc/custom_content/#multiqc-specific-data-file

    Binned values like ``"50+"`` are converted to their numeric prefix
    (e.g. 50) so they sort and plot correctly.
    """
    logger.info("Parsing read collapser cluster size distributions: %s", input_file)
    distributions: dict[str, dict[int | float, int | float]] = {}
    with input_file.open() as f:
        for row in _dict_reader(f):
            cluster_size = _parse_cluster_size(row.get("cluster_sizes", "NA"))
            if cluster_size is None:
                continue
            for col, val in row.items():
                if col == "cluster_sizes":
                    continue
                numeric_val = _num(val)
                if numeric_val == "NA":
                    continue
                if col not in distributions:
                    distributions[col] = {}
                distributions[col][cluster_size] = numeric_val
    return {
        name: [[x, points[x]] for x in sorted(points)]
        for name, points in distributions.items()
    }
