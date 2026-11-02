"""Unit tests for accuracy metric parsers.

Covers ``parse_accuracy_metrics`` (base_level_accuracy_summary.tsv) and
``parse_hp_accuracy`` (hp_errors.tsv), including the no-data regression:
when the base-quality filter excludes every base, ``alignment_metrics``
emits a Phred value capped at 93 even though the error rate is undefined
(``percentage = NA``). That capped 93 must not be reported as genuine
accuracy; it is surfaced as the ``NA`` no-data state instead.
"""

import textwrap
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from parsers import parse_accuracy_metrics, parse_hp_accuracy


# ---------------------------------------------------------------------------
# Fixtures - base_level_accuracy_summary.tsv
# ---------------------------------------------------------------------------

# Normal case: bases pass the filter, error rate is defined.
POPULATED_ACCURACY = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    type\tcount\tdenominator\tpercentage\tphred
    all_errors_with_indel_bases\t1000\ttotal_aligned_bases\t0.0100\t40
    indel_bases\t100\ttotal_aligned_bases\t0.0010\t50
    substitutions\t900\ttotal_aligned_bases_minus_indel_bases\t0.0090\t45
""")

# No-data case: every base excluded by the quality filter. count=0 and the
# error rate is undefined (percentage=NA), but the tool still caps phred at 93.
ZERO_OBSERVATION_ACCURACY = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    type\tcount\tdenominator\tpercentage\tphred
    all_errors_with_indel_bases\t0\ttotal_aligned_bases\tNA\t93
    indel_bases\t0\ttotal_aligned_bases\tNA\t93
    substitutions\t0\ttotal_aligned_bases_minus_indel_bases\tNA\t93
""")

# A genuinely error-free sample: bases were evaluated (percentage=0), so the
# high Phred is real and must be preserved - this is NOT the no-data state.
ERROR_FREE_ACCURACY = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    type\tcount\tdenominator\tpercentage\tphred
    all_errors_with_indel_bases\t0\ttotal_aligned_bases\t0.0\t93
    indel_bases\t0\ttotal_aligned_bases\t0.0\t93
    substitutions\t0\ttotal_aligned_bases_minus_indel_bases\t0.0\t93
""")


# ---------------------------------------------------------------------------
# Fixtures - hp_errors.tsv
# ---------------------------------------------------------------------------

POPULATED_HP = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    hp_base\thp_length\tcount\tpercentage\tphred
    any\t8\t500\t0.05\t33
    any\t10\t250\t0.02\t37
""")

ZERO_OBSERVATION_HP = textwrap.dedent("""\
    ##RocheCommandLine=<ID=alignment_metrics,Version="1.0.0">
    hp_base\thp_length\tcount\tpercentage\tphred
    any\t8\t0\tNA\t93
    any\t10\t0\tNA\t93
""")


# ---------------------------------------------------------------------------
# parse_accuracy_metrics
# ---------------------------------------------------------------------------

class TestParseAccuracyMetrics:
    def test_populated_extracts_phred(self, tmp_path):
        f = tmp_path / "base_level_accuracy_summary.tsv"
        f.write_text(POPULATED_ACCURACY)
        result = parse_accuracy_metrics(f)

        assert result == {
            "Overall Phred-scale Accuracy": 40,
            "Indel Phred-scale Accuracy": 50,
            "Substitution Phred-scale Accuracy": 45,
        }

    def test_zero_observation_reports_na_not_q93(self, tmp_path):
        """Regression: no evaluated bases must yield NA, not 93."""
        f = tmp_path / "base_level_accuracy_summary.tsv"
        f.write_text(ZERO_OBSERVATION_ACCURACY)
        result = parse_accuracy_metrics(f)

        assert result == {
            "Overall Phred-scale Accuracy": "NA",
            "Indel Phred-scale Accuracy": "NA",
            "Substitution Phred-scale Accuracy": "NA",
        }
        # The meaningless capped Phred must never surface as a real value.
        assert 93 not in result.values()

    def test_error_free_sample_preserves_high_phred(self, tmp_path):
        """A defined 0% error rate is real accuracy and must be kept."""
        f = tmp_path / "base_level_accuracy_summary.tsv"
        f.write_text(ERROR_FREE_ACCURACY)
        result = parse_accuracy_metrics(f)

        assert result == {
            "Overall Phred-scale Accuracy": 93,
            "Indel Phred-scale Accuracy": 93,
            "Substitution Phred-scale Accuracy": 93,
        }


# ---------------------------------------------------------------------------
# parse_hp_accuracy
# ---------------------------------------------------------------------------

class TestParseHpAccuracy:
    def test_populated_extracts_phred(self, tmp_path):
        f = tmp_path / "hp_errors.tsv"
        f.write_text(POPULATED_HP)
        result = parse_hp_accuracy(f)

        assert result == {
            "8-base HP Phred-scale Accuracy": 33,
            "10-base HP Phred-scale Accuracy": 37,
        }

    def test_zero_observation_reports_na_not_q93(self, tmp_path):
        """Regression for HP accuracy: no data must yield NA."""
        f = tmp_path / "hp_errors.tsv"
        f.write_text(ZERO_OBSERVATION_HP)
        result = parse_hp_accuracy(f)

        assert result == {
            "8-base HP Phred-scale Accuracy": "NA",
            "10-base HP Phred-scale Accuracy": "NA",
        }
        assert 93 not in result.values()
