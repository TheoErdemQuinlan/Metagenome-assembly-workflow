from __future__ import annotations

import pytest

from gh13_cbm_hmm.validation import (
    calibration_report,
    confusion_at_threshold,
    confusion_from_curve,
    roc_auc,
    roc_curve,
    score_distribution_summary,
    select_zero_observed_false_positive_cutoff,
    summarize_threshold,
    wilson_interval,
)


def test_confusion_matrix_uses_inclusive_bit_score_boundary() -> None:
    result = confusion_at_threshold(
        positive_scores=[9.0, 7.0, 5.0],
        negative_scores=[8.0, 6.0, 2.0],
        threshold_bits=6.0,
    )

    assert result == {
        "bit_score_threshold_bits": 6.0,
        "correctly_detected_positives": 2,
        "missed_positives": 1,
        "incorrectly_accepted_negatives": 2,
        "correctly_rejected_negatives": 1,
        "sensitivity": pytest.approx(2 / 3),
        "specificity": pytest.approx(1 / 3),
        "false_positive_rate": pytest.approx(2 / 3),
        "false_negative_rate": pytest.approx(1 / 3),
    }


def test_curve_lookup_uses_step_function_not_nearest_observed_score() -> None:
    curve = roc_curve(
        positive_scores=[10.0, 20.0],
        negative_scores=[5.0, 7.0],
    )

    result = confusion_from_curve(curve, threshold_bits=7.6)

    assert result["bit_score_threshold_bits"] == 10.0
    assert result["correctly_detected_positives"] == 2
    assert result["incorrectly_accepted_negatives"] == 0


def test_curve_lookup_agrees_with_direct_calculation_between_scores() -> None:
    positives = [3.0, 11.0, 19.0]
    negatives = [2.0, 8.0, 14.0]
    requested_threshold = 8.25

    direct = confusion_at_threshold(positives, negatives, requested_threshold)
    from_curve = confusion_from_curve(roc_curve(positives, negatives), requested_threshold)

    for field in (
        "correctly_detected_positives",
        "missed_positives",
        "incorrectly_accepted_negatives",
        "correctly_rejected_negatives",
        "sensitivity",
        "specificity",
    ):
        assert from_curve[field] == direct[field]


def test_zero_observed_false_positive_rule_selects_lowest_eligible_score() -> None:
    curve = roc_curve(
        positive_scores=[5.0, 9.0, 12.0],
        negative_scores=[1.0, 4.0, 8.0],
    )

    selected = select_zero_observed_false_positive_cutoff(curve)

    assert selected["bit_score_threshold_bits"] == 9.0
    assert selected["incorrectly_accepted_negatives"] == 0
    assert selected["correctly_detected_positives"] == 2


def test_wilson_interval_for_zero_observed_events_is_not_zero_width() -> None:
    lower, upper = wilson_interval(successes=0, trials=453)

    assert lower == 0.0
    assert upper == pytest.approx(0.008408, abs=1e-6)


def test_threshold_summary_reports_counts_rates_units_and_intervals() -> None:
    summary = summarize_threshold(
        positive_scores=[12.0, 10.0, 4.0],
        negative_scores=[7.0, 3.0],
        threshold_bits=8.0,
    )

    assert summary["bit_score_threshold_bits"] == 8.0
    assert summary["number_of_positive_references"] == 3
    assert summary["number_of_negative_references"] == 2
    assert summary["correctly_detected_positives"] == 2
    assert summary["correctly_rejected_negatives"] == 2
    assert len(summary["sensitivity_wilson_95_interval"]) == 2
    assert len(summary["specificity_wilson_95_interval"]) == 2


def test_calibration_report_labels_internal_scope_and_robustness_points() -> None:
    report = calibration_report(
        positive_scores=[12.0, 10.0, 7.0],
        negative_scores=[6.0, 3.0, 1.0],
        comparison_threshold_bits=8.0,
        robustness_delta_bits=2.0,
    )

    cutoff = report["selected_cutoff"]["bit_score_threshold_bits"]
    evaluated = [row["bit_score_threshold_bits"] for row in report["robustness"]]
    assert report["scope"] == "internal calibration; not independent external validation"
    assert evaluated == [cutoff - 2.0, cutoff, cutoff + 2.0]
    assert report["comparison_threshold"]["bit_score_threshold_bits"] == 8.0
    assert 0.0 <= roc_auc(report["curve"]) <= 1.0


def test_even_length_score_summary_uses_the_standard_median() -> None:
    summary = score_distribution_summary([1.0, 4.0, 8.0, 20.0])

    assert summary["median_bit_score_bits"] == 6.0
    assert summary["mean_bit_score_bits"] == pytest.approx(8.25)


@pytest.mark.parametrize(
    ("successes", "trials"),
    [(-1, 10), (11, 10), (0, 0)],
)
def test_wilson_interval_rejects_invalid_counts(successes: int, trials: int) -> None:
    with pytest.raises(ValueError):
        wilson_interval(successes, trials)


def test_validation_rejects_empty_or_nonfinite_score_collections() -> None:
    with pytest.raises(ValueError, match="positive"):
        confusion_at_threshold([], [1.0], 1.0)
    with pytest.raises(ValueError, match="negative"):
        confusion_at_threshold([1.0], [float("nan")], 1.0)
