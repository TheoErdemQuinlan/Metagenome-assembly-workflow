"""Cutoff calibration and internal robustness statistics."""

from __future__ import annotations

import math
import statistics
from collections.abc import Iterable

Z_95 = 1.959963984540054


def _finite_scores(values: Iterable[float], label: str) -> list[float]:
    scores = [float(value) for value in values]
    if not scores:
        raise ValueError(f"At least one {label} score is required")
    if not all(math.isfinite(value) for value in scores):
        raise ValueError(f"All {label} scores must be finite")
    return scores


def confusion_at_threshold(
    positive_scores: Iterable[float],
    negative_scores: Iterable[float],
    threshold_bits: float,
) -> dict[str, float | int]:
    """Return the exact confusion matrix for ``score >= threshold_bits``."""

    positives = _finite_scores(positive_scores, "positive")
    negatives = _finite_scores(negative_scores, "negative")
    threshold = float(threshold_bits)
    if not math.isfinite(threshold):
        raise ValueError("threshold_bits must be finite")

    correctly_detected_positives = sum(score >= threshold for score in positives)
    incorrectly_accepted_negatives = sum(score >= threshold for score in negatives)
    missed_positives = len(positives) - correctly_detected_positives
    correctly_rejected_negatives = len(negatives) - incorrectly_accepted_negatives
    return {
        "bit_score_threshold_bits": threshold,
        "correctly_detected_positives": correctly_detected_positives,
        "missed_positives": missed_positives,
        "incorrectly_accepted_negatives": incorrectly_accepted_negatives,
        "correctly_rejected_negatives": correctly_rejected_negatives,
        "sensitivity": correctly_detected_positives / len(positives),
        "specificity": correctly_rejected_negatives / len(negatives),
        "false_positive_rate": incorrectly_accepted_negatives / len(negatives),
        "false_negative_rate": missed_positives / len(positives),
    }


def roc_curve(
    positive_scores: Iterable[float],
    negative_scores: Iterable[float],
) -> list[dict[str, float | int]]:
    """Evaluate the classifier at every observed score and two sentinels."""

    positives = _finite_scores(positive_scores, "positive")
    negatives = _finite_scores(negative_scores, "negative")
    observed = sorted(set(positives) | set(negatives))
    lower = math.nextafter(observed[0], -math.inf)
    upper = math.nextafter(observed[-1], math.inf)
    thresholds = [lower, *observed, upper]
    return [confusion_at_threshold(positives, negatives, value) for value in thresholds]


def confusion_from_curve(
    curve: list[dict[str, float | int]],
    threshold_bits: float,
) -> dict[str, float | int]:
    """Read an exact ``>=`` decision from a complete step-function curve.

    Between observed scores, the correct row is the smallest stored threshold
    greater than or equal to the requested threshold.  Choosing the nearest
    stored threshold is incorrect and can change the confusion matrix.
    """

    if not curve:
        raise ValueError("curve must contain at least one row")
    threshold = float(threshold_bits)
    if not math.isfinite(threshold):
        raise ValueError("threshold_bits must be finite")
    ordered = sorted(curve, key=lambda row: float(row["bit_score_threshold_bits"]))
    for row in ordered:
        stored_threshold = float(row["bit_score_threshold_bits"])
        if not math.isfinite(stored_threshold):
            raise ValueError("curve thresholds must be finite")
        if stored_threshold >= threshold:
            return row
    return ordered[-1]


def roc_auc(curve: list[dict[str, float | int]]) -> float:
    """Calculate trapezoidal area under a receiver-operating curve."""

    if not curve:
        raise ValueError("curve must contain at least one row")
    points = sorted(
        {(float(row["false_positive_rate"]), float(row["sensitivity"])) for row in curve}
    )
    if any(
        not math.isfinite(x) or not math.isfinite(y) or not 0 <= x <= 1 or not 0 <= y <= 1
        for x, y in points
    ):
        raise ValueError("curve rates must be finite proportions between zero and one")
    area = 0.0
    for (x0, y0), (x1, y1) in zip(points, points[1:], strict=False):
        area += (x1 - x0) * (y0 + y1) / 2
    return area


def select_zero_observed_false_positive_cutoff(
    curve: list[dict[str, float | int]],
) -> dict[str, float | int]:
    """Select the lowest evaluated score admitting no observed negatives."""

    eligible = [row for row in curve if row["incorrectly_accepted_negatives"] == 0]
    if not eligible:
        raise ValueError("No evaluated threshold excludes every negative")
    return min(eligible, key=lambda row: float(row["bit_score_threshold_bits"]))


def wilson_interval(
    successes: int,
    trials: int,
    z: float = Z_95,
) -> tuple[float, float]:
    """Return a two-sided Wilson score interval as proportions."""

    if isinstance(successes, bool) or not isinstance(successes, int):
        raise TypeError("successes must be an integer")
    if isinstance(trials, bool) or not isinstance(trials, int):
        raise TypeError("trials must be an integer")
    if trials <= 0:
        raise ValueError("trials must be greater than zero")
    if successes < 0 or successes > trials:
        raise ValueError("successes must be between zero and trials")
    if not math.isfinite(z) or z <= 0:
        raise ValueError("z must be finite and greater than zero")
    estimate = successes / trials
    denominator = 1 + z**2 / trials
    centre = (estimate + z**2 / (2 * trials)) / denominator
    margin = (
        z * math.sqrt(estimate * (1 - estimate) / trials + z**2 / (4 * trials**2)) / denominator
    )
    return max(0.0, centre - margin), min(1.0, centre + margin)


def score_distribution_summary(scores: Iterable[float]) -> dict[str, float | int]:
    """Summarize a score distribution in bits using the standard median.

    :func:`statistics.median` averages the two central observations for an
    even-sized set. Selecting the upper central observation would bias the
    reported median upward.
    """

    values = _finite_scores(scores, "score")
    return {
        "number_of_scores": len(values),
        "minimum_bit_score_bits": min(values),
        "median_bit_score_bits": float(statistics.median(values)),
        "mean_bit_score_bits": float(statistics.fmean(values)),
        "maximum_bit_score_bits": max(values),
    }


def summarize_threshold(
    positive_scores: Iterable[float],
    negative_scores: Iterable[float],
    threshold_bits: float,
) -> dict[str, float | int | list[float]]:
    """Return counts, rates and Wilson 95% intervals at one threshold."""

    positives = _finite_scores(positive_scores, "positive")
    negatives = _finite_scores(negative_scores, "negative")
    result = confusion_at_threshold(positives, negatives, threshold_bits)
    sensitivity_interval = wilson_interval(
        int(result["correctly_detected_positives"]), len(positives)
    )
    specificity_interval = wilson_interval(
        int(result["correctly_rejected_negatives"]), len(negatives)
    )
    result["sensitivity_wilson_95_interval"] = list(sensitivity_interval)
    result["specificity_wilson_95_interval"] = list(specificity_interval)
    result["number_of_positive_references"] = len(positives)
    result["number_of_negative_references"] = len(negatives)
    return result


def calibration_report(
    positive_scores: Iterable[float],
    negative_scores: Iterable[float],
    *,
    comparison_threshold_bits: float | None = None,
    robustness_delta_bits: float = 5.0,
) -> dict:
    """Build a complete, sequence-free cutoff calibration report."""

    if not math.isfinite(robustness_delta_bits) or robustness_delta_bits < 0:
        raise ValueError("robustness_delta_bits must be finite and non-negative")
    positives = _finite_scores(positive_scores, "positive")
    negatives = _finite_scores(negative_scores, "negative")
    curve = roc_curve(positives, negatives)
    selected = select_zero_observed_false_positive_cutoff(curve)
    cutoff = float(selected["bit_score_threshold_bits"])
    report = {
        "selection_rule": (
            "lowest evaluated bit-score threshold with zero observed false "
            "positives in the calibration negatives"
        ),
        "scope": "internal calibration; not independent external validation",
        "receiver_operating_characteristic_auc": roc_auc(curve),
        "positive_score_distribution": score_distribution_summary(positives),
        "negative_score_distribution": score_distribution_summary(negatives),
        "selected_cutoff": summarize_threshold(positives, negatives, cutoff),
        "robustness": [
            summarize_threshold(positives, negatives, cutoff - robustness_delta_bits),
            summarize_threshold(positives, negatives, cutoff),
            summarize_threshold(positives, negatives, cutoff + robustness_delta_bits),
        ],
        "curve": curve,
    }
    if comparison_threshold_bits is not None:
        report["comparison_threshold"] = summarize_threshold(
            positives,
            negatives,
            comparison_threshold_bits,
        )
    return report
