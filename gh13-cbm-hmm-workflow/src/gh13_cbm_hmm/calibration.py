"""Score labelled references and construct sequence-free calibration reports."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from .references import (
    ReferenceRecord,
    assert_disjoint_reference_classes,
    reference_set_manifest,
)
from .scan import DomainHit, profile_name, scan_records
from .validation import calibration_report, score_distribution_summary


def best_profile_scores(
    hits: Iterable[DomainHit],
    expected_target_ids: Iterable[str],
) -> dict[str, float]:
    """Return the highest domain bit score per target, assigning zero to no hit.

    Zero represents absence of a domain reported under the calibration search
    settings. It is not a measured alignment score.
    """

    target_ids = [str(target_id) for target_id in expected_target_ids]
    if len(target_ids) != len(set(target_ids)):
        raise ValueError("expected_target_ids must be unique")
    expected = set(target_ids)
    best = {target_id: 0.0 for target_id in target_ids}
    for hit in hits:
        if hit.target_id not in expected:
            raise ValueError(f"Unexpected target identifier in HMM result: {hit.target_id}")
        best[hit.target_id] = max(best[hit.target_id], hit.domain_bit_score_bits)
    return best


def _labelled_records(
    positives: list[ReferenceRecord],
    hard_negatives: list[ReferenceRecord],
    easy_negatives: list[ReferenceRecord],
) -> list[tuple[str, ReferenceRecord, str]]:
    labelled: list[tuple[str, ReferenceRecord, str]] = []
    for label, records in (
        ("positive", positives),
        ("hard_negative", hard_negatives),
        ("easy_negative", easy_negatives),
    ):
        for index, record in enumerate(records):
            target_id = f"{label}:{index:08d}"
            labelled.append((target_id, record, label))
    return labelled


def calibrate_profile(
    family: str,
    hmm: Any,
    positives: list[ReferenceRecord],
    hard_negatives: list[ReferenceRecord],
    easy_negatives: list[ReferenceRecord],
    *,
    evalue_search_space: float,
    comparison_threshold_bits: float | None = None,
    robustness_delta_bits: float = 5.0,
    calibration_target_reporting_evalue_threshold: float = 10.0,
    calibration_domain_reporting_evalue_threshold: float = 10.0,
    calibration_independent_domain_evalue_threshold: float = 10.0,
    cpus: int = 0,
) -> dict[str, Any]:
    """Calibrate one family profile against labelled reference classes.

    The calibration search deliberately uses permissive reporting thresholds
    and no HMM-coverage filter so that score discrimination can be evaluated.
    ``evalue_search_space`` is still explicit and is recorded in the result.
    The selected boundary is an internal operational cutoff, not independent
    biological validation.
    """

    normalized_family = str(family).strip().removesuffix(".hmm")
    if not normalized_family:
        raise ValueError("family must not be empty")
    hmm_name = profile_name(hmm)
    if hmm_name != normalized_family:
        raise ValueError(
            f"Requested family {normalized_family!r} does not match HMM profile {hmm_name!r}"
        )
    if not positives:
        raise ValueError("At least one positive reference is required")
    if not hard_negatives and not easy_negatives:
        raise ValueError("At least one negative reference is required")
    assert_disjoint_reference_classes(positives, hard_negatives, easy_negatives)

    labelled = _labelled_records(positives, hard_negatives, easy_negatives)
    records = [(target_id, record["sequence"]) for target_id, record, _label in labelled]
    target_ids = [target_id for target_id, _record, _label in labelled]
    hits = scan_records(
        records,
        [hmm],
        evalue_search_space=evalue_search_space,
        target_reporting_evalue_threshold=calibration_target_reporting_evalue_threshold,
        domain_reporting_evalue_threshold=calibration_domain_reporting_evalue_threshold,
        independent_domain_evalue_threshold=(calibration_independent_domain_evalue_threshold),
        minimum_hmm_coverage=0.0,
        cpus=cpus,
    )
    scores_by_target = best_profile_scores(hits, target_ids)
    positive_scores = [
        scores_by_target[target_id] for target_id, _record, label in labelled if label == "positive"
    ]
    hard_negative_scores = [
        scores_by_target[target_id]
        for target_id, _record, label in labelled
        if label == "hard_negative"
    ]
    easy_negative_scores = [
        scores_by_target[target_id]
        for target_id, _record, label in labelled
        if label == "easy_negative"
    ]
    negative_scores = hard_negative_scores + easy_negative_scores
    report = calibration_report(
        positive_scores,
        negative_scores,
        comparison_threshold_bits=comparison_threshold_bits,
        robustness_delta_bits=robustness_delta_bits,
    )
    report.update(
        {
            "family": normalized_family,
            "profile_length_match_states": int(hmm.M),
            "reference_set_manifest": reference_set_manifest(
                positives, hard_negatives, easy_negatives
            ),
            "reference_class_counts": {
                "positive": len(positive_scores),
                "hard_negative": len(hard_negative_scores),
                "easy_negative": len(easy_negative_scores),
            },
            "hard_negative_score_distribution": score_distribution_summary(hard_negative_scores)
            if hard_negative_scores
            else None,
            "easy_negative_score_distribution": score_distribution_summary(easy_negative_scores)
            if easy_negative_scores
            else None,
            "number_of_positive_references_without_reported_domain": sum(
                score == 0.0 for score in positive_scores
            ),
            "number_of_negative_references_without_reported_domain": sum(
                score == 0.0 for score in negative_scores
            ),
            "calibration_search_parameters": {
                "target_reporting_evalue_threshold": (
                    calibration_target_reporting_evalue_threshold
                ),
                "domain_reporting_evalue_threshold": (
                    calibration_domain_reporting_evalue_threshold
                ),
                "independent_domain_evalue_threshold": (
                    calibration_independent_domain_evalue_threshold
                ),
                "minimum_hmm_coverage_fraction": 0.0,
                "evalue_search_space_Z": evalue_search_space,
                "cpu_worker_count": cpus,
                "missing_reported_domain_score_bits": 0.0,
            },
        }
    )
    return report
