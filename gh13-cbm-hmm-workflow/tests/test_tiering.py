from __future__ import annotations

import pytest

from gh13_cbm_hmm.tiering import (
    BELOW_OPERATIONAL_CUTOFF,
    PASSES_OPERATIONAL_CUTOFF,
    UNCALIBRATED_FAMILY,
    classify_score,
    distinct_domains,
    overlap_fraction_of_shorter,
)


def _hit(start: int, end: int, independent_evalue: float) -> dict[str, float | int]:
    return {
        "target_start_amino_acid_position": start,
        "target_end_amino_acid_position": end,
        "independent_domain_evalue": independent_evalue,
    }


def test_score_equal_to_cutoff_passes() -> None:
    assert classify_score("CBM20", 31.5, {"CBM20": 31.5}) == PASSES_OPERATIONAL_CUTOFF


def test_score_below_cutoff_and_unknown_family_have_distinct_labels() -> None:
    cutoffs = {"CBM20": 31.5}

    assert classify_score("CBM20", 31.49, cutoffs) == BELOW_OPERATIONAL_CUTOFF
    assert classify_score("CBM41", 100.0, cutoffs) == UNCALIBRATED_FAMILY


def test_overlap_fraction_uses_one_based_inclusive_coordinates() -> None:
    first = _hit(10, 20, 1e-20)
    second = _hit(20, 29, 1e-15)

    assert overlap_fraction_of_shorter(first, second) == pytest.approx(0.1)


def test_distinct_domains_retains_best_evalue_in_overlapping_region() -> None:
    weaker = _hit(10, 30, 1e-10)
    stronger = _hit(12, 32, 1e-30)
    separate = _hit(60, 80, 1e-12)

    retained = distinct_domains(
        [weaker, separate, stronger],
        maximum_overlap_fraction=0.5,
    )

    assert retained == [stronger, separate]


def test_invalid_coordinates_and_overlap_parameter_are_rejected() -> None:
    with pytest.raises(ValueError, match="start"):
        overlap_fraction_of_shorter(_hit(20, 10, 1e-10), _hit(1, 5, 1e-10))
    with pytest.raises(ValueError, match="between zero and one"):
        distinct_domains([], maximum_overlap_fraction=1.1)
