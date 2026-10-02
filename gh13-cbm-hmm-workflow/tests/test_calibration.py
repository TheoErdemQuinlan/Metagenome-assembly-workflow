from __future__ import annotations

from types import SimpleNamespace

import pytest

from gh13_cbm_hmm import calibration
from gh13_cbm_hmm.calibration import best_profile_scores, calibrate_profile
from gh13_cbm_hmm.scan import DomainHit


def _hit(target_id: str, score_bits: float) -> DomainHit:
    return DomainHit(
        profile_name="profile_alpha",
        target_id=target_id,
        target_full_sequence_evalue=1e-20,
        independent_domain_evalue=1e-18,
        domain_bit_score_bits=score_bits,
        hmm_start_match_state=1,
        hmm_end_match_state=80,
        hmm_length_match_states=100,
        hmm_alignment_coverage_fraction=0.8,
        target_start_amino_acid_position=2,
        target_end_amino_acid_position=81,
        target_length_amino_acids=120,
    )


def test_best_profile_scores_uses_maximum_and_assigns_zero_to_no_hit() -> None:
    result = best_profile_scores(
        [_hit("target_a", 14.0), _hit("target_a", 17.5)],
        ["target_a", "target_b"],
    )

    assert result == {"target_a": 17.5, "target_b": 0.0}


def test_best_profile_scores_rejects_ambiguous_or_unexpected_targets() -> None:
    with pytest.raises(ValueError, match="unique"):
        best_profile_scores([], ["target_a", "target_a"])
    with pytest.raises(ValueError, match="Unexpected"):
        best_profile_scores([_hit("target_b", 5.0)], ["target_a"])


def test_calibrate_profile_records_search_semantics_and_class_counts(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_scan(records, hmms, **parameters):
        captured["records"] = records
        captured["hmms"] = hmms
        captured["parameters"] = parameters
        return [
            _hit("positive:00000000", 12.0),
            _hit("hard_negative:00000000", 8.0),
        ]

    monkeypatch.setattr(calibration, "scan_records", fake_scan)
    monkeypatch.setattr(
        calibration,
        "assert_disjoint_reference_classes",
        lambda *_classes: None,
    )
    monkeypatch.setattr(
        calibration,
        "reference_set_manifest",
        lambda *_classes: {"cross_class_exact_sequence_overlap_count": 0},
    )
    positive = {
        "accession": "positive_record",
        "tags": ["CBM20"],
        "sequence": "opaque_positive",
        "sequence_sha256": "positive_digest",
    }
    hard_negative = {
        "accession": "hard_record",
        "tags": ["CBM41"],
        "sequence": "opaque_hard",
        "sequence_sha256": "hard_digest",
    }
    easy_negative = {
        "accession": "easy_record",
        "tags": ["GH5"],
        "sequence": "opaque_easy",
        "sequence_sha256": "easy_digest",
    }
    hmm = SimpleNamespace(name=b"profile_alpha", M=100)

    report = calibrate_profile(
        "profile_alpha",
        hmm,
        [positive],
        [hard_negative],
        [easy_negative],
        evalue_search_space=875.0,
        cpus=2,
    )

    assert report["family"] == "profile_alpha"
    assert report["reference_class_counts"] == {
        "positive": 1,
        "hard_negative": 1,
        "easy_negative": 1,
    }
    assert report["number_of_negative_references_without_reported_domain"] == 1
    assert report["calibration_search_parameters"] == {
        "target_reporting_evalue_threshold": 10.0,
        "domain_reporting_evalue_threshold": 10.0,
        "independent_domain_evalue_threshold": 10.0,
        "minimum_hmm_coverage_fraction": 0.0,
        "evalue_search_space_Z": 875.0,
        "cpu_worker_count": 2,
        "missing_reported_domain_score_bits": 0.0,
    }
    assert captured["records"] == [
        ("positive:00000000", "opaque_positive"),
        ("hard_negative:00000000", "opaque_hard"),
        ("easy_negative:00000000", "opaque_easy"),
    ]
