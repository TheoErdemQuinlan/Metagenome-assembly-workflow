from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import Mock

import pytest

import gh13_cbm_hmm.scan as scan_module
from gh13_cbm_hmm.scan import (
    DomainHit,
    _collect_domain_hits,
    _pipeline_options,
    _validate_thresholds,
    iter_hit_rows,
    profile_name,
    scan_fasta,
)


def _domain(
    independent_evalue: float,
    score_bits: float,
    hmm_from: int,
    hmm_to: int,
    hmm_length: int,
    target_from: int = 11,
    target_to: int = 40,
) -> SimpleNamespace:
    alignment = SimpleNamespace(
        hmm_from=hmm_from,
        hmm_to=hmm_to,
        hmm_length=hmm_length,
        target_from=target_from,
        target_to=target_to,
    )
    return SimpleNamespace(
        i_evalue=independent_evalue,
        score=score_bits,
        alignment=alignment,
    )


class _ProfileHits:
    def __init__(self, profile: str, targets: list[SimpleNamespace]):
        self.query = SimpleNamespace(name=profile.encode())
        self._targets = targets

    def __iter__(self):
        return iter(self._targets)


def test_domain_collection_applies_independent_evalue_and_hmm_coverage() -> None:
    retained = _domain(1e-20, 73.25, 1, 40, 100)
    weak_evalue = _domain(1e-4, 90.0, 1, 90, 100)
    short_alignment = _domain(1e-30, 50.0, 1, 20, 100)
    target = SimpleNamespace(
        name=b"abstract_target",
        evalue=1e-22,
        length=100,
        domains=[retained, weak_evalue, short_alignment],
    )
    for domain in target.domains:
        domain.alignment.target_length = target.length
    results = [_ProfileHits("CBM20.hmm", [target])]

    hits = _collect_domain_hits(
        results,
        independent_domain_evalue_threshold=1e-15,
        minimum_hmm_coverage=0.35,
    )

    assert hits == [
        DomainHit(
            profile_name="CBM20",
            target_id="abstract_target",
            target_full_sequence_evalue=1e-22,
            independent_domain_evalue=1e-20,
            domain_bit_score_bits=73.25,
            hmm_start_match_state=1,
            hmm_end_match_state=40,
            hmm_length_match_states=100,
            hmm_alignment_coverage_fraction=0.4,
            target_start_amino_acid_position=11,
            target_end_amino_acid_position=40,
            target_length_amino_acids=100,
        )
    ]


def test_domain_collection_does_not_reapply_target_reporting_threshold() -> None:
    domain = _domain(1e-20, 60.0, 1, 60, 100)
    domain.alignment.target_length = 100
    target = SimpleNamespace(
        name=b"abstract_target",
        evalue=1e6,
        length=100,
        domains=[domain],
    )

    hits = _collect_domain_hits(
        [_ProfileHits("CBM20.hmm", [target])],
        independent_domain_evalue_threshold=1e-15,
        minimum_hmm_coverage=0.35,
    )

    assert len(hits) == 1
    assert hits[0].target_full_sequence_evalue == 1e6
    assert hits[0].independent_domain_evalue == 1e-20


def test_streamed_scan_rejects_terminal_stop_during_fasta_preflight(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta_path = tmp_path / "abstract.input"
    fasta_path.write_text("placeholder\n", encoding="utf-8")
    monkeypatch.setattr(
        scan_module,
        "_iter_text_fasta",
        lambda _path: iter([("abstract_target", "opaque_payload*")]),
    )
    sequence_file = Mock()
    monkeypatch.setattr(scan_module.pyhmmer.easel, "SequenceFile", sequence_file)

    with pytest.raises(ValueError, match="terminal stop marker"):
        scan_fasta(
            fasta_path,
            [SimpleNamespace(name=b"CBM20")],
            evalue_search_space=875.0,
        )

    sequence_file.assert_not_called()


def test_hit_rows_express_hmm_coverage_as_percent() -> None:
    hit = DomainHit(
        profile_name="CBM20",
        target_id="abstract_target",
        target_full_sequence_evalue=1e-22,
        independent_domain_evalue=1e-20,
        domain_bit_score_bits=70.0,
        hmm_start_match_state=1,
        hmm_end_match_state=80,
        hmm_length_match_states=100,
        hmm_alignment_coverage_fraction=0.8,
        target_start_amino_acid_position=5,
        target_end_amino_acid_position=90,
        target_length_amino_acids=100,
    )

    row = next(iter_hit_rows([hit]))

    assert "hmm_alignment_coverage_fraction" not in row
    assert row["hmm_alignment_coverage_percent"] == pytest.approx(80.0)


def test_pipeline_options_include_explicit_search_space() -> None:
    explicit = _pipeline_options(
        target_reporting_evalue_threshold=1e-3,
        domain_reporting_evalue_threshold=10.0,
        evalue_search_space=875,
    )

    assert explicit == {"E": 1e-3, "domE": 10.0, "Z": 875.0}


@pytest.mark.parametrize(
    "arguments",
    [
        (0.0, 10.0, 1e-15, 0.35, 1, 875.0),
        (10.0, 0.0, 1e-15, 0.35, 1, 875.0),
        (10.0, 10.0, 0.0, 0.35, 1, 875.0),
        (10.0, 10.0, 1e-15, -0.1, 1, 875.0),
        (10.0, 10.0, 1e-15, 1.1, 1, 875.0),
        (10.0, 10.0, 1e-15, 0.35, -1, 875.0),
        (10.0, 10.0, 1e-15, 0.35, 1, 0.0),
    ],
)
def test_invalid_scan_parameters_are_rejected(arguments: tuple[float, ...]) -> None:
    with pytest.raises(ValueError):
        _validate_thresholds(*arguments)


def test_profile_name_removes_only_a_trailing_hmm_suffix() -> None:
    assert profile_name(SimpleNamespace(name=b"GH13.hmm")) == "GH13"
    assert profile_name(SimpleNamespace(name=b"hmm_family")) == "hmm_family"
