from __future__ import annotations

import hashlib
from unittest.mock import Mock

import pytest

from gh13_cbm_hmm import references
from gh13_cbm_hmm.references import (
    ReferencePools,
    build_pools,
    fetch_uniprot_reference_table,
    labelled_set_for_family,
    parse_cazy_header,
)


def _opaque_records() -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    for index in range(40):
        records.append((f">entry_a_{index}|CBM20", f"payload_a_{index}"))
        records.append((f">entry_b_{index}|CBM41", f"payload_b_{index}"))
        records.append((f">entry_h_{index}|GH13", f"payload_h_{index}"))
        records.append((f">entry_e_{index}|GH5", f"payload_e_{index}"))
    return records


def test_parse_cazy_header_supports_accession_first_and_tag_first_layouts() -> None:
    accession_first = parse_cazy_header(">entry_one|CBM20|GH13")
    tag_first = parse_cazy_header(">CBM20|entry_two|GH13_5")

    assert accession_first == ("entry_one", {"CBM20", "GH13"})
    assert tag_first == ("entry_two", {"CBM20", "GH13_5"})


def test_reference_sampling_is_independent_of_requested_family_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = _opaque_records()
    monkeypatch.setattr(references, "iter_fasta", lambda _path: iter(records))

    forward = build_pools(
        "unused-input",
        ["CBM20", "CBM41"],
        positives_per_family=9,
        hard_pool_size=11,
        easy_pool_size=7,
        seed=42,
    )
    reverse = build_pools(
        "unused-input",
        ["CBM41", "CBM20"],
        positives_per_family=9,
        hard_pool_size=11,
        easy_pool_size=7,
        seed=42,
    )

    assert forward == reverse


def test_invalid_reference_sequence_policy_errors_or_records_exclusion(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def erroring_records(_path):
        raise ValueError("abstract invalid record")

    monkeypatch.setattr(references, "iter_fasta", erroring_records)

    with pytest.raises(ValueError, match="abstract invalid record"):
        build_pools(
            "unused-input",
            ["CBM20"],
            positives_per_family=2,
            hard_pool_size=2,
            easy_pool_size=2,
        )

    raw_records = [
        (">invalid_entry|CBM20", "invalid_payload"),
        (">retained_entry|CBM20", "retained_payload"),
    ]
    monkeypatch.setattr(
        references,
        "_iter_raw_fasta",
        lambda _path: iter(raw_records),
    )

    def normalize_or_reject(value, *, identifier):
        if value == "invalid_payload":
            raise ValueError(f"invalid abstract payload for {identifier}")
        return value

    monkeypatch.setattr(references, "normalise_protein_sequence", normalize_or_reject)

    pools = build_pools(
        "unused-input",
        ["CBM20"],
        positives_per_family=2,
        hard_pool_size=2,
        easy_pool_size=2,
        invalid_sequence_policy="exclude",
    )

    assert pools.number_of_source_entries == 2
    assert pools.number_of_entries_excluded_for_invalid_sequence == 1
    assert pools.number_of_unique_source_sequences == 1
    assert pools.positive_totals == {"CBM20": 1}
    assert pools.positive_unique_sequence_totals == {"CBM20": 1}


def test_labelled_sets_are_mutually_exclusive_and_ignore_ambiguous_tags(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        references,
        "normalise_protein_sequence",
        lambda value, *, identifier: value,
    )
    shared = {"accession": "shared", "tags": ["CBM20"], "sequence": "payload_shared"}
    pools = ReferencePools(
        positives={"CBM20": [shared, shared]},
        positive_totals={"CBM20": 2},
        hard_negative_pool=[
            shared,
            {"accession": "other", "tags": ["CBM41"], "sequence": "payload_other"},
            {"accession": "ambiguous", "tags": ["CBM0"], "sequence": "payload_ambiguous"},
        ],
        easy_negative_pool=[{"accession": "easy", "tags": ["GH5"], "sequence": "payload_easy"}],
    )

    positives, hard_negatives, easy_negatives = labelled_set_for_family(
        pools,
        "CBM20",
        number_of_positives=10,
        number_of_hard_negatives=10,
        number_of_easy_negatives=10,
        seed=17,
    )

    assert positives == [shared]
    assert [record["accession"] for record in hard_negatives] == ["other"]
    assert [record["accession"] for record in easy_negatives] == ["easy"]


def test_uniprot_helper_returns_tags_schema(monkeypatch: pytest.MonkeyPatch) -> None:
    response = Mock()
    response.text = "Entry\tCAZy\tSequence\nentry_one\tCBM20; GH13\topaque_payload\n"
    response.raise_for_status = Mock()
    request = Mock(return_value=response)
    monkeypatch.setattr(references.requests, "get", request)
    monkeypatch.setattr(
        references,
        "normalise_protein_sequence",
        lambda value, *, identifier: value,
    )

    rows = fetch_uniprot_reference_table(
        "reviewed:true",
        maximum_records=1,
        delay_seconds=0,
        timeout_seconds=3,
    )

    assert rows == [
        {
            "accession": "entry_one",
            "tags": ["CBM20", "GH13"],
            "sequence": "opaque_payload",
            "sequence_sha256": hashlib.sha256(b"opaque_payload").hexdigest(),
        }
    ]
    assert "cazy" not in rows[0]
    request.assert_called_once()
