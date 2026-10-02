"""Deterministic construction of labelled CAZy reference samples.

No reference sequences are distributed with this package. Callers supply a
locally obtained CAZy FASTA and remain responsible for its licence, release
provenance and biological interpretation.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
import random
import re
import time
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TypedDict

import requests

from .scan import normalise_protein_sequence

_TAG_PATTERN = re.compile(r"^(GH|GT|PL|CE|AA|CBM)\d+(_\d+)?$")
_BASE_FAMILY_PATTERN = re.compile(r"^(GH|GT|PL|CE|AA|CBM)\d+$")

AMBIGUOUS_TAGS = frozenset({"CBM0"})

DEFAULT_HARD_NEGATIVE_FAMILIES = frozenset(
    {
        "CBM1",
        "CBM2",
        "CBM3",
        "CBM4",
        "CBM5",
        "CBM6",
        "CBM9",
        "CBM13",
        "CBM20",
        "CBM32",
        "CBM34",
        "CBM37",
        "CBM41",
        "CBM48",
        "CBM50",
        "CBM69",
        "GH13",
    }
)


class ReferenceRecord(TypedDict):
    """One in-memory labelled-reference record."""

    accession: str
    tags: list[str]
    sequence: str
    sequence_sha256: str


def parse_cazy_header(header: str) -> tuple[str, set[str]]:
    """Parse accession-first and family-tag-first CAZy FASTA headers.

    In the JGI-style layout the first field is itself a family tag. Treating
    it as an accession would silently discard a true label, so every field is
    examined before the accession is selected.
    """

    if not isinstance(header, str):
        raise TypeError("CAZy FASTA headers must be text")
    text = header.strip()
    if text.startswith(">"):
        text = text[1:].strip()
    if not text:
        raise ValueError("CAZy FASTA header is empty")
    parts = [part.strip() for part in text.split("|")]
    if any(not part for part in parts):
        raise ValueError(f"CAZy FASTA header contains an empty field: {text!r}")
    tags = {part for part in parts if _TAG_PATTERN.fullmatch(part)}
    if _TAG_PATTERN.fullmatch(parts[0]):
        non_tag_fields = [part for part in parts if not _TAG_PATTERN.fullmatch(part)]
        if not non_tag_fields:
            raise ValueError(f"CAZy FASTA header has no accession field: {text!r}")
        accession = "|".join(non_tag_fields)
    else:
        accession = parts[0]
    if not accession:
        raise ValueError(f"CAZy FASTA header has no accession: {text!r}")
    return accession, tags


def base_family(tag: str) -> str:
    """Return the family component of a CAZy family or subfamily tag."""

    normalized = str(tag).strip()
    if not _TAG_PATTERN.fullmatch(normalized):
        raise ValueError(f"Invalid CAZy family tag: {tag!r}")
    return normalized.split("_", maxsplit=1)[0]


def has_family(tags: set[str] | frozenset[str], family: str) -> bool:
    """Return whether any tag belongs to ``family`` or one of its subfamilies."""

    normalized_family = _validate_base_family(family)
    return any(base_family(tag) == normalized_family for tag in tags)


def _validate_base_family(family: str) -> str:
    normalized = str(family).strip()
    if not _BASE_FAMILY_PATTERN.fullmatch(normalized):
        raise ValueError(f"Invalid base CAZy family: {family!r}")
    return normalized


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield validated FASTA header and normalized protein-sequence pairs."""

    for header, raw_sequence in _iter_raw_fasta(path):
        yield header, normalise_protein_sequence(raw_sequence, identifier=header)


def _iter_raw_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield raw FASTA records while validating only the container format."""

    fasta_path = Path(path)
    if not fasta_path.is_file():
        raise FileNotFoundError(fasta_path)
    header: str | None = None
    chunks: list[str] = []
    number_of_records = 0
    with fasta_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                    number_of_records += 1
                header = stripped
                chunks = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence text precedes the first FASTA header at line {line_number}"
                    )
                chunks.append(stripped)
    if header is not None:
        yield header, "".join(chunks)
        number_of_records += 1
    if number_of_records == 0:
        raise ValueError(f"Reference FASTA contains no records: {fasta_path}")


class Reservoir:
    """Fixed-memory uniform reservoir sample."""

    def __init__(self, size: int, random_generator: random.Random):
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            raise ValueError("Reservoir size must be a positive integer")
        self.size = size
        self.random_generator = random_generator
        self.number_seen = 0
        self.items: list[ReferenceRecord] = []

    def add(self, item: ReferenceRecord) -> None:
        self.number_seen += 1
        if len(self.items) < self.size:
            self.items.append(item)
            return
        index = self.random_generator.randrange(self.number_seen)
        if index < self.size:
            self.items[index] = item


@dataclass(slots=True)
class ReferencePools:
    """Sampled reference pools and explicit source-population counts."""

    positives: dict[str, list[ReferenceRecord]] = field(default_factory=dict)
    positive_totals: dict[str, int] = field(default_factory=dict)
    positive_unique_sequence_totals: dict[str, int] = field(default_factory=dict)
    hard_negative_pool: list[ReferenceRecord] = field(default_factory=list)
    easy_negative_pool: list[ReferenceRecord] = field(default_factory=list)
    number_of_source_entries: int = 0
    number_of_unique_source_sequences: int = 0
    number_of_duplicate_sequence_entries: int = 0
    number_of_duplicate_sequences_with_differing_tags: int = 0
    number_of_entries_excluded_for_invalid_sequence: int = 0
    random_seed: int = 0


def _sequence_digest(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii", errors="strict")).hexdigest()


def _validated_positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def build_pools(
    fasta_path: str | Path,
    families: list[str],
    *,
    positives_per_family: int = 1500,
    hard_pool_size: int = 12000,
    easy_pool_size: int = 3000,
    seed: int = 20260923,
    hard_negative_families: set[str] | frozenset[str] | None = None,
    invalid_sequence_policy: str = "error",
) -> ReferencePools:
    """Build deterministic pools invariant to requested family-list order.

    Exact sequences are deduplicated before reservoir sampling. The function
    makes two streaming passes: first it unions CAZy tags across exact
    duplicate sequences, then it samples each unique sequence once. This
    prevents file order from assigning a duplicated sequence to contradictory
    classes. The selected reservoirs remain dependent on source-record order,
    so the source checksum is part of the reproducibility record. By default,
    malformed protein sequences stop the run. The explicit ``"exclude"``
    policy instead omits them and records their count in
    :class:`ReferencePools`.
    """

    _validated_positive_integer(positives_per_family, "positives_per_family")
    _validated_positive_integer(hard_pool_size, "hard_pool_size")
    _validated_positive_integer(easy_pool_size, "easy_pool_size")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("seed must be an integer")
    if invalid_sequence_policy not in {"error", "exclude"}:
        raise ValueError("invalid_sequence_policy must be 'error' or 'exclude'")
    normalized_families = [_validate_base_family(family) for family in families]
    if len(normalized_families) != len(set(normalized_families)):
        raise ValueError("families must not contain duplicates")
    ordered_families = sorted(normalized_families)
    if not ordered_families:
        raise ValueError("At least one family is required")
    supplied_hard_families = (
        DEFAULT_HARD_NEGATIVE_FAMILIES if hard_negative_families is None else hard_negative_families
    )
    hard_families = {_validate_base_family(family) for family in supplied_hard_families}

    positive_reservoirs = {
        family: Reservoir(
            positives_per_family,
            random.Random(f"{seed}:positive:{family}"),
        )
        for family in ordered_families
    }
    hard_reservoir = Reservoir(hard_pool_size, random.Random(f"{seed}:hard"))
    easy_reservoir = Reservoir(easy_pool_size, random.Random(f"{seed}:easy"))
    positive_entry_totals = {family: 0 for family in ordered_families}
    positive_unique_totals = {family: 0 for family in ordered_families}

    tags_by_sequence_digest: dict[str, set[str]] = {}
    seen_accessions: dict[str, str] = {}
    number_of_entries = 0
    number_of_invalid_entries = 0
    number_of_duplicates = 0
    differing_tag_duplicate_digests: set[str] = set()

    def reference_records(*, count_source_entries: bool) -> Iterator[tuple[str, str]]:
        nonlocal number_of_entries, number_of_invalid_entries
        if invalid_sequence_policy == "error":
            for candidate_header, normalized_sequence in iter_fasta(fasta_path):
                if count_source_entries:
                    number_of_entries += 1
                yield candidate_header, normalized_sequence
            return
        for candidate_header, raw_sequence in _iter_raw_fasta(fasta_path):
            if count_source_entries:
                number_of_entries += 1
            try:
                normalized_sequence = normalise_protein_sequence(
                    raw_sequence, identifier=candidate_header
                )
            except (TypeError, ValueError):
                if count_source_entries:
                    number_of_invalid_entries += 1
                continue
            yield candidate_header, normalized_sequence

    for header, sequence in reference_records(count_source_entries=True):
        accession, tags = parse_cazy_header(header)
        for family in ordered_families:
            if has_family(tags, family):
                positive_entry_totals[family] += 1

        digest = _sequence_digest(sequence)
        previous_digest = seen_accessions.get(accession)
        if previous_digest is not None and previous_digest != digest:
            raise ValueError(
                f"Reference accession {accession!r} is associated with multiple sequences"
            )
        seen_accessions[accession] = digest
        previous_tags = tags_by_sequence_digest.get(digest)
        if previous_tags is not None:
            if previous_tags != tags:
                differing_tag_duplicate_digests.add(digest)
                previous_tags.update(tags)
            number_of_duplicates += 1
            continue
        tags_by_sequence_digest[digest] = set(tags)

    sampled_digests: set[str] = set()
    for header, sequence in reference_records(count_source_entries=False):
        accession, _entry_tags = parse_cazy_header(header)
        digest = _sequence_digest(sequence)
        if digest in sampled_digests:
            continue
        sampled_digests.add(digest)
        tags = tags_by_sequence_digest[digest]
        ordered_tags = tuple(sorted(tags))

        record: ReferenceRecord = {
            "accession": accession,
            "tags": list(ordered_tags),
            "sequence": sequence,
            "sequence_sha256": digest,
        }
        for family in ordered_families:
            if has_family(tags, family):
                positive_reservoirs[family].add(record)
                positive_unique_totals[family] += 1
        if any(has_family(tags, family) for family in hard_families):
            hard_reservoir.add(record)
        elif tags and not any(tag.startswith("CBM") for tag in tags):
            easy_reservoir.add(record)

    return ReferencePools(
        positives={family: positive_reservoirs[family].items for family in ordered_families},
        positive_totals=positive_entry_totals,
        positive_unique_sequence_totals=positive_unique_totals,
        hard_negative_pool=hard_reservoir.items,
        easy_negative_pool=easy_reservoir.items,
        number_of_source_entries=number_of_entries,
        number_of_unique_source_sequences=len(tags_by_sequence_digest),
        number_of_duplicate_sequence_entries=number_of_duplicates,
        number_of_duplicate_sequences_with_differing_tags=(len(differing_tag_duplicate_digests)),
        number_of_entries_excluded_for_invalid_sequence=number_of_invalid_entries,
        random_seed=seed,
    )


def _record_digest(record: ReferenceRecord) -> str:
    sequence = normalise_protein_sequence(
        record["sequence"], identifier=record.get("accession", "reference")
    )
    calculated = _sequence_digest(sequence)
    recorded = record.get("sequence_sha256")
    if recorded is not None and recorded != calculated:
        raise ValueError(
            f"Sequence checksum does not match for reference {record.get('accession', 'unknown')!r}"
        )
    return calculated


def _deduplicate_by_sequence(records: list[ReferenceRecord]) -> list[ReferenceRecord]:
    seen: set[str] = set()
    retained: list[ReferenceRecord] = []
    for record in records:
        digest = _record_digest(record)
        if digest not in seen:
            seen.add(digest)
            retained.append(record)
    return retained


def _sample_records(
    records: list[ReferenceRecord],
    number_requested: int,
    *,
    seed_text: str,
) -> list[ReferenceRecord]:
    if isinstance(number_requested, bool) or not isinstance(number_requested, int):
        raise TypeError("Requested sample sizes must be integers")
    if number_requested < 0:
        raise ValueError("Requested sample sizes must be non-negative")
    ordered = sorted(records, key=_record_digest)
    random.Random(seed_text).shuffle(ordered)
    return ordered[:number_requested]


def _class_digest_set(records: list[ReferenceRecord]) -> set[str]:
    digests = [_record_digest(record) for record in records]
    if len(digests) != len(set(digests)):
        raise ValueError("A labelled reference class contains duplicate sequences")
    return set(digests)


def assert_disjoint_reference_classes(
    positives: list[ReferenceRecord],
    hard_negatives: list[ReferenceRecord],
    easy_negatives: list[ReferenceRecord],
) -> None:
    """Fail if an exact sequence occurs in more than one label class."""

    class_digests = {
        "positive": _class_digest_set(positives),
        "hard_negative": _class_digest_set(hard_negatives),
        "easy_negative": _class_digest_set(easy_negatives),
    }
    pairs = (
        ("positive", "hard_negative"),
        ("positive", "easy_negative"),
        ("hard_negative", "easy_negative"),
    )
    for first, second in pairs:
        overlap = class_digests[first] & class_digests[second]
        if overlap:
            raise ValueError(
                f"Reference classes {first!r} and {second!r} share {len(overlap)} exact sequence(s)"
            )


def reference_set_manifest(
    positives: list[ReferenceRecord],
    hard_negatives: list[ReferenceRecord],
    easy_negatives: list[ReferenceRecord],
) -> dict[str, object]:
    """Return a sequence-free checksum manifest for three label classes."""

    assert_disjoint_reference_classes(positives, hard_negatives, easy_negatives)

    def class_record(records: list[ReferenceRecord], class_label: str) -> dict[str, object]:
        digests = sorted(_class_digest_set(records))
        aggregate = hashlib.sha256(("\n".join(digests) + "\n").encode("ascii")).hexdigest()
        return {
            "number_of_unique_sequences": len(digests),
            "aggregate_sequence_sha256": aggregate,
            "members": sorted(
                (
                    {
                        "accession": record["accession"],
                        "class_label": class_label,
                        "sequence_sha256": _record_digest(record),
                    }
                    for record in records
                ),
                key=lambda member: (
                    member["sequence_sha256"],
                    member["accession"],
                ),
            ),
        }

    return {
        "schema_version": "1.0",
        "hash_algorithm": "SHA-256",
        "positive": class_record(positives, "positive"),
        "hard_negative": class_record(hard_negatives, "hard_negative"),
        "easy_negative": class_record(easy_negatives, "easy_negative"),
        "cross_class_exact_sequence_overlap_count": 0,
    }


def labelled_set_for_family(
    pools: ReferencePools,
    family: str,
    *,
    number_of_positives: int = 500,
    number_of_hard_negatives: int = 1000,
    number_of_easy_negatives: int = 500,
    seed: int = 20260923,
) -> tuple[list[ReferenceRecord], list[ReferenceRecord], list[ReferenceRecord]]:
    """Return mutually exclusive, sequence-deduplicated reference classes."""

    normalized_family = _validate_base_family(family)
    if normalized_family not in pools.positives:
        raise KeyError(normalized_family)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("seed must be an integer")

    positives = _deduplicate_by_sequence(pools.positives[normalized_family])

    def is_clean_negative(record: ReferenceRecord) -> bool:
        tags = set(record["tags"])
        return not has_family(tags, normalized_family) and not (tags & AMBIGUOUS_TAGS)

    hard_negatives = _deduplicate_by_sequence(
        [record for record in pools.hard_negative_pool if is_clean_negative(record)]
    )
    easy_negatives = _deduplicate_by_sequence(
        [record for record in pools.easy_negative_pool if is_clean_negative(record)]
    )
    selected_positives = _sample_records(
        positives,
        number_of_positives,
        seed_text=f"{seed}:labelled:{normalized_family}:positive",
    )
    selected_hard = _sample_records(
        hard_negatives,
        number_of_hard_negatives,
        seed_text=f"{seed}:labelled:{normalized_family}:hard_negative",
    )
    selected_easy = _sample_records(
        easy_negatives,
        number_of_easy_negatives,
        seed_text=f"{seed}:labelled:{normalized_family}:easy_negative",
    )
    assert_disjoint_reference_classes(selected_positives, selected_hard, selected_easy)
    return selected_positives, selected_hard, selected_easy


def fetch_uniprot_reference_table(
    query: str,
    *,
    maximum_records: int = 500,
    delay_seconds: float = 1.0,
    timeout_seconds: float = 60.0,
) -> list[ReferenceRecord]:
    """Fetch a UniProt table using the same ``tags`` schema as CAZy records.

    Live queries are not immutable. A scientific run should save its query,
    retrieval time, accession list and response checksum outside this code
    repository.
    """

    normalized_query = str(query).strip()
    if not normalized_query:
        raise ValueError("query must not be empty")
    _validated_positive_integer(maximum_records, "maximum_records")
    if maximum_records > 500:
        raise ValueError("maximum_records must not exceed the UniProt page limit of 500")
    for value, name, allow_zero in (
        (delay_seconds, "delay_seconds", True),
        (timeout_seconds, "timeout_seconds", False),
    ):
        if not math.isfinite(value) or value < 0 or (not allow_zero and value == 0):
            qualifier = "non-negative" if allow_zero else "greater than zero"
            raise ValueError(f"{name} must be finite and {qualifier}")

    time.sleep(delay_seconds)
    response = requests.get(
        "https://rest.uniprot.org/uniprotkb/search",
        params={
            "query": normalized_query,
            "format": "tsv",
            "fields": "accession,xref_cazy,sequence",
            "size": maximum_records,
        },
        headers={"User-Agent": "gh13-cbm-hmm-workflow/0.1"},
        timeout=timeout_seconds,
    )
    response.raise_for_status()
    reader = csv.DictReader(io.StringIO(response.text), delimiter="\t")
    expected_columns = {"Entry", "CAZy", "Sequence"}
    if reader.fieldnames is None or not expected_columns.issubset(reader.fieldnames):
        raise ValueError(
            "UniProt response is missing one or more required columns: "
            + ", ".join(sorted(expected_columns))
        )

    rows: list[ReferenceRecord] = []
    seen_accessions: set[str] = set()
    for row_number, row in enumerate(reader, start=2):
        accession = (row.get("Entry") or "").strip()
        raw_sequence = row.get("Sequence") or ""
        if not accession or not raw_sequence:
            raise ValueError(f"UniProt response row {row_number} is incomplete")
        if accession in seen_accessions:
            raise ValueError(f"UniProt response contains duplicate accession {accession!r}")
        seen_accessions.add(accession)
        sequence = normalise_protein_sequence(raw_sequence, identifier=accession)
        raw_tags = [tag.strip() for tag in (row.get("CAZy") or "").split(";") if tag.strip()]
        invalid_tags = [tag for tag in raw_tags if not _TAG_PATTERN.fullmatch(tag)]
        if invalid_tags:
            raise ValueError(
                f"UniProt response row {row_number} contains unrecognized CAZy tags: "
                + ", ".join(invalid_tags)
            )
        rows.append(
            {
                "accession": accession,
                "tags": sorted(set(raw_tags)),
                "sequence": sequence,
                "sequence_sha256": _sequence_digest(sequence),
            }
        )
    return rows
