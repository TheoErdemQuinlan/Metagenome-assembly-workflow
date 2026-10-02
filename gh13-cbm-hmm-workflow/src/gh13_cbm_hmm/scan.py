"""Profile-HMM searches with explicit statistical and coordinate semantics.

The implementation uses :func:`pyhmmer.hmmer.hmmsearch`. Target and domain
reporting thresholds are separate parameters. HMMER's effective search space
(``Z``) is mandatory because an implicit, batch-dependent value makes
otherwise identical E-value decisions non-reproducible. Coordinates in the
returned records are one-based and end-inclusive, matching HMMER output.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pyhmmer

_ALLOWED_PROTEIN_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWYBJZOUX")


@dataclass(frozen=True, slots=True)
class FastaSummary:
    """Sequence-free summary of a validated protein FASTA file."""

    number_of_protein_records: int
    total_amino_acid_residues: int
    minimum_length_amino_acids: int
    maximum_length_amino_acids: int


@dataclass(frozen=True, slots=True)
class DomainHit:
    """One retained profile-HMM domain alignment.

    All positions are one-based and end-inclusive. Bit scores are measured in
    bits; E-values and coverage fractions are dimensionless.
    """

    profile_name: str
    target_id: str
    target_full_sequence_evalue: float
    independent_domain_evalue: float
    domain_bit_score_bits: float
    hmm_start_match_state: int
    hmm_end_match_state: int
    hmm_length_match_states: int
    hmm_alignment_coverage_fraction: float
    target_start_amino_acid_position: int
    target_end_amino_acid_position: int
    target_length_amino_acids: int

    def as_dict(self) -> dict[str, Any]:
        """Return a field-name-preserving representation of the hit."""

        return asdict(self)


def _decode(value: bytes | str | None) -> str:
    if value is None:
        return ""
    return value.decode("utf-8", errors="strict") if isinstance(value, bytes) else value


def profile_name(hmm: Any) -> str:
    """Return a validated profile name without a trailing ``.hmm``."""

    name = _decode(hmm.name).strip().removesuffix(".hmm")
    if not name:
        raise ValueError("Every HMM profile must have a non-empty name")
    return name


def load_hmms(
    hmm_path: str | Path,
    selected_profiles: Iterable[str] | None = None,
) -> list[Any]:
    """Load HMMER profiles, optionally retaining an explicit name subset.

    Duplicate model names and missing requested names are errors. Either
    condition could otherwise make a nominally complete search ambiguous.
    """

    path = Path(hmm_path)
    if not path.is_file():
        raise FileNotFoundError(path)

    requested: set[str] | None = None
    if selected_profiles is not None:
        normalized = [str(name).strip().removesuffix(".hmm") for name in selected_profiles]
        if any(not name for name in normalized):
            raise ValueError("Requested HMM profile names must not be empty")
        if len(normalized) != len(set(normalized)):
            raise ValueError("Requested HMM profile names must be unique")
        requested = set(normalized)
        if not requested:
            raise ValueError("selected_profiles must contain at least one profile name")

    retained: list[Any] = []
    names_in_file: set[str] = set()
    with pyhmmer.plan7.HMMFile(str(path)) as handle:
        for hmm in handle:
            name = profile_name(hmm)
            if name in names_in_file:
                raise ValueError(f"Duplicate HMM profile name in {path}: {name}")
            names_in_file.add(name)
            if requested is None or name in requested:
                retained.append(hmm)

    if requested is not None:
        missing = sorted(requested - names_in_file)
        if missing:
            raise ValueError(f"Requested HMM profiles were not found: {', '.join(missing)}")
    if not retained:
        raise ValueError(f"No HMM profiles were loaded from {path}")
    return retained


def _validate_thresholds(
    target_reporting_evalue_threshold: float,
    domain_reporting_evalue_threshold: float,
    independent_domain_evalue_threshold: float,
    minimum_hmm_coverage: float,
    cpus: int,
    evalue_search_space: float,
) -> None:
    for value, name in (
        (target_reporting_evalue_threshold, "target_reporting_evalue_threshold"),
        (domain_reporting_evalue_threshold, "domain_reporting_evalue_threshold"),
        (independent_domain_evalue_threshold, "independent_domain_evalue_threshold"),
        (evalue_search_space, "evalue_search_space"),
    ):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and greater than zero")
    if not math.isfinite(minimum_hmm_coverage) or not 0 <= minimum_hmm_coverage <= 1:
        raise ValueError("minimum_hmm_coverage must be finite and between zero and one")
    if isinstance(cpus, bool) or not isinstance(cpus, int) or cpus < 0:
        raise ValueError("cpus must be zero or a positive integer")


def _pipeline_options(
    *,
    target_reporting_evalue_threshold: float,
    domain_reporting_evalue_threshold: float,
    evalue_search_space: float,
) -> dict[str, float]:
    return {
        "E": float(target_reporting_evalue_threshold),
        "domE": float(domain_reporting_evalue_threshold),
        "Z": float(evalue_search_space),
    }


def _collect_domain_hits(
    search_results: Iterable[Any],
    *,
    independent_domain_evalue_threshold: float,
    minimum_hmm_coverage: float,
) -> list[DomainHit]:
    retained: list[DomainHit] = []
    for profile_hits in search_results:
        query_name = _decode(profile_hits.query.name).strip().removesuffix(".hmm")
        if not query_name:
            raise ValueError("A search result has an empty HMM profile name")
        for target_hit in profile_hits:
            target_id = _decode(target_hit.name).strip()
            if not target_id:
                raise ValueError(
                    f"Profile {query_name} produced a hit with an empty target identifier"
                )
            target_evalue = float(target_hit.evalue)
            if not math.isfinite(target_evalue) or target_evalue < 0:
                raise ValueError(f"Invalid target E-value for {query_name} against {target_id}")
            target_length = int(target_hit.length)
            if target_length <= 0:
                raise ValueError(f"Invalid target length for {query_name} against {target_id}")

            for domain in target_hit.domains:
                independent_evalue = float(domain.i_evalue)
                score_bits = float(domain.score)
                if (
                    not math.isfinite(independent_evalue)
                    or independent_evalue < 0
                    or not math.isfinite(score_bits)
                ):
                    raise ValueError(
                        f"Invalid domain statistic for {query_name} against {target_id}"
                    )
                if independent_evalue > independent_domain_evalue_threshold:
                    continue

                alignment = domain.alignment
                hmm_start = int(alignment.hmm_from)
                hmm_end = int(alignment.hmm_to)
                hmm_length = int(alignment.hmm_length)
                target_start = int(alignment.target_from)
                target_end = int(alignment.target_to)
                alignment_target_length = int(alignment.target_length)
                if hmm_length <= 0 or not 1 <= hmm_start <= hmm_end <= hmm_length:
                    raise ValueError(
                        f"Invalid HMM coordinates for {query_name} against {target_id}"
                    )
                if (
                    alignment_target_length != target_length
                    or not 1 <= target_start <= target_end <= target_length
                ):
                    raise ValueError(
                        f"Invalid target coordinates for {query_name} against {target_id}"
                    )

                aligned_match_states = hmm_end - hmm_start + 1
                coverage = aligned_match_states / hmm_length
                if coverage < minimum_hmm_coverage:
                    continue
                retained.append(
                    DomainHit(
                        profile_name=query_name,
                        target_id=target_id,
                        target_full_sequence_evalue=target_evalue,
                        independent_domain_evalue=independent_evalue,
                        domain_bit_score_bits=score_bits,
                        hmm_start_match_state=hmm_start,
                        hmm_end_match_state=hmm_end,
                        hmm_length_match_states=hmm_length,
                        hmm_alignment_coverage_fraction=coverage,
                        target_start_amino_acid_position=target_start,
                        target_end_amino_acid_position=target_end,
                        target_length_amino_acids=target_length,
                    )
                )
    return sorted(
        retained,
        key=lambda hit: (
            hit.profile_name,
            hit.target_id,
            hit.target_start_amino_acid_position,
            hit.target_end_amino_acid_position,
            -hit.domain_bit_score_bits,
        ),
    )


def normalise_protein_sequence(sequence: str, *, identifier: str = "protein") -> str:
    """Normalize and validate one unaligned amino-acid sequence.

    ASCII whitespace is removed and one terminal stop marker is accepted and
    removed. Internal stops, gaps and non-IUPAC protein symbols are rejected.
    """

    if not isinstance(sequence, str):
        raise TypeError(f"Protein sequence for {identifier!r} must be text")
    cleaned = "".join(sequence.split()).upper()
    if cleaned.endswith("*"):
        cleaned = cleaned[:-1]
    if not cleaned:
        raise ValueError(f"Protein sequence is empty for identifier {identifier!r}")
    invalid = sorted(set(cleaned) - _ALLOWED_PROTEIN_RESIDUES)
    if invalid:
        rendered = ", ".join(repr(symbol) for symbol in invalid)
        raise ValueError(
            f"Protein sequence for {identifier!r} contains invalid symbol(s): {rendered}"
        )
    return cleaned


def _normalise_record(record_id: str, sequence: str) -> tuple[str, str]:
    if not isinstance(record_id, str):
        raise TypeError("Protein identifiers must be text")
    identifier = record_id.strip()
    if not identifier:
        raise ValueError("Protein identifiers must not be empty")
    if any(character.isspace() for character in identifier):
        raise ValueError(f"Protein identifier contains whitespace: {identifier!r}")
    return identifier, normalise_protein_sequence(sequence, identifier=identifier)


def _iter_text_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header_identifier: str | None = None
    sequence_chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if header_identifier is not None:
                    yield header_identifier, "".join(sequence_chunks)
                header_text = stripped[1:].strip()
                if not header_text:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                header_identifier = header_text.split(maxsplit=1)[0]
                sequence_chunks = []
            else:
                if header_identifier is None:
                    raise ValueError(
                        f"Sequence text precedes the first FASTA header at line {line_number}"
                    )
                sequence_chunks.append(stripped)
    if header_identifier is not None:
        yield header_identifier, "".join(sequence_chunks)


def inspect_protein_fasta(protein_fasta: str | Path) -> FastaSummary:
    """Validate a protein FASTA and reject duplicate identifiers.

    The function is streaming with respect to sequence content. It retains
    identifiers in memory to guarantee uniqueness before a potentially long
    search begins.
    """

    path = Path(protein_fasta)
    if not path.is_file():
        raise FileNotFoundError(path)
    seen_identifiers: set[str] = set()
    number_of_records = 0
    total_residues = 0
    minimum_length: int | None = None
    maximum_length = 0
    for record_id, sequence in _iter_text_fasta(path):
        if "".join(sequence.split()).upper().endswith("*"):
            raise ValueError(
                "Protein FASTA contains a terminal stop marker for identifier "
                f"{record_id!r}; remove terminal stop markers before a streamed search"
            )
        identifier, cleaned = _normalise_record(record_id, sequence)
        if identifier in seen_identifiers:
            raise ValueError(f"Duplicate protein identifier in FASTA: {identifier}")
        seen_identifiers.add(identifier)
        length = len(cleaned)
        number_of_records += 1
        total_residues += length
        minimum_length = length if minimum_length is None else min(minimum_length, length)
        maximum_length = max(maximum_length, length)
    if number_of_records == 0 or minimum_length is None:
        raise ValueError(f"Protein FASTA contains no records: {path}")
    return FastaSummary(
        number_of_protein_records=number_of_records,
        total_amino_acid_residues=total_residues,
        minimum_length_amino_acids=minimum_length,
        maximum_length_amino_acids=maximum_length,
    )


def _digital_block(records: Iterable[tuple[str, str]]) -> pyhmmer.easel.DigitalSequenceBlock:
    alphabet = pyhmmer.easel.Alphabet.amino()
    seen: set[str] = set()
    digital_sequences = []
    for record_id, sequence in records:
        identifier, cleaned = _normalise_record(record_id, sequence)
        if identifier in seen:
            raise ValueError(f"Duplicate protein identifier: {identifier}")
        seen.add(identifier)
        digital_sequences.append(
            pyhmmer.easel.TextSequence(
                name=identifier.encode("utf-8"),
                sequence=cleaned,
            ).digitize(alphabet)
        )
    if not digital_sequences:
        raise ValueError("At least one protein record is required")
    return pyhmmer.easel.DigitalSequenceBlock(alphabet, digital_sequences)


def scan_records(
    protein_records: Iterable[tuple[str, str]],
    hmms: Sequence[Any],
    *,
    evalue_search_space: float,
    target_reporting_evalue_threshold: float = 10.0,
    domain_reporting_evalue_threshold: float = 10.0,
    independent_domain_evalue_threshold: float = 1e-15,
    minimum_hmm_coverage: float = 0.35,
    cpus: int = 0,
) -> list[DomainHit]:
    """Search an in-memory protein collection.

    ``evalue_search_space`` is HMMER's ``Z`` value and must be supplied
    explicitly. For dbCAN runs intended to reproduce hmmscan-style
    statistics, use the number of profiles in the complete HMM library, even
    when searching a subset. Raw bit scores do not depend on ``Z``.
    """

    _validate_thresholds(
        target_reporting_evalue_threshold,
        domain_reporting_evalue_threshold,
        independent_domain_evalue_threshold,
        minimum_hmm_coverage,
        cpus,
        evalue_search_space,
    )
    if not hmms:
        raise ValueError("At least one HMM profile is required")
    block = _digital_block(protein_records)
    results = pyhmmer.hmmer.hmmsearch(
        hmms,
        block,
        cpus=cpus,
        **_pipeline_options(
            target_reporting_evalue_threshold=target_reporting_evalue_threshold,
            domain_reporting_evalue_threshold=domain_reporting_evalue_threshold,
            evalue_search_space=evalue_search_space,
        ),
    )
    return _collect_domain_hits(
        results,
        independent_domain_evalue_threshold=independent_domain_evalue_threshold,
        minimum_hmm_coverage=minimum_hmm_coverage,
    )


def scan_fasta(
    protein_fasta: str | Path,
    hmms: Sequence[Any],
    *,
    evalue_search_space: float,
    target_reporting_evalue_threshold: float = 10.0,
    domain_reporting_evalue_threshold: float = 10.0,
    independent_domain_evalue_threshold: float = 1e-15,
    minimum_hmm_coverage: float = 0.35,
    cpus: int = 0,
    preflight_validation: bool = True,
) -> list[DomainHit]:
    """Search a protein FASTA without loading all sequences into Python.

    By default the FASTA is first validated and checked for duplicate target
    identifiers. PyHMMER then reopens and streams the file for the search.
    Disabling preflight validation saves one read pass but delegates format
    validation to PyHMMER and cannot guarantee identifier uniqueness.
    """

    _validate_thresholds(
        target_reporting_evalue_threshold,
        domain_reporting_evalue_threshold,
        independent_domain_evalue_threshold,
        minimum_hmm_coverage,
        cpus,
        evalue_search_space,
    )
    if not hmms:
        raise ValueError("At least one HMM profile is required")
    path = Path(protein_fasta)
    if not path.is_file():
        raise FileNotFoundError(path)
    if preflight_validation:
        inspect_protein_fasta(path)

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(
        str(path),
        digital=True,
        alphabet=alphabet,
    ) as sequences:
        results = pyhmmer.hmmer.hmmsearch(
            hmms,
            sequences,
            cpus=cpus,
            **_pipeline_options(
                target_reporting_evalue_threshold=target_reporting_evalue_threshold,
                domain_reporting_evalue_threshold=domain_reporting_evalue_threshold,
                evalue_search_space=evalue_search_space,
            ),
        )
        return _collect_domain_hits(
            results,
            independent_domain_evalue_threshold=independent_domain_evalue_threshold,
            minimum_hmm_coverage=minimum_hmm_coverage,
        )


def iter_hit_rows(hits: Iterable[DomainHit]) -> Iterator[dict[str, Any]]:
    """Yield tabular rows with HMM coverage expressed in percent."""

    for hit in hits:
        row = hit.as_dict()
        fraction = row.pop("hmm_alignment_coverage_fraction")
        row["hmm_alignment_coverage_percent"] = fraction * 100.0
        yield row
