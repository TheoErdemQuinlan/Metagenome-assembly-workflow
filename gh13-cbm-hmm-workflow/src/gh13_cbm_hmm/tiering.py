"""Apply operational score boundaries and resolve overlapping domain hits."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

PASSES_OPERATIONAL_CUTOFF = "passes_operational_cutoff"
BELOW_OPERATIONAL_CUTOFF = "below_operational_cutoff"
UNCALIBRATED_FAMILY = "uncalibrated_family"


def classify_score(
    family: str,
    domain_bit_score_bits: float,
    cutoffs_by_family: Mapping[str, float],
) -> str:
    """Classify one score without implying independent biological validation."""

    normalized_family = str(family).strip()
    if not normalized_family:
        raise ValueError("family must not be empty")
    score = float(domain_bit_score_bits)
    if not math.isfinite(score):
        raise ValueError("domain_bit_score_bits must be finite")
    if normalized_family not in cutoffs_by_family:
        return UNCALIBRATED_FAMILY
    cutoff = float(cutoffs_by_family[normalized_family])
    if not math.isfinite(cutoff):
        raise ValueError(f"Operational cutoff for {normalized_family} must be finite")
    if score >= cutoff:
        return PASSES_OPERATIONAL_CUTOFF
    return BELOW_OPERATIONAL_CUTOFF


def overlap_fraction_of_shorter(first: Mapping[str, Any], second: Mapping[str, Any]) -> float:
    """Return inclusive coordinate overlap divided by the shorter interval.

    Input positions are interpreted as one-based and end-inclusive amino-acid
    coordinates. The same arithmetic also applies to zero-based inclusive
    positions, but mixing coordinate systems within one call is invalid.
    """

    first_start = int(first["target_start_amino_acid_position"])
    first_end = int(first["target_end_amino_acid_position"])
    second_start = int(second["target_start_amino_acid_position"])
    second_end = int(second["target_end_amino_acid_position"])
    if first_start < 1 or second_start < 1:
        raise ValueError("Domain positions must be positive integers")
    if first_start > first_end or second_start > second_end:
        raise ValueError("Domain start positions must not exceed end positions")
    intersection = max(0, min(first_end, second_end) - max(first_start, second_start) + 1)
    shorter = min(first_end - first_start + 1, second_end - second_start + 1)
    return intersection / shorter


def distinct_domains(
    hits: list[Mapping[str, Any]],
    *,
    maximum_overlap_fraction: float = 0.5,
) -> list[Mapping[str, Any]]:
    """Collapse overlapping matches while retaining the strongest evidence.

    Domains on different targets are never compared. Within a target, lower
    independent-domain E-value wins; ties are resolved by higher bit score,
    greater HMM coverage and then coordinates. ``maximum_overlap_fraction``
    is the largest overlap fraction that two retained domains may share.
    """

    if not math.isfinite(maximum_overlap_fraction) or not 0 <= maximum_overlap_fraction <= 1:
        raise ValueError("maximum_overlap_fraction must be between zero and one")

    def target_key(hit: Mapping[str, Any]) -> str:
        value = str(hit.get("target_id", ""))
        if "target_id" in hit and not value.strip():
            raise ValueError("target_id must not be empty when provided")
        return value

    def evidence_key(hit: Mapping[str, Any]) -> tuple[float, float, float, str, int, int]:
        evalue = float(hit["independent_domain_evalue"])
        score = float(hit.get("domain_bit_score_bits", 0.0))
        coverage = float(hit.get("hmm_alignment_coverage_fraction", 0.0))
        if not math.isfinite(evalue) or evalue < 0:
            raise ValueError("independent_domain_evalue must be finite and non-negative")
        if not math.isfinite(score):
            raise ValueError("domain_bit_score_bits must be finite when provided")
        if not math.isfinite(coverage) or not 0 <= coverage <= 1:
            raise ValueError(
                "hmm_alignment_coverage_fraction must be between zero and one when provided"
            )
        overlap_fraction_of_shorter(hit, hit)
        return (
            evalue,
            -score,
            -coverage,
            target_key(hit),
            int(hit["target_start_amino_acid_position"]),
            int(hit["target_end_amino_acid_position"]),
        )

    accepted: list[Mapping[str, Any]] = []
    ordered = sorted(hits, key=evidence_key)
    for hit in ordered:
        if all(
            target_key(hit) != target_key(retained)
            or overlap_fraction_of_shorter(hit, retained) <= maximum_overlap_fraction
            for retained in accepted
        ):
            accepted.append(hit)
    return sorted(
        accepted,
        key=lambda hit: (
            target_key(hit),
            int(hit["target_start_amino_acid_position"]),
            int(hit["target_end_amino_acid_position"]),
        ),
    )
