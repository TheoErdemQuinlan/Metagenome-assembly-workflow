"""Command-line interface for the standalone HMM workflow."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict
from pathlib import Path
from typing import Any

from . import __version__
from .provenance import build_scan_manifest, sha256_file, write_json
from .scan import (
    inspect_protein_fasta,
    iter_hit_rows,
    load_hmms,
    profile_name,
    scan_fasta,
)
from .tiering import UNCALIBRATED_FAMILY, classify_score
from .validation import calibration_report

HIT_COLUMNS = [
    "profile_name",
    "target_id",
    "target_full_sequence_evalue",
    "independent_domain_evalue",
    "domain_bit_score_bits",
    "hmm_start_match_state",
    "hmm_end_match_state",
    "hmm_length_match_states",
    "target_start_amino_acid_position",
    "target_end_amino_acid_position",
    "target_length_amino_acids",
    "hmm_alignment_coverage_percent",
]

ACCEPTED_NEGATIVE_LABELS = {"negative", "hard_negative", "easy_negative"}


class CliInputError(ValueError):
    """A concise, user-facing input validation error."""


def _positive_float(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a number greater than zero") from error
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("must be a finite number greater than zero")
    return parsed


def _percentage(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a percentage from 0 to 100") from error
    if not math.isfinite(parsed) or not 0 <= parsed <= 100:
        raise argparse.ArgumentTypeError("must be a finite percentage from 0 to 100")
    return parsed


def _nonnegative_float(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be zero or a positive number") from error
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("must be a finite, non-negative number")
    return parsed


def _nonnegative_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be zero or a positive integer") from error
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be zero or a positive integer")
    return parsed


def _flatten_profiles(profile_groups: list[list[str]] | None) -> list[str] | None:
    if not profile_groups:
        return None
    names: list[str] = []
    for group in profile_groups:
        for value in group:
            names.extend(item.strip() for item in value.split(",") if item.strip())
    return list(dict.fromkeys(names)) or None


def _atomic_text_writer(path: Path):
    """Return a temporary path next to ``path`` for an atomic replacement."""

    path.parent.mkdir(parents=True, exist_ok=True)
    return path.with_name(f".{path.name}.{os.getpid()}.tmp")


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    temporary = _atomic_text_writer(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(columns),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    temporary = _atomic_text_writer(path)
    try:
        write_json(temporary, payload)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise CliInputError(f"Input table does not exist: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CliInputError(f"Input table has no header: {path}")
        return list(reader.fieldnames), list(reader)


def _require_columns(fieldnames: Sequence[str], required: set[str], path: Path) -> None:
    missing = sorted(required - set(fieldnames))
    if missing:
        raise CliInputError(f"Input table {path} is missing required columns: {', '.join(missing)}")


def _require_output_separate(
    output_path: Path,
    input_paths: Iterable[Path],
    *,
    label: str,
) -> None:
    resolved_output = output_path.expanduser().resolve()
    for input_path in input_paths:
        if resolved_output == input_path.expanduser().resolve():
            raise CliInputError(f"{label} must not overwrite input file {input_path}")


def _run_scan(arguments: argparse.Namespace) -> int:
    output_path = Path(arguments.output)
    manifest_path = (
        Path(arguments.manifest)
        if arguments.manifest
        else output_path.with_suffix(".manifest.json")
    )
    _require_output_separate(
        output_path,
        [Path(arguments.proteins), Path(arguments.hmms)],
        label="Domain-hit output",
    )
    _require_output_separate(
        manifest_path,
        [Path(arguments.proteins), Path(arguments.hmms), output_path],
        label="Manifest output",
    )
    selected_profiles = _flatten_profiles(arguments.profile)
    fasta_summary = (
        None if arguments.skip_fasta_preflight else inspect_protein_fasta(arguments.proteins)
    )
    hmms = load_hmms(arguments.hmms, selected_profiles)
    selected_names = [profile_name(hmm) for hmm in hmms]
    minimum_coverage_fraction = arguments.minimum_hmm_coverage_percent / 100.0

    hits = scan_fasta(
        arguments.proteins,
        hmms,
        target_reporting_evalue_threshold=arguments.target_reporting_evalue_threshold,
        domain_reporting_evalue_threshold=arguments.domain_reporting_evalue_threshold,
        independent_domain_evalue_threshold=arguments.independent_domain_evalue_threshold,
        minimum_hmm_coverage=minimum_coverage_fraction,
        cpus=arguments.cpus,
        evalue_search_space=arguments.evalue_search_space,
        preflight_validation=False,
    )
    hits.sort(
        key=lambda hit: (
            hit.target_id,
            hit.target_start_amino_acid_position,
            hit.profile_name,
            hit.independent_domain_evalue,
        )
    )
    _write_tsv(output_path, HIT_COLUMNS, iter_hit_rows(hits))

    parameters = {
        "target_reporting_evalue_threshold": arguments.target_reporting_evalue_threshold,
        "domain_reporting_evalue_threshold": arguments.domain_reporting_evalue_threshold,
        "independent_domain_evalue_threshold": (arguments.independent_domain_evalue_threshold),
        "minimum_hmm_alignment_coverage_percent": (arguments.minimum_hmm_coverage_percent),
        "cpus": arguments.cpus,
        "evalue_search_space_Z": arguments.evalue_search_space,
        "protein_fasta_preflight_validation": not arguments.skip_fasta_preflight,
        "decision_rule": (
            "independent domain E-value <= independent-domain threshold and "
            "HMM alignment coverage >= coverage threshold"
        ),
    }
    manifest = build_scan_manifest(
        protein_fasta=arguments.proteins,
        hmm_file=arguments.hmms,
        profile_names=selected_names,
        parameters=parameters,
        output_file=output_path,
    )
    manifest["inputs"]["protein_fasta"]["sequence_free_summary"] = (
        asdict(fasta_summary) if fasta_summary is not None else None
    )
    manifest["output"].update(
        {
            "number_of_retained_domain_hits": len(hits),
            "sha256": sha256_file(output_path),
        }
    )
    _write_json_atomic(manifest_path, manifest)
    print(
        f"Wrote {len(hits)} retained domain hits to {output_path} "
        f"and the run manifest to {manifest_path}."
    )
    return 0


def _score_rows(
    path: Path,
    family: str,
) -> tuple[list[float], list[float], dict[str, int]]:
    fieldnames, rows = _read_tsv(path)
    _require_columns(
        fieldnames,
        {"profile_name", "class_label", "domain_bit_score_bits"},
        path,
    )
    positives: list[float] = []
    negatives: list[float] = []
    class_counts = {
        "positive": 0,
        "hard_negative": 0,
        "easy_negative": 0,
        "negative": 0,
    }
    for line_number, row in enumerate(rows, start=2):
        if row["profile_name"].strip() != family:
            continue
        label = row["class_label"].strip().casefold()
        if label != "positive" and label not in ACCEPTED_NEGATIVE_LABELS:
            raise CliInputError(
                f"Unrecognized class_label {row['class_label']!r} at {path}:{line_number}"
            )
        try:
            score = float(row["domain_bit_score_bits"])
        except ValueError as error:
            raise CliInputError(f"Invalid domain_bit_score_bits at {path}:{line_number}") from error
        if not math.isfinite(score):
            raise CliInputError(f"Non-finite domain_bit_score_bits at {path}:{line_number}")
        class_counts[label] += 1
        if label == "positive":
            positives.append(score)
        else:
            negatives.append(score)
    if not positives:
        raise CliInputError(f"No positive scores were found for profile {family!r}")
    if not negatives:
        raise CliInputError(f"No negative scores were found for profile {family!r}")
    return positives, negatives, class_counts


def _run_calibrate_scores(arguments: argparse.Namespace) -> int:
    input_path = Path(arguments.scores)
    output_path = Path(arguments.output)
    _require_output_separate(
        output_path,
        [input_path],
        label="Calibration output",
    )
    positives, negatives, class_counts = _score_rows(input_path, arguments.family)
    report = calibration_report(
        positives,
        negatives,
        comparison_threshold_bits=arguments.comparison_threshold_bits,
        robustness_delta_bits=arguments.robustness_delta_bits,
    )
    cutoff = float(report["selected_cutoff"]["bit_score_threshold_bits"])
    payload = {
        "schema_version": "1.0",
        "software": {
            "package": "gh13-cbm-hmm-workflow",
            "package_version": __version__,
        },
        "family": arguments.family,
        "score_units": "bits",
        "classification_rule": "domain_bit_score_bits >= bit_score_threshold_bits",
        "input": {
            "file_name": input_path.name,
            "sha256": sha256_file(input_path),
            "class_counts": class_counts,
        },
        "cutoffs": {
            arguments.family: {
                "bit_score_threshold_bits": cutoff,
                "scope": "internal operational calibration",
            }
        },
        "calibration": report,
    }
    _write_json_atomic(output_path, payload)
    print(f"Wrote the internal calibration report for {arguments.family} to {output_path}.")
    return 0


def _extract_cutoffs(payload: Any, source: Path) -> dict[str, float]:
    if not isinstance(payload, dict) or not isinstance(payload.get("cutoffs"), dict):
        raise CliInputError(f"Cutoff file {source} must contain a JSON object named 'cutoffs'")
    cutoffs: dict[str, float] = {}
    for family, record in payload["cutoffs"].items():
        raw_value = record.get("bit_score_threshold_bits") if isinstance(record, dict) else record
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as error:
            raise CliInputError(f"Cutoff for {family!r} in {source} is not a number") from error
        if not math.isfinite(value):
            raise CliInputError(f"Cutoff for {family!r} in {source} is not finite")
        cutoffs[str(family)] = value
    if not cutoffs:
        raise CliInputError(f"Cutoff file {source} contains no family cutoffs")
    return cutoffs


def _load_cutoffs(paths: Sequence[Path]) -> dict[str, float]:
    combined: dict[str, float] = {}
    for path in paths:
        if not path.is_file():
            raise CliInputError(f"Cutoff file does not exist: {path}")
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            raise CliInputError(f"Cutoff file is not valid JSON: {path}") from error
        for family, value in _extract_cutoffs(payload, path).items():
            if family in combined and combined[family] != value:
                raise CliInputError(f"Conflicting cutoffs were supplied for family {family!r}")
            combined[family] = value
    return combined


def _run_apply_cutoff(arguments: argparse.Namespace) -> int:
    input_path = Path(arguments.hits)
    cutoff_paths = [Path(value) for value in arguments.cutoffs]
    output_path = Path(arguments.output)
    _require_output_separate(
        output_path,
        [input_path, *cutoff_paths],
        label="Classified output",
    )
    fieldnames, rows = _read_tsv(input_path)
    _require_columns(fieldnames, {"profile_name", "domain_bit_score_bits"}, input_path)
    cutoffs = _load_cutoffs(cutoff_paths)

    output_rows: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    for line_number, row in enumerate(rows, start=2):
        family = row["profile_name"].strip()
        try:
            score = float(row["domain_bit_score_bits"])
        except ValueError as error:
            raise CliInputError(
                f"Invalid domain_bit_score_bits at {input_path}:{line_number}"
            ) from error
        if not math.isfinite(score):
            raise CliInputError(f"Non-finite domain_bit_score_bits at {input_path}:{line_number}")
        classification = classify_score(family, score, cutoffs)
        output_row: dict[str, Any] = dict(row)
        output_row["operational_cutoff_bit_score_bits"] = cutoffs.get(family, "")
        output_row["operational_cutoff_classification"] = classification
        output_rows.append(output_row)
        counts[classification] = counts.get(classification, 0) + 1

    if arguments.require_calibrated and counts.get(UNCALIBRATED_FAMILY, 0):
        raise CliInputError(
            f"{counts[UNCALIBRATED_FAMILY]} hit rows use profiles without supplied cutoffs"
        )
    output_columns = [
        *fieldnames,
        "operational_cutoff_bit_score_bits",
        "operational_cutoff_classification",
    ]
    _write_tsv(output_path, output_columns, output_rows)
    print(f"Wrote {len(output_rows)} classified domain-hit rows to {output_path}.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gh13-cbm-hmm",
        description=(
            "Profile-HMM screening and internally calibrated operational "
            "cutoffs for GH13 and associated CBM families."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    scan = subparsers.add_parser(
        "scan",
        help="Search a protein FASTA with one or more profile HMMs.",
    )
    scan.add_argument("--proteins", type=Path, required=True, help="Input protein FASTA.")
    scan.add_argument("--hmms", type=Path, required=True, help="Input HMMER profile file.")
    scan.add_argument(
        "--profile",
        nargs="+",
        action="append",
        help="Profile name(s) to retain; repeat or provide comma-separated names.",
    )
    scan.add_argument("--output", type=Path, required=True, help="Output domain-hit TSV.")
    scan.add_argument(
        "--manifest",
        type=Path,
        help="Output manifest JSON; defaults to <output stem>.manifest.json.",
    )
    scan.add_argument(
        "--target-reporting-evalue-threshold",
        type=_positive_float,
        default=10.0,
        help="HMMER per-target reporting E-value threshold (default: 10).",
    )
    scan.add_argument(
        "--domain-reporting-evalue-threshold",
        type=_positive_float,
        default=10.0,
        help="HMMER conditional-domain reporting E-value threshold (default: 10).",
    )
    scan.add_argument(
        "--independent-domain-evalue-threshold",
        type=_positive_float,
        default=1e-15,
        help="Maximum independent-domain E-value (default: 1e-15).",
    )
    scan.add_argument(
        "--minimum-hmm-coverage-percent",
        type=_percentage,
        default=35.0,
        help="Minimum aligned percentage of HMM match states (default: 35).",
    )
    scan.add_argument(
        "--cpus",
        type=_nonnegative_integer,
        default=0,
        help="Worker CPU count; zero uses PyHMMER's default (default: 0).",
    )
    scan.add_argument(
        "--evalue-search-space",
        type=_positive_float,
        required=True,
        metavar="Z",
        help=(
            "Explicit HMMER Z value. For dbCAN hmmscan-equivalent statistics, "
            "use the profile count in the complete HMM library."
        ),
    )
    scan.add_argument(
        "--skip-fasta-preflight",
        action="store_true",
        help=(
            "Skip the preliminary format, residue and duplicate-identifier check. "
            "The default preflight provides stronger input validation."
        ),
    )
    scan.set_defaults(handler=_run_scan)

    calibrate = subparsers.add_parser(
        "calibrate-score",
        aliases=["calibrate-scores"],
        help="Calibrate one operational bit-score cutoff from labelled scores.",
    )
    calibrate.add_argument("--scores", type=Path, required=True, help="Labelled score TSV.")
    calibrate.add_argument("--family", required=True, help="Exact profile name to calibrate.")
    calibrate.add_argument("--output", type=Path, required=True, help="Output report JSON.")
    calibrate.add_argument(
        "--comparison-threshold-bits",
        type=float,
        help="Optional pre-existing bit-score threshold for comparison, in bits.",
    )
    calibrate.add_argument(
        "--robustness-delta-bits",
        type=_nonnegative_float,
        default=5.0,
        help="Offset around the selected cutoff for the robustness table, in bits (default: 5).",
    )
    calibrate.set_defaults(handler=_run_calibrate_scores)

    apply_cutoff = subparsers.add_parser(
        "apply-cutoff",
        aliases=["apply-cutoffs"],
        help="Apply one or more operational bit-score cutoffs to domain-hit rows.",
    )
    apply_cutoff.add_argument("--hits", type=Path, required=True, help="Input domain-hit TSV.")
    apply_cutoff.add_argument(
        "--cutoffs",
        type=Path,
        required=True,
        action="append",
        help="Calibration JSON; repeat to combine non-conflicting family cutoffs.",
    )
    apply_cutoff.add_argument("--output", type=Path, required=True, help="Output classified TSV.")
    apply_cutoff.add_argument(
        "--require-calibrated",
        action="store_true",
        help="Fail if any hit belongs to a family without a supplied cutoff.",
    )
    apply_cutoff.set_defaults(handler=_run_apply_cutoff)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        return int(arguments.handler(arguments))
    except (CliInputError, FileNotFoundError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    sys.exit(main())
