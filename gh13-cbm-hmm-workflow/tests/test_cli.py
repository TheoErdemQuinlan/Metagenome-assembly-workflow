from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from gh13_cbm_hmm.cli import main
from gh13_cbm_hmm.tiering import (
    BELOW_OPERATIONAL_CUTOFF,
    PASSES_OPERATIONAL_CUTOFF,
    UNCALIBRATED_FAMILY,
)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_calibrate_score_writes_sequence_free_exact_threshold_report(
    tmp_path: Path,
) -> None:
    scores = tmp_path / "scores.tsv"
    scores.write_text(
        "profile_name\tclass_label\tdomain_bit_score_bits\n"
        "profile_alpha\tpositive\t12\n"
        "profile_alpha\tpositive\t10\n"
        "profile_alpha\thard_negative\t8\n"
        "profile_alpha\teasy_negative\t3\n",
        encoding="utf-8",
    )
    output = tmp_path / "calibration.json"

    exit_code = main(
        [
            "calibrate-score",
            "--scores",
            str(scores),
            "--family",
            "profile_alpha",
            "--output",
            str(output),
            "--comparison-threshold-bits",
            "8.5",
        ]
    )

    report = json.loads(output.read_text(encoding="utf-8"))
    assert exit_code == 0
    assert report["score_units"] == "bits"
    assert report["classification_rule"] == ("domain_bit_score_bits >= bit_score_threshold_bits")
    assert report["cutoffs"]["profile_alpha"]["bit_score_threshold_bits"] == 10.0
    assert report["input"]["class_counts"] == {
        "easy_negative": 1,
        "hard_negative": 1,
        "negative": 0,
        "positive": 2,
    }
    assert report["calibration"]["comparison_threshold"]["incorrectly_accepted_negatives"] == 0


def test_apply_cutoff_preserves_rows_and_adds_operational_classification(
    tmp_path: Path,
) -> None:
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "profile_name\tdomain_bit_score_bits\tnote\n"
        "profile_alpha\t10\tone\n"
        "profile_alpha\t9.9\ttwo\n"
        "profile_beta\t90\tthree\n",
        encoding="utf-8",
    )
    cutoffs = tmp_path / "cutoffs.json"
    cutoffs.write_text(
        json.dumps({"cutoffs": {"profile_alpha": {"bit_score_threshold_bits": 10.0}}}),
        encoding="utf-8",
    )
    output = tmp_path / "classified.tsv"

    exit_code = main(
        [
            "apply-cutoff",
            "--hits",
            str(hits),
            "--cutoffs",
            str(cutoffs),
            "--output",
            str(output),
        ]
    )

    rows = _read_tsv(output)
    assert exit_code == 0
    assert [row["note"] for row in rows] == ["one", "two", "three"]
    assert [row["operational_cutoff_classification"] for row in rows] == [
        PASSES_OPERATIONAL_CUTOFF,
        BELOW_OPERATIONAL_CUTOFF,
        UNCALIBRATED_FAMILY,
    ]
    assert rows[0]["operational_cutoff_bit_score_bits"] == "10.0"
    assert rows[2]["operational_cutoff_bit_score_bits"] == ""


def test_require_calibrated_fails_before_writing_output(tmp_path: Path) -> None:
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "profile_name\tdomain_bit_score_bits\nprofile_beta\t90\n",
        encoding="utf-8",
    )
    cutoffs = tmp_path / "cutoffs.json"
    cutoffs.write_text(
        json.dumps({"cutoffs": {"profile_alpha": 10.0}}),
        encoding="utf-8",
    )
    output = tmp_path / "classified.tsv"

    with pytest.raises(SystemExit) as error:
        main(
            [
                "apply-cutoff",
                "--hits",
                str(hits),
                "--cutoffs",
                str(cutoffs),
                "--output",
                str(output),
                "--require-calibrated",
            ]
        )

    assert error.value.code == 2
    assert not output.exists()
