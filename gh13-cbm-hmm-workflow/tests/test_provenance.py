from __future__ import annotations

import hashlib
from pathlib import Path

from gh13_cbm_hmm.provenance import build_scan_manifest, file_record, sha256_file


def _all_strings(value: object):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from _all_strings(key)
            yield from _all_strings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _all_strings(item)


def test_file_hash_and_record_are_content_addressed(tmp_path: Path) -> None:
    path = tmp_path / "input.bin"
    payload = b"abstract test payload\n"
    path.write_bytes(payload)

    expected = hashlib.sha256(payload).hexdigest()
    assert sha256_file(path) == expected
    assert file_record(path) == {
        "file_name": "input.bin",
        "size_bytes": len(payload),
        "sha256": expected,
    }


def test_manifest_records_checksums_parameters_units_and_relative_names(
    tmp_path: Path,
) -> None:
    proteins = tmp_path / "proteins.input"
    profiles = tmp_path / "profiles.input"
    proteins.write_bytes(b"input collection")
    profiles.write_bytes(b"profile collection")

    manifest = build_scan_manifest(
        protein_fasta=proteins,
        hmm_file=profiles,
        profile_names=["GH13", "CBM20"],
        parameters={
            "target_reporting_evalue_threshold": 10.0,
            "domain_reporting_evalue_threshold": 10.0,
            "independent_domain_evalue_threshold": 1e-15,
            "minimum_hmm_alignment_coverage_percent": 35.0,
            "evalue_search_space_Z": 875.0,
        },
        output_file=tmp_path / "hits.tsv",
    )

    coordinate_convention = manifest["coordinate_convention"]
    assert "one-based" in coordinate_convention
    assert "end-inclusive" in coordinate_convention
    assert manifest["profiles"] == {
        "number_searched": 2,
        "names": ["GH13", "CBM20"],
    }
    assert manifest["parameters"]["minimum_hmm_alignment_coverage_percent"] == 35.0
    assert manifest["inputs"]["protein_fasta"]["file_name"] == "proteins.input"
    assert manifest["output"] == {"file_name": "hits.tsv"}
    assert all(not value.startswith(str(tmp_path)) for value in _all_strings(manifest))
