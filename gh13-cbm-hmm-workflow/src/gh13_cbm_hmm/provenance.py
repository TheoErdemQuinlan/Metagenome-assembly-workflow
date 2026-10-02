"""Run manifests and content checksums."""

from __future__ import annotations

import hashlib
import json
import math
import platform
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pyhmmer

from . import __version__


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    """Return the lowercase SHA-256 digest of a regular file."""

    if isinstance(block_size, bool) or not isinstance(block_size, int) or block_size <= 0:
        raise ValueError("block_size must be a positive integer")
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(file_path)
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: str | Path) -> dict[str, Any]:
    """Return path-independent size and content identity for one file."""

    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(file_path)
    return {
        "file_name": file_path.name,
        "size_bytes": file_path.stat().st_size,
        "sha256": sha256_file(file_path),
    }


def build_scan_manifest(
    *,
    protein_fasta: str | Path,
    hmm_file: str | Path,
    profile_names: list[str],
    parameters: dict[str, Any],
    output_file: str | Path,
) -> dict[str, Any]:
    """Build a path-independent, sequence-free profile-search manifest.

    The output file is content-addressed when it already exists. HMMER's
    effective search space ``Z`` is required in ``parameters`` so an
    unrecorded batch-dependent default cannot enter a reproducibility record.
    """

    required_parameters = {
        "target_reporting_evalue_threshold",
        "domain_reporting_evalue_threshold",
        "independent_domain_evalue_threshold",
        "evalue_search_space_Z",
    }
    missing_parameters = sorted(required_parameters - parameters.keys())
    if missing_parameters:
        raise ValueError(
            "Scan manifest parameters are missing required values: " + ", ".join(missing_parameters)
        )
    for key in required_parameters:
        value = parameters[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise TypeError(f"Manifest parameter {key!r} must be numeric")
        if not math.isfinite(float(value)) or float(value) <= 0:
            raise ValueError(f"Manifest parameter {key!r} must be finite and positive")

    normalized_profiles = [str(name).strip().removesuffix(".hmm") for name in profile_names]
    if not normalized_profiles or any(not name for name in normalized_profiles):
        raise ValueError("profile_names must contain non-empty profile names")
    if len(normalized_profiles) != len(set(normalized_profiles)):
        raise ValueError("profile_names must be unique")

    output_path = Path(output_file)
    output_record: dict[str, Any] = {"file_name": output_path.name}
    if output_path.exists():
        output_record = file_record(output_path)

    manifest = {
        "schema_version": "1.0",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "software": {
            "package": "gh13-cbm-hmm-workflow",
            "package_version": __version__,
            "python_version": platform.python_version(),
            "pyhmmer_version": pyhmmer.__version__,
        },
        "inputs": {
            "protein_fasta": file_record(protein_fasta),
            "hmm_file": file_record(hmm_file),
        },
        "profiles": {
            "number_searched": len(normalized_profiles),
            "names": normalized_profiles,
        },
        "parameters": dict(parameters),
        "coordinate_convention": ("protein and HMM positions are one-based and end-inclusive"),
        "output": output_record,
    }
    # Reject NaN, infinity and values that the strict JSON format cannot encode.
    json.dumps(manifest, allow_nan=False)
    return manifest


def write_json(
    path: str | Path,
    payload: dict[str, Any],
    *,
    overwrite: bool = False,
) -> None:
    """Write strict UTF-8 JSON without silently replacing an existing file."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if overwrite else "x"
    serialized = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with output_path.open(mode, encoding="utf-8", newline="\n") as handle:
        handle.write(serialized)
