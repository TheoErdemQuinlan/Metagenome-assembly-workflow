#!/usr/bin/env python3
"""Fail when a proposed public repository contains data or sensitive material.

This check is deliberately conservative.  It inspects the working tree rather
than only Git-tracked files so that an accidental data copy is detected before
it can be staged.  Cache, build and version-control directories are excluded.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

MAXIMUM_PUBLIC_FILE_SIZE_BYTES = 5 * 1024 * 1024

EXCLUDED_DIRECTORY_NAMES = {
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    ".tox",
    ".venv",
    "__pycache__",
    "build",
    "dist",
}

PROHIBITED_SUFFIXES = {
    ".2bit",
    ".a2m",
    ".aln",
    ".bam",
    ".bz2",
    ".cram",
    ".db",
    ".dmnd",
    ".fa",
    ".faa",
    ".fast5",
    ".fasta",
    ".fastq",
    ".fna",
    ".fq",
    ".ffn",
    ".gb",
    ".gbk",
    ".gz",
    ".h3f",
    ".h3i",
    ".h3m",
    ".h3p",
    ".hdf5",
    ".hmm",
    ".mmi",
    ".npy",
    ".npz",
    ".pickle",
    ".pkl",
    ".sam",
    ".sra",
    ".sqlite",
    ".sqlite3",
    ".sto",
    ".tar",
    ".vcf",
    ".xz",
    ".zip",
}

PROHIBITED_FILE_NAMES = {
    ".env",
    "credentials.json",
    "id_dsa",
    "id_ed25519",
    "id_rsa",
    "secrets.json",
}

TEXT_SUFFIXES = {
    "",
    ".cfg",
    ".cff",
    ".csv",
    ".ini",
    ".json",
    ".lock",
    ".md",
    ".py",
    ".rst",
    ".toml",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}


@dataclass(frozen=True)
class Violation:
    path: Path
    reason: str


def _content_patterns() -> tuple[tuple[re.Pattern[str], str], ...]:
    home_prefixes = (
        re.compile(re.escape("/" + "Users" + "/")),
        re.compile(re.escape("/" + "home" + "/")),
        re.compile(r"[A-Za-z]:\\\\Users\\\\"),
    )
    patterns: list[tuple[re.Pattern[str], str]] = [
        *((pattern, "contains an absolute user-home path") for pattern in home_prefixes),
        (
            re.compile(r"\b" + "Con" + "sett" + r"\b", re.IGNORECASE),
            "contains the study-site name",
        ),
        (
            re.compile(r"\bPS-[A-Z0-9]{8}\b"),
            "contains a project-specific protein identifier",
        ),
        (
            re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{60,}"),
            "contains a long protein-like character string",
        ),
        (
            re.compile(r"[ACGTUN]{80,}"),
            "contains a long nucleotide-like character string",
        ),
        (
            re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
            "contains private-key material",
        ),
        (
            re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
            "contains an AWS access-key identifier",
        ),
        (
            re.compile(r"\bgh[pousr]_[A-Za-z0-9]{30,}\b"),
            "contains a GitHub credential-like value",
        ),
        (
            re.compile(r"\bsk-[A-Za-z0-9]{20,}\b"),
            "contains an API credential-like value",
        ),
        (
            re.compile(
                r"\b(?:postgres(?:ql)?|mysql|mariadb|mongodb)://"
                r"[^\s/:@]+:[^\s/@]+@",
                re.IGNORECASE,
            ),
            "contains a database URL with embedded credentials",
        ),
    ]
    return tuple(patterns)


CONTENT_PATTERNS = _content_patterns()


def iter_repository_files(root: Path) -> Iterable[Path]:
    """Yield regular files while pruning generated and local tool directories."""

    for path in sorted(root.rglob("*")):
        try:
            relative_parts = path.relative_to(root).parts
        except ValueError:
            continue
        if any(part in EXCLUDED_DIRECTORY_NAMES for part in relative_parts):
            continue
        if path.is_symlink() or path.is_file():
            yield path


def audit_tree(
    root: str | Path,
    *,
    maximum_file_size_bytes: int = MAXIMUM_PUBLIC_FILE_SIZE_BYTES,
) -> list[Violation]:
    """Return all policy violations found below ``root``."""

    root_path = Path(root).resolve()
    if not root_path.is_dir():
        raise NotADirectoryError(root_path)
    violations: list[Violation] = []

    for path in iter_repository_files(root_path):
        relative_path = path.relative_to(root_path)
        if path.is_symlink():
            violations.append(Violation(relative_path, "symbolic links are not allowed"))
            continue
        lowered_name = path.name.lower()
        suffix = path.suffix.lower()
        is_environment_example = lowered_name == ".env.example"
        is_environment_file = lowered_name.startswith(".env.") and not is_environment_example
        if lowered_name in PROHIBITED_FILE_NAMES or is_environment_file:
            violations.append(Violation(relative_path, "prohibited secret-bearing file name"))
        if suffix in PROHIBITED_SUFFIXES:
            violations.append(
                Violation(relative_path, "prohibited data, model or archive extension")
            )

        file_size = path.stat().st_size
        if file_size > maximum_file_size_bytes:
            violations.append(
                Violation(
                    relative_path,
                    f"file exceeds {maximum_file_size_bytes} bytes",
                )
            )
            continue

        if suffix not in TEXT_SUFFIXES and not is_environment_example:
            if suffix in PROHIBITED_SUFFIXES:
                continue
            violations.append(
                Violation(
                    relative_path,
                    "file type is not allowed in the code-only public tree",
                )
            )
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            violations.append(Violation(relative_path, "expected text file is not valid UTF-8"))
            continue
        for pattern, reason in CONTENT_PATTERNS:
            if pattern.search(text):
                violations.append(Violation(relative_path, reason))

    return violations


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit a code-only public repository for data and sensitive material."
    )
    parser.add_argument("root", nargs="?", default=".", help="repository root")
    arguments = parser.parse_args(argv)

    violations = audit_tree(arguments.root)
    if violations:
        print("Public-tree audit failed:", file=sys.stderr)
        for violation in violations:
            print(f"  {violation.path}: {violation.reason}", file=sys.stderr)
        return 1
    print("Public-tree audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
