# Contributing

Contributions should preserve the scientific scope and the code-only nature of this
repository.

## Before submitting a change

1. Create a branch from the current default branch.
2. Install the development environment with `python -m pip install -e '.[dev]'`.
3. Run `ruff check .`, `ruff format --check .` and `pytest`.
4. Describe any change to a statistical definition, threshold rule, coordinate convention
   or output schema in both the pull request and `CHANGELOG.md`.

Tests should use small abstract fixtures. Do not commit real or synthetic amino-acid or
nucleotide sequences, metagenomic identifiers, sample names, accession lists, database
dumps or derived candidate tables. Integration tests that require biological data should
read paths from environment variables and skip cleanly when those variables are absent.

Third-party profile HMMs and reference sequences remain under the terms imposed by their
providers. Do not add those files to this repository. Documentation may record a resource
name, release, checksum and retrieval date when required for reproducibility.

Reports of incorrect scientific assumptions are especially useful. A proposed correction
should state which result could change, supply a minimal non-biological regression case and
distinguish code validation from experimental validation.
