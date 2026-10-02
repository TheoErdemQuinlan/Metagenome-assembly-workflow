# Changelog

All notable changes to this repository are recorded here. Version numbers follow
[Semantic Versioning](https://semver.org/).

## 0.1.0 — 2026-10-02

Initial audited, code-only release.

- Adapted the profile-HMM workflow into a standalone package independent of the project
  database and web application.
- Added FASTA scanning with separate HMMER target/domain reporting parameters, an explicit
  independent-domain E-value decision boundary, HMM coverage and HMMER `Z`.
- Added deterministic reference sampling with separate random streams for each requested
  family and negative pool.
- Added exact threshold evaluation, complete confusion matrices, receiver-operating
  characteristic curves and Wilson 95% confidence intervals.
- Added operational cutoff application without describing calibration as independent
  biological validation.
- Added manifests containing input checksums, software versions, parameter values and
  coordinate conventions.
- Pinned Requests 2.34.2 for the optional UniProt retrieval helper.
- Excluded protein sequences, candidate identifiers, database exports, derived research
  results and third-party HMM files from the repository.
