# Resource and software provenance

This repository contains code only. The records below identify the environment and HMM
resource examined during extraction of the standalone implementation; the referenced HMM
file is not distributed here.

## Audited software environment

| Component | Version |
| --- | --- |
| Python | 3.12.14 |
| PyHMMER | 0.11.0 |
| Requests | 2.34.2 |
| pytest | 8.4.2 |
| Standalone package | 0.1.0 |

Runtime dependencies are pinned in `pyproject.toml`. The exact environment used for a new
analysis should be reported from that run rather than inferred from this table.

## Audited dbCAN HMM resource

| Field | Value |
| --- | --- |
| Local filename | `dbCAN.hmm` |
| Associated dbCAN software release | 5.2.9 |
| Resource snapshot label | `db_v5-2-9_5-5-2026` |
| Local retrieval date | 2026-09-18 |
| File size | 129,842,960 bytes |
| Number of profiles | 875 |
| MD5 | `8b277c30b9c3ac5c36c28944241d33f1` |
| SHA-256 | `396c791fb13defab152864d8046687ef03b492937ec6f9bbab008bc77cbaa7a5` |

The local file's byte size and MD5 match the object in the official
[`db_v5-2-9_5-5-2026` S3 snapshot](https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/dbCAN.hmm),
checked on 2026-10-02. The SHA-256 above identifies the local file more strongly. Users
should obtain resources through the official
[run_dbcan distribution](https://github.com/bcb-unl/run_dbcan), verify their own download and
retain provider licensing or access information. A matching filename alone is not adequate
provenance because databases can be replaced without a filename change.

## Reference data

The applied calibration used locally obtained CAZy-labelled references, including related
CAZy families as hard negatives. No reference sequence, accession list or database export is
present in this repository. A reproducible analysis archive should separately retain:

- provider and release or retrieval date;
- original query or extraction command;
- checksum of the immutable source file;
- complete accession list assigned to each label;
- requested and realized class sizes;
- deduplication rule;
- sampling seed; and
- checksum of the resulting sequence-free score table.

Live UniProt requests are mutable. When UniProt is used, preserve the query, retrieval date,
returned accession list and response checksum. The source code returns UniProt records using
the same `tags` field expected by the CAZy reference-processing functions.
