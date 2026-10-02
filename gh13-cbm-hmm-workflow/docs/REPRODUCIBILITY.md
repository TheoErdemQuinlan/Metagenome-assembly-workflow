# Reproducibility checklist

Use this checklist for each reported analysis.

## Before execution

- Freeze the source-code commit and record the package version.
- Archive a resolved environment file or `pip freeze` output when exact reconstruction of
  transitive dependencies is required.
- Place biological inputs and third-party HMM files outside the repository.
- Verify the SHA-256 checksum of every immutable input.
- Keep input files unchanged for the duration of the run.
- Record the database release, retrieval date and provider.
- Confirm that protein identifiers are unique within the input FASTA.
- Define positive, hard-negative and easy-negative labels before inspecting score outcomes.
- State whether the goal is exploratory calibration or independent validation.

## Search parameters

- Record selected profile names.
- Record the target-reporting E-value threshold.
- Record the conditional-domain reporting E-value threshold.
- Record the independent-domain E-value threshold.
- Record minimum HMM coverage in percent.
- Record the explicit HMMER `Z` value and justify how it was defined.
- Record the Python, PyHMMER and package versions.
- If `cpus=0` is used, record that it requests PyHMMER's automatic worker selection rather
  than a fixed resolved worker count.
- Retain the generated manifest and exact command.

## Calibration parameters

- Record family name and score unit.
- Record the score assigned to references with no reported domain.
- Record all reference-class counts.
- Record the threshold selection rule before calculation.
- Retain unrounded bit-score thresholds.
- Report all four confusion-matrix counts with sensitivity and specificity.
- Report 95% confidence intervals and their method.
- Label reuse of calibration references as internal robustness analysis.
- Reserve “independent validation” for data excluded from every calibration decision.

## Reporting

- Distinguish domain-hit counts from unique-protein counts.
- Use one-based, end-inclusive coordinates unless a conversion is stated.
- Describe HMM evidence as family similarity rather than confirmed activity.
- State that biochemical properties require experimental evidence.
- Report exclusions, missing data and uncertain profiles explicitly.
- Archive outputs in a controlled research-data location, not this source repository.
