# Input and output specification

## Protein FASTA input

The `scan` command accepts an amino-acid FASTA. Its default preflight requires non-empty,
unique identifiers and non-empty sequences containing accepted protein symbols. It removes
whitespace and normalizes case for validation. Terminal stop markers must be removed before
a streamed FASTA search; accepting and normalizing them in preflight would disagree with the
unchanged file subsequently read by PyHMMER. The separate in-memory API accepts and removes
one terminal stop marker because it passes its normalized sequence directly to PyHMMER.
Protein length and sequence coordinates are measured in amino-acid residues.

The duplicate-identifier check retains identifiers in memory. `--skip-fasta-preflight`
avoids that preliminary read and memory cost for very large files, but delegates residue and
format checks to PyHMMER and cannot guarantee identifier uniqueness.

Protein inputs are research data and must remain outside this repository.

## HMM input

`--hmms` accepts a HMMER profile file readable by `pyhmmer.plan7.HMMFile`. Profile names are
read from the HMM `NAME` field and normalized by removing a terminal `.hmm` suffix. Requested
profile names that are absent cause the command to stop.

## Domain-hit table

The `scan` command writes a tab-separated table with one row per retained domain alignment.

| Column | Definition | Unit |
| --- | --- | --- |
| `profile_name` | Normalized HMM profile name | none |
| `target_id` | Protein FASTA identifier | none |
| `target_full_sequence_evalue` | HMMER full-sequence target expected count | dimensionless |
| `independent_domain_evalue` | HMMER independent-domain expected count | dimensionless |
| `domain_bit_score_bits` | Domain alignment bit score | bits |
| `hmm_start_match_state` | First aligned HMM match state | one-based position |
| `hmm_end_match_state` | Last aligned HMM match state | one-based position, inclusive |
| `hmm_length_match_states` | Full model length | match states |
| `target_start_amino_acid_position` | First aligned target residue | one-based amino-acid position |
| `target_end_amino_acid_position` | Last aligned target residue | one-based amino-acid position, inclusive |
| `target_length_amino_acids` | Full target sequence length | amino-acid residues |
| `hmm_alignment_coverage_percent` | Aligned HMM span divided by model length | percent |

The table is sorted by target identifier, target start position, profile name and independent-
domain E-value. Several rows can refer to one protein, either because different profiles
match or because a repeated domain is present.

Conditional-domain E-values are not emitted. They depend on HMMER's separate effective
domain search-space convention (`domZ`) and are not used by this workflow's retention rule.

## Scan manifest

The JSON manifest contains:

- package, Python and PyHMMER versions;
- SHA-256 checksum and byte size for each input;
- a sequence-free protein-input summary containing record count, total amino-acid residues
  and minimum and maximum protein lengths when FASTA preflight is enabled;
- searched profile names and profile count;
- HMMER target- and conditional-domain reporting E-value thresholds;
- the independent-domain E-value decision threshold;
- minimum HMM alignment coverage in percent;
- the explicit HMMER `Z` value;
- processor count;
- coordinate convention;
- output filename, row count and SHA-256 checksum; and
- UTC creation time.

The manifest contains input filenames and hashes, but no biological sequences.

## Labelled score table

The `calibrate-score` command reads a tab-separated table with three required columns.

| Column | Allowed content | Unit |
| --- | --- | --- |
| `profile_name` | Exact HMM profile name | none |
| `class_label` | `positive`, `hard_negative`, `easy_negative` or `negative` | none |
| `domain_bit_score_bits` | Best profile score; zero if no domain was reported | bits |

The table must contain one row per unique reference sequence. It deliberately omits
reference identifiers and sequences, so the command cannot detect duplicate records from
this table alone. Deduplicate the labelled references before producing the score table and
archive the sequence-free reference manifest separately.

Rows for other profiles are ignored when `--family` selects one profile. Unknown labels,
non-numeric scores and non-finite scores cause the command to stop. At least one positive and
one negative score are required.

## Calibration report

The JSON calibration report records:

- selected family and score units;
- input table filename and SHA-256 checksum;
- count of each reference class;
- the `score >= threshold` decision rule;
- selected operational cutoff in bits;
- confusion counts, sensitivity, specificity and Wilson 95% confidence intervals;
- receiver-operating characteristic area;
- threshold-robustness results; and
- the complete evaluated receiver-operating characteristic curve.

Rates and confidence-interval limits are proportions from zero to one. Confusion counts are
numbers of reference records. A calibration report also exposes its selected value under the
`cutoffs` object so it can be supplied directly to `apply-cutoff`.

## Classified domain-hit table

`apply-cutoff` retains every input column and appends:

| Column | Definition | Unit |
| --- | --- | --- |
| `operational_cutoff_bit_score_bits` | Family cutoff, blank if unavailable | bits |
| `operational_cutoff_classification` | Score-rule result | none |

Allowed classifications are `passes_operational_cutoff`, `below_operational_cutoff` and
`uncalibrated_family`. These values do not encode a biological activity claim.
