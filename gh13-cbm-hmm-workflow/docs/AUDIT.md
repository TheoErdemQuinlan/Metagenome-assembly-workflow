# Scientific and software audit

## Scope

The standalone release was reviewed for separation from the application database,
reproducibility of the profile-HMM calculation, correctness of threshold arithmetic,
terminology, provenance and accidental disclosure of biological data. The audit covers the
source code in this repository. It does not validate the biological labels assigned by an
external database or establish enzyme function.

## Corrections incorporated into version 0.1.0

### Exact comparison-threshold evaluation

An earlier project helper returned the receiver-operating characteristic row with the score
nearest to a requested threshold. That is incorrect when the threshold lies between observed
scores and can alter the confusion matrix by one or more records. The standalone
implementation evaluates the exact `score >= threshold` rule. Its curve lookup, when used,
selects the step-function state appropriate to that rule.

### Common reference-record schema

An earlier UniProt retrieval helper stored CAZy annotations under a field name that did not
match the downstream validation function. The standalone helper uses `tags` consistently for
both CAZy-derived and UniProt-derived records.

### Family-list-order-independent sampling

An earlier reference builder shared one pseudorandom stream across all requested families.
Changing family order could therefore alter the sampled records. The standalone builder
derives a separate deterministic stream for each family and negative pool. The family list
is normalized before iteration. Reservoir membership still depends on the order of records
in the source FASTA, so exact reproduction requires the same source-file content and order.

### Explicit HMMER search-space semantics

Earlier scanning code did not expose HMMER's effective E-value search space `Z`. The current
interface requires an explicit numeric value and records it in the manifest. It also
distinguishes HMMER's `E` and `domE` reporting parameters from the manually applied
independent-domain E-value decision boundary. Reporting thresholds default to a permissive
value and all three values are recorded. These changes make the statistical convention
inspectable; they do not retroactively change archived runs.

### Scientific terminology

The public interface uses `passes_operational_cutoff` rather than language implying that an
HMM score independently validates a domain. Documentation distinguishes family similarity,
threshold classification and biochemical function. Counts and coordinates include their
units, and confusion-matrix fields use complete descriptive names.

### Data isolation

The public tree excludes sequences, sample names, candidate identifiers, reference-accession
lists, database exports, result tables, cutoffs and third-party model files. Ignore rules
cover common sequence, read, HMM and DIAMOND formats. Release checks should additionally scan
the staged Git tree because ignore rules do not remove files already tracked by Git.

## Consequences for earlier outputs

Version 0.1.0 does not claim that earlier project results have been regenerated. Results that
depend on a between-score comparison threshold should be recalculated with exact threshold
semantics. Calibration samples should also be rebuilt if invariance to requested family order
is required. A scan should be repeated when its `Z` value or batching context was not recorded
and E-value comparability matters.

Bit-score rankings are unaffected by the choice of `Z`, but E-value-based retention can be
affected. Existing operational cutoffs should remain frozen until the consequences of any
changed reference sampling are assessed and documented.

## Remaining scientific limitations

- CAZy family labels can incorporate profile-based annotation, creating partial circularity
  when they are used to assess related profile HMMs.
- “Zero observed false positives” applies only to the sampled calibration negatives. Its
  uncertainty depends on the negative sample size.
- Selecting and evaluating a threshold on the same records does not estimate external
  predictive performance.
- Exact-sequence deduplication does not remove dependence among close homologues. Wilson
  intervals treat reference records as binomial observations and do not account for
  phylogenetic dependence; homology-aware held-out evaluation is preferable for external
  performance claims.
- Short profiles and homologous neighbouring families can be difficult to separate by bit
  score alone.
- Partial open reading frames can contain a real domain while failing to represent a complete
  protein or complete catalytic architecture.
- One protein can produce several domain hits. Domain-hit and unique-protein counts answer
  different questions.
- Family-profile evidence cannot establish alpha-amylase activity, substrate preference,
  temperature optimum, pH tolerance or detergent performance.

These limitations should accompany interpretation of downstream candidate panels.

## Release acceptance checks

A release is acceptable only when:

1. the unit-test suite passes in the pinned Python environment;
2. static checks pass without suppressing scientific logic errors;
3. an installation from `pyproject.toml` exposes the documented commands;
4. the staged Git tree contains no sequence or candidate data;
5. command help and output schemas agree with the documentation; and
6. the package build contains only the intended package source, metadata and licence.

Passing these checks validates implementation behaviour. It does not validate the underlying
biological labels or any candidate's experimental activity.
