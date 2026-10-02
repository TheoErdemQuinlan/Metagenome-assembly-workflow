# Computational methods

## Scope

This workflow searches protein sequences against family-level profile hidden Markov models
and evaluates score boundaries for specified profiles. It is suitable for reproducible
screening of predicted metagenomic proteins when the input proteins and third-party HMM
resources are maintained separately from the code.

The method assigns sequence-similarity evidence at the family or subfamily level represented
by the selected model. It does not assign an experimentally established enzyme activity.
Within GH13, multiple hydrolytic, branching, debranching and transferase activities can share
homologous catalytic machinery. Functional interpretation therefore requires further
evidence beyond a family HMM match.

## Inputs

The search requires:

1. a protein FASTA containing unique, non-empty identifiers;
2. a HMMER-format file containing one or more profile HMMs;
3. an explicit list of profiles when the complete HMM file should not be searched; and
4. the statistical and coverage parameters recorded in the run manifest.

The code does not distribute dbCAN or CAZy files. Profile HMMs and labelled references are
external scientific inputs whose versions, checksums and usage terms must be managed by the
researcher.

## Profile-HMM search

Profiles are loaded with `pyhmmer.plan7.HMMFile` and searched against amino-acid sequences
using `pyhmmer.hmmer.hmmsearch`. HMMER's `E` and `domE` parameters are target and
conditional-domain reporting thresholds, respectively. They are set independently and
permissively. Neither parameter is labelled as the independent-domain decision threshold.
The workflow reads `domain.i_evalue` and applies that boundary directly after the search.

For an aligned domain, HMM coverage is calculated as:

\[
\mathrm{HMM\ coverage} =
\frac{\mathrm{hmm\_to} - \mathrm{hmm\_from} + 1}
     {\mathrm{HMM\ length\ in\ match\ states}}.
\]

A domain is retained when both conditions hold:

\[
E_{\mathrm{independent\ domain}} \leq E_{\mathrm{threshold}}
\]

and

\[
\mathrm{HMM\ coverage} \geq C_{\mathrm{minimum}}.
\]

The applied project screen used an independent-domain E-value threshold of
`1 × 10⁻¹⁵` and minimum HMM coverage of `35%`. Its historical implementation also supplied
`1 × 10⁻¹⁵` as HMMER's target-reporting parameter. This standalone release uses permissive
defaults of `10` for both reporting parameters and records them separately from the actual
decision rule. These are operational screening settings rather than universal biological
boundaries. Any observed bit score associated with the weakest retained hit is a consequence
of that run; it is not the definition of the original screen.

The output records bit score in bits, independent-domain E-value as a dimensionless expected
count under the HMMER statistical model and stated search space, HMM coverage in percent, and
HMM and target coordinates as one-based, end-inclusive positions. Conditional-domain
E-values are not emitted because they use the separate `domZ` convention and do not enter the
retention rule. Each output row is a domain hit. Counts of rows must not be reported as counts
of unique proteins without deduplication by target identifier.

## E-value search space

Target and independent-domain E-values depend on the effective search space `Z`. Automatic
`Z` is derived by HMMER
from the target search context. Consequently, E-values may differ when the same sequences are
searched in batches or when a `hmmsearch` arrangement is compared with an `hmmscan`
arrangement. Raw bit scores are not rescaled by `Z`.

Conditional-domain E-values use a separate effective search space, `domZ`. This workflow
does not report them or use them for retention, so an automatic `domZ` cannot alter the final
decision rule. The permissive `domE` value is nevertheless recorded as a pipeline reporting
parameter.

The command-line interface therefore requires a positive numeric
`--evalue-search-space` value and records it in the manifest. For a transposed search intended
to reproduce the statistical convention of scanning proteins against an HMM database, the
appropriate database-wide profile count should be recorded rather than inferred from an
arbitrary selected subset. The audited dbCAN resource contained 875 profiles, so that run
uses `Z = 875` even when only selected profiles are loaded.

## Reference-set construction

The reference utilities accept a locally supplied CAZy FASTA and parse both accession-first
and family-tag-first headers. Family tags are normalized to their base family where required.
The unclassified `CBM0` tag is treated as ambiguous and is not used as clean negative
evidence for a named CBM family.

Sample reservoirs have fixed maximum sizes. Each requested positive family and each negative
pool has a separate deterministic random-number stream derived from the stated seed and
class name. This makes a family's sampled records invariant to the order in which families
are requested. Reservoir membership remains dependent on source-record order, so exact
reproduction requires the same source-file content and order. Exact duplicate amino-acid
sequences are removed before reservoir sampling; when duplicate records carry different
tags, their tags are united before class assignment. The preliminary digest-to-tag map
scales with the number of unique reference sequences and must therefore be included in
memory planning for a large reference collection.

Malformed reference sequences stop construction by default
(`invalid_sequence_policy="error"`). The explicit `"exclude"` policy omits them and records
the number excluded in `ReferencePools`. A scientific run that uses exclusion must report
that count. The generated sequence-free reference manifest records the accession, class and
SHA-256 sequence digest for every retained member, allowing the labelled set to be audited
without placing its sequences in this source repository.

For one family, the classes are:

- **positive reference:** annotated with the target family;
- **hard negative reference:** drawn from related CAZy families but without the target-family
  tag or an ambiguous tag;
- **easy negative reference:** drawn from other annotated CAZy families and likewise lacking
  the target-family and ambiguous tags.

The composition, sample sizes, seed, release, retrieval date and immutable accession list
must be archived outside this public code repository. CAZy annotations may themselves use
profile-based evidence, so agreement with those labels is partly circular and may overstate
performance on independently established biology.

## Score calibration

Each labelled reference is represented by its best domain bit score for the profile under
evaluation. A reference with no reported alignment is assigned zero bits for the calibration
table. Positive and negative labels are defined before threshold selection.

For a bit-score threshold \(t\), a score is called positive when:

\[
s \geq t.
\]

The implementation evaluates every unique observed score, with lower and upper sentinel
values. It calculates:

- correctly detected positives;
- missed positives;
- incorrectly accepted negatives;
- correctly rejected negatives;
- sensitivity;
- specificity;
- false-positive rate; and
- false-negative rate.

If a requested comparison threshold falls between observed scores, its confusion matrix is
calculated from the exact `score >= threshold` rule. It is not approximated by selecting the
nearest observed threshold.

The implemented operational cutoff is the lowest evaluated bit score at which no calibration
negative is observed above or equal to the boundary. This maximizes sensitivity subject to
the stated zero-observed-false-positive constraint. The word “observed” must be retained:
finite negative samples cannot demonstrate that the population false-positive rate is zero.

Sensitivity and specificity are reported with two-sided Wilson 95% confidence intervals.
The receiver-operating characteristic area is calculated by trapezoidal integration over
false-positive rate and sensitivity. Robustness is summarized at the selected threshold and
at a stated offset below and above it, in bits. Because the same references are reused, this
is a boundary-sensitivity analysis rather than out-of-sample validation.

## Applying a cutoff

A domain score greater than or equal to the family-specific threshold is labelled
`passes_operational_cutoff`; a lower score is labelled `below_operational_cutoff`. A profile
without a supplied cutoff is labelled `uncalibrated_family` or causes the command to stop
when strict mode is requested.

The label describes a numerical decision only. It should not be translated into “validated
domain,” “confirmed alpha-amylase” or a claim about assay performance. Domain boundaries,
active-site conservation, open-reading-frame completeness, architecture, structural
similarity and experimental measurements remain separate evidence.

## Software environment

The audited release targets Python 3.12 and pins PyHMMER 0.11.0 and Requests 2.34.2. Run
manifests record the Python, package and PyHMMER versions. Reproducibility also requires the
external HMM and protein-input checksums, because software versioning alone cannot identify
the scientific inputs.
