# Changelog

## [1.0.3] - 2026-XX-XX

In development

- Seed samples listed in `include_samples` are now matched separately from the
  daily samples, with recombination effectively disallowed. These samples are
  usually far diverged from the current ARG, which previously led to spurious
  recombinations.

- Add basic support for non-SARS-CoV-2 genomes via an optional reference FASTA.
  Supply `--reference` to `import-alignments` and a `reference_fasta` key in the
  inference config; both default to the built-in SARS-CoV-2 reference, so
  existing workflows are unchanged. When `reference_fasta` is set, a
  `reference_date` config key is required (the date assigned to the reference).
  The redundant `genbank_id` reference-sequence metadata field has been dropped.

## [1.0.2] - 2026-03-05

Maintenance release.

- Require Python >= 3.11

## [1.0.1] - 2025-11-28

Maintenance release.

- Minor packaging update to track newly released tskit 1.0

## [1.0.0] - 2025-11-23

Initial stable release of sc2ts.
