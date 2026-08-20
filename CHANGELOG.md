# Changelog

## [1.0.3] - 2026-XX-XX

In development

- Add `seed_groups`, replacing the previous `include_samples` option. It is a
  list of seed *groups*, each a list of strain IDs that are inserted into the
  ARG together as a single local tree. A tree is inferred over the group's
  haplotypes, the haplotype of the group's inferred ancestor is matched against
  the ARG, and the whole group is then attached at that placement. A group is
  inserted on the minimum date over its members, and members with later dates
  keep their real dates via negative ("in the future") node times. The previous
  `include_samples` formats (bare strain IDs, or `(strain, match_date)` tuples)
  are no longer supported.

- Fix retrospective matching on a day with no new samples. The table of samples
  already in the ARG was only rebuilt when there were samples to match, so on
  such a day the retrospective query ran against the previous day's table and
  could add a second copy of a sample that had just been inserted.

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
