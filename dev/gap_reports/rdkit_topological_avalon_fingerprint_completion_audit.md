# RDKFingerprint and Avalon Completion Audit

This audit records the final state of the active source-port plan. Avalon is
closed for its selected explicit-bit API. The plan's broader workspace gate is
not green because pre-existing repository checks outside these families still
fail; those failures are listed rather than hidden.

## Closed scope

- Avalon source-backed Rust implementation covers the selected RDKit adapter
  boundary (`n_bits`, `is_query`, `bit_flags`) with source markers across
  conversion, preprocessing, hashing, traversal, rings, aromaticity, and all
  active flag families.
- RDKit Python oracle generated 5,000 rows × 23 profiles = 115,000 exact
  comparisons. Release strict Rust test passed with zero mismatches.
- Focused Rust boundary test passed for Python/C++ defaults, query mode,
  non-SSS flags, and 8/9/31/32/33/511/512/513-bit boundaries.
- Rust combined interaction tests passed for Avalon, RDK/topological, Morgan,
  and MACCS ordering, repeated calls, parallel determinism, and molecule
  non-mutation.
- Python focused parity/API tests passed (`64 passed`); complete Python suite
  passed (`574 passed, 37 skipped`). Stub surface passed (`2 passed`).
- `cargo fmt --all -- --check`, workspace feature checks, Rustdoc, Sphinx
  warnings-as-errors, and the fingerprint example passed.
- Support metadata, parity scope, inventory, README, Python docs, and generated
  stub now describe the validated Avalon explicit-bit surface and its explicit
  out-of-scope overloads.

## Reproducible artifacts

- Generator: `tools/testdata/rdkit/_generate_avalon_fingerprint_golden.py`
- Profile: `tools/testdata/rdkit/avalon_fingerprint_profile.json`
- Corpus golden SHA-256:
  `235b66b30ebb73e966c540e4d89f74736baaa72b4ef2a25e0d592e5d77ba7c2f`
- Focused golden SHA-256:
  `12b131ed3e1722652f6d9a274f2193727fe9c001fb6f47d52e06426ad983adfb`
- Generation command, engine identity, and counts:
  `tmp/parity-audit/avalon_fingerprint_generation.log`
- Detailed family validation: [`avalon_fingerprint_full_port_validation.md`](avalon_fingerprint_full_port_validation.md)

## Remaining blockers

The first workspace release test run reached 2,738 passing tests but exposed
five stale expected-data checksum guards. The repository-standard
`smiles_small --suite all --jobs 4` preparation was then run, and the complete
workspace release strict matrix passed; the refreshed run includes all
workspace test and doctest targets with zero failures.

`cargo clippy -p cosmolkit-core --tests --features op-contracts-strict -- -D warnings`
remains non-clean because the existing source-port crates emit 1,035 naming
and related warnings, predominantly in `cosmolkit-inchi`. Suppressing or
renaming those copied source identifiers would violate the source-reproduction
marker and provenance policy, so this is retained as a repository-level
follow-up rather than misreported as an Avalon defect.

The active plan therefore has a completed Avalon family validation and one
remaining repository-wide hygiene issue in the legacy clippy baseline. No
unsupported Avalon overload is represented as implemented.
