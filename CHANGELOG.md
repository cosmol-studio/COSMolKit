# Changelog

<!-- release-header:start -->
**COSMolKit** is a pure-Rust cheminformatics toolkit.

[Source repository](https://github.com/cosmol-studio/COSMolKit) ·
[Documentation](https://kit.cosmol.org/) ·
[Web tools](https://tools.cosmol.org/) ·
[Rust crate](https://crates.io/crates/cosmolkit) ·
[Python package](https://pypi.org/project/cosmolkit/)
<!-- release-header:end -->

All notable changes to COSMolKit are recorded in this file.

The text between the `release-header` markers is prepended to every GitHub
Release so the project, documentation, and web tools remain easy to discover.
For each release, move the completed entries from `Unreleased` into a section
named `## [x.y.z] - YYYY-MM-DD`. The release workflow extracts the section whose
version matches the pushed `v*` tag and fails when that section is missing or
empty.

## [Unreleased]

### Changed

- Replaced the RDKit-shaped Python InChI surface with COSMolKit's object-oriented
  naming: `Molecule.to_inchi()`, `Molecule.to_inchi_key()`,
  `Molecule.from_inchi()`, and the top-level `inchi_to_key()` function. The
  former `Chem.Mol*` namespace and `InchiToInchiKey()` name were removed.

## [0.2.11] - 2026-08-11

This release substantially strengthens source-level and behavioral parity with
the pinned RDKit and official InChI implementations. It also establishes
explicit support boundaries for fingerprint APIs that are not yet complete.

### Added

- Added source-backed InChI canonicalization, stereochemistry, normalization,
  reverse conversion, and supporting source types, with line-level provenance
  markers and parity coverage.
- Added direct `cosmolkit-inchi` examples for InChIKey generation, conversion
  from InChI to the supported neutral-graph representation, and conversion from
  that representation back to InChI.
- Added broader MMFF and UFF parity fixtures, parameter coverage, coordinate and
  minimization comparisons, golden-data generators, and regression reports.
- Added a source-port plan for RDKit topological and Avalon fingerprints.

### Changed

- Aligned MMFF and UFF atom typing, parameter lookup, force-field construction,
  torsion handling, minimization behavior, and topology-operation contracts
  more closely with RDKit.
- Improved source-aligned InChI hot paths while preserving the official
  algorithm and output behavior.
- Updated Rust crates and the Python package to version 0.2.11.
- Hardened the 5000-row InChI coverage job for the larger stack requirement of
  the intentionally unoptimized coverage build.

### Fixed

- Registered the existing borrowed optional-block read operation as an allowed
  sanitize capability, restoring the strict operation-capability test without
  weakening access validation.
- Corrected smaller stereochemistry, SMARTS documentation, oracle-script, and
  source-marker issues found during the release audit.

### Upgrade notes and known limitations

- `Molecule.topological_fingerprint()` and
  `Molecule.avalon_fingerprint()` remain unsupported. They fail explicitly
  instead of returning approximate or substituted fingerprints; their
  source-exact ports are planned for a later release.
- The standalone neutral-graph InChI examples intentionally expose only their
  documented supported boundary. Unsupported graph chemistry returns a
  structured error rather than a heuristic conversion.
- Existing Morgan and MACCS fingerprint users do not need to migrate for this
  release.

Full comparison: [v0.2.10...v0.2.11](https://github.com/cosmol-studio/COSMolKit/compare/v0.2.10...v0.2.11)
