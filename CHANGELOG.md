# Changelog

<!-- release-header:start -->
**COSMolKit** is a pure-Rust cheminformatics toolkit.

[Source repository](https://github.com/cosmol-studio/COSMolKit) ·
[Documentation](https://kit.cosmol.org/) ·
[Web tools](https://tools.cosmol.org/) ·
[Rust crate](https://crates.io/crates/cosmolkit) ·
[Python package](https://pypi.org/project/cosmolkit/).
<!-- release-header:end -->

All notable changes to COSMolKit are recorded in this file.

The text between the `release-header` markers is prepended to every GitHub
Release so the project, documentation, and web tools remain easy to discover.
For each release, move the completed entries from `Unreleased` into a section
named `## [x.y.z] - YYYY-MM-DD`. The release workflow extracts the section whose
version matches the pushed `v*` tag and fails when that section is missing or
empty.

## [Unreleased]

### Added

- Added source-backed AtomPair sparse-count, folded-count, sparse-bit, and
  explicit-bit fingerprints in Rust and Python, including 2D and 3D distances,
  chirality, count simulation, selectors, custom atom invariants, complete
  provenance, metadata/JSON restoration, and ordered batch execution through
  one generator core.

### Changed

- Reproduced RDKit's logically read-only 2D and resolved-conformer 3D distance
  matrix caching in shared molecule computed state, with topology/coordinate
  invalidation and clone-local synchronization. AtomPair and distance-geometry
  consumers reuse this single matrix implementation.

### Fixed

- Corrected `Molecule.from_rdkit()` so its default preserves the imported RDKit
  graph while copying hybridization and computing the explicit-valence cache.
  `sanitize=True` remains the explicit full-sanitize path and `sanitize=False`
  remains the raw unprepared path.

### Validation

- The complete AtomPair audit processed all 2,897,819 ChEMBL 37 source records
  over 128 shards. All 2,897,804 mutually parseable records matched pinned
  RDKit across ten profiles, four result forms, and one complete provenance
  output, totaling 118,809,964 exact comparisons with zero mismatch or failed
  task.
- After the source-defined distance-cache port, all 24 deterministic scalar
  AtomPair benchmarks and all eight 3,072-molecule batch benchmarks were faster
  than the pinned RDKit build. The former 80-atom 2D long tail fell from
  3.224–4.768 times RDKit to 0.458–0.493 times without changing output.
- The maintained 5,000-molecule gate runs ten AtomPair profiles as a committed
  exhaustive regression. Neither the full-corpus audit nor the maintained gate
  uses sampling, tolerance, similarity thresholds, or fallback fingerprints.

## [0.2.12] - 2026-08-19

This release adds source-backed RDKit topological and Avalon fingerprints,
moves Python InChI conversion onto `Molecule`, strengthens operation contracts,
and incorporates the source-level fixes found by the complete ChEMBL 37 parity
profile.

### Added

- Added source-backed RDKit topological fingerprints to Rust and Python,
  including path-length limits, branched-path enumeration, hydrogen and bond
  controls, custom atom invariants, source-compatible random-bit expansion,
  density folding, and from-atoms selection.
- Added `topological_fingerprint_with_output()` and the typed
  `TopologicalFingerprintResult`, exposing the source `atomBits` and `bitInfo`
  provenance branches without changing the ordinary fingerprint return type.
- Added source-backed Avalon/REACCS explicit-bit fingerprints with configurable
  vector length, query preprocessing, typed feature flags, aromaticity passes,
  and source-compatible byte rounding.
- Added the repeatable ChEMBL 37 parity and composition profile: 2,897,819
  source records across 22 sharded phases, exhaustive option matrices,
  operation-order permutations, scalar/batch comparisons, and shared-object
  concurrent reads.

### Changed

- Replaced the RDKit-shaped Python InChI namespace with the object-oriented
  `Molecule.to_inchi()`, `Molecule.to_inchi_key()`, and
  `Molecule.from_inchi()` methods plus the top-level `inchi_to_key()` function.
  The former `Chem.Mol*` entry points and `InchiToInchiKey()` were removed.
- Removed self-reported operation `NoOp` outcomes and `allows_noop` metadata.
  Operation obligations now follow actual block lifecycle, mappings, writes,
  and derived-state traces.
- Narrowed operation read views to their declared capabilities and added scoped
  mutation helpers that always return checked-out blocks on success or error.
  In-place operations retain their efficient basic failure guarantee: they do
  not roll back partial source-defined changes, but internal molecule storage
  remains complete after a returned error.
- Made source-compatible computed descriptor caching independent from
  user-visible molecule state and preserved the correct cache lifecycle across
  clones, operation invalidation, and repeated Crippen calls.
- Updated all Rust crates and the Python package to version 0.2.12. GitHub
  Release notes are now assembled from the matching version section in this
  changelog.

### Fixed

- Corrected Morgan `num_bits_per_feature` propagation and its interaction with
  every requested `AdditionalOutput` branch.
- Aligned reverse InChI conversion with the source's final stereochemistry
  cleanup, double-bond reconstruction, tautomer endpoint allocation, and
  sanitize/remove-H property-cache behavior.
- Corrected explicit-hydrogen state, distance-geometry bounds, canonical and
  parameterized SMILES writing, SVG output, and MolBlock/SDF atom, bond,
  coordinate, and stereo-direction behavior found by the large-corpus profile.
- Aligned conformer outcomes, UFF/MMFF preparation, aromatic-ring atom typing,
  parameter-unavailable classification, optimizer rejection, energy, and final
  coordinate paths with pinned RDKit.
- Hardened operation-contract enforcement for uncommitted writable blocks,
  undeclared reads, cache effects, topology mappings, and in-place error paths.

### Validation

- The complete ChEMBL 37 audit uses pinned RDKit 2026.03.1 and evaluates
  2,897,804 mutually parseable records across 22 sharded chemistry,
  composition, batch, and concurrency phases.
- Ten InChI option profiles produce 28,543,620 matching InChI strings and
  28,543,620 matching InChIKeys. Each complete 768-profile SMILES writer matrix
  performs 2,192,150,016 comparisons; three complete matrices were executed.
- Descriptor coverage performs 84,036,316 comparisons over 29 branches;
  distance-geometry coverage compares 2,757,910,995 matrix entries.
- Scalar-versus-batch validation covers 1,027,660 records at both one and eight
  jobs, totaling 12,331,920 output comparisons. Repeated concurrent-read runs
  each compare 1,032,384 values and 1,032,384 graph states.
- Every retained ChEMBL finding was replayed after its source-aligned fix across
  the complete affected branch set, with no unexplained mismatch. This includes
  descriptor, Morgan, explicit-H, bounds, InChI, SMILES, SVG, MolBlock/SDF,
  conformer, and force-field regressions; it is not represented as a post-fix
  rerun of the entire ChEMBL corpus.
- The complete ChEMBL 37 fingerprint audit processed all 2,897,819 source
  records across 128 shards: 2,897,804 were mutually parseable and 15 were
  rejected by both libraries. RDKit topological fingerprints matched exactly
  across 14 vector profiles (40,569,256 comparisons) and two complete
  `atomBits`/`bitInfo` provenance profiles (5,795,608 comparisons); Avalon
  matched exactly across 23 explicit-bit profiles (66,649,492 comparisons).
  All 113,014,356 comparisons completed with zero mismatch or failed task.
- The maintained 5,000-molecule gates continue to run the same 14 topological
  profiles (70,000 comparisons) and 23 Avalon profiles (115,000 comparisons)
  as committed exhaustive regressions. Neither the full-corpus audit nor the
  maintained gates use sampling, tolerance, similarity thresholds, or fallback
  fingerprints.
- The 3,480 archived ChEMBL InChIs and 3,480 archived InChIKeys that differ from
  current output are reference-corpus version differences: COSMolKit matches
  the current pinned RDKit output on those records.

### Upgrade notes

- Python callers must replace `Chem.MolToInchi(mol)` with `mol.to_inchi()`,
  `Chem.MolToInchiKey(mol)` with `mol.to_inchi_key()`,
  `Chem.MolFromInchi(text)` with `Molecule.from_inchi(text)`, and
  `InchiToInchiKey(text)` with `inchi_to_key(text)`.
- The former placeholder fingerprint signatures have been replaced by their
  source APIs. `topological_fingerprint()` now uses `min_path`, `max_path`,
  `fp_size`, `num_bits_per_feature`, `use_hs`, `target_density`, `min_size`,
  `branched_paths`, `use_bond_order`, `atom_invariants`, and `from_atoms`;
  `avalon_fingerprint()` uses `n_bits`, `is_query`, and `bit_flags`.
- Code that requires rollback on an in-place failure should use the value-style
  operation instead. In-place methods may retain valid partial changes after an
  error, while guaranteeing complete internal storage.

Full comparison: [v0.2.11...v0.2.12](https://github.com/cosmol-studio/COSMolKit/compare/v0.2.11...v0.2.12)

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
