# Changelog

<!-- release-header:start -->
**COSMolKit** is a pure-Rust cheminformatics toolkit.

[Source repository](https://github.com/cosmol-studio/COSMolKit) ·
[Documentation](https://kit.cosmol.org/) ·
[Web tools](https://tools.cosmol.org/) ·
[Rust crate](https://crates.io/crates/cosmolkit) ·
[Python package](https://pypi.org/project/cosmolkit/).
<!-- release-header:end -->

<!-- release-footer:start -->
COSMolKit advances its minor compatibility series through successive prime
numbers, just for fun. 😈
<!-- release-footer:end -->

All notable changes to COSMolKit are recorded in this file.

The text between the `release-header` markers is prepended to every GitHub
Release, and the text between the `release-footer` markers is appended. For
each release, move the completed entries from `Unreleased` into a section named
`## [x.y.z] - YYYY-MM-DD`. The release workflow extracts the section whose
version matches the pushed `v*` tag and fails when that section is missing or
empty.

## [Unreleased]

This release adds the public modern RDKit CIPLabeler, ordinary molecular
alignment and RMSD, source-backed Pattern and Layered fingerprints, Rust
tautomer enumeration, and typed potential-stereo and lazy stereoisomer
enumeration. It also completes the ordinary-molecule SMARTS parser, query,
writer, and substructure-matching path through source-backed Rust
implementations with explicit state lifecycle and pinned parity validation.

### Added

- Added registered Rust and Python CIP-label assignment for complete molecules
  and selected atoms or bonds, including source-width recursion limits and
  typed `R/S/r/s`, `E/Z`, and `M/P/m/p` descriptor queries on atoms and bonds.
  Value-style and explicit trailing-underscore in-place forms share the same
  molecular-context implementation; neither delegates to the legacy CIP-rank
  helper.
- Added focused RDKit oracle fixtures, Rust and Python examples, generated
  expected-data support, a maintained exhaustive 5,000-row CI gate, and a
  reproducible ChEMBL 37 parity phase for the complete observable CIP state.
- Added source-backed Rust SMILES parser options so sanitization and hydrogen
  removal can be selected independently when constructing exact intermediate
  states.
- Added the canonical query-bearing `Molecule` SMARTS parser/compiler and
  source-backed atom, bond, molecule, fragment, and CXSMARTS writers in Rust
  and Python.
- Added ordinary-molecule substructure matching with recursive queries,
  query-query and generic matching, callbacks, source ordering, uniqueness and
  match limits, and typed matching parameters.
- Added the source-backed RDKit Layered fingerprint 0.7.0 to Rust and Python,
  including all six active layers, linear and branched rooted paths, masks,
  seeded atom counts, structured errors, complete result state, and ordered
  batch execution through the sole shared fingerprint core.
- Added source-backed ordinary-molecule Pattern fingerprints in Rust and
  Python with the 13 pinned built-in SMARTS patterns, ordinary and
  tautomer-aware hashing, query-bearing molecule behavior, configurable
  nonzero widths, and ordered parallel batch execution through one
  compile-once table and one shared core.
- Added deterministic Pattern oracle generation, focused/small/5,000-row exact
  regression gates, screening and architecture guards, composition tests, a
  reproducible ChEMBL 37 phase, and representative cross-engine benchmarks.
- Added source-backed Rust tautomer catalogs, enumeration options, ordered
  multi-molecule results, canonical selection and scoring, callbacks, current
  and V1 transform sets, structured errors, and focused, PCS, 5,000-row, and
  ChEMBL 37 validation workflows.
- Added source-backed ordinary RDKit MolAlign in Rust and Python with typed
  explicit/automatic atom maps, weighted and reflected transforms, read-only
  best and coordinate-frame RMSD, all-conformer pair output, and explicit
  value-style or trailing-underscore coordinate mutation. O3A remains a
  separately scoped future capability.
- Added project-native Rust and Python potential-stereo analysis with ordered,
  typed atom/bond centers, specified state, descriptors, permutations, and
  controlling atoms, returning an isolated analyzed molecule without mutating
  the source.
- Added a single lazy Rust and Python stereoisomer iterator over the pinned
  RDKit Python behavior, including exact option defaults, arbitrary-width
  counts, atom, double-bond, and enhanced-group flippers, exhaustive and finite
  unique-random configurations, Python-compatible integer and custom random
  sources, canonical uniqueness, optional embedding, and structured lazy
  errors.
- Added focused source fixtures, committed 152-row and 5,000-row exact parity
  profiles, a complete reproducible ChEMBL 37 phase, architecture guards,
  composition tests, public examples, and deterministic stage benchmarks for
  potential stereo and stereoisomer enumeration.

### Changed

- Modeled RDKit computed-property membership explicitly for molecule, atom,
  and bond properties. Sanitization and hydrogen operations now reproduce
  `clearComputedProps()` by membership, preserving same-named user properties
  that were not registered as computed.
- Added a generated CIP-state policy to every molecule operation. Modern CIP
  assignment is the only assignment owner; other operations explicitly
  preserve or clear computed CIP state according to their source transition.
- Advanced molecule binary state encoding while retaining reads of earlier
  archive and legacy payloads, so CIP values and computed-property
  classification survive clone and serialization roundtrips.
- Consolidated every descriptor, fingerprint, force-field, coordinate, SDF,
  and public SMARTS consumer onto one typed query model, parser, writer, and
  matcher, removing the parallel `SmartsMolecule` representation and
  consumer-local conversion paths.
- Extended the operation-contract system with validated ordered multi-molecule
  output while retaining the existing value-style and in-place separation for
  single-molecule operations.
- Preserved RDKit's experimental classification for Layered fingerprints.
  COSMolKit returns the source-documented bond-path result for unrooted linear
  calls instead of reproducing the pinned implementation's input-dependent
  process crash.
- Consolidated depiction, distance geometry, MolAlign, and molecular
  coordinate transforms onto one private pure-Rust quaternion/Jacobi alignment
  and 4x4 transform implementation; pinned source behavior required no
  numerical FFI backend.

### Fixed

- Corrected the shared SMILES stereo dispatcher so writer-side reassignment
  with stereo cleanup disabled stops before RDKit's `cleanIt`-only ring and
  terminal cleanup, preserving canonical disconnected-fragment, bridgehead,
  and directional double-bond stereochemistry.
- Corrected MMFF force-field construction for sparse or explicitly named
  conformer IDs by resolving the source conformer ID exactly once and reading
  coordinates from that selected storage row.
- Completed the modern CIP digraph, configuration, auxiliary-label, and
  sequence-rule call paths, including source-prefix sorter ownership in Rule4b
  and Rule5New, exact selected-call dispatch, shared recursion budgeting,
  atropisomer carrier ordering, and source-observable partial state on errors.
- Removed duplicate atropisomer carrier extraction from the CIP path and made
  it reuse the source-backed shared topology helper.

### Validation

- The complete modern CIPLabeler audit processed all 2,897,819 ChEMBL 37
  source records over 128 shards. Fifteen records were rejected by both
  parsers and 43,442 accepted records exceeded the configured 80-atom phase
  boundary. Every one of the remaining 2,854,362 records matched pinned RDKit
  after full, selected-atom, selected-bond, and empty-selection assignment,
  totaling 11,417,448 exact complete-state comparisons with zero mismatch or
  retained finding.
- Comparisons include success and exact error state, `_CIPComputed`, atom and
  bond `_CIPCode` and `_CIPNeighborOrder`, atom `_CIPRank`, computed-property
  membership, chiral tags, bond stereo, and stereo atoms. Focused fixtures also
  cover all ten emitted descriptor spellings, repeated calls, recursion
  limits, malformed state, and the declared dispatcher boundary.
- Source-equivalent process-global cancellation and non-tetrahedral
  configuration construction outside pinned `findConfigs` remain separately
  identified capabilities; no mismatch inside the supported dispatcher is
  classified as unsupported.
- The pinned RDKit 2026.03.1 SMARTS corpus covers 162 source rows: 91 accepted
  parser inputs, nine expected rejects, and 62 ordered matcher inputs. Exact
  parse acceptance, SMARTS output, graph size, and atom mappings match the
  reference; strict Rust, Python, consumer-composition, and architecture tests
  cover parameters, callbacks, recursive queries, fragments, and CXSMARTS.
- Reaction SMARTS and database/container SMARTS are separate source closures
  and remain outside the ordinary-molecule boundary. The Bison diagnostic
  stream requested by `debug_parse` is explicitly unsupported rather than
  silently ignored.
- The complete Layered audit processed all 2,897,819 ChEMBL 37 source records
  over 128 shards. All 2,897,804 mutually parseable records matched pinned
  RDKit across 18 profiles and 52,160,472 exact complete-result comparisons,
  with zero mismatch, invalid profile, failed task, or retained finding.
- The committed focused, 152-row, and 5,000-row Layered gates cover 5,168 input
  rows, 92,808 branch executions, and 278,424 exact output-field comparisons.
  They cover every active layer, path/root mode, mask, seeded-count state,
  query branch, structured error, scalar/batch form, and composition path.
- The complete Pattern audit processed all 2,897,819 ChEMBL 37 records over
  128/128 deterministic shards. All 2,897,804 mutually parseable records
  matched RDKit `2026.03.1` revision
  `351f8f378f8ad6bbd517980c38896e66bf907af8` across ten complete profiles,
  totaling 28,978,040 exact vectors with zero mismatch; 15 records were
  rejected by both parsers.
- Pattern validation covers default and tautomeric output, query-bearing
  graphs, boundary and collision-heavy widths, source-inert count/mask
  arguments, exact errors, screening containment, immutability, ordering,
  concurrency, and interleaving with every implemented fingerprint family.
  RDKit's separate `MolBundle` intersection overload remains outside the
  ordinary molecule and ordered batch model. The behavioral parity claim does
  not include performance parity: representative warm benchmarks retain a
  measured `1.508-3.104x` constant-factor runtime gap.
- The complete ChEMBL 37 tautomer validation processed all 2,897,819 inputs
  across 128 shards and covered 211,539,707 exact output and molecular-state
  observations with zero mismatch. The matrix covers parse behavior plus
  default, V1, one-tautomer, one-transform, retained atom stereo, retained bond
  stereo, retained isotopic hydrogen, and disabled stereo-reassignment
  profiles.
- Every row identified during source alignment was exhaustively replayed after
  correction, and a deterministic 524,288-record regression drawn across all
  128 ChEMBL shards re-executed the complete branch matrix. Focused fixtures
  permanently lock the distinct legacy/modern stereo-center predicates, ring
  stereo cleanup, canonical traversal, and fragment computed-state lifecycle.
- The complete ordinary MolAlign audit processed 2,854,362 eligible ChEMBL 37
  records across 128/128 shards and performed 11,417,207 comparisons of RMSD,
  complete transforms, selected maps, conformer-pair order, changed
  coordinates, and source/reference immutability with zero mismatch or
  retained finding. Fifteen rows were rejected by both parsers and 43,442
  accepted rows exceeded the configured 80-atom boundary.
- The complete potential-stereo and stereoisomer-enumeration audit processed
  all 2,897,819 ChEMBL 37 records over 128/128 deterministic shards. It fully
  compared 2,897,804 rows accepted by both parsers across four potential-stereo
  and four bounded enumeration profiles. The remaining 15 rows were rejected
  by both parsers at parse entry. It completed 89,831,939 exact parse, state,
  output, and source-preservation observations with zero blocking,
  informational, or other mismatch.
- The committed 152-row and 5,000-row profiles contribute 89,047 and 4,680,490
  exact leaf comparisons. Focused tests additionally lock option defaults,
  source-record order, cleanup state, enhanced stereo groups, arbitrary-width
  counts, exact exhaustive and random sequences, uniqueness, lazy errors,
  embedding outcomes, repeated calls, parallel isolation, and mixed hydrogen,
  conformer, descriptor, and serialization calls.
- Exact-output benchmark preflights passed for candidate discovery, count and
  setup, one-configuration finalization, lazy-prefix, exhaustive, bounded
  random, embedding, and parallel profiles. No enumerator-specific complexity
  regression was found; the measured optional-embedding path retains the
  existing distance-geometry composition cost. The parity boundary is RDKit's
  Python enumerator and does not claim the separate newer C++ enumerator or
  atropisomer enumeration.

## [0.2.13] - 2026-08-24

This release completes source-backed AtomPair and Topological Torsion
fingerprints, expands complete structural input into `BioStructure`, and adds
Gemmi-aligned mmCIF serialization through consistent Rust and Python APIs.

### Added

- Added source-backed AtomPair sparse-count, folded-count, sparse-bit, and
  explicit-bit fingerprints in Rust and Python, including 2D and 3D distances,
  chirality, count simulation, selectors, custom atom invariants, complete
  provenance, metadata/JSON restoration, and ordered batch execution through
  one generator core.
- Added source-backed Topological Torsion sparse-count, folded-count,
  sparse-bit, and explicit-bit fingerprints in Rust and Python, including
  modern and legacy APIs, live options and JSON restoration, chirality,
  selectors, custom atom invariants, complete provenance, and ordered batch
  execution through the shared fingerprint generator core.
- Added Gemmi-aligned mmCIF serialization and file writing to the complete
  `BioStructure` model in Rust and Python. The single writer path covers all 34
  source `MmcifOutputGroups`, canonical CIF document construction and quoting,
  structural hierarchy, entities, crystallographic metadata, assemblies,
  secondary structure, connections, cis-peptides, modified residues, TLS, NCS,
  and the other categories represented by the typed model.
- Added public `MmcifOutputGroups` and `MmcifWriteOptions` controls in Rust and
  Python, together with Rust's structured `BioWriteError`; Python maps writer
  state and IO failures to `ValueError` and `OSError`, respectively. Structural
  format conversion now composes explicitly through `BioStructure`, including
  PDB-to-mmCIF workflows, without adding parallel writers to the lossy
  `Protein` projection or the chemical-graph `Molecule` model.

### Changed

- Updated all Rust crates and the Python package to version 0.2.13.
- Reproduced RDKit's logically read-only 2D and resolved-conformer 3D distance
  matrix caching in shared molecule computed state, with topology/coordinate
  invalidation and clone-local synchronization. AtomPair and distance-geometry
  consumers reuse this single matrix implementation.
- Reproduced the source-defined AddHs explicit-valence cache transition through
  the operation contract so molecules returned by `with_hydrogens()` can be
  consumed directly by modern and legacy Topological Torsion APIs.
- Completed the remaining Gemmi-backed PDB, mmCIF, and mmJSON reader paths for
  state representable by `BioStructure`, including hybrid-36 identifiers,
  modified-residue identity, all source connection kinds, entity and sequence
  metadata, crystallographic transforms, assemblies, secondary structure,
  anisotropic displacement values, and indexed chain/residue construction.
- Clarified the structural ownership boundary throughout the public
  documentation and examples: `BioStructure` owns complete modeled structural
  IO, `Protein` is an explicit amino-acid projection, and `Molecule` retains the
  separate RDKit-compatible chemical-graph and PDB-block boundary.
- Renamed the canonical structural-read support constants to
  `BIO_PDB_READ_FEATURE` and `BIO_MMCIF_READ_FEATURE`; the former subset-named
  constants remain deprecated compatibility aliases and are not separate
  registry entries.

### Fixed

- Corrected `Molecule.from_rdkit()` so its default preserves the imported RDKit
  graph while copying hybridization and computing the explicit-valence cache.
  `sanitize=True` remains the explicit full-sanitize path and `sanitize=False`
  remains the raw unprepared path.
- Corrected BioStructure source fidelity across optional sequence identifiers,
  residue-qualified structural addresses, connection and cis-peptide numeric
  precision, TLS residue ranges, space-group transforms, PDB fixed-field
  parsing, and full model/chain/residue/atom metadata propagation.

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
- The complete Topological Torsion audit processed all 2,897,819 ChEMBL 37
  source records over 128 shards. All 2,897,804 mutually parseable records
  matched pinned RDKit across nine profiles, four vector forms, and complete
  requested provenance, totaling 127,503,376 exact comparisons with zero
  mismatch, invalid profile task, or retained finding.
- The maintained 5,000-molecule Topological Torsion gate runs all nine profiles
  as an exhaustive CI regression; focused source, public-API, atom-invariant,
  batch, and Python tests cover live states, edge graphs, errors, deterministic
  parallel execution, and source immutability.
- Pinned Gemmi exact-output tests compare canonical mmCIF bytes and reparsed
  semantic state for complete PDB and mmCIF fixtures. Focused Rust and Python
  integration tests additionally cover non-mutating conversion, hierarchy and
  metadata preservation, multi-model coordinates, output-group suppression,
  file output, structured IO errors, and the public runtime/stub surface.
- The completed Bio reader/writer line passes the release strict workspace test
  gate, Python static and runtime API checks, and the Sphinx documentation
  build. The declared boundary is the state represented by `BioStructure`;
  arbitrary private CIF categories and source-document layout are not claimed
  as byte-preserving round trips.

Full comparison: [v0.2.12...v0.2.13](https://github.com/cosmol-studio/COSMolKit/compare/v0.2.12...v0.2.13)

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
