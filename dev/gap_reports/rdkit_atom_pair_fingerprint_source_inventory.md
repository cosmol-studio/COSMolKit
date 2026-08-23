# RDKit AtomPair Fingerprint Source Inventory

Status: completed pre-port source inventory retained as an audit baseline. This
report fixes the selected call graph, ownership boundaries, source ranges, and
line-coverage obligations as they existed before the AtomPair production port;
its historical disposition columns are intentionally not a current support
matrix. Current implementation status and exact focused, maintained-corpus,
and complete ChEMBL 37 evidence are recorded in
[`rdkit_atom_pair_fingerprint_full_port_validation.md`](rdkit_atom_pair_fingerprint_full_port_validation.md).

## Pinned source and selected boundary

- RDKit source tree: `third_party/rdkit`.
- Pinned release/revision: RDKit `2026.03.1`,
  `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- Normative entry point: `RDKit::AtomPair::getAtomPairGenerator()` in
  `Code/GraphMol/Fingerprints/AtomPairGenerator.cpp:251-283` feeding the common
  `FingerprintGenerator` pipeline.
- Compatibility sources: `Code/GraphMol/Fingerprints/AtomPairs.cpp`,
  `Code/GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp`, and
  `rdkit/Chem/AtomPairs/Pairs.py`. They define defaults and legacy projection
  semantics, but do not authorize a second pair-enumeration or encoding engine.
- Output boundary: sparse count, folded count, sparse bit, explicit bit, count
  simulation, all AtomPair-supported `AdditionalOutput`, custom atom
  invariants, chirality, 2D/3D distances, conformer selection, rooted and
  ignored atoms, metadata/JSON, and ordered batch projection.
- Out of scope: topological-torsion fingerprint generation and the diagnostic
  Python helpers `pyScorePair`, `ExplainPairScore`, `ExplainAtomCode`, and
  `BitsInCommon`. The torsion correction flag in the shared atom-invariant
  generator remains in scope because it is part of the selected C++ type.

## Active source call graph

1. Both `getAtomPairGenerator()` overloads construct `AtomPairArguments`, an
   `AtomPairEnvGenerator`, and the default `AtomPairAtomInvGenerator`, then
   instantiate the common `FingerprintGenerator`
   (`AtomPairGenerator.cpp:251-283`).
2. A selected output method enters
   `FingerprintGenerator::getFingerprintHelper()`, prepares stereochemistry on
   a clone when requested, resets requested provenance containers, obtains atom
   invariants, enumerates environments, projects their identifiers, adds
   deterministic extra bits, and records provenance
   (`FingerprintGenerator.cpp:323-435`).
3. `AtomPairAtomInvGenerator::getAtomInvariants()` visits atoms in index order
   and calls `AtomPairs::getAtomCode()`
   (`AtomPairGenerator.cpp:31-43`; `FingerprintUtil.cpp:45-97`). The latter
   calls `RDKit::numPiElectrons()` (`Atom.cpp:934-953`) and consumes CIP labels
   only when chirality is enabled.
4. `AtomPairEnvGenerator::getEnvironments()` selects the 2D topological or 3D
   Euclidean distance matrix, visits pairs in `i < j` order, applies ignored
   atoms before rooted selection, floors distances, applies inclusive bounds,
   and emits immutable two-atom environments
   (`AtomPairGenerator.cpp:177-237`; `Matrices.cpp:167-236,392-434`).
5. `AtomPairAtomEnv::getBitId()` either packs the canonical exact atom-pair
   code or combines endpoint invariants and distance through three ordered
   Boost hash-combine calls (`AtomPairGenerator.cpp:130-167`;
   `FingerprintUtil.cpp:99-107`).
6. `AtomPairAtomEnv::updateAdditionalOutput()` updates `bitInfoMap`, duplicate-
   preserving `atomToBits`, both endpoint counts, and ordered `atomsPerBit`
   (`AtomPairGenerator.cpp:108-128`). Count simulation duplicates provenance
   only through the common projector (`FingerprintGenerator.cpp:437-500`).
7. The four scalar projections and the shared multi-thread batch helper remain
   common-generator responsibilities (`FingerprintGenerator.cpp:503-756`).
8. Deprecated count and bit-vector functions construct the modern generator;
   the special legacy `nBitsPerEntry` thresholds are a post-count compatibility
   projection, not a second pair loop (`AtomPairs.cpp:51-157`).

## Function and ownership ledger

| Source function or unit | Exact pinned range | Reachability | Current COSMolKit ownership/state | Planned sole owner |
| --- | --- | --- | --- | --- |
| AtomPair bit-width constants and element table | `FingerprintUtil.h:27-41` | Direct dependency of atom/pair code | Missing | `properties/fingerprint/atom_pair.rs` |
| `RDKit::numPiElectrons` | `Atom.cpp:934-953` | Called by `getAtomCode` | No shared source-equivalent helper | Shared atom-property module |
| `AtomPairs::getAtomCode` | `FingerprintUtil.cpp:45-97` | Called for every default atom invariant | Missing | `properties/fingerprint/atom_pair.rs` |
| `AtomPairs::getAtomPairCode` | `FingerprintUtil.cpp:99-107` | Exact sparse-count identifiers | Missing | `properties/fingerprint/atom_pair.rs` |
| `AtomPairAtomInvGenerator` constructor | `AtomPairGenerator.cpp:26-29` | Default and public invariant factories | Missing | `properties/fingerprint/atom_pair.rs` |
| `AtomPairAtomInvGenerator::getAtomInvariants` | `AtomPairGenerator.cpp:31-43` | Default generator preparation | Missing | `properties/fingerprint/atom_pair.rs` |
| invariant `infoString` | `AtomPairGenerator.cpp:45-48` | Metadata/JSON inspection | Missing | `properties/fingerprint/atom_pair.rs` |
| invariant `toJSON` | `AtomPairGenerator.cpp:50-55` | Generator serialization | Missing | `properties/fingerprint/atom_pair.rs` plus common JSON envelope |
| invariant `fromJSON` | `AtomPairGenerator.cpp:56-61` | Generator restoration | Missing | `properties/fingerprint/atom_pair.rs` plus common JSON envelope |
| invariant `clone` | `AtomPairGenerator.cpp:63-66` | Owned factory state | Missing | Rust `Clone` on the sole invariant type |
| `AtomPairEnvGenerator::getResultSize` | `AtomPairGenerator.cpp:68-75` | Sparse-result allocation/range | Missing | `properties/fingerprint/atom_pair.rs` |
| `AtomPairArguments` constructor | `AtomPairGenerator.cpp:77-87` | Every modern factory | Missing | `properties/fingerprint/atom_pair.rs`, layered on common args |
| argument `infoString` | `AtomPairGenerator.cpp:89-93` | Generator metadata | Missing | AtomPair args plus common metadata route |
| argument `toJSON` | `AtomPairGenerator.cpp:94-100` | Generator serialization | Missing | AtomPair args plus common JSON route |
| argument `fromJSON` | `AtomPairGenerator.cpp:101-106` | Generator restoration | Missing | AtomPair args plus common JSON route |
| `AtomPairAtomEnv::updateAdditionalOutput` | `AtomPairGenerator.cpp:108-128` | Every requested provenance output | Missing; common containers/projector exist | AtomPair environment using existing `AdditionalOutput` |
| `AtomPairAtomEnv::getBitId` | `AtomPairGenerator.cpp:130-167` | Every emitted environment | Missing; `hash_combine` exists | `properties/fingerprint/atom_pair.rs` |
| `AtomPairAtomEnv` constructor | `AtomPairGenerator.cpp:169-175` | Every accepted atom pair | Missing | Sole immutable AtomPair environment type |
| `AtomPairEnvGenerator::getEnvironments` | `AtomPairGenerator.cpp:177-237` | Every scalar/batch output | Missing | `properties/fingerprint/atom_pair.rs` using shared matrices |
| environment `infoString` | `AtomPairGenerator.cpp:239-242` | Generator metadata | Missing | AtomPair environment generator |
| environment `toJSON` | `AtomPairGenerator.cpp:244-249` | Generator serialization | Missing | AtomPair environment generator plus common envelope |
| factory taking owned arguments/invariant generator | `AtomPairGenerator.cpp:251-268` | Modern root | Missing | Thin AtomPair factory into shared core |
| parameter factory | `AtomPairGenerator.cpp:270-283` | Public/default construction | Missing | Thin forwarding AtomPair factory |
| explicit C++ template instantiations | `AtomPairGenerator.cpp:285-308` | Link-time C++ artifact only | No Rust behavior to reproduce | Covered by all four concrete Rust output projections |
| common argument constructor/info/JSON | `FingerprintGenerator.cpp:36-87` | Every fingerprint family | Already represented by `FingerprintArguments` | `properties/fingerprint/generator.rs` |
| common generator construction/destruction | `FingerprintGenerator.cpp:89-131` | Every fingerprint family | Lifetime/state currently Morgan-bound | Shared Rust-owned generator core |
| reinitialize `AdditionalOutput` | `FingerprintGenerator.cpp:133-150` | Every fingerprint call | Existing implementation | Shared generator core only |
| generator info/JSON dispatch | `FingerprintGenerator.cpp:161-317` | Metadata/restoration | Morgan-oriented route exists | Shared generator core with family discriminators |
| `getFingerprintHelper` | `FingerprintGenerator.cpp:323-435` | All four scalar projections | Implemented as Morgan methods | `properties/fingerprint/generator.rs` only |
| `duplicateAdditionalOutputBit` | `FingerprintGenerator.cpp:437-479` | Count-simulation provenance | Existing free helper | Shared generator core only |
| `setupTempAdditionalOutput` | `FingerprintGenerator.cpp:481-500` | Count simulation | Existing free helper | Shared generator core only |
| sparse-count projection | `FingerprintGenerator.cpp:503-508` | Modern output | Morgan method | Shared generator core only |
| sparse-bit projection | `FingerprintGenerator.cpp:510-573` | Modern output | Morgan method | Shared generator core only |
| folded-count projection | `FingerprintGenerator.cpp:575-589` | Modern output | Morgan method | Shared generator core only |
| explicit-bit projection | `FingerprintGenerator.cpp:591-652` | Modern output | Morgan method | Shared generator core only |
| MT helper/four batch methods | `FingerprintGenerator.cpp:654-756` | Modern batch outputs | Existing COSMolKit scheduler and Morgan callbacks | Shared scheduler plus family callbacks |
| Floyd-Warshall selected implementation | `Matrices.cpp:39-100` | 2D topological matrix | Private reproduction inside `distgeom.rs:6106-6223` | Shared chemistry matrix module |
| `MolOps::getDistanceMat` | `Matrices.cpp:167-236` | AtomPair `use2D=true` | Private distgeom implementation covers selected unweighted mode | Shared chemistry matrix module |
| `MolOps::get3DDistanceMat` | `Matrices.cpp:392-434` | AtomPair `use2D=false` | Missing; conformer storage exists | Shared chemistry matrix module |
| legacy `getAtomPairFingerprintInternal` | `AtomPairs.cpp:51-78` | Deprecated facade | Missing | Thin adapter to modern generator |
| legacy sparse overloads | `AtomPairs.cpp:81-103` | Deprecated facade | Missing | Thin adapter to modern generator |
| legacy hashed count | `AtomPairs.cpp:105-120` | Deprecated facade | Missing | Thin adapter to modern generator |
| legacy hashed bit vector | `AtomPairs.cpp:122-157` | Deprecated facade | Missing | Thin threshold projection over common count output |
| `updateElement` / `setAtomPairBit` | `AtomPairs.cpp:27-48` | Unreachable from pinned modern path and current deprecated facade | Not implemented | Explicitly excluded; source-history only |
| Python generator/invariant wrappers | `Wrap/AtomPairWrapper.cpp:23-106` | Public Python behavior/defaults | Missing | Project-native Python facade to same Rust generator |
| descriptor wrappers/defaults | `rdMolDescriptors.cpp:213-330,969-1039` | Legacy Python semantics | Missing | Oracle/tests; no duplicate naming surface required |
| Python diagnostic helpers | `Pairs.py:27-150`; `Utils.py:14-109` | Outside production generation graph except legacy bit-vector facade | Missing | Lower-level encoding tests; no production duplicate |

## Existing infrastructure that must remain singular

The current `properties/fingerprint.rs` already owns `AdditionalOutput`,
`FingerprintArguments`, `FingerprintFuncArguments`, `Fingerprint`,
`SparseBitFingerprint`, `SparseCountFingerprint`, `hash_combine`, the RDKit MT
random stream, count-simulation provenance duplication, and temporary
provenance setup. These facilities are common source behavior and must be moved
or reused, never copied into the AtomPair module.

The common scalar pipeline was extracted from `MorganFingerprintGenerator`
into the family-neutral `properties/fingerprint/generator.rs` core. Morgan is
now a thin family adapter, and AtomPair can become the second family without
duplicating projection logic.

The caller audit used exact symbol searches across `crates/cosmolkit-core/src`.
Before removal, `build_fingerprint()` and `compute_initial_invariants()` each
had only their own definition; `fold_invariant()` and `atom_is_excluded()` were
called only by `build_fingerprint()`, and `morgan_bond_invariant()` likewise had
only its definition. Those five functions therefore formed an unreachable
alternative Morgan engine. They were removed together. Shared helpers with
independent callers, including `compute_feature_invariants()`,
`compute_connectivity_invariants()`, and
`morgan_additional_output_from_rdkit_output()`, were retained.

The only 2D topological matrix reproduction is currently private to
`chemistry/distgeom.rs`. It must be relocated rather than copied. The existing
`Molecule::conformers_3d()`/`Conformer3D::coordinates()` storage is sufficient
for the source 3D computation; only the shared matrix behavior is missing.
Existing CIP labeling and batch scheduling are also dependencies, not
AtomPair-owned algorithms.

## Exact behavioral boundaries discovered by the audit

- The atom-type table contains 16 entries but only 15 explicit atomic numbers:
  `5,6,7,8,9,14,15,16,17,33,34,35,51,52,53`; the implicit final zero and
  fallback index are source-significant.
- Widths are exactly: type 4, pi 2, branch 3, chirality 2, non-chiral atom code
  9, path 5, and unfolded non-chiral pair fingerprint 23 bits.
- `numPiElectrons()` returns one for aromatic atoms, zero for SP3 atoms, and
  otherwise explicit valence minus physical bonds, where physical bonds include
  explicit hydrogens and only incident bonds with nonzero valence contribution.
- Custom invariant reduction uses modulo `(1 << atom_code_width) - 1` (511 or
  2047), not a bit mask or modulo the full code-space size.
- Exact pair encoding requires distance below 31. Folded/hash output can hash a
  larger caller-selected maximum distance because the C++ argument constructor
  validates only `minDistance <= maxDistance`.
- Pair iteration is `i < j`. Ignored endpoints are rejected before rooted-pair
  inclusion; a rooted pair is retained when either endpoint is selected.
- 3D distances are converted with `floor()`, not rounding. Minimum and maximum
  distance comparisons are inclusive.
- Default AtomPair options are count simulation enabled, bounds `[1,2,4,8]`,
  one bit per feature, fingerprint size 2048, minimum distance 1, maximum
  distance 30, chirality disabled, and 2D mode enabled.
- `AdditionalOutput` supports all four existing containers for AtomPair,
  including `atomsPerBit`; endpoint order and duplicate insertion are observable.
- The legacy bit-vector facade uses thresholds `[1,2,4,8]` when
  `nBitsPerEntry == 4`; for other values it sets slot `i` only when count is
  strictly greater than `i`.
- Public Rust validation must require a custom-invariant entry for every atom.
  The C++ endpoint precondition is weaker than its subsequent indexing, so
  reproducing successful behavior requires a structured bounds error rather
  than reproducing undefined out-of-bounds access.
- Chirality requires the common helper's clone-and-assign stereo branch. The
  current Rust helper returns unsupported when `_StereochemDone` is absent;
  that is a dependency gap, not an AtomPair-specific fallback opportunity.

## Source-line coverage ledger

Legend: `inventory` means the source behavior and target are fixed here but no
port credit is claimed; `reuse/refactor` means a source-backed implementation
exists but ownership or family generality must still be corrected; `excluded`
means reachability was proven absent from the selected graph.

| Source range | Behavior | Target artifact | Audit status |
| --- | --- | --- | --- |
| `FingerprintUtil.h:27-41` | widths, limits, atom table | AtomPair constants | inventory |
| `Atom.cpp:934-953` | pi-electron count | shared atom helper | inventory |
| `FingerprintUtil.cpp:45-107` | atom and pair codes | AtomPair encoding | inventory |
| `AtomPairGenerator.cpp:26-66` | invariant generator state/behavior/JSON | AtomPair invariant generator | inventory |
| `AtomPairGenerator.cpp:68-106` | result size and arguments | AtomPair argument/environment metadata | inventory |
| `AtomPairGenerator.cpp:108-175` | provenance, bit ID, environment state | AtomPair environment | inventory |
| `AtomPairGenerator.cpp:177-249` | ordered enumeration and environment metadata | AtomPair environment generator | inventory |
| `AtomPairGenerator.cpp:251-283` | both factories | thin shared-core constructors | inventory |
| `AtomPairGenerator.cpp:285-308` | concrete C++ instantiations | four output tests | inventory |
| `FingerprintGenerator.cpp:36-150` | common args/lifetime/output reset | shared generator | reuse/refactor |
| `FingerprintGenerator.cpp:161-317` | metadata/JSON dispatch | shared generator | reuse/refactor |
| `FingerprintGenerator.cpp:323-435` | environment-to-count pipeline | shared generator | reuse/refactor |
| `FingerprintGenerator.cpp:437-500` | count-simulation provenance | shared generator | reuse/refactor |
| `FingerprintGenerator.cpp:503-652` | four scalar outputs | shared generator | reuse/refactor |
| `FingerprintGenerator.cpp:654-756` | ordered batch outputs | shared scheduler/adapters | reuse/refactor |
| `Matrices.cpp:39-100,167-236` | unweighted 2D distance matrix | shared matrix module | reuse/refactor |
| `Matrices.cpp:392-434` | conformer 3D distance matrix | shared matrix module | inventory |
| `AtomPairs.cpp:51-157` | legacy adapters and thresholds | compatibility delegates | inventory |
| `AtomPairs.cpp:27-48` | dead historical bit updater | none | excluded |
| `AtomPairWrapper.cpp:23-106` | modern Python defaults/factories | project-native Python facade | inventory |
| `rdMolDescriptors.cpp:213-330,969-1039` | legacy Python defaults/errors | adapters/tests | inventory |
| `Pairs.py:27-150`; `Utils.py:14-109` | diagnostics/legacy presentation | encoding and facade tests | excluded from production implementation |

Every implemented source function must later carry the project's line-level
two-axis source markers. Moving existing shared code moves its authoritative
source frame; it must not leave a second marked copy behind.

## Required source tests and oracle evidence

The focused source oracle includes AtomPair sections in
`testFingerprintGenerators.cpp`, `catch_tests.cpp`,
`fpgen_catch_tests.cpp`, `Wrap/testGenerators.py`, and
`Descriptors/Wrap/testMolDescriptors.py`. Their cases cover defaults, old/new
API equality, chiral atom codes, dative-bond pi behavior, rooted/ignored atoms,
custom invariants, all output families, `AdditionalOutput`, JSON, and bulk
generation. COSMolKit evidence must additionally exercise checked Rust error
boundaries for invalid atom indices, invariant lengths, missing conformers,
coordinate-count mismatch, and non-finite coordinates.

Completion requires exact focused fixtures followed by every active small and
5,000-row profile and every valid ChEMBL 37 row for which the source operation
is defined. No mismatch allowlist, row filtering, similarity-only comparison,
or corpus-conditioned production branch can close this inventory.
