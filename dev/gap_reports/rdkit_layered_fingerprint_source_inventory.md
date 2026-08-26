# RDKit Layered Fingerprint Source Inventory

## Scope And Provenance

This inventory fixes the source boundary for the complete COSMolKit port of
RDKit's Layered fingerprint. The normative source is RDKit `2026.03.1` at
revision `351f8f378f8ad6bbd517980c38896e66bf907af8`, vendored at
`third_party/rdkit/`. Line numbers below refer to that revision.

The port boundary is `RDKit::LayeredFingerprintMol()` and every reachable
helper that affects its observable output, option semantics, or failure
behavior. Pattern Fingerprint is a separate algorithm and is not part of this
port. Existing COSMolKit fingerprint infrastructure is a dependency to reuse,
not permission to preserve an older or parallel Layered implementation.

## Normative Call Graph

```text
Python/Rust public Layered API
  -> wrapper width and sequence conversion
  -> LayeredFingerprintMol
     -> preconditions
     -> MolOps::findSSSR when RingInfo is uninitialized
     -> QueryOps classifiers
        -> isComplexQuery(Bond)
        -> isComplexQuery(Atom)
           -> _complexQueryHelper
        -> isAtomAromatic
        -> queryIsBondInRing
        -> queryBondMinRingSize
     -> path selection
        -> findAllSubgraphsOfLengthsMtoN (branched)
        -> findAllPathsOfLengthsMtoN (linear)
     -> six active layer encoders
     -> sorted layer payload plus path-atom count and layer number
     -> gboost::hash_range<uint32_t>
        -> gboost::hash_combine<uint32_t>
     -> ExplicitBitVect::setBit
     -> optional setOnlyBits gate and atomCounts update
```

## Source-Line Ledger

| Source unit | Exact pinned range | Behavior owned by the range |
|---|---:|---|
| `LayeredFingerprintMol` declaration and documentation | `Code/GraphMol/Fingerprints/Fingerprints.h:73-114` | Parameters, defaults, six documented layers, ten declared layer slots, version `0.7.0`, and substructure mask `0x07`. |
| Complete `LayeredFingerprintMol` body | `Code/GraphMol/Fingerprints/Fingerprints.cpp:252-518` | Preconditions, ring preparation, caches, enumeration, all encoders, projection, masking, counts, and result. |
| Preconditions and ring preparation | `Fingerprints.cpp:257-267` | Reject zero `minPath`, reversed range, zero width, short counts, and width-mismatched mask; call exact SSSR only when ring state is uninitialized. |
| Query/bond cache | `Fingerprints.cpp:269-288` | Cache bonds by bond index and form the three-bit path query mask. |
| Atom caches | `Fingerprints.cpp:290-301` | Cache query-aware aromaticity and atomic number by atom index. |
| Path selection and root aggregation | `Fingerprints.cpp:303-328` | Branched/linear choice, null versus present roots, duplicates, root order, and prepend aggregation. |
| Path state and neighbor counts | `Fingerprints.cpp:330-388` | Atom/bond membership, path-local endpoint degree, pairwise bond adjacency, and query-mask union. |
| Pure topology layer `0x01` | `Fingerprints.cpp:389-401` | Endpoint degree canonicalization, three-bit modulo fields, and packing. |
| Bond-order layer `0x02` | `Fingerprints.cpp:402-424` | Complex-bond suppression, aromatic/single normalization, degree canonicalization, and packing. |
| Atom-type layer `0x04` | `Fingerprints.cpp:425-447` | Complex-endpoint suppression, atomic-number modulo, endpoint canonicalization, and packing. |
| Ring-presence layer `0x08` | `Fingerprints.cpp:448-453` | Complex-endpoint suppression and sparse emission only for a ring bond. |
| Minimum-ring-size layer `0x10` | `Fingerprints.cpp:454-458` | Complex-endpoint suppression and minimum ring size modulo eight, including zero. |
| Aromaticity layer `0x20` | `Fingerprints.cpp:459-473` | Complex-endpoint suppression, query-aware endpoint aromaticity ordering, neighbor count, and sparse packing. |
| Layer projection | `Fingerprints.cpp:475-514` | Empty-layer omission, sorting, distinct-atom suffix, one-based layer suffix, hash, modulo, mask, bit insertion, and once-per-accepted-path counts. |
| Bond query complexity | `Code/GraphMol/QueryOps.cpp:728-766` | No-query, negation, simple order, Single-or-Aromatic, And/Xor, and exact two-child Or classification. |
| Recursive atom query helper | `QueryOps.cpp:769-795` | Negation, atomic-number/type discovery, Or/Xor rejection, and recursive And traversal with shared `hasAtNum`. |
| Atom query complexity | `QueryOps.cpp:892-919` | No-query, simple queries, compound queries, and missing atomic-number behavior. |
| Query-aware atom aromaticity | `QueryOps.cpp:920-959` | Concrete atoms and AtomicNum, IsAromatic, IsAliphatic, AtomType, and ordered AtomAnd query branches. |
| Ring query accessors | `Code/GraphMol/QueryOps.h:303-311` | Boolean bond-in-ring and minimum bond-ring-size access through `RingInfo`. |
| Exact SSSR entry | `Code/GraphMol/FindRings.cpp:769-1003` | RingInfo reset/initialization, eligible-bond graph, SSSR construction, and ring-state population. |
| Branched subgraph range | `Code/GraphMol/Subgraphs/Subgraphs.cpp:347-401` | Root filtering, bond-index iteration, forbidden state, recursion, and per-length result ordering. |
| Linear path range | `Subgraphs.cpp:443-549` | Atom adjacency, bond-length conversion, rooted traversal, bond-path conversion, and duplicate elimination. |
| Hash result width | `Code/RDGeneral/hash/hash_fwd.hpp:17-20` | `std::hash_result_t` is pinned to `std::uint32_t`. |
| `hash_combine` and `hash_range` | `Code/RDGeneral/hash/hash.hpp:209-230` | Ordered seed-zero 32-bit wrapping hash accumulation. |
| Bit insertion | `Code/DataStructs/ExplicitBitVect.cpp:84-94` | In-range bit setting and idempotent repeated writes. |
| Python wrapper | `Code/GraphMol/Wrap/MolOps.cpp:633-665` | Root conversion, seeded mutable counts, short-count error, core call, and copy-back. |
| Python registration | `MolOps.cpp:2638-2651` | Source-facing defaults and exported version/substructure-layer constants. |

## Exact Defaults And Public Metadata

The source defaults are `layerFlags=0xffffffff`, `minPath=1`, `maxPath=7`,
`fpSize=2048`, no atom-count vector, no bit mask, `branchedPaths=true`, and no
root selection. `maxFingerprintLayers` is ten, but only layer indices zero
through five have encoders. Consequently flags above `0x20` are accepted and
produce no components; they are not an unsupported error and must not be
invented as additional features.

The upstream algorithm version is `0.7.0`. `substructLayers=0x07` is metadata
for the topology, bond-order, and atom-type prefix. The upstream header labels
the API experimental, so COSMolKit documentation must preserve that source
metadata without presenting the Rust API as unstable by accident.

## Exact Preconditions And Errors

The C++ core uses preconditions with these source messages:

| Condition | Source message |
|---|---|
| `minPath == 0` | `minPath==0` |
| `maxPath < minPath` | `maxPath<minPath` |
| `fpSize == 0` | `fpSize==0` |
| `atomCounts.size() < mol.getNumAtoms()` | `bad atomCounts size` |
| `setOnlyBits.getNumBits() != fpSize` | `bad setOnlyBits size` |

The Python wrapper checks the count-vector length before entering the core and
raises `ValueError("atomCounts shorter than the number of atoms")`. Root
conversion uses `pythonObjectToVect(fromAtoms, mol.getNumAtoms())`; its index
validation is part of the wrapper boundary and must be represented by the
project-native typed error surface. Rust must not reproduce C++ assertion
termination, but it must preserve each rejection boundary as a structured
error and must not clamp, filter, or silently repair invalid options.

## Enumeration And Ordering Semantics

- A null root pointer enumerates the entire molecule. A present empty root
  vector enumerates no paths. These states are observably different.
- Branched mode calls `findAllSubgraphsOfLengthsMtoN(..., useHs=false)`;
  linear unrooted mode calls `findAllPathsOfLengthsMtoN` with the header
  defaults equivalent to bond paths, no hydrogens, no root, and no
  shortest-path restriction.
- Rooted branched mode passes each input root to the subgraph enumerator.
  Rooted linear mode explicitly requests bond paths with `useHs=false` and
  the selected root.
- Each root's per-length path list is inserted at the beginning of the
  accumulated list. Later root groups therefore precede earlier groups.
  Duplicate and reordered roots remain observable and are not deduplicated by
  `LayeredFingerprintMol`.
- Linear path duplicate elimination is performed by the pinned subgraph
  helper on bond membership. Layered orchestration must not add another
  uniqueness or sorting pass.

## Cache And Encoder Semantics

The query mask for a bond is bit `0x01` for a complex bond query, `0x02` for a
complex begin-atom query, and `0x04` for a complex end-atom query. The mask is
ORed over the complete path. Layer `0x02` is suppressed by any complex bond;
layers `0x04`, `0x08`, `0x10`, and `0x20` are suppressed by a complex endpoint;
layer `0x01` is never query-suppressed.

Layer payloads are multisets: their per-bond components are sorted, not
deduplicated. The distinct atom count and one-based layer number are appended
after sorting. Hashing is ordered `uint32_t` gboost hashing even though the
result is assigned to C++ `unsigned long`; native pointer width does not alter
the seed. The bit id is `seed % fpSize`.

`setOnlyBits` gates every projected bit. Atom counts are incremented for every
atom in a path at most once for that path, and only when at least one active
layer produces a bit accepted by the mask. Counts are seeded/incremental, not
cleared on entry. A repeated bit collision still represents an accepted path
and therefore still affects counts.

## Upstream Tests To Reproduce

| Upstream test | Exact range | Covered source behavior |
|---|---:|---|
| `test1Layers` | `Code/GraphMol/Fingerprints/test1.cpp:389-602` | Determinism and individual/cumulative layer distinctions across atom, bond, ring, and aromatic examples. |
| `test2Layers` | `test1.cpp:604-672` | Seeded atom counts, repeat accumulation, bit masks, and mask-dependent counts. |
| `test3Layers` | `test1.cpp:674-865` | SMARTS/query suppression and substructure-layer behavior. |
| rooted fingerprint regression block | `test1.cpp:2182-2299`, Layered cases `2260-2297` | Single/multiple roots and branched/linear rooted results. |
| `test55LayeredFingerprint` | `Code/GraphMol/Wrap/rough_test.py:2685-2715` | Python defaults, counts copy-back, seeded repeats, empty/full masks, and mask-dependent counts. |

These tests are necessary source fixtures but are not sufficient closure. The
COSMolKit matrix must additionally cover explicit empty roots, duplicate and
reordered roots, inactive high flags, width/range failures, non-power-of-two
widths, query-mask intermediates, hash overflow, path collisions, all mask
shapes, and exact complete-corpus output/count vectors.

## Dependency And Reuse Ownership

The port must have one owner for each dependency:

- typed SMARTS/query representation and evaluation remain owned by
  `search/query.rs`; Layered receives narrow crate-private source-backed
  classifiers rather than a second query tree;
- ring perception remains owned by the existing exact SSSR implementation;
  Layered may trigger it at the source boundary but must not implement a local
  cycle algorithm;
- linear and branched enumeration are shared with RDKFingerprint through one
  crate-private source-backed interface, with Layered-specific argument and
  aggregation behavior expressed explicitly;
- 32-bit gboost hash combine/range have a single implementation next to the
  current fingerprint hashing infrastructure;
- `Fingerprint` remains the sole public bit-vector value and owns bit
  insertion and similarity behavior;
- scalar, batch, Rust, and Python surfaces delegate to one Layered core;
- oracle generation and ChEMBL execution extend the repository-owned
  testdata/parity workflows instead of adding an independent corpus runner.

No operation-registry entry is required: Layered fingerprinting is read-only
and does not mutate molecule state. If SSSR preparation requires internal
derived state, it must use the existing read-operation/cache lifecycle rather
than presenting Layered as a public topology mutation.

## Required Artifacts

Completion requires source-marked Rust helpers and one Layered core; focused
intermediate fixtures; exact small and 5,000-row goldens; deterministic full
ChEMBL 37 run artifacts and identities; Rust and Python scalar/batch tests;
benchmarks and allocation/complexity evidence; architecture guards against
duplicate implementations; generated stubs; public documentation/examples;
feature/support metadata; parity-scope evidence; and a validation closure
report. Machine-local corpora, binaries, caches, and benchmark output remain
outside tracked source unless the repository organization policy explicitly
designates them as committed expected data.

## Current COSMolKit Reuse And Duplication Gap Audit

This section records the pre-port repository audit. It distinguishes reusable
source-backed infrastructure from behavior that must be completed for
Layered. A reusable helper is not considered Layered coverage until the exact
Layered argument combination and observable ordering have dedicated tests.

| Area | Current owner/evidence | Reuse decision | Concrete gap or required consolidation |
|---|---|---|---|
| Public bit vector | `properties/fingerprint.rs::Fingerprint` | Reuse as the only explicit bit-vector result. | Add a crate-private checked bit insertion path or construct from the exact on-bit set; do not add a Layered vector type. Test collisions and width preservation. |
| Fingerprint error | `properties/fingerprint.rs::FingerprintError` | Reuse the single structured fingerprint error surface. | Map every Layered precondition/root/ring failure precisely; do not clamp options or create string-only Python behavior. |
| 32-bit hash combine | `properties/fingerprint.rs::hash_combine` | Reuse unchanged as the sole gboost combine. | Its source body is currently summarized instead of copied line-for-line. Bring it into protocol compliance when hash range is consolidated. |
| 32-bit hash range | `properties/fingerprint.rs::rdkit_fp_hash_range` | Move/consolidate beside `hash_combine` and reuse from RDKFingerprint and Layered. | It is private, RDKFingerprint-named, and its copied source signature incorrectly says `std::size_t`; pinned RDKit uses `std::hash_result_t == uint32_t`. A second Layered hash implementation is forbidden. |
| Linear path enumeration | private `find_all_paths_of_lengths_m_to_n` plus `enumerate_rdkit_fp_paths` in `properties/fingerprint.rs` | Reuse through one crate-private shared interface. | Exact rooted Layered arguments need tests. Rooted prepend order is present, but explicit empty roots must remain distinct from no roots. |
| Branched subgraph enumeration | RDKFingerprint helpers in `properties/fingerprint.rs` | Reuse through the same shared interface. | `enumerate_rdkit_fp_paths` currently uses `roots = from_atoms.unwrap_or(&[])` and treats `roots.is_empty()` as unrooted in the branched branch. This is a confirmed semantic gap: `Some([])` must yield no paths. Root validation, duplicates, reordered roots, and group-prepend ordering need exact tests. |
| Query tree/model | `search/query.rs::{QueryNode, AtomQueryPredicate, BondQueryPredicate}` and SMARTS parser | Reuse as the sole typed query representation. | No Layered-local query model may be introduced. |
| Recursive atom complexity | private `complex_atom_query_helper` in `search/query.rs` | Reuse after narrow crate-private exposure. | Existing source anchors and complexity review are strong, but typed compressed list-query branches and every sibling-carried `hasAtNum` case need Layered-focused verification. |
| Atom query complexity | private `is_complex_atom_query` in `search/query.rs` | Reuse after narrow crate-private exposure. | Visibility prevents Layered use. Add exact source-derived local tests rather than copying the classifier. |
| Bond query complexity | private `is_complex_bond_query` in `search/query.rs` | Reuse after narrow crate-private exposure. | Visibility prevents Layered use. Complete tests must distinguish simple two-order Or from every complex Or/And/Xor/negation shape. |
| Query-aware aromaticity | no reusable `is_atom_aromatic` equivalent was found | Add one narrow typed helper in `search/query.rs`. | This is the only missing QueryOps behavior in the Layered call graph. It must reproduce the exact ordered `AtomAnd` branch; ordinary `Atom::is_aromatic` is not a valid substitute for query atoms. |
| Ring accessors | private `query_is_bond_in_ring` and `query_bond_min_ring_size` in `search/query.rs` | Reuse after narrow crate-private exposure. | Visibility prevents Layered use; add direct initialized/uninitialized and fused-ring accessor tests. |
| SSSR | `chemistry/rings.rs::find_sssr` and `RingInfo` | Reuse as the sole ring engine. | Layered must select an already initialized RingInfo when present or compute exact SSSR locally when absent. The public fingerprint remains read-only; no operation-registry mutation or fast-ring substitution is allowed. |
| Fingerprint API style | `TopologicalFingerprintParams`, typed output request/result, free functions, and `Molecule` methods | Follow this project-native pattern. | Add typed Layered flags/params/output/result and method facade. Do not expose RDKit's mutable pointer API or create a generator abstraction not present in the source boundary. |
| Rust batch engine | `properties/batch.rs::collect_optional_values_with_options` and ordered Rayon collection | Reuse for scalar-core delegation, stable ordering, errors, job counts, and progress. | There is no Layered batch surface. Add one without copying the core or parallelizing internal path enumeration per molecule. |
| Python scalar facade | PyO3 fingerprint methods in `python/src/lib.rs` | Follow the existing project naming, typed `Fingerprint`, and `fingerprint_pyerr`. | No Layered method/result types exist. Preserve absent versus explicitly empty roots; unlike the older helper `normalize_fingerprint_indices`, an empty Layered root list must not be normalized to `None`. |
| Python ordered batch facade | existing Morgan/AtomPair batch methods | Reuse the Rust ordered batch path. | No Layered batch method exists. Do not implement a Python loop or a second scalar algorithm. |
| Support metadata | `support.rs` fingerprint `FeatureSpec` values and `PUBLIC_FEATURES` | Add one separate Layered feature specification after validation. | Layered is absent. It must not be folded into Morgan, AtomPair, Topological Torsion, Avalon, or the currently broad Morgan feature text. |
| Golden schema/runner | `tools/testdata/rdkit/generate_all.py`, `_expected_schema.py`, and fingerprint generators/profiles | Extend the existing generator registry and committed expected-data layout. | No Layered profile, generator, schema, focused fixture, or expected file exists. Use the pinned reference identity and existing small/5,000 profiles. |
| Focused and maintained parity | fingerprint integration targets under `crates/cosmolkit-core/tests/` | Add Layered-specific exact tests. | No Layered tests exist. Existing RDKFingerprint tests prove dependencies only, not Layered encoders, masks, counts, or query suppression. |
| ChEMBL 37 workflow | `dev/tools/chembl_parity/audit_fingerprints.py` plus checksummed/resumable `run.py` | Add a Layered mode/profile to this repository-owned workflow. | Layered is absent from `KNOWN_SCRIPTS` modes and comparison functions. Preserve exact attempted/parsed/compared/error/mismatch counts and external artifacts. |
| Documentation/stubs/examples | public README, Python docs/examples, generated `python/cosmolkit.pyi` | Extend only after the validated API is fixed. | There is no Layered surface to document. Generated stubs must come from `stub_gen`, never hand edits. |
| Architecture guards | existing operation/fingerprint architecture tests | Add Layered ownership assertions. | There is currently no guard preventing a duplicate Layered core, local query classifier, local hash range, or separate scalar/batch implementation. |

### Existing Facilities That Are Not Reuse Proofs

Topological/RDKFingerprint shares enumeration and hash primitives, but its
feature encoding, multiple-bit projection, optional folding, `useHs`, bond
order, and provenance semantics are different. It is not a Layered
implementation. Morgan, MACCS, AtomPair, Topological Torsion, Avalon, and
Pattern Fingerprint likewise own different algorithms and must not be called
as a fallback or used to infer missing Layered results.

The repository currently has no symbol, public method, Python method, feature
entry, golden schema, profile, parity test, or ChEMBL comparison named for the
Layered fingerprint. There is therefore no historical Layered core to retain.
The implementation mainline begins with shared-dependency consolidation and
then adds exactly one source-backed Layered core.

### Confirmed Pre-Port Defects Or Risks

1. Branched `enumerate_rdkit_fp_paths` conflates absent and present-empty root
   selections. This must be fixed at the shared enumerator boundary before the
   Layered orchestration is added.
2. `rdkit_fp_hash_range` embeds the right wrapping computation but cites the
   wrong pinned result type and is unnecessarily tied to RDKFingerprint by
   name and visibility.
3. Query complexity/ring helpers already exist with source anchors but are
   private; copying them into the Layered module would create two cores.
4. Query-aware `isAtomAromatic` has no shared implementation and cannot be
   replaced by the concrete aromatic flag on query atoms.
5. Existing Python index normalization deliberately collapses empty lists for
   other fingerprint generators. Layered must bypass that normalization
   because the pinned API observes an explicitly empty `fromAtoms` selection.
6. No existing batch method covers RDKFingerprint/Avalon, so Layered batch
   work must deliberately reuse the generic ordered collector rather than
   assuming a matching topological batch facade exists.
