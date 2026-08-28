# RDKit Python EnumerateStereoisomers Source Audit
## Reference Boundary

This audit fixes the source boundary for the COSMolKit stereoisomer
enumeration port. The pinned reference is:

| Field | Value |
|---|---|
| RDKit version | `2026.03.1` |
| RDKit revision | `351f8f378f8ad6bbd517980c38896e66bf907af8c` |
| Python runtime | CPython `3.13.12` |
| Public parity entry | `rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers` |
| Count entry | `rdkit.Chem.EnumerateStereoisomers.GetStereoisomerCount` |
| Potential-stereo entry | `Chem.FindPotentialStereo` |

The public contract is the Python implementation in
`third_party/rdkit/rdkit/Chem/EnumerateStereoisomers.py`. The newer C++
`StereoisomerEnumerator` is a different API and is not a substitute for that
contract. In particular, it has different option defaults, cleanup behavior,
random generation, uniqueness behavior, conformer handling, and atropisomer
coverage.

## Source Call Graph

The active Python call graph is:

```text
GetStereoisomerCount(molecule, options)
  clone molecule
  _getFlippers(clone, options)
    Chem.FindPotentialStereo(clone, cleanIt=false, flagPossible=true)
    _AtomFlipper / _BondFlipper / _StereoGroupFlipper
  return 2 ** flipper_count

EnumerateStereoisomers(molecule, options, verbose)
  clone molecule
  clear atom _CIPCode properties
  clear EITHERDOUBLE and UNKNOWN bond directions
  _getFlippers(workspace, options)
  choose _RangeBitsGenerator or _UniqueRandomBitsGenerator
  apply bits to flippers in source order, least-significant bit first
  remove stereo groups from each output
  SetDoubleBondNeighborDirections
  ClearComputedProps(includeRings=false)
  AssignStereochemistry(cleanIt=true, force=true,
                        flagPossibleStereoCenters=true)
  optionally deduplicate by canonical isomeric SMILES
  optionally AddHs, EmbedMolecule, and copy heavy-atom coordinates
  yield lazily until maxIsomers successful outputs or source exhaustion
```

`Chem.FindPotentialStereo` binds the mutable overload of
`Chirality::findPotentialStereo(mol, cleanIt=false, flagPossible=true)` in
`ChiralityOps.cpp`. Its implementation in `FindStereo.cpp` has the following
single pipeline:

```text
symmetrizeSSSR when required
updatePropertyCache(false) when required
runCleanup
  initAtomInfo
    isAtomPotentialTetrahedralCenter
    isAtomPotentialNontetrahedralCenter
    getStereoInfo(atom)
  initBondInfo
    isBondPotentialStereoBond
    getStereoInfo(bond)
  flagRingStereo
  rankFragmentAtoms
  updateAtoms and updateBonds until stable
  cleanMolStereo when cleanIt=true
  restore possible candidates and repeat when flagPossible requires it
  findChiralAtomSpecialCases
  store _ChiralAtomRank on every atom
store _potentialStereo on the molecule
```

The transitive behavioral dependencies include symmetrized ring order,
canonical fragment ranks, `countSwapsToInterconvert`, protium-neighbor tests,
bridgehead queries, conjugated-bond tests, atom property-cache semantics,
unknown-stereo properties, non-tetrahedral chiral permutations, stereo group
state, and ring-transmitted atom/bond distinctions. These are source behavior,
not optional post-processing.

## Exact Python Semantics

The source defaults are:

| Option | Default |
|---|---|
| `tryEmbedding` | `false` |
| `onlyUnassigned` | `true` |
| `onlyStereoGroups` | `false` |
| `maxIsomers` | `1024` |
| `rand` | `None` |
| `unique` | `true` |

When no flipper exists, enumeration yields one cloned molecule and does not
set `_MolFileChiralFlag`. Otherwise the workspace receives a chiral flag of
`1`. Exhaustive generation visits integer bit patterns from zero upward.
Random generation uses Python `random.Random.getrandbits(nCenters)`, retains a
set of configurations already drawn, and terminates after every possible bit
pattern has been seen. The default random seed is the Python hash of the
sorted `(degree, atomic_number)` tuple, which is invariant to atom order.

`maxIsomers` limits successfully yielded results, not attempted
configurations. Embedding failures therefore do not consume the output limit.
The finite unique-random generator is required to prevent an infinite loop
when embedding rejects configurations. `tryEmbedding` is an upstream-defined
heuristic and must be reproduced exactly; no additional COSMolKit feasibility
filter belongs in this path.

Atropisomer enumeration is absent from the Python `_getFlippers` body. Existing
atropisomer state remains part of molecule preservation and potential-stereo
records, but the Python-parity enumerator must not add an atropisomer flipper.

## Pre-Port COSMolKit Baseline

`crates/cosmolkit-core/src/chemistry/stereo_enumerate.rs` is publicly reachable
as `cosmolkit_core::stereo_enumerate`. Before this port it exported:

- `EnumerationError`, `EnumerationStrategy`, and `EnumerationParams`;
- separate tetrahedral, double-bond, and combined eager enumeration functions;
- three separate count functions.

There was no project-native `Molecule` method, Python binding, typed
potential-stereo public result, lazy iterator, or focused RDKit oracle profile.
Coverage consisted of local unit tests with broad assertions such as non-empty
output or an upper bound on result count.

Those public functions cloned their input and therefore did not visibly mutate
the source molecule. Their internal output construction directly edited the
private topology block. Enumeration is a multi-result derivation rather than
an in-place public mutation, so the replacement does not require a separate
registered mutation operation. The completed implementation confines mutable
state to one owned private workspace, exposes no public mutable storage, and
runs the same molecule invariant checks used by project-native derived outputs.

## Pre-Port Gaps Closed By This Work

The following findings describe the implementation replaced by this port. The
current closure and exact evidence are recorded in
[`rdkit_python_enumerate_stereoisomers_full_port_validation.md`](./rdkit_python_enumerate_stereoisomers_full_port_validation.md).

### Wrong mainline and options

The module comments and source blocks primarily claimed the newer C++
enumerator while its public behavior resembled neither that API nor the Python
API. The Rust defaults used `max_isomers = 0` instead of `1024`; the options
omitted `tryEmbedding`, `onlyUnassigned`, and `onlyStereoGroups`; and they introduced
`Default`, `Random`, `Symmetry`, and `max_tries` controls that do not exist in
the Python contract. Several C++ source blocks were marked behaviorally exact
despite these observable differences.

### Candidate discovery

`find_tetrahedral_centers` combined CIP ranks with a broad potential-center
helper. `find_stereo_bonds` and `is_valid_double_bond_config` independently
guessed E/Z eligibility and small-ring feasibility. They did not reproduce
`StereoInfo`, iterative symmetry refinement, ring transmission, unknown state,
cleanup, controlling-atom order, missing controlling atoms, descriptors,
permutations, or specified state. These alternate candidate detectors were
removed.

### Flippers and enhanced stereo

The replaced implementation had no typed flipper sequence or enhanced stereo
group flipper. It did not implement source ordering, source skip rules,
relative group inversion, `onlyUnassigned`, `onlyStereoGroups`, missing
double-bond controller handling, or the Python decision not to enumerate
atropisomers.

### Configuration generation

`generate_combinations` eagerly allocated every Boolean vector. Public paths
rejected more than 20 atom or bond centers and used machine-word shifts for
counts. The Python generator is lazy and intentionally supports much larger
center counts when `maxIsomers` bounds sampling. The former xorshift generator,
time-dependent default seed, SMILES-derived hash, retry cap, and per-bit random
draws were all observably different from Python `random.Random.getrandbits`.

### Output finalization

The replaced builders set atom or bond state but omitted the Python
preprocessing and finalization sequence. Missing behavior included CIP
clearing, bond-direction cleanup, `_MolFileChiralFlag`, stereo-group removal, double-bond
neighbor directions, computed-property clearing with ring preservation, forced
stereochemistry assignment, canonical-isomeric uniqueness timing, and the
precise output property/cache state.

### Embedding

`tryEmbedding` was not implemented. The required source composition includes
hydrogen addition, a 31-bit configuration-derived seed, exactly one distance
geometry conformer attempt, rejection on negative conformer id, heavy-atom
coordinate copying, and appending a conformer to the heavy-atom result.

### Errors and no-center behavior

The former `NoStereoCenters` result conflicted with Python, which yields one
clone. `TooManyCenters` and `TooManyBonds` encoded local limitations absent from
the source. `IsomerGenerationFailed(0)` lost the actual operation context. The
replacement provides structured errors for real model or operation failures,
including errors delivered while consuming a lazy iterator.

### Tests and evidence

No committed oracle previously proved exact potential-stereo record order,
intermediate molecule state, option behavior, random sequence, emitted order,
full output state, or embedding coordinates. The old unit tests did not cover
the source cases in `UnitTestMol3D.py`, the relevant `catch_chirality.cpp`
sections, enhanced stereo fixtures, unknown/either bonds, large center counts,
or composition and concurrency.

## Reusable Source-Backed Infrastructure

The port must review and reuse these existing implementations where their
source boundary matches, rather than copy them into the enumerator:

| Capability | Existing owner | Required action |
|---|---|---|
| Swap parity | `source_port_helpers::count_swaps_to_interconvert` | Reuse directly after type adaptation. |
| Canonical ranks | `notation::canon_rank::rank_fragment_atoms` | Expose narrowly inside the crate and reuse the exact option set. |
| Symmetrized rings | `chemistry::rings::symmetrize_sssr` and read-parts variant | Reuse and preserve source ring order/cache behavior. |
| Tetrahedral eligibility helpers | `chemistry::stereo` | Consolidate only source-identical branches into the new engine. |
| Bridgehead and ring stereo helpers | `chemistry::stereo` and ring infrastructure | Reuse after line-level source comparison; do not fork them. |
| Canonical isomeric SMILES | notation writer | Reuse with the exact Python finalization state. |
| Hydrogen addition | registered hydrogen operation | Invoke on the owned isomer through the established operation path. |
| Distance geometry | `chemistry::distgeom` | Reuse the pinned one-conformer seeded path. |
| Stereo state | typed `ChiralTag`, `BondStereo`, `BondDirection`, `StereoGroup` | Extend with typed potential-stereo records, not integer sentinels. |

Existing helpers are not automatically proven suitable merely because their
names correspond to RDKit. Each call site must be checked against the exact
source overload, options, mutation order, property-cache state, and error
semantics before reuse.

## Required Convergence

Implementation must leave exactly four coordinated layers:

1. One typed potential-stereo engine owns all `FindPotentialStereo` chemistry.
2. One private typed flipper core consumes only potential-stereo records and
   stereo groups.
3. One lazy configuration and finalization engine implements the Python
   generator and count behavior.
4. Rust and Python public surfaces are thin project-native adapters over that
   single core.

The old candidate detectors, eager combination generator, xorshift state,
center caps, local feasibility checks, and three enumeration loops must be
removed or reduced to deliberate thin compatibility delegates. No compatibility
wrapper may retain independent chemistry or configuration logic.

## Observable Parity Requirements

Potential-stereo oracle rows must compare ordered `StereoInfo` records and the
complete affected molecule state, including typed center kind, specified
state, descriptor, permutation, ordered optional controlling atoms, atom and
bond stereo/directions, stereo atoms, unknown state, computed properties,
chiral ranks, ring state, and source preservation.

Enumeration oracle rows must compare option identity, flipper kind and order,
configuration sequence, theoretical count, emitted count and order, canonical
isomeric SMILES, complete atom/bond/stereo-group/property state, conformer
count, embedding outcome, and successful coordinates. Discrete state is exact;
only distance-geometry coordinates may use the repository's already justified
pinned-source tolerance.

The focused source fixtures, repository small corpus, maintained 5,000-row
corpus, and bounded complete ChEMBL 37 phase must all finish with zero
unexplained mismatch before support metadata can claim RDKit parity.
