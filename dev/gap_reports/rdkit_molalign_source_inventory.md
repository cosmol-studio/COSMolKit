# RDKit MolAlign Source Inventory

## Source Pin And Boundary

The behavior source is pinned RDKit `2026.03.1`, revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`, vendored under
`third_party/rdkit/`.

The ordinary alignment boundary is:

| Source | Lines | Owned behavior |
|---|---:|---|
| `Code/Numerics/Alignment/AlignPoints.cpp` | 16-340 | Quaternion point alignment, optional weights/reflection, Jacobi eigensolver, SSR, and transform construction. |
| `Code/Numerics/Alignment/AlignPoints.h` | 22-44 | Public defaults and the SSD-not-RMS result contract. |
| `Code/GraphMol/MolAlign/AlignMolecules.cpp` | 29-470 | Molecule/conformer mapping, symmetry, best-map selection, RMS calculations, and coordinate application. |
| `Code/GraphMol/MolAlign/AlignMolecules.h` | 45-404 | Parameter defaults, public overloads, result ordering, and mutation descriptions. |
| `Code/GraphMol/MolAlign/Wrap/rdMolAlign.cpp` | 35-320, 669-945 | Python conversion, overload defaults, list/map validation, and Python-visible result shapes. |
| `Code/GraphMol/MolTransforms/MolTransforms.cpp` | 445-451 | Application of a 4x4 transform to every conformer position. |

`RandomTransform`, `O3AAlignMolecules.{h,cpp}`, MMFF O3A, Crippen O3A,
constraint O3A, and the `PyO3A` wrapper are not reachable from the functions
above and are outside this port. They remain separately named future
capabilities rather than unsupported branches inside ordinary MolAlign.

## Numerics Call Graph

```text
AlignPoints (272-340)
  -> _sumOfWeights (57-65)
  -> _weightedSumOfPoints (22-38)
  -> _weightedSumOfLenSq (40-55)
  -> _computeCovarianceMat (67-95)
  -> reflectCovMat (263-270)
  -> _covertCovMatToQuad (97-132)
  -> jacobi (163-261)
  -> Transform3D::SetRotationFromQuaternion
  -> Transform3D::Reflect
  -> Transform3D::TransformPoint
  -> Transform3D::SetTranslation
```

`AlignPoints` requires equal point counts. When weights are supplied it also
requires equal weight count and every weight is strictly greater than zero.
It deliberately uses the source loop and arithmetic order, a 4x4 Jacobi
solver with tolerance `1e-6`, and the minimum eigenvector as quaternion. The
return value is SSR; molecule adapters divide by the number of mapped points
and take the square root. A negative SSR is clamped only when its magnitude is
below `1e-6`.

The pinned implementation has no BLAS, LAPACK, Eigen, or other mathematical
FFI dependency. The exact C++ implementation can be built behind a test-only C
ABI shim as an oracle, but production must not invent an FFI dependency that
is absent from the source call graph.

## MolAlign Call Graph And Observable Behavior

```text
getAlignmentTransform (215-239)
  -> Conformer lookup by ID
  -> optional one-result SubstructMatch(reference, probe)
  -> alignConfsOnAtomMap (57-70)
  -> AlignPoints

alignMol (241-252)
  -> getAlignmentTransform
  -> MolTransforms::transformConformer(probe)

getBestAlignmentTransform (254-270)
getBestRMS (273-305)
getAllConformerBestRMS (307-372)
CalcRMS (375-396)
  -> getAllMatchesPrbRef (72-110), unless maps are supplied
  -> getBestRMSInternal (139-211)
     -> alignConfsOnAtomMap, or calcMSDInternal (112-137)

alignMolConformers (418-470)
  -> _fillAtomPositions (398-416)
  -> AlignPoints
  -> MolTransforms::transformConformer
```

### Mapping

- Atom-map pairs are `(probe atom index, reference atom index)`.
- `getAlignmentTransform` without a map requests one match with recursive
  queries enabled, chirality disabled, and query-query matching enabled.
- Best-map APIs match the probe as query into the reference with
  `uniquify=false`, recursive queries enabled, chirality disabled,
  query-query matching disabled, and default `maxMatches=1_000_000`.
- Optional symmetry creates a probe clone, neutralizes matched terminal O/N,
  and replaces the terminal bond with the exact single-or-double query.
- `ignoreHs` removes map entries whose probe atom is hydrogen after matching.
- No match throws `MolAlignException`; it is not an empty successful result.

### RMS And Mutation

- `getAlignmentTransform` returns aligned RMSD and a transform without
  changing either molecule.
- `alignMol` applies that transform to the selected probe conformer.
- `getBestAlignmentTransform` returns best RMSD, transform, and atom map
  without mutation.
- RDKit `getBestRMS` computes the same best transform and then mutates the
  selected probe conformer. COSMolKit will separate the read-only measurement
  from explicitly named coordinate transformation APIs.
- `CalcRMS` enumerates maps but does not align coordinates; it measures the
  minimum RMS in the current coordinate frames and does not mutate graph or
  coordinates despite accepting a mutable C++ molecule reference.
- `getAllConformerBestRMS` is read-only and returns pairs in triangular order
  `(1,0), (2,0), (2,1), ...` using stored conformer IDs.
- `alignMolConformers` uses the first selected conformer, or the first stored
  conformer, as reference and transforms every subsequent selected conformer.

### Selection, Weights, And Threads

- A conformer ID of `-1` selects the first stored conformer; nonnegative IDs
  select by `Conformer::getId()`, not container index.
- Explicit maps are evaluated in caller order.
- A supplied weight vector must match each evaluated map. All weights must be
  positive in aligned-RMS paths; coordinate-frame RMS multiplies weights but
  relies on the wrapper/map checks for length.
- Candidate replacement uses strict `msd < msdBest`; equal values retain the
  first result in the source's effective merge order.
- Threaded best-map evaluation distributes match indices by
  `index % numThreads`, then scans thread result buckets by thread index.
- Threaded all-conformer calculation stores each result by its original pair
  index, so output order remains triangular and deterministic.

## Python Wrapper Boundary

The wrapper converts Python lists into maps and weight vectors, returns
transforms as NumPy-compatible matrices, and exposes defaults of conformer ID
`-1`, `maxMatches=1_000_000`, `symmetrizeConjugatedTerminalGroups=true`,
`reflect=false`, `maxIters=50`, and `numThreads=1`. Empty Python map/list
arguments conventionally mean automatic mapping rather than an explicit
zero-point alignment.

COSMolKit will retain these chemical defaults but use project-native method
names. No API named only as an RMSD calculation will reproduce RDKit's hidden
probe mutation.

## Required Dependency Reuse

- Reuse the canonical substructure matcher and its query recursion/order.
- Reuse one real conformer-ID resolver across MolAlign, MolTransforms, matrix
  helpers, and coordinate-dependent matching.
- Reuse one terminal-group symmetry helper across distgeom pruning and
  MolAlign mapping.
- Reuse one 4x4 transform type/application path across depiction, distgeom,
  MolTransforms, and MolAlign.
- Reuse the operation registry for every public coordinate mutation.
- Use Rayon only where its scheduling and merge order explicitly reproduce
  the source result order.

## Source-Level Acceptance Boundary

Ordinary MolAlign is complete only when every function above is represented by
one source-marked Rust implementation or an exact reuse proof, its public
observable values and errors are compared with pinned RDKit, and no separate
alignment math or match-enumeration implementation remains in a caller.

## Initial COSMolKit Audit

The repository already contains the following source-backed pieces:

| Existing location | Finding | Required disposition |
|---|---|---|
| `chemistry/coordinates.rs:5415-5715` | A fixed-array quaternion/Jacobi implementation returns a 4x4 transform and is used by constrained 2D depiction. It currently accepts only unweighted points and owns its transform helpers locally. | Move the behavior into the shared alignment module, retain the 2D caller, and add the weighted/reflection/error paths from `AlignPoints.cpp` without changing validated operation order. |
| `chemistry/distgeom.rs:2816-3125` | A second fixed-array implementation uses `ForceFieldVec3`, duplicates the covariance/quad/Jacobi code, and exposes an unweighted SSR-only helper to conformer pruning and tests. | Replace the duplicate implementation with a conversion-free call to the shared core; keep the SSR-only caller contract as an adapter. |
| `chemistry/coordinates.rs:5740-5815` | Existing private 2D best-alignment code enumerates maps and selects strict lower RMS, but uses a local `max_matches` path and 2D-specific conformer selection. | Reuse only its proven map-selection behavior; do not expose it as a third MolAlign implementation. |
| `chemistry/mol_transforms.rs:62-91` | `conf_coords` and `conf_coords_mut` resolve a nonnegative value as a vector index. | Add one ID-based resolver and migrate coordinate-dependent public paths before claiming MolAlign conformer parity; preserve explicit index APIs only where their names document index semantics. |
| `search/substruct.rs:41-52,188-285,4399-4430` | Public matcher returns `query -> target` mappings and supports recursive queries, non-unique results, max matches, chirality and query-query controls. | Use this as the sole mapping engine; adapt probe-query/reference-target orientation at one private boundary. |
| `chemistry/distgeom.rs:3296-3410` | Terminal O/N conjugated-group symmetrization is already source-marked but is private to distgeom and returns a cloned molecule. | Promote the helper into a shared private MolAlign support module; preserve clone-only query edits and cache lifecycle. |
| `operations/ops.rs:883-1004` | Coordinate writes are already registered for coordinate replacement/add/remove operations. | Add only explicitly named alignment mutation operations with coordinate write access and drawing invalidation; no registry entry for read-only RMSD. |
| `chemistry/mol_transforms.rs:1009-1045` | Transform application exists as a value-style function but resolves conformers by index and is not the MolAlign public operation. | Reuse the point-application primitive after ID resolver consolidation; route public alignment mutation through registered operations. |
| `properties/batch.rs` and `python/src/lib.rs` | Batch execution and generated stub patterns already exist, but no MolAlign methods are present. | Add scalar/batch wrappers only after the core API is stable; preserve order and error indexes without copying alignment logic. |

### Duplication And Numerical Findings

The two existing numerical implementations are structurally close but are not
yet proven identical for all source behavior. The coordinates path has a
different empty-input error type, no weights, and a separate 2D conformer
resolver. The distgeom path computes only unweighted SSR and does not build or
apply a transform. Both currently use Rust `sqrt` in some paths while the
coordinates path uses an unsafe C `libm` wrapper. This is a numerical-policy
decision to resolve in the shared module and test, not a reason to add a second
algorithm.

`ForceFieldVec3` is a force-field domain type and should not become the public
MolAlign coordinate type. The shared core should operate on a narrow fixed
`[f64; 3]` view; distgeom may adapt its existing values at the boundary without
retaining a second covariance/Jacobi implementation.

The current matcher has a default `max_matches` of 1000, while MolAlign's
automatic mapping default is 1,000,000. MolAlign must construct explicit
`SubstructMatchParams` rather than rely on matcher defaults. Existing
`rdkit_substruct_matches_unordered` also enumerates with a 1000 limit and is
not suitable for BestRMS without changing its call contract.

The existing terminal symmetrization implementation is behaviorally useful,
but its `DgBoundsError` return type is domain-specific. Promotion requires a
narrow shared error or a private adapter that preserves the source invariant
failure without making distgeom's error type part of MolAlign's API.

No current Python or Rust public API mutates a molecule under an RMSD-named
calculation. The new surface can therefore preserve the project policy by
making all calculations read-only and introducing mutation only under explicit
`align_*_` operations.

## Final Closure

The initial reuse and duplication findings above have been resolved. The sole
numerical implementation is `chemistry/numerics/alignment.rs`; depiction,
distance geometry, MolAlign, and transform application delegate to it. The
shared conformer-ID resolver, matcher, and terminal-symmetry helper now own
their respective paths. Rust and Python expose the read-only measurements and
explicit coordinate mutations defined in the API design, and both coordinate
mutation families are registered under `MOLALIGN_FEATURE`.

Focused, project-small, maintained 5,000-row, and complete ChEMBL 37 parity
gates pass. The final audit and exact counters are recorded in
`dev/gap_reports/rdkit_molalign_full_port_validation.md`. O3A remains the only
named future boundary identified by this inventory.
