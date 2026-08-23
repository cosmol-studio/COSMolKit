# RDKit Topological Torsion Fingerprint Full-Port Validation

## Audit identity

- RDKit version: `2026.03.1`
- RDKit source revision: `351f8f378f8ad6bbd517980c38896e66bf907af8`
- Reference identity source: `testdata/reference/rdkit.json`
- Audited COSMolKit feature: `fingerprint.topological_torsion`
- Distinct non-target family: `TopologicalFingerprint*` / RDKit
  `RDKFingerprintMol`
- Runtime-oracle policy: production code has no RDKit, Python-oracle, or
  fixture-table dependency

This report audits the source-port boundary implemented by
`dev/plans/rdkit_topological_torsion_fingerprint_full_port_plan.md`. It does
not broaden the feature to Atom Pair fingerprints, RDKFingerprint, or other
RDKit fingerprint families.

## Result

The documented Topological Torsion boundary is implemented and exact for the
modeled input state. Modern and legacy public adapters terminate in one atom
invariant, path enumeration, torsion code/hash, and shared fingerprint-vector
driver. The focused and maintained-corpus gates have zero mismatches and do
not use tolerance, sampling, expected failure, skipped row, or accepted
mismatch policies.

Two intentionally visible reproduction qualifications remain:

1. The missing-`_StereochemDone` branch uses the repository's shared
   `assign_stereochemistry_cleanup_subset()` implementation. Its private-copy
   ownership and original-versus-copy selection match RDKit, and every
   chirality profile in this port passes exactly, but the shared stereo helper
   is not a claim for every upstream `MolOps::assignStereochemistry` state.
2. RDKit emits runtime deprecation log messages from the three legacy C++
   functions. Rust provides `#[deprecated]` compile-time diagnostics instead;
   vector and error behavior is exact, while warning timing/channel is an
   intentional language-level difference.

These lines remain `RDKit❗✔️` and are not silently promoted. There is no
first-axis `RDKit❌` line in an implemented Topological Torsion chemistry or
vector path. The visible `RDKit❌❌` cases in fingerprint type/JSON dispatch
belong to separately scoped Atom Pair and RDKFingerprint generator families.

## Source-line closure

The audited calling path contains inline verbatim source anchors and two-axis
markers for:

- AtomPairs constants, `numPiElectrons`, `getAtomCode`, and
  `AtomPairAtomInvGenerator` construction/invariants/metadata/clone;
- `getTopologicalTorsionCode` and `getTopologicalTorsionHash`;
- `extendPaths`, `pathFinderHelper`, `findAllPathsOfLengthsMtoN`, and
  `findAllPathsOfLengthN`;
- shared `MolOps::getDistanceMat` and `FloydWarshall` behavior;
- common fingerprint arguments, `AdditionalOutput`, generator ownership,
  JSON, information strings, `getFingerprintHelper`, the four scalar output
  forms, count simulation, provenance duplication, and multithread bulk
  collection;
- Topological Torsion arguments, environment result size, environment
  generation, environment provenance, metadata, and both factories;
- legacy unfolded count, hashed count, `TorsionFpCalc`, and hashed-bit
  threshold/block sizing;
- Python generator/options/scalar/bulk/JSON bindings and the legacy/helper
  surface.

Source-undefined empty paths, shifts at or beyond 64 bits, zero divisors, and
impossible size arithmetic are deliberately represented by structured Rust
errors. They are documented beside the source anchors and are not fabricated
as chemistry results.

## Single-core architecture

The single-core contract is satisfied:

- `properties::fingerprint::generator::FingerprintGenerator` is the sole
  implementation of the RDKit
  `FingerprintGenerator.cpp:323-652` environment-to-count, folding,
  count-simulation, random-extra-bit, explicit-bit, and provenance driver.
- `MorganFingerprintGenerator` and
  `TopologicalTorsionFingerprintGenerator` only supply family arguments,
  invariant generation, environment generation, and result width to that
  driver.
- No historical `build_fingerprint`, `fold_invariant`, or
  `atom_is_excluded` implementation remains.
- Modern scalar conveniences, typed dispatch, generator bulk calls, Python
  calls, and Rust batch conveniences delegate to the same scalar/shared
  generator code.
- Legacy hashed count and bit adapters call the modern generator core.
  Legacy unfolded count switches only the source-required endpoint/internal
  atom-code transformation, then uses the same environment/path/code/vector
  machinery.
- Topological Torsion hashing calls the existing Boost-compatible
  `hash_combine`; no family-local hash formula exists.
- `only_shortest_paths` calls the shared
  `chemistry::matrices::topological_distance_matrix`, also used by
  distance geometry. There is no second Topological Torsion Floyd-Warshall
  implementation.
- RDKFingerprint's branched path/subgraph enumeration remains separate and is
  not used as a shortcut for torsion paths.

## Exact parity evidence

### Focused matrix

- Rows: 152, consumed contiguously and completely.
- Profiles per row: 12.
- Profiles: all provenance, chirality, count-simulation provenance, custom
  atom invariants, default, extra bits, folding collisions, first-atom root,
  first-atom ignore, shortest paths, custom three-atom count bounds, and
  unfolded without count simulation.
- The same rows also carry legacy outputs, helper outputs, JSON metadata, and
  structured error expectations where applicable.
- Golden SHA-256:
  `a47c321f6556a10092e8d44c4b17eb0109f1432b00c0e155cc1e680e92f62158`.
- Exact focused Rust parity target: passed.

### Maintained 5,000-row matrix

- Rows: exactly 5,000, with contiguous indexes `0..4999`.
- Profiles per row: 9.
- Profiles: all provenance, chirality, count-simulation provenance, custom
  atom invariants, default, extra-bit folding, five-atom shortest paths,
  custom three-atom count bounds, and unfolded without count simulation.
- Golden SHA-256:
  `28d91a8ca475ffedf43ac0e50621be2b0b6ae394225dffbe76b6fff061affe7d`.
- The manifest records the same checksum and 5,000 records.
- Release exact-parity command passed `1 passed; 0 failed`; the test body ran
  in 17.23 seconds on the audit host after compilation.

The parity reader streams JSONL, evaluates rows in deterministic indexed
parallel execution, retains only row/failure summaries, sorts by input row,
and asserts the full contiguous corpus and exact profile set. It cannot pass
by ignoring an output field or a failed row.

The 5,000-row target is also a required coverage-workflow step. CI restores or
regenerates the pinned expected-data artifact, validates its manifest before
the test starts, and propagates a parity failure after coverage artifacts have
been collected. The matrix is therefore an active gate rather than a manual
audit command.

### Complete ChEMBL 37 matrix

- Source rows: 2,897,819 over 128 immutable shards.
- Mutually parseable rows: 2,897,804; both libraries rejected the remaining
  15 rows.
- Profiles per mutually parseable row: 9.
- Every profile compares sparse-count, folded-count, sparse-bit, and
  explicit-bit vectors. The two provenance profiles also compare atom counts,
  atom-to-bits, bit-info maps, bit paths, and atoms-per-bit for all four vector
  forms.
- Exact comparisons: 127,503,376.
- Result: 128/128 successful shards, zero mismatch, zero invalid profile task,
  and zero retained finding.
- Aggregate SHA-256:
  `e0c753ba945b0ce9cde98f4df26282ae0ea94ac5eb89398ba345774ad907ffe4`.

This is the same reproducible ChEMBL 37 source and pinned RDKit boundary used
by the complete AtomPair audit. The profile definition, runner, aggregation,
and zero-mismatch acceptance rules are version controlled; corpus shards and
run artifacts remain outside Git.

## Live options and molecule-state boundaries

The generator reads live options at each call, matching the mutable RDKit
argument object rather than caching construction-time derived state. In
particular, changing `torsionAtomCount` or `includeChirality` changes the
unfolded sparse-vector width on the next call. Width mappings that would shift
by 64 bits or more return a structured error instead of invoking undefined
integer behavior.

The constructor retains the source preconditions for non-empty count bounds
when count simulation is enabled, nonzero folded-vector size, and nonzero
`numBitsPerFeature`. Direct live-option mutation and JSON restoration can
create states that bypass those construction checks in RDKit, so generation
reproduces the source call path instead of reapplying constructor policy:

- `numBitsPerFeature = 0` emits the primary bit for each feature and no extra
  random bits;
- empty count bounds remain valid for sparse-count generation, which does not
  consume them, but return a structured error for count-simulated bit output;
- zero `fpSize` remains valid for unfolded sparse-count and sparse-bit output
  when count simulation is disabled, but returns a structured error for folded
  count and explicit-bit output.

Atom-pair invariants read the typed explicit-valence cache, matching
`Atom::getValence(EXPLICIT)`. An unsanitized molecule without that cache now
returns `ExplicitValenceCacheNotInitialized`; calling the public valence-cache
assignment operation makes the same molecule eligible for generation. The
implementation does not silently recompute a different state inside the
fingerprint call.

Focused full-API regressions additionally lock query atoms and bonds, directed
dative bonds, a hypervalent sulfur graph, and an explicit-hydrogen graph. Both
modern unfolded sparse-count output and the legacy unfolded adapter are
compared against exact RDKit identifiers. Boron endpoint cases (`BCCC` and
`BCCCC`) separately lock the source-required legacy atom-code correction
boundary. A dedicated operation-composition regression verifies that
`with_hydrogens()` immediately supports both modern and legacy Topological
Torsion calls. Hydrogen addition commits the source-defined explicit-valence
cache transition through the operation contract; the fingerprint path does not
recompute it or require an extra public operation.

## Public and batch behavior

- Rust exposes distinct modern parameter, generator, four-vector,
  output-request/result, legacy parameter/result, and helper types.
- Python exposes generator construction, live options, four scalar and four
  ordered bulk vector forms, shared `AdditionalOutput`, JSON roundtrip, three
  legacy functions, AtomPairs constants/code explanation, and torsion score,
  explanation, and id helpers.
- Rust `MoleculeBatch` exposes four vector conveniences and a typed
  vector-plus-provenance convenience. Existing invalid records remain `None`
  at their original positions. New computation failures carry
  `BatchRecordError.index`, operation, and message in `BatchValidationError`.
- Focused strict batch tests compare every output with scalar calls for
  default, 1-, 2-, 3-configured-, and 4-thread execution, repeat parallel
  provenance generation eight times, verify invalid-record alignment and new
  error indexes, and assert batch/source molecule immutability. The target
  passes `2 passed; 0 failed`.

## Errors and fail-closed behavior

The active public boundary returns `FingerprintError` for invalid constructor
arguments, empty count bounds where source output construction consumes them,
zero folded or explicit-bit sizes, source-undefined result-width shifts,
missing explicit-valence cache state, invalid atom selection indexes,
insufficient custom invariant arrays, malformed or incompatible JSON,
arithmetic overflow, and parallel pool creation failure. Live
`num_bits_per_feature = 0` is intentionally not in this list because the
source generation path accepts it as described above. Legacy zero
sizes/divisors and compatibility-boundary overflow also fail closed. Python
maps these to typed Python exceptions. No error branch selects Morgan,
RDKFingerprint, an approximate torsion, or an empty placeholder.

## Determinism, allocation shape, and complexity

- Scalar computation is deterministic. Extra bits use the single
  RDKit-compatible fingerprint MT RNG seeded per source environment.
- Generator bulk collection uses Rayon indexed parallel iterators, preserving
  input order. Each worker owns its call arguments; no shared mutable scratch
  output exists.
- Path enumeration preserves source atom-index scan order. Its dense adjacency
  matrix is `O(V^2)`, optional Floyd-Warshall is `O(V^3)`, and path expansion
  remains output-sensitive/exponential in the same way as the source.
- Bond-set deduplication uses ordered vectors of bit membership and a linear
  search, matching the source's vector-plus-`std::find` complexity. `BTreeMap`
  length buckets have the same logarithmic map class as the C++ ordered map.
- Each emitted path owns its atom sequence as in the source. Requested
  `bitPaths` and `atomsPerBit` each own their recorded sequences; no provenance
  is allocated when it is not requested.
- Chirality preparation clones a molecule only when chirality is requested
  and `_StereochemDone` is absent, matching the source private-copy branch and
  preserving the caller's molecule.
- The legacy hashed-count adapter returns the shared core's owned sparse
  vector directly, eliminating the C++ wrapper's second allocation and
  element-copy pass; the nearby `RDKit✔️🔝` rationale records this safely.
- Rust `MoleculeBatch` convenience methods intentionally construct the scalar
  call per record to guarantee a single public scalar behavior path. This adds
  small per-record generator setup versus direct generator bulk calls without
  changing asymptotic chemistry work. Callers requiring the lower-overhead
  bulk shape can use the generator's four `get*Fingerprints` methods, which
  reuse one generator.

No feature-specific Criterion benchmark exists. The 5,000-row release parity
runtime above is a validation observation, not a stable benchmark. Local code
inspection finds no algorithmic regression in the shared generator, path,
hash, distance, or vector cores. A future benchmark may compare scalar-loop,
generator-bulk, and `MoleculeBatch` convenience throughput without changing
the parity claim or marker status.

## Validation gates

The final sequential release gates completed as follows:

- `cargo fmt --all`: passed;
- `cargo check -p cosmolkit-core --features op-contracts-strict`: passed;
- `cargo test -p cosmolkit-core --release --features op-contracts-strict`:
  passed, including 2,900 core library tests with 46 ignored tests and every
  non-ignored integration/doc test;
- `cargo test --workspace --release --features
  cosmolkit-core/op-contracts-strict`: passed across the workspace;
- `uv sync --group dev` and `.venv`-targeted `maturin develop`: passed;
- the post-integration focused Python fingerprint selection passed 23 tests
  with 463 tests deselected; the Topological Torsion target contributes 12
  focused tests;
- pinned expected-data validation regenerated and published all 31
  `smiles_small/all` outputs plus the dedicated 5,000-row Topological Torsion
  output under the current generator identity; each output and manifest was
  schema- and checksum-validated before publication;
- Sphinx HTML documentation build: passed with zero warnings;
- `.venv/bin/basedpyright python/tests python/examples`: zero errors and 444
  repository-wide warning-only diagnostics. The command exits nonzero under
  the repository's current warning policy; the new Topological Torsion test
  and example files independently report zero errors and zero warnings;
- generated Python stub surface and its focused surface checks: passed.

The focused exact Rust parity target, exact 5,000-row release parity target,
complete ChEMBL 37 audit, and focused strict Rust batch target also passed as
described above. The feature-specific suite contains 64 source-function and
boundary Rust tests, six public-API tests, five atom-invariant tests, two
concentrated batch tests, and 12 focused Python tests.

## Final gap classification

There is no known output mismatch, hidden fallback, duplicate chemistry core,
or unconsumed maintained-corpus profile inside the declared feature. The
remaining explicit boundaries are:

- broader `assignStereochemistry` states outside the repository's modeled
  cleanup subset;
- runtime versus compile-time legacy deprecation notification;
- independently implemented AtomPair and RDKFingerprint families, which remain
  outside this feature-specific audit;
- the absence of a dedicated throughput benchmark.

These are documented boundaries or measurement work, not accepted mismatches
inside the exact parity matrix.
