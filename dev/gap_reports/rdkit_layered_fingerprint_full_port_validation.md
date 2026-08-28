# RDKit Layered Fingerprint Full-Port Validation

## Boundary And Source Pin

This report closes the implementation audit for the Layered fingerprint
boundary defined in
[`rdkit_layered_fingerprint_full_port_plan.md`](../plans/rdkit_layered_fingerprint_full_port_plan.md).
The oracle is RDKit `2026.03.1`, source revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`, vendored under
`third_party/rdkit/`. The upstream algorithm version is `0.7.0` and upstream
continues to describe this fingerprint family as experimental.

The audited boundary is `LayeredFingerprintMol()` and the reachable query,
ring, path/subgraph, hash, and bit-vector behavior listed in
[`rdkit_layered_fingerprint_source_inventory.md`](rdkit_layered_fingerprint_source_inventory.md).
Pattern Fingerprint and other fingerprint families are outside this report.

## Source-Line Closure

| Source-owned behavior | Sole COSMolKit owner | Closure evidence |
|---|---|---|
| constants, flags, defaults, parameters, result | `properties/fingerprint.rs:4315-4422` | ten source slots, version `0.7.0`, `0x07` substructure mask, six active flags, retained high flags, exact defaults, seeded counts, mask, and absent/present roots |
| atom query complexity and aromaticity | `search/query.rs:1001-1169` | copied source bodies, sibling-carried atomic-number state, ordered compound-query handling, and focused branch tests |
| bond query complexity | `search/query.rs:3178-3247` | copied source body and exact simple/complex query-shape tests |
| ring query accessors | `search/query.rs:4269-4316` | typed access through the sole `RingInfo`, with initialized and fused-ring tests |
| 32-bit `gboost::hash_range` | `properties/fingerprint.rs:5063` | the shared sole hash implementation, empty/order/wrapping/non-native-width tests |
| linear and branched enumeration | `properties/fingerprint.rs:6144-6277` | shared sole enumerator, absent/empty/duplicate/reordered root and prepend-order tests |
| preconditions, SSSR, and caches | `properties/fingerprint.rs:6287-6413` | exact source errors, initialized-ring reuse or exact SSSR, bond-index cache, query masks, aromaticity, and atomic numbers |
| six layer encoders | `properties/fingerprint.rs:6416-6650` | copied source blocks plus packing, modulo, ordering, query suppression, omission, ring, and fused-ring tests |
| projection | `properties/fingerprint.rs:6653-6738` | copied source block plus sorting, suffix, 32-bit hash, mask, collision, repeated-path, and once-per-path count tests |
| complete orchestration | `properties/fingerprint.rs:6740-6991` | one scalar core over the shared dependencies; no Layered-local ring, query, hash, enumerator, or vector implementation |

No `RDKit❌` marker remains in the Layered implementation or its shared
QueryOps helpers. The one `RDKit❗✔️` line is the pinned unrooted-linear source
defect documented below; it marks the intentional process-safety difference
instead of misrepresenting the upstream crash as reproduced. Static
architecture tests inventory the production owners and reject a second
Layered encoder, projection core, query helper, hash range, path enumerator,
vector type, scalar core, or Python/batch algorithm copy.

## Branch And Option Closure

The focused tests cover all six active layers independently and in cumulative
combinations; `0x07`, `0x3f`, `0xffffffff`, zero, and inactive high-only flags;
default and non-power-of-two widths; path bounds; branched and linear paths;
absent, explicitly empty, single, duplicate, and reordered roots; exact root
group prepend order; all-zero, all-one, even, modulo-three, and sparse masks;
unseeded and seeded counts; repeated hashes and collisions; query masks
`0..7`; initialized and uninitialized ring state; ordinary, aromatic,
query-bearing, disconnected, fused, and unusual-bond inputs.

The five committed source-derived query fixtures compare the complete
intermediate bond query-mask vector, query-aware aromatic atom vector, and
final fingerprint output. The complete maintained profile has 18 branches,
including every individual layer, active/default/substructure flags, two
width/path configurations, branched and rooted-linear paths, two rooted
branched configurations, two mask/count configurations, no active layer, and
inactive high source bits.

## Output And Error Closure

Successful comparisons include all observable Layered outputs exposed by the
source boundary:

- exact fingerprint width;
- complete ordered on-bit set;
- complete optional atom-count vector, including caller-provided seed values;
- one count increment per atom per accepted path, independent of the number
  of accepted layers;
- no mutation of the input molecule or caller-owned mask.

The Rust surface returns structured `FingerprintError` values for `minPath==0`,
`maxPath<minPath`, `fpSize==0`, `bad atomCounts size`, `bad setOnlyBits size`,
out-of-range roots, ring-preparation failure, and an invalid enumerated bond
index. Python maps the argument failures to `ValueError` without clamping,
filtering, or repairing input. `None` roots and an explicitly empty root list
remain observably distinct.

### Upstream Unrooted-Linear Defect

Pinned `Fingerprints.cpp:312` calls
`findAllPathsOfLengthsMtoN(mol, minPath, maxPath, false)` for an unrooted
linear request. In the pinned signature, that positional `false` selects atom
indices instead of the bond indices consumed by the following Layered loop.
Isolated oracle processes confirmed `SIGSEGV` for `CC`, `CCO`, `CCCC`, and
`CC.CCCC`; a ring input happened to remain in range and returned a value.

COSMolKit deliberately does not reproduce this input-dependent process crash.
It follows the function's documented bond-path result and the source's rooted
linear call, producing a deterministic bond-path fingerprint. This is a
process-safety compatibility exception, not an inferred chemistry fallback.
It must remain documented wherever exact failure behavior is described. All
maintained oracle comparisons use the valid rooted linear source path; no
zero-mismatch total below counts the crashing upstream branch as a match.

## Allocation And Complexity Closure

Local source review found the same dominant work as RDKit: one ordered atom
and bond cache pass, the shared source-backed path/subgraph enumeration, the
same `O(P^2)` bond-neighbor loop for a path of `P` bonds, fixed `O(1)` layer
encoders, one sort and one 32-bit hash per emitted layer, and at most one atom
mask scan per accepted path. The implementation does not clone the molecule,
does not enumerate paths twice, and uses the existing dense `Fingerprint`
words.

The preserved benchmark artifact is
`tmp/layered-fingerprint-benchmark/raw.json`, schema
`cosmolkit-layered-fingerprint-benchmark-v1`, profile SHA-256
`e24e8ffc182f6bff89aee40412a3025c6df5219165176b714de9842e1a687051`.
Across path-rich, ring-rich, query-bearing, disconnected, and large cases,
the COSMolKit/RDKit median ratios were respectively `1.0773`, `0.8430`,
`1.0666`, `0.8460`, and `1.1324`. COSMolKit peak working-set deltas were
`0-20 KiB`. The measurements show no order-of-magnitude regression,
unexpected allocation growth, duplicate enumeration, or worse asymptotic
shape.

## Rust And Python API Closure

Rust exposes typed `LayeredFingerprintLayers`, `LayeredFingerprintParams`, and
`LayeredFingerprintResult`, free-function delegates, method-oriented
`Molecule` calls, and ordered `MoleculeBatch` calls. Python exposes scalar and
ordered batch methods plus the typed result object. All surfaces delegate to
the sole Rust core. Composition tests interleave Layered with Morgan, MACCS,
Avalon, RDKFingerprint, AtomPair, and Topological Torsion on shared molecules,
including repeated and concurrent calls, and verify stable output, input
immutability, property preservation, order, and option isolation.

Generated stub tests lock method names, parameter order, defaults, return
annotations, uniqueness, and runtime signatures. The focused Python Layered
suite passed `5` tests. `basedpyright python/tests python/examples` reported
`0 errors` and `528 warnings`; its warning-only nonzero exit remains a
repository-wide baseline and is not represented here as a clean warning gate.

## Maintained Corpus Evidence

| Corpus | Input rows | Mutually parsed | Mutually rejected | Branch executions | Exact output-field comparisons | Mismatches |
|---|---:|---:|---:|---:|---:|---:|
| focused | 16 | 16 | 0 | 288 | 864 | 0 |
| `smiles_small` | 152 | 140 | 12 | 2,520 | 7,560 | 0 |
| `smiles_5000` | 5,000 | 5,000 | 0 | 90,000 | 270,000 | 0 |
| total | 5,168 | 5,156 | 12 | 92,808 | 278,424 | 0 |

Each successful branch compares width, the complete ordered bit set, and the
complete optional atom-count vector. Rejected rows are retained and require
COSMolKit to reject parsing as well; there is no row filtering, tolerance,
allowlist, or fallback.

## ChEMBL 37 Evidence

The complete external run is preserved under `tmp/chembl37-layered-full/`.
Its manifest records:

| Measure | Value |
|---|---:|
| source rows processed | 2,897,819 |
| rows accepted by both parsers | 2,897,804 |
| rows rejected by both parsers at parse entry | 15 |
| deterministic shards completed | 128 / 128 |
| active profiles per parsed row | 18 |
| exact Layered comparisons | 52,160,472 |
| blocking mismatches | 0 |
| informational mismatches | 0 |
| retained mismatch examples | 0 |

The corpus source SHA-256 is
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`,
the corpus manifest SHA-256 is
`f8b1a516a7794c8d4428a8309c80c049b5c4f34c6c9c54381f0aebccd9ccb976`,
and the effective run-profile SHA-256 is
`eac2f19382bb454a1ba9e0387da3b46328ef6edc199d15f9d57ebadce87adf1d`.
The Layered profile SHA-256 matches the benchmark and maintained generator
profile above.

## Validation Commands And Results

The plan executed and passed the focused query, hash, enumeration,
preparation, encoder, projection, end-to-end, query-fixture, workflow,
composition, architecture, maintained-corpus, Python API, and stub-surface
targets immediately after their artifacts were added. Final pre-report checks
also produced these results:

```text
cargo check -p cosmolkit-core --features op-contracts-strict
  passed

cargo test -p cosmolkit-core --release --features op-contracts-strict
  core library: 3349 passed; 0 failed; 46 ignored
  all default integration targets passed
  doctests: 7 passed; 0 failed

.venv/bin/pytest -q python/tests/test_layered_fingerprint.py
  5 passed

.venv/bin/pytest -q python/tests/test_stub_surface.py -k layered_fingerprint_methods
  1 passed

.venv/bin/basedpyright python/tests python/examples
  0 errors; 528 warnings; nonzero warning exit

cargo test --workspace --release --features cosmolkit-core/op-contracts-strict
  passed across every workspace crate, integration target, and doctest
```

The focused and complete Layered parity targets were rerun after documentation
and generated-artifact changes, and the final strict workspace release suite
completed successfully after all Layered artifacts were in place.

## Closure Decision

The implemented Layered chemical-output boundary, its six active encoders,
options, maintained error surface, Rust/Python APIs, sole-core architecture,
composition behavior, complexity shape, and all non-crashing oracle profiles
are closed by source anchors and zero-mismatch evidence. No Layered-specific
behavioral, allocation, API, or corpus gap remains inside that validated
boundary.

Exact emulation of the pinned upstream unrooted-linear process crash is
intentionally excluded. COSMolKit provides the source-documented bond-path
operation instead. Documentation and support claims must retain this narrow
compatibility note and must not describe upstream crash behavior as matched.
Repository-wide Python warning cleanup is also outside the Layered port and
remains separate environmental debt.
