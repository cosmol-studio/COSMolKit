# RDKit AtomPair Fingerprint Full-Port Validation

## Scope And Source Pin

This report validates the project-native AtomPair fingerprint surface against
RDKit `2026.03.1`, revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`. The selected boundary is the
modern `FingerprintGenerator` call graph rooted at
`RDKit::AtomPair::getAtomPairGenerator()`, including sparse-count,
folded-count, sparse-bit, explicit-bit, 2D and 3D distances, chirality,
selectors, custom atom invariants, count simulation, provenance, metadata,
JSON restoration, legacy-semantic adapters, and ordered batch execution.

## Validation Status

Behavioral and selected-source validation are complete with zero known
mismatch or unexplained implementation gap. The shared molecular
distance-matrix cache now reproduces RDKit's logically read-only reuse,
resolved-conformer keying, and topology/coordinate invalidation behavior.

| Boundary | Result |
|---|---|
| AtomPair source functions and branches | Implemented through the single modern generator core; no unresolved marker remains in `properties/fingerprint/atom_pair.rs` |
| Shared vector projection | One implementation in `properties/fingerprint/generator.rs` for all four output forms and count-simulation provenance |
| Pair enumeration | One production `i < j` loop in `properties/fingerprint/atom_pair.rs` |
| Hashing | Existing single shared `hash_combine()` implementation reused |
| Mutation behavior | Public fingerprint calls preserve logical molecule state; stereochemical preparation uses a working value |
| Structured errors | Invalid arguments, missing conformers, invalid invariant widths, missing valence state, and sparse-index overflow remain typed errors |
| Determinism and concurrency | Exact committed goldens, ordered batch checks, and the 128-shard ChEMBL run are deterministic and complete |
| Distance-matrix caching | Source-aligned 2D and resolved-conformer 3D entries are reused through the shared computed-property cache; focused tests cover cold/warm calls, clone isolation, invalidation, keying, and parallel initialization |

## Exact Parity Evidence

The committed golden matrix covers 12 focused records, 152 project-small
records, and all 5,000 maintained-strict records. Ten source-meaningful
profiles exercise default behavior, chirality, count-simulation disablement,
custom collision bounds, distance bounds, `fromAtoms`, `ignoreAtoms`, custom
atom invariants, `numBitsPerFeature`, and complete `AdditionalOutput`.
Every profile compares exact vector length and exact sparse elements or on-bit
sets for all four result forms; the provenance profile additionally compares
ordered atom counts, atom-to-bit rows, bit-info pairs, and atoms-per-bit rows.

The complete ChEMBL 37 audit processed all 2,897,819 source records over 128
shards. COSMolKit and RDKit both parsed 2,897,804 records and both rejected 15.
All 41 compared outputs for every mutually parseable record matched exactly:

- 40 exact fingerprint vectors: four output forms across ten profiles;
- one complete provenance output for the dedicated additional-output profile;
- 118,809,964 exact comparisons in total;
- zero blocking or informational mismatch;
- zero invalid profile task, missing shard, failed shard, or retained finding.

The ignored aggregate artifact is
[`tmp/atom-pair-chembl37/full/aggregate.json`](../../tmp/atom-pair-chembl37/full/aggregate.json),
SHA-256
`3a83c40b1ceeb63eab829121861907aa19a913800939f0ca9b66a02123415f15`.
The repository-owned workflow rejects partial shard sets, profile drift,
identity drift, and mismatched expected counts before accepting an aggregate.

## Single-Core And Source-Line Audit

`AtomPairFingerprintGenerator` is the sole production generator. Project-native
Rust methods, Python scalar methods, Python batch methods, JSON restoration,
and retained legacy-semantic adapters all delegate to it. AtomPair defines
only its invariant generation, pair environment, source code packing, and
family arguments. It reuses the shared fingerprint arguments, vector types,
random extra-bit generator, count simulation, provenance duplication, batch
scheduler, CIP labeler, hash combine, and distance-matrix implementation.

The source marker scan finds no unresolved AtomPair-family line. Remaining
`RDKit❌❌` matrix lines are the explicitly unselected bond-order-weighted and
atom-weighted matrix branches. The selected unweighted 2D and 3D computation,
cache lookup/store, resolved conformer identity, and operation-driven cache
invalidation are source-aligned. RDKit's internal path-matrix property is not
retained because no selected or modeled API observes it; the same path buffer
is computed for the distance result without retaining an otherwise
unobservable second quadratic allocation.

## Allocation And Complexity Audit

Pair enumeration is quadratic in atom count, matching RDKit, with one
environment allocation vector and the shared result projector. Hashing,
selectors, invariant access, and provenance writes have the same asymptotic
shape as the source. The Floyd-Warshall matrix construction is cubic for a
cold molecule and the resulting matrix is reused thereafter, matching RDKit's
repeated-call complexity. Cache entries are shared by logically read-only
clones but the lock containers remain clone-local; topology changes clear 2D
and 3D entries, while coordinate-only changes clear only 3D entries.

The pre-fix deterministic benchmark artifact is
[`tmp/atom-pair-bench/benchmark.json`](../../tmp/atom-pair-bench/benchmark.json),
SHA-256
`0d788c98bdc842931c270191faeaac0a3addbc6768a84176bfd8f0d4bbdb2a20`.
It records three-round medians for small, medium, and 80-atom molecules; 2D
and fixed-coordinate 3D modes; all four result forms; and a 3,072-molecule
mixed batch at one and 32 workers.

| Measurement | COSMolKit / RDKit |
|---|---:|
| Scalar 2D, all cases | 0.435–4.768; median 0.716 |
| Scalar 2D, 80 atoms | 3.224–4.768 |
| Scalar 3D, all cases | 0.409–0.517; median 0.486 |
| Mixed batch, one worker | 2.619–4.078; median 3.664 |
| Mixed batch, 32 workers | 1.046–1.744; median 1.181 |

The size-dependent 2D regression directly identified the missing cache source
lines: repeated Floyd-Warshall work dominated the 80-atom case. After the
source-aligned cache port, the same deterministic benchmark definition was
rerun. The post-port artifact is
[`tmp/atom-pair-bench/benchmark_post_cache.json`](../../tmp/atom-pair-bench/benchmark_post_cache.json),
SHA-256
`e99ef48ac6650a4133653f671b8eb25f47c9412b8b6ad6b4066aaa549e67aa68`.

| Post-port measurement | COSMolKit / RDKit |
|---|---:|
| Scalar 2D, all cases | 0.363–0.493; median 0.414 |
| Scalar 2D, 80 atoms | 0.458–0.493 |
| Scalar 3D, all cases | 0.360–0.478; median 0.390 |
| Mixed batch, one worker | 0.271–0.595; median 0.438 |
| Mixed batch, 32 workers | 0.241–0.336; median 0.311 |

All 24 scalar and eight batch measurements are now faster than the pinned
RDKit build. More importantly, the former atom-count-dependent long tail is
absent: the 80-atom 2D cases fell from 3.224–4.768 times RDKit to
0.458–0.493 times RDKit without changing any fingerprint result.

## Current Verdict

All observed AtomPair chemistry and public-output behavior is exact across the
committed corpora and the complete ChEMBL 37 corpus. The selected RDKit source
path has one production implementation, the distance-cache behavior and
invalidation rules are exercised under release strict contracts, focused
goldens remain exact after caching, and the post-port benchmark has no retained
performance gap. The final verdict is zero unexplained mismatch and zero
selected-boundary implementation gap; the support status can be promoted to
`supported_with_rdkit_parity`.
