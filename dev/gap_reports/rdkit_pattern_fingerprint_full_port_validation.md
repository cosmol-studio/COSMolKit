# RDKit Pattern Fingerprint Full-Port Validation

## Validated Boundary

This report closes the source-reproduction audit for the ordinary-molecule
Pattern fingerprint in RDKit `2026.03.1`, revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`. The validated boundary contains:

```text
PatternFingerprintMol(const ROMol &)
all 13 active built-in Pattern SMARTS queries, in source order
ordinary and query-bearing molecules
tautomericFingerprint=false and true
all valid nonzero fingerprint widths
the observable validation and inert-value behavior of atomCounts/setOnlyBits
Rust scalar and ordered batch APIs
Python scalar and ordered batch APIs
```

`PatternFingerprintMol(const MolBundle &)` is not part of this claim. A bundle
fingerprint is the intersection of member fingerprints and is not equivalent
to the project's ordered batch API. COSMolKit has no public bundle model, so
the overload remains an explicit model-boundary exclusion instead of a hidden
fallback or a duplicate slice API.

## Source-Line Closure

The exact call graph and line ownership are fixed in
[`rdkit_pattern_fingerprint_source_inventory.md`](./rdkit_pattern_fingerprint_source_inventory.md).
Every production unit in the selected boundary has one disposition:

| Source unit | Closure |
|---|---|
| `ss_matcher`, `pqs`, `pattern_flyweight` | One immutable, lazy, thread-safe cache compiles the exact 13 SMARTS strings once through the canonical SMARTS parser. |
| `isComplexQuery(Atom/Bond)` | Shared typed query classifiers are source-marked in the canonical query module; Pattern adds only its source-specific bond classifiers. |
| `isPatternComplexQuery`, `isTautomerBondQuery` | Ported as typed root-query tests with exact negation and description-equivalent behavior. |
| `updatePatternFingerprint` | Ported line-by-line, including validation, query masks, match parameters, occurrence hashes, structural hashes, tautomer hashes, suppression branches, modulo, and writes. |
| `PatternFingerprintMol(ROMol)` | One project-native facade validates, allocates the result, and delegates to the sole core. |
| `fastFindRings`, `SmartsToMol`, `SubstructMatch`, `hash_combine`, `ExplicitBitVect::setBit` | Reuse the repository's sole ring, parser, VF2, 32-bit hash, and packed fingerprint implementations. |
| `detail::getAtomNumbers` | Proven unreachable from both pinned Pattern overloads; intentionally not copied into production. |
| `PatternFingerprintMol(MolBundle)` | Explicitly excluded pending a deliberate public bundle model. |

The production implementation carries verbatim C++ review anchors and
two-axis markers in `properties/fingerprint/pattern.rs`. Architecture guards
prove one Pattern table/cache, one compilation site, one core, and thin
scalar/batch delegates; they also reject Pattern-local replacements for the
shared query classifier, matcher, hash, and `Fingerprint` representation.

## Branch, Option, Output, And Error Closure

Focused unit and integration tests exercise absent and nested atom/bond query
trees, negation, Pattern's `AtomNull` override, both recognized tautomer bond
queries, every built-in pattern, zero/repeated/overlapping/symmetric matches,
the exact non-unique match order, occurrence-count seed evolution, atom and
bond hash evolution, aromatic normalization, query suppression,
`u32::MAX` tautomer hashing, collisions, and single-atom patterns.

The public and golden matrices cover default, tautomeric, widths `1`, `7`,
`127`, and `2048`, zeroed and pre-populated atom counts, and all-zero, all-one,
and sparse set-only masks. They reproduce the pinned ordinary overload's
observable behavior: valid count and mask values are inert after their size
checks. They do not implement the stale header comment that claims mutation or
masking that the source implementation does not perform.

The following error boundaries are exact and executable:

```text
fpSize == 0                         -> structured empty-fingerprint error
atomCounts shorter than atom count -> structured invalid-argument error
setOnlyBits width != fpSize         -> structured invalid-argument error
repository-owned SMARTS failure     -> internal Pattern invariant error
substructure matching failure       -> structured Pattern fingerprint error
```

Complete vectors compare exact width and the complete ordered on-bit set. No
tolerance, row allowlist, output-field omission, fallback, or accepted
alternative vector is used. Screening regressions additionally require exact
vectors and all probe bits to be contained in every reference-positive target.

## Corpus Evidence

| Evidence set | Rows | Profiles or branches | Exact result |
|---|---:|---:|---|
| Focused Pattern fixtures | 18 | 11 branches per row | All output and error branches pass. |
| `smiles_small` | 152 | 11 branches per row | Every row and branch passes. |
| `smiles_5000` | 5,000 | 11 branches per row | Every row and branch passes. |
| ChEMBL 37 | 2,897,819 processed; 2,897,804 mutually parsed; 15 mutually rejected | 10 complete profiles | 28,978,040 exact fingerprint comparisons, zero mismatch. |

The ChEMBL execution completed all `128/128` deterministic shards. Its corpus
source SHA-256 is
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.
The external run manifest and aggregate remain outside tracked repository
content at:

```text
/home/datahouse-raid5/wjt/chembl37-pattern-full/manifest.json
/home/datahouse-raid5/wjt/chembl37-pattern-full/pattern/aggregate.json
```

The committed focused/small/5,000 expected data are generated by the pinned
repository workflow and resolved through checksummed manifests; committed
tests do not depend on `tmp/`.

## API, Composition, And Tooling Evidence

Rust and Python expose method-oriented scalar and ordered batch APIs that
delegate to the same core. Tests cover exact defaults, invalid arguments,
source immutability, stable order, repeated calls, shared values, worker-count
validation, concurrent calls, and interleaving with Avalon, Morgan, AtomPair,
RDKFingerprint, MACCS, and Topological Torsion fingerprints. No Pattern call
mutates molecule state or leaks options into another call.

The generated Python stub contains the scalar and batch signatures. Focused
runtime/stub tests pass (`8 passed` for the focused Pattern Python suite). The
repository-wide `basedpyright python/tests python/examples` run reports
`0 errors, 510 warnings`; its nonzero status comes from pre-existing warnings,
so this is not a claim that the warning gate is clean.

The strict core release suite passed with `3338 passed`, `46 ignored`, and all
integration targets passing; doctests passed `7/7`. Pattern golden, screening,
architecture, public API, and composition targets are included in that result.

## Allocation And Complexity Closure

The implementation has the same asymptotic outer structure as the pinned
source: one process-lifetime compilation of 13 queries, one target match
context per fingerprint call, 13 non-unique VF2 searches, one pass over each
returned mapping, and constant-time wrapping hashes and packed-bit writes.
There is no per-molecule SMARTS reparse, query clone, target clone, duplicate
ring engine, alternate matcher, or second bit-vector conversion path.

The representative benchmark covers path-rich, ring-rich, query-bearing,
disconnected, and large inputs. All five 2,048-bit vectors match RDKit exactly,
and warm calls observe cache reuse. The measured median warm-call
COSMolKit/RDKit ratios are:

| Case | Ratio |
|---|---:|
| Path-rich | 2.120 |
| Ring-rich | 3.038 |
| Query-bearing | 1.508 |
| Disconnected | 2.525 |
| Large | 3.104 |

Therefore this validation does **not** claim performance parity. The current
implementation retains a measured constant-factor matcher/allocation gap while
preserving the same algorithmic structure. The external benchmark artifact is
`tmp/pattern-fingerprint-bench/benchmark.json`, SHA-256
`c370d6f7bbc81bbd50f4b9c36ba6c4dd770b2427ebba46d15a317a04e5b4ac6c`;
machine artifacts are intentionally not tracked.

## Validation Conclusion

The ordinary-molecule Pattern fingerprint behavioral boundary is closed for
the pinned source revision: source lines, reachable branches, valid options,
validation errors, exact outputs, query behavior, Rust/Python APIs,
composition, architecture, and complete ChEMBL 37 corpus evidence have no
unexplained mismatch. The claim excludes `MolBundle` and excludes performance
parity. The only measured implementation gap in this report is the explicit
constant-factor runtime difference above; it does not conceal a behavioral or
asymptotic divergence.
