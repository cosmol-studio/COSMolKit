# RDKit Pattern Fingerprint Source Inventory

> Status: this inventory records the pre-implementation Step 3/5 source and
> reuse audit. Its `Missing` and `must add` rows are audit-time obligations,
> all now resolved by the implementation and evidence summarized in the
> [full-port validation report](rdkit_pattern_fingerprint_full_port_validation.md).

## Audit Boundary

This inventory fixes the source boundary for the COSMolKit Pattern fingerprint
port to RDKit `2026.03.1`, revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`. The selected production boundary
is the ordinary-molecule overload rooted at:

```text
Code/GraphMol/Fingerprints/PatternFingerprints.cpp:339-352
```

The port includes ordinary molecules and query-bearing molecules, both values
of `tautomericFingerprint`, all 13 active built-in SMARTS patterns, all valid
nonzero fingerprint widths, and the observable validation behavior of
`atomCounts` and `setOnlyBits`. The `MolBundle` overload is outside the selected
public boundary because COSMolKit has no public bundle model.

## Reachable Call Graph

```text
MolOps.cpp:666-695 wrapPatternFingerprint
  -> PatternFingerprints.cpp:339-352 PatternFingerprintMol(ROMol)
     -> PatternFingerprints.cpp:161-334 updatePatternFingerprint
        -> PatternFingerprints.cpp:57-80 pqs / pattern_flyweight
           -> PatternFingerprints.cpp:42-46 ss_matcher(string)
              -> SmilesParse.h:189-203 SmartsToMol
        -> FindRings.cpp:1131-1158 fastFindRings, conditionally
        -> QueryOps.cpp:892-919 isComplexQuery(Atom)
           -> QueryOps.cpp:769-795 _complexQueryHelper
        -> PatternFingerprints.cpp:140-152 isPatternComplexQuery
        -> PatternFingerprints.cpp:154-159 isTautomerBondQuery
        -> SubstructMatch.cpp:481-525 SubstructMatch vector overload
        -> hash/hash.hpp:211-219 gboost::hash_combine
        -> ExplicitBitVect.cpp:84-95 ExplicitBitVect::setBit
```

`QueryOps.cpp:728-766 isComplexQuery(Bond)` is not called directly by Pattern,
but it owns the general bond-query classification behavior that the crate must
consolidate into its sole query model. Pattern intentionally uses its narrower
`isPatternComplexQuery()` classifier instead.

## Exact Source Ledger

| Unit | Exact pinned range | Reachability and required behavior |
|---|---|---|
| `ss_matcher()` | `PatternFingerprints.cpp:41` | Default constructor is not used by the selected call graph; Rust needs no meaningless equivalent. |
| `ss_matcher(string)` | `PatternFingerprints.cpp:42-46` | Reachable. Compile a SMARTS string and assert the repository-owned pattern is valid. |
| `ss_matcher::getMatcher` | `PatternFingerprints.cpp:48-49` | Reachable. Return the immutable cached query without cloning it per molecule. |
| `pqs` | `PatternFingerprints.cpp:57-77` | Reachable. Preserve the 13 active strings, their order, and the empty sentinel semantics. Commented disabled patterns are source history, not active behavior. |
| `pattern_flyweight` | `PatternFingerprints.cpp:78-80` | Reachable. Preserve lazy, thread-safe, no-eviction sharing for the fixed table. |
| `detail::getAtomNumbers` | `PatternFingerprints.cpp:83-135` | Dead from both Pattern overloads in the pinned revision. Record this proof; do not add dead production code. |
| `isPatternComplexQuery` | `PatternFingerprints.cpp:140-152` | Reachable for every molecule bond. No query is simple; negation is complex; only description `BondOrder` is simple. |
| `isTautomerBondQuery` | `PatternFingerprints.cpp:154-159` | Reachable only after Pattern complexity is true. Recognizes exactly `SingleOrDoubleOrAromaticBond` and `SingleOrAromaticBond`. |
| `updatePatternFingerprint` preconditions | `PatternFingerprints.cpp:161-170` | Reachable. Reject zero width, short atom-count storage, and a mask whose width differs from `fpSize`. |
| Pattern cache materialization | `PatternFingerprints.cpp:172-184` | Reachable. Iterate in fixed order to the sentinel and collect immutable compiled queries. |
| Ring preparation | `PatternFingerprints.cpp:186-188` | Reachable. Run fast ring finding only if current ring information is below fast-or-better. |
| Query masks | `PatternFingerprints.cpp:190-209` | Reachable. Pattern overrides general `AtomNull` simplicity and marks it complex; classify Pattern bonds and conditional tautomer bonds. |
| Pattern matching | `PatternFingerprints.cpp:211-225` | Reachable. For each pattern use all pinned `SubstructMatchParameters` defaults except `uniquify=false` and `maxMatches=100000000`. |
| Occurrence-count bit | `PatternFingerprints.cpp:226-233` | Reachable once per returned match. Start from pattern index plus query atom and bond counts; combine `0xBEEF` cumulatively and set modulo width. |
| Atom hash and map | `PatternFingerprints.cpp:238-258` | Reachable. Iterate match pairs in matcher order, stop on a complex target atom, hash atomic numbers, and map query atom index to target atom index. |
| Bond and tautomer hash | `PatternFingerprints.cpp:259-315` | Reachable. Iterate query edges in source order, suppress complex target bonds except accepted tautomer query bonds, hash `u32::MAX` into the tautomer stream for tautomer/single/double/aromatic bonds, and hash normalized target bond order into the ordinary stream unless a tautomer query was encountered. |
| Final writes | `PatternFingerprints.cpp:317-331` | Reachable. Suppress all structural bits for rejected query features, suppress the ordinary structural bit after a tautomer query, and write the tautomer structural bit whenever tautomer mode remains admissible. |
| Ordinary facade | `PatternFingerprints.cpp:339-352` | Selected public root. Repeat validation, allocate a zeroed vector of exactly `fpSize`, call the sole core, and return it. |
| Bundle facade | `PatternFingerprints.cpp:355-374` | Excluded. It fingerprints each bundle member and intersects all member vectors. This is not a slice API and remains unavailable until a deliberate bundle model exists. |
| Public declaration/version | `Fingerprints.h:116-153` | Defaults are width 2048 and non-tautomeric mode. Preserve version metadata `1.0.0`. |
| Python wrapper | `MolOps.cpp:666-695,2653-2678` | Oracle for list conversion, defaults, errors, mutable list copy-back, version exposure, and overload boundary. COSMolKit retains project-native method names. |

## Active Built-In Pattern Order

The exact active sequence from `PatternFingerprints.cpp:57-77` is:

```text
 1 [*]~[*]
 2 [*]~[*]~[*]
 3 [R]~1~[R]~[R]~1
 4 [*]~[*](~[*])~[*]
 5 [R]~1[R]~[R]~[R]~1
 6 [*]~[*]~[*](~[*])~[*]
 7 [R]~1~[R]~[R]~[R]~[R]~1
 8 [R]~1~[R]~[R]~[R]~[R]~[R]~1
 9 [R](@[R])(@[R])~[R]~[R](@[R])(@[R])
10 [R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])
11 [*]~[R](@[R])@[R](@[R])~[*]
12 [*]~[R](@[R])@[R]@[R](@[R])~[*]
13 [*]
```

The one-based `pIdx` is behaviorally significant and participates directly in
both occurrence and structural hashes.

## Parameter And Error Semantics

The C++ facade and core both enforce:

```text
fpSize != 0
atomCounts == null || atomCounts.size >= molecule atom count
setOnlyBits == null || setOnlyBits.width == fpSize
```

The Python wrapper rejects a short `atomCounts` list before calling the core.
It copies valid list values into a C++ vector and copies them back afterward.
In the pinned ordinary overload, neither `atomCounts` values nor `setOnlyBits`
bits are read or modified after validation. Therefore valid arguments are
observably inert, despite the stale header documentation at
`Fingerprints.h:126-137` claiming count updates and result masking. Exact
source reproduction must preserve implementation behavior, not implement that
unrealized comment.

`fpSize=0`, a short count vector, and a wrong-width mask are errors. SMARTS
failure for the repository-owned constants and a null cached matcher are
internal invariants. Matcher errors retain their existing structured search
error boundary. No molecule-specific fallback or row exclusion is allowed.

## Dependency Ownership

| Dependency | Pinned source range | Ownership decision |
|---|---|---|
| General bond-query complexity | `QueryOps.cpp:728-766` | Consolidate into the existing query module; Pattern uses its own narrower wrapper. |
| Atom-query recursion | `QueryOps.cpp:769-795` | Consolidate into the existing query module with sibling-carried `hasAtNum`. |
| General atom-query complexity | `QueryOps.cpp:892-919` | Consolidate into the existing query module; Pattern adds the explicit `AtomNull` override. |
| SMARTS compilation | `SmilesParse.h:189-203` plus reachable parser implementation | Reuse the existing sole SMARTS parser. |
| Ring preparation | `FindRings.cpp:1131-1158` and `_DFS` at `1095-1129` | Reuse the existing sole fast-ring implementation. |
| Substructure matching | `SubstructMatch.cpp:481-525`; defaults in `SubstructMatch.h:42-86` | Reuse the existing sole matcher with only the two Pattern overrides. |
| Hashing | `hash/hash.hpp:211-219` | Reuse the existing wrapping 32-bit gboost hash combine. The `-1` source argument becomes `u32::MAX`. |
| Bit insertion | `ExplicitBitVect.cpp:84-95` | Reuse `Fingerprint`; repeated and colliding writes remain idempotent. |

The selected call graph does not reach subgraph enumeration despite the
includes at `PatternFingerprints.cpp:17-18`, nor does it reach random-number or
sorting behavior from the unused includes. They do not justify new local
implementations.

## Upstream Regression Inventory

The source port must retain the following upstream evidence:

| Upstream test | Exact range | Covered behavior |
|---|---|---|
| GitHub 151 | `Fingerprints/test1.cpp:2706-2778` | Four positive substructure-screening containment cases, including fused ring systems. |
| GitHub 258 | `Fingerprints/test1.cpp:2861-2962` | Five query-molecule screening cases with aromatic, atomic-number, negation, degree, and explicit-H queries. |
| Multithreaded Pattern FP | `Fingerprints/test1.cpp:2964-3049` | Concurrent initial compilation and repeated equality with serial reference vectors. |
| GitHub 879 | `Fingerprints/test1.cpp:3427-3457` | Single-atom and disconnected single-atom fragments set Pattern bits. |
| GitHub 1496 | `Fingerprints/test1.cpp:3459-3509` | Degree-zero atom query screening containment. |
| GitHub 2051 | `Fingerprints/catch_tests.cpp:33-51` | Wildcard one- and two-atom query containment. |
| GitHub 2614 | `Fingerprints/catch_tests.cpp:53-66` | Match count exceeds the old 1,000-result default without losing screening bits. |
| Bundle intersection | `Fingerprints/catch_tests.cpp:112-132` | Source evidence for the excluded bundle overload only; not a production obligation for this boundary. |
| Java query/bundle tests | `JavaWrappers/gmwrapper/src-test/org/RDKit/FingerprintsTests.java:146-165` | Additional ordinary/query screening evidence and excluded bundle evidence. |

The TautomerQuery consumers at `TautomerQuery.cpp:290-299` prove that both the
template and target paths select `tautomericFingerprint=true`. Generalized
substructure, substructure-library, PostgreSQL, MinimalLib, reaction, Pandas,
and synthon consumers are downstream uses; they do not add another Pattern
algorithm.

## Reachability And Exclusion Proofs

`detail::getAtomNumbers()` has no caller in the pinned source tree. The only
occurrences are its definition and declaration context in
`PatternFingerprints.cpp`; neither `updatePatternFingerprint()` nor either
facade references it. Adding it to production would create dead duplicate
query interpretation.

The `MolBundle` overload is a distinct semantic operation: it computes one
Pattern vector per member and intersects them. COSMolKit currently lacks the
bundle ownership and query semantics required to expose that behavior. A slice
or batch is not equivalent because batch APIs preserve one result per input.
The overload therefore remains an explicit documented exclusion, not an
unsupported branch hidden in the ordinary facade.

## Step-3 Gap Result

At the Step 3 audit point, no production Pattern core existed. The source
boundary was nevertheless closed and implementable with existing SMARTS,
matching, ring, hash, and bit-vector facilities, subject to the separate
repository-reuse audit. Every
production line in `PatternFingerprints.cpp:41-80,140-352` is assigned either a
reachable port obligation or a precise Rust-language omission. Lines `83-135`
are proven dead, and lines `355-374` are explicitly outside the public model.

## Repository Reuse And Duplication Audit

This section records the Step-5 audit of the current COSMolKit tree. The reuse
decisions below are binding for the Pattern implementation: a missing adapter
may be added at the narrow owner boundary, but the underlying facility must not
be copied into a Pattern-local substitute.

### Query Model And Classifiers

The sole typed query tree is `search/query.rs::QueryNode<T>`, with atom and bond
leaves represented by `AtomQueryPredicate` and `BondQueryPredicate`. Atom and
bond values expose immutable `query()` accessors. There is no need for a
description-string query model.

`complex_atom_query_helper()` and `is_complex_atom_query()` already reproduce
`QueryOps.cpp:769-795,892-919` with verbatim source anchors and two-axis
markers. They are currently private and are not exhaustively tested as a
standalone classifier boundary. Step 7 must make the narrow classifier callable
inside the crate and complete focused tests; it must not create a Pattern copy.

The general `is_complex_bond_query()` also exists in the same module and must
remain the only general bond classifier. Pattern's source-specific narrower
classifier is a separate semantic unit because it considers
`SingleOrAromaticBond` complex while the general source function does not.

Gap result:

```text
reuse QueryNode and typed predicates
consolidate visibility/tests for existing general classifiers
add only the two Pattern-specific bond classifiers
do not introduce description strings or duplicate recursion
```

### SMARTS Compilation

`search/smarts_parse.rs::mol_from_smarts()` is the sole parser and accepts
`SmartsParseParams`. It is re-exported by `search/mod.rs`; the public parse
boundary is already used for source-owned fingerprint pattern tables.

`properties/fingerprint.rs::SsMatcher` is the existing fingerprint-owned
compiled SMARTS wrapper. It stores one typed `Molecule` and returns an immutable
reference. `cached_default_feature_matchers()` and
`cached_rdkit_maccs_pattern_matchers()` demonstrate the repository's
`OnceLock<Result<Vec<_>, FingerprintError>>` compile-once convention.

Pattern must reuse `SsMatcher` and this cache convention. It needs one distinct
fixed table/cache because the source table is distinct, but no second wrapper,
parser, or generic SMARTS cache abstraction is justified.

### Substructure Matching

`search/substruct.rs` owns the only VF2 matcher. `SubstructMatchParams` contains
all pinned defaults from RDKit plus the required Pattern overrides.
`try_get_substruct_matches_with_params()` returns the ordered
`SubstructMatchResult` list with structured errors; each result already contains
query-to-target atom and bond mappings.

Pattern must use:

```text
SubstructMatchParams::default()
params.uniquify = false
params.max_matches = 100_000_000
try_get_substruct_matches_with_params()
```

The convenience `get_substruct_matches_with_params()` erases matcher errors and
is therefore not appropriate for a parity-sensitive fingerprint core. No local
match loop, VF2 adapter, atom map reconstruction, deduplication, sorting, or
match cap is allowed.

### Ring Information

`chemistry/rings.rs::fast_find_rings()` and its parts/context helpers are the
sole port of `FindRings.cpp:1095-1158`. It returns typed `RingInfo` and carries
complete line-level source anchors. Query matching already builds its ring
context through the same ring facilities.

Pattern must not add ring DFS or ring predicates. COSMolKit's public
fingerprinting contract is non-mutating, so a required fast-ring result may be
computed as local derived state instead of reproducing RDKit's const-molecule
cache side effect. This changes storage only, not fingerprint or matcher
behavior; the benchmark step must verify that it does not create an unexplained
hot-path regression.

### Hash And Vector Infrastructure

`properties/fingerprint.rs::hash_combine()` is the sole wrapping `u32`
implementation of `gboost::hash_combine` and has the correct shift/add/XOR
semantics. Pattern must call it directly, including `u32::MAX` for source `-1`.

`properties/fingerprint.rs::Fingerprint` is the sole explicit bit-vector type.
It stores packed `Vec<u64>`, reports ordered on-bits, and has idempotent
construction from on-bit indices. It currently lacks an internal mutable
`setBit` equivalent. Building a temporary set/vector and reconstructing the
fingerprint would duplicate vector work and worsen the source hot path. The
narrow owner-level gap is therefore a crate-private zeroed constructor and
idempotent checked/known-in-range bit insertion on `Fingerprint`; Pattern must
not define another bit vector.

### Fingerprint Module And Public Rust API

The existing main module owns `Fingerprint`, `FingerprintError`, shared SMARTS
matchers, Morgan, MACCS, and topological fingerprints, with narrow submodules
for larger families such as AtomPair. `Molecule` exposes method-oriented
fingerprint methods that delegate to the sole core.

Pattern should live in one narrow `properties/fingerprint/pattern.rs` module,
re-export its parameter/facade types through `fingerprint.rs`, and add one
`Molecule::pattern_fingerprint()` method. There is no existing Pattern
production implementation to preserve or wrap.

`FingerprintError` already represents zero width, invalid arguments, SMARTS
invariants, unsupported options, and parallel execution. Matcher errors need a
lossless Pattern mapping; a new structured variant is preferable to a string
fallback if the existing error cannot be converted directly.

### Batch And Python Facilities

`properties/batch.rs` owns deterministic Rayon-backed ordered mapping, progress,
worker-count validation, and indexed error preservation. Existing fingerprint
families expose scalar and list methods through this facility. Pattern needs a
thin list facade that calls `Molecule::pattern_fingerprint()` for each item; no
batch core or alternate Pattern algorithm is allowed.

`python/src/lib.rs` owns the `Molecule`, `MoleculeBatch`, `Fingerprint`, and
fingerprint error conversion boundaries. Pattern scalar and ordered batch
methods must follow those current names and return types, delegate to Rust, and
leave stub generation to the repository generator. They must not expose
`SsMatcher`, query classifiers, or mutable packed storage.

### Golden And Corpus Workflow

The canonical golden entrypoint is `tools/testdata/rdkit/generate_all.py`.
`GeneratorSpec` registers generator script, output schema, domain/suites,
profile inputs, generator dependencies, deterministic sharding, and manifest
identity. `_expected_schema.py` owns top-level JSONL validation.

Existing committed corpus identities are:

```text
testdata/smiles/corpus/smiles_small.smi
testdata/smiles/corpus/smiles_5000.smi
```

Focused fingerprint inputs live under
`testdata/fingerprint/fixtures/rdkit/`. Generated oracle data belongs under
`testdata/fingerprint/expected/rdkit/<profile>/` and is resolved and checksum-
validated by `cosmolkit-test-support`; temporary parity data must not be loaded
from `tmp/` by committed tests.

The Pattern workflow needs one profile JSON, one generator, one expected schema,
focused/small/5000 registrations, and exact Rust consumers. It can follow the
AtomPair/topological registration shape but must not edit another family's
profile or overload its schema.

### ChEMBL 37 Workflow

`dev/tools/chembl_parity/` already provides pinned-version checks, deterministic
record iteration, checksummed resumable sharding, mismatch artifacts, profile
validation, and complete-run accounting. `audit_fingerprints.py` is the sole
fingerprint audit runner and currently handles topological/Avalon, AtomPair,
and Topological Torsion modes. `profiles/complete.json` is the complete-run
registry.

Pattern must be added as another exact mode/profile in that runner, comparing:

```text
exact n_bits
complete ordered on_bits
reference and COSMolKit error outcomes
processed, parse, compared, error, and mismatch counts
```

No sampling, row allowlist, silent parse filtering, or second ChEMBL runner is
needed.

### Concrete Gap Matrix

| Facility | Current status | Pattern action | Duplicate implementation forbidden |
|---|---|---|---|
| Query tree | Complete typed owner | Reuse | Yes |
| General query complexity | Implemented, narrow visibility/tests incomplete | Consolidate and test | Yes |
| Pattern bond classifiers | Missing | Add in Pattern module | Not applicable; source-specific |
| SMARTS parser | Complete | Reuse | Yes |
| Compiled matcher wrapper | Complete | Reuse | Yes |
| Pattern table/cache | Missing | Add one fixed `OnceLock` cache | Yes, exactly one table/cache |
| VF2 matcher | Complete | Reuse exact params/result order | Yes |
| Fast ring finding | Complete | Reuse | Yes |
| `hash_combine` | Complete | Reuse | Yes |
| Packed `Fingerprint` | Complete | Add narrow internal mutation primitives | Yes |
| Pattern core/facade | Missing | Add one implementation | Yes |
| Molecule method | Missing | Add delegate | Yes |
| Batch method | Generic engine complete | Add thin delegate | Yes |
| Python scalar/batch | Missing | Add thin bindings | Yes |
| Golden workflow | Complete framework | Register Pattern artifacts | Yes |
| ChEMBL workflow | Complete framework | Extend fingerprint runner/profile | Yes |

The audit found no historical Pattern implementation that needs migration. The
main risk is accidental duplication of mature infrastructure, not reconciliation
with an old Pattern branch.
