# COSMolKit Testing Contract

This document maps COSMolKit policy invariants to concrete testing requirements.

[`policy_invariants.md`](./policy_invariants.md) defines what COSMolKit promises. This document defines how those promises should be tested, which checks already exist, and which checks still need implementation.

This document is intentionally mechanical. It should help future contributors and coding agents decide what tests are required before merging a feature.

---

## 1. Status Labels

Use these labels when tracking test coverage.

| Status | Meaning |
|---|---|
| `implemented` | Test exists and is expected to run in CI |
| `implemented_failing` | Test exists and intentionally exposes an unresolved failure |
| `implemented_with_xfail` | Test exists and unresolved cases are tracked as executable expected failures |
| `partial` | Some related tests exist, but coverage is incomplete |
| `planned` | Required by policy, but not implemented yet |
| `not_applicable` | Policy does not apply to this feature or operation |
| `blocked` | Cannot be implemented until supporting API exists |
| `manual_review` | Requires human review or source-level investigation |

---

## 2. Test Categories

COSMolKit tests should be organized around these categories.

| Category | Purpose |
|---|---|
| Structural invariant tests | Ensure molecule objects remain internally valid |
| Policy invariant tests | Ensure COSMolKit design promises are preserved |
| RDKit parity tests | Compare supported behavior against pinned RDKit |
| Golden tests | Compare deterministic generated outputs |
| Known-failure tests | Track expected mismatches without hiding them |
| Batch tests | Ensure ordered, traceable, parallel-safe batch behavior |
| I/O roundtrip tests | Ensure parse/write cycles remain valid |
| Error tests | Ensure unsupported or invalid behavior fails visibly |
| Performance smoke tests | Catch severe regressions without overfitting benchmarks |
| Unsafe soundness review | Ensure unsafe code has local justification |

---

## 3. Policy-to-Test Mapping

| Policy area | Required tests | Current status | Notes |
|---|---|---|---|
| Support status | Public features classified as supported / parity / preserved-only / experimental / unsupported | partial | `unknown` is internal triage only |
| Support matrix | `SUPPORT_MATRIX` connects `FeatureSpec` entries to registered operations | planned | Matrix should be generated from `molecule_ops!`, not hand-maintained |
| Unsupported behavior | Unsupported cases return explicit errors | partial | Include tests that no public unsupported path panics or returns placeholder chemistry |
| Value-style API | Public transforms do not mutate source molecule | implemented | Covered by Python value-semantics tests and Rust molecule tests |
| Explicit mutation | Mutation only through clearly marked APIs | partial | Requires API naming/documentation checks |
| COW value semantics | Clone + transform does not alter original | implemented | Internal sharing checks are optimization-level tests |
| Strong topology ops | Atom/bond changes keep graph, coords, stereo, props, caches valid | partial | Unified runner still needed |
| Weak topology ops | Chemical graph state changes invalidate affected caches and preserve stable indices | partial | `sanitize` / `kekulize` need explicit classification |
| Coordinate alignment | Coordinate rows match atom count after atom remapping | partial | Remove-H covered; broader matrix needed |
| Atom/bond remapping | Atom/bond props and dependent state remapped/dropped/error explicitly | planned | Needs property-level corpus |
| Cache invalidation | Build cache, mutate, re-access, compare fresh recompute | planned | Needs test helper APIs |
| Molecule-level props | Single-source transforms preserve molecule props | partial | Extend across topology operations |
| Fragment/combine props | Multi-source/splitting ops define operation-specific prop policy | planned | Only required when these APIs exist |
| Stereo validity | Stereo indices remain valid or stereo is recomputed/dropped/error | partial | Needs unified stereo invariant runner |
| RDKit parity | Supported parity features match pinned RDKit | implemented_failing | Existing `rdkit_*_parity.rs` and golden generation scripts exist; unresolved parity failures must stay visible or become executable xfails |
| Parity policy | Registry uses `ParityPolicy`, not `rdkit_parity: bool` | planned | `RequiredWhenSupported` may exist before implementation; `RequiredNow` requires full parity assets |
| Known failures | Known failures run as executable expected failures | planned | Must not become skips; topology xfail schema has a first target record |
| Golden reproducibility | Golden outputs generated from pinned scripts/env | implemented | Metadata may be script-level, README-level, or CI-level |
| Batch order | Batch outputs preserve input order | implemented | Python batch tests cover this |
| Batch errors | Error index maps back to original input | implemented | Covered by `errors="keep"/"skip"/"raise"` style tests |
| Error visibility | Errors expose kind + context | partial | Python wrappers may still be string-heavy |
| No silent chemistry guessing | Unclear behavior marked unsupported / known mismatch / TODO | partial | Requires review discipline and tests |
| Performance | Optimizations do not change public behavior | partial | Mostly functional tests currently |
| Unsafe code | Unsafe blocks/wrappers have soundness comments | partial | Existing FFI wrappers should be reviewed |
| Parallelism | Parallel batch behavior is deterministic and order-preserving | partial | Batch tests exist; stress tests needed |
| I/O roundtrip | Supported formats roundtrip at documented strict/lossy level | partial | Existing parity/roundtrip tests; topology-op roundtrip matrix needed |
| AI-native export | Export schema documented and versioned | not_applicable | Applies when stable AI export APIs exist |
| Documentation | Public behavior documented at support level | partial | Needs review during API changes |
| RDKit differences | Intentional differences documented | partial | Should be linked from mismatch/known-failure records |
| InChI scalar APIs | Four public scalar APIs have exact source-defined official-C/RDKit tests, Python tests, structured diagnostics, and fail-closed unsupported-state coverage | implemented | `INCHI_API_PARITY_MATRIX` is the structured non-operation parity source. The official-C undefined `NormalizeAndCompare` allocation path maps to deterministic `allocation_failed` and is not an exact-parity claim. |

---

## 4. Topology Operation Classification

Every topology-related operation must be classified before tests are written.

The source of truth for operation classification is the generated operation
registry, not a hand-written test list. Test discovery should use:

```text
MOLECULE_OPS
SUPPORT_MATRIX
OPERATION_INVARIANT_MATRIX
PARITY_MATRIX
```

Every registered topology or coordinate operation must appear in
`OPERATION_INVARIANT_MATRIX`, including operations whose current support status
is `unsupported`. For unsupported operations, the invariant test should verify:

```text
operation returns structured UnsupportedFeature
source molecule remains unchanged
no partial molecule with placeholder chemistry is produced
operation metadata still declares kind/domain/support/parity
```

### Strong topology-changing operations

These change atom count, bond count, atom order, bond-table order, or atom/bond index mapping.

Examples:

```text
add_atom
remove_atom
add_bond
remove_bond
with_hydrogens
without_hydrogens
renumber_atoms
fragment
combine
```

Required checks:

```text
graph index validity
coordinate row alignment
atom mapping if available
bond mapping if available
stereo index validity
cache invalidation
molecule-level prop preservation unless documented otherwise
atom/bond-level prop remap/drop/error
source molecule unchanged
I/O roundtrip if applicable
RDKit parity if feature is parity-supported
```

### Weak topology-state operations

These keep atom and bond indices stable, but change chemical graph state.

Examples:

```text
kekulize
sanitize
set_bond_type
set_formal_charge
set_aromaticity
update_implicit_hydrogens
update_valence_state
```

Required checks:

```text
graph index validity
coordinate preservation unless documented otherwise
stereo validity if affected
cache invalidation
molecule-level prop preservation unless documented otherwise
source molecule unchanged
I/O roundtrip if applicable
RDKit parity if feature is parity-supported
```

---

## 5. Minimal Invariant Runner Contract

A unified invariant runner should eventually provide these helpers.

```rust
assert_graph_invariants(&mol);
assert_conformer_invariants(&mol);
assert_stereo_indices_valid(&mol);
assert_cache_consistency(&mol);
assert_property_consistency(&before, &after, mapping);
assert_source_unchanged(&before_snapshot, &before_after_op);
assert_io_roundtrip_if_applicable(&mol);
```

The runner should expose structured `check_*` helpers and assertion wrappers.

Preferred internal shape:

```rust
check_molecule_invariants(&mol, options) -> Result<(), InvariantFailure>;
check_topology_operation(&before, &after, &report, &spec) -> Result<(), InvariantFailure>;
```

Test-facing assertion wrappers may then be thin panic-on-error adapters:

```rust
assert_molecule_invariants(&mol, options);
assert_topology_operation(&before, &after, &report, &spec);
```

The runner may also expose a report-based API:

```rust
let report = check_molecule_invariants_with_report(&mol, options);
assert!(report.is_ok(), "{report:#?}");
```

or:

```rust
assert_molecule_invariants(&mol, options);
```

### Required behavior

The invariant runner should return structured failures, not only panic strings.

A useful failure should include:

```text
invariant name
operation name if available
molecule identifier if available
atom index or bond index if relevant
expected condition
actual condition
```

---

## 6. Strong Operation Test Template

Every strong topology-changing operation should have tests following this structure.

```rust
#[test]
fn invariant_without_hydrogens() {
    for case in topology_corpus() {
        let mol0 = Molecule::from_smiles(case.smiles).unwrap();

        let snapshot = MoleculeSnapshot::capture(&mol0);

        force_all_relevant_caches(&mol0);

        let mol1 = mol0.with_hydrogens().unwrap();
        let mol2 = mol1.without_hydrogens().unwrap();

        assert_source_unchanged(&snapshot, &mol0);

        assert_graph_invariants(&mol1);
        assert_graph_invariants(&mol2);

        assert_conformer_invariants(&mol1);
        assert_conformer_invariants(&mol2);

        assert_stereo_indices_valid(&mol1);
        assert_stereo_indices_valid(&mol2);

        assert_cache_consistency(&mol2);

        assert_io_roundtrip_if_applicable(&mol2);
    }
}
```

The exact APIs may differ. The required logic should remain the same.

---

## 7. Weak Operation Test Template

Every weak topology-state operation should have tests following this structure.

```rust
#[test]
fn invariant_kekulize() {
    for case in aromatic_corpus() {
        let mol0 = Molecule::from_smiles(case.smiles).unwrap();

        let snapshot = MoleculeSnapshot::capture(&mol0);
        let coords_before = capture_coords(&mol0);

        force_all_relevant_caches(&mol0);

        let mol1 = mol0.with_kekulized_bonds().unwrap();

        assert_source_unchanged(&snapshot, &mol0);
        assert_graph_invariants(&mol1);
        assert_conformer_invariants(&mol1);

        assert_coords_unchanged_if_expected(&coords_before, &mol1);
        assert_stereo_indices_valid(&mol1);
        assert_cache_consistency(&mol1);

        assert_io_roundtrip_if_applicable(&mol1);
    }
}
```

---

## 8. Cache Invalidation Test Contract

Cache invalidation tests must intentionally build the cache before mutation.

Bad test pattern:

```text
perform operation
access cache
```

This may miss stale-cache bugs.

Required test pattern:

```text
build/access cache
perform operation
access cache again
compare with fresh recomputation
```

Example:

```rust
#[test]
fn remove_atom_invalidates_adjacency_cache() {
    let mol0 = Molecule::from_smiles("CCO").unwrap();

    let _ = mol0.adjacency();

    let mol1 = mol0.remove_atom(1).unwrap();

    assert_adjacency_equals_fresh_recompute(&mol1);
}
```

If public cache APIs do not exist, provide internal `#[cfg(test)]` hooks.

---

## 9. Coordinate Remapping Test Contract

For strong operations that remove or reorder atoms, tests should verify coordinate mapping.

Required cases:

```text
remove hydrogens from molecule with 2D coords
remove hydrogens from molecule with 3D coords
renumber atoms with coords
remove atom with coords
fragment molecule with coords
```

Expected checks:

```text
coords.len() == atom_count
remaining atom coordinates match old coordinates via mapping
no stale coordinate rows remain
no NaN or infinite coordinate values are introduced
```

If an operation cannot preserve coordinates safely, it must explicitly drop coordinates or return an error.

Silent stale coordinates are forbidden.

---

## 10. Stereo Test Contract

Stereo-sensitive operations require both index-validity tests and parity tests when applicable.

Required structural checks:

```text
stereo center index < atom_count
all stereo ligand indices < atom_count
no duplicated ligand index in a tetrahedral ligand set unless representation explicitly allows it
bond stereo references valid atom/bond indices
removed atoms are not referenced by stereo state
```

Required corpus examples:

```text
F[C@](Cl)(Br)I
C[C@H](O)Cl
N[C@@H](C)C(=O)O
C/C=C/C
C/C=C\C
[H][C@](F)(Cl)Br
```

Required parity checks when supported:

```text
operation -> isomeric SMILES parity with RDKit
operation -> MolBlock stereo parity with RDKit if applicable
```

If stereo preservation is not supported for a case, the test should assert explicit unsupported behavior or executable known failure.

---

## 11. Property Test Contract

### Molecule-level props

Single-source transformations should preserve molecule-level props by default.

Required tests:

```text
with_hydrogens preserves molecule props
without_hydrogens preserves molecule props
kekulize preserves molecule props
sanitize preserves molecule props unless documented otherwise
with_2d_coords preserves molecule props
```

### Atom/bond-level props

Atom-level and bond-level props must not silently misalign.

Required outcomes after remapping operations:

```text
remap
drop
recompute
explicit error
```

Required tests once atom/bond props are exposed:

```text
remove atom with atom props
renumber atoms with atom props
remove bond with bond props
renumber bonds with bond props
add/remove hydrogens with atom props
```

---

## 12. RDKit Parity Test Contract

Any feature marked `supported_with_rdkit_parity` must have parity tests against the pinned RDKit version.

Operation-level parity obligations are declared through `ParityPolicy`:

```text
NotApplicable
RequiredWhenSupported
RequiredNow
```

`RequiredWhenSupported` means the operation should have a parity matrix profile
and may currently assert explicit `UnsupportedFeature` behavior. When support
moves beyond `unsupported`, the same matrix entry becomes the driver for real
RDKit comparison tests.

`RequiredNow` means parity corpus, golden generation, comparison schema, and
known-failure handling are required immediately.

Required parity metadata:

```text
pinned RDKit version
input corpus
generation script
expected output format
feature under test
```

Allowed mismatch responses:

```text
fix COSMolKit
update expected behavior with RDKit source evidence
mark unsupported
add executable known failure
```

Forbidden mismatch responses:

```text
skip failing field silently
relax assertion without explanation
molecule-specific hack
accept multiple incompatible outputs without reason
```

---

## 13. Known-Failure Test Contract

Known failures must be executable expected failures.

Known-failure records should use common fields where possible and add domain-specific fields where needed.

Common fields:

```text
case_id
feature
expected_failure_kind
reason
date_added or created_at
input or input_reference
link_to_issue_or_comment if available
```

Domain-specific fields are allowed when they make the failure precise.

Examples:

```text
operation
invariant
expected_error_kind
rdkit_version if parity-related
branch
```

Required behavior:

```text
case still runs
expected failure passes the xfail test
unexpected pass is reported
different failure kind is reported
known failure can be promoted to normal regression test after fix
```

Known failures must not be implemented as silent skips.

Suggested file:

```text
testdata/known_failures/rdkit_parity_known_failures.jsonl
testdata/known_failures/topology_invariants.jsonl
```

Example:

```json
{"case_id":"molblock_1065","feature":"molblock_parity","expected_failure_kind":"MolBlockMismatch","reason":"Known RDKit parity gap under kekulize=False","date_added":"2026-05-07"}
```

---

## 14. Golden Test Contract

Golden files should be generated, not hand-edited.

A golden test should be traceable to:

```text
generator script
input corpus
RDKit version or COSMolKit version
feature under test
output schema
```

Allowed metadata locations:

```text
README in golden directory
generator version assertion
metadata sidecar file
CI environment pin
file-level metadata
```

Required checks:

```text
golden files can be regenerated
regeneration is deterministic
CI uses pinned dependency versions
unexpected golden diff requires explanation
```

---

## 15. Batch Test Contract

Batch APIs must preserve order and traceability.

Required tests:

```text
input order == output order
errors keep original input index
errors="keep" preserves length and records failures
errors="skip" drops failures only by explicit user request
errors="raise" fails visibly
n_jobs produces same output as single-threaded path
progress bar does not alter results
```

Required failure metadata:

```text
input index
input string or molecule id when available
error kind or message
operation name when available
```

---

## 16. Error Visibility Test Contract

Unsupported or invalid behavior should fail visibly.

Required tests:

```text
unsupported feature returns explicit error
invalid molecule state returns explicit error
batch error records original index
parse failure contains input context
operation failure contains operation name when practical
```

Error messages may contain words like `failed`, but must not erase useful context.

Preferred structure:

```text
error kind + operation + context + reason
```

Python bindings may temporarily wrap Rust errors as standard Python exceptions, but tests should move toward preserving coarse error categories. At minimum, Python-facing failures should remain distinguishable for:

```text
unsupported feature
parse failure
sanitize failure
coordinate failure
batch validation failure
```

---

## 17. I/O Roundtrip Test Contract

Roundtrip tests must define their strictness.

Examples:

```text
strict topology roundtrip
lossy SMILES roundtrip
coordinate-preserving SDF roundtrip
property-preserving SDF roundtrip
depiction output golden comparison
```

Required tests for supported I/O:

```text
SMILES -> Molecule -> SMILES
SDF -> Molecule -> SDF
MolBlock -> Molecule -> MolBlock
Molecule -> SDF -> Molecule
```

Topology-related operations should be followed by roundtrip tests when the affected state is representable by the target format.

Example:

```text
from SMILES -> with_hydrogens -> without_hydrogens -> to_smiles -> from_smiles
from SDF -> without_hydrogens -> to_sdf -> from_sdf
```

---

## 18. Parallelism and Performance Test Contract

Functional correctness is required before performance.

Parallel tests should check:

```text
same input -> same output
same input order -> same output order
same failures -> same failure indices
n_jobs=1 and n_jobs>1 produce equivalent results
```

Performance tests should avoid brittle microbenchmark thresholds in normal CI unless the threshold is intentionally coarse.

Recommended levels:

```text
smoke test: catches extreme slowdown or accidental O(n^3)
benchmark: run manually or in scheduled CI
profile: used during optimization work
```

Performance optimizations must be accompanied by correctness tests.

---

## 19. Unsafe Code Test and Review Contract

Every unsafe block or unsafe wrapper should have a soundness explanation.

Required review items:

```text
why unsafe is needed
what assumptions must hold
which invariants protect those assumptions
whether inputs can violate the assumptions
what tests cover the wrapper
```

For centralized FFI wrappers, a module-level soundness comment is acceptable if it covers all unsafe calls in the module.

Suggested check:

```text
grep for unsafe blocks during review
confirm each unsafe block or wrapper has a SAFETY comment
```

---

## 20. AI-Native Export Test Contract

This section applies when AI-native exports are introduced or stabilized.

Required schema documentation:

```text
schema name
schema version
node ordering
edge ordering
feature dimensions
feature meanings
coordinate convention
chirality convention
missing-value convention
dtype/layout if tensor output
```

Required tests:

```text
schema shape stability
node order determinism
edge order determinism
feature dimension checks
roundtrip or reconstruction test if applicable
backward compatibility test for stable schema versions
```

Breaking schema changes must update the schema version.

---

## 21. Required Corpus Layout

Recommended corpus organization:

```text
testdata/
  topology/corpus/
    core.csv
    cow_small.csv
  sdf/corpus/
    sdf_v2000_cases.sdf
    sdf_v3000_cases.sdf
  smiles/corpus/
    known_edge_cases.smi
```

Topology CSV corpora should include stable case IDs and tags:

```text
case_id,smiles,tags
basic_ethanol_001,CCO,basic hydrogens coords props
stereo_tetra_h_001,[H][C@](F)(Cl)Br,stereo explicit_h remove_h
```

Use tags to select operation-specific subsets instead of relying on row numbers.

Recommended minimum cases:

```text
C
CCO
CC(=O)O
c1ccccc1
C1CCCCC1
[NH4+]
C[N+](C)(C)C
O=C([O-])C
[N+](=O)[O-]
c1ccncc1
c1[nH]cccc1
F[C@](Cl)(Br)I
C[C@H](O)Cl
N[C@@H](C)C(=O)O
C/C=C/C
C/C=C\C
[H][C@](F)(Cl)Br
```

Large external corpora may be used for stress tests, but small curated corpora should remain fast enough for normal CI.

---

## 22. Pull Request Checklist

Before merging a change, answer these questions.

```text
[ ] Does this change public behavior?
[ ] Does this change molecule topology or topology state?
[ ] Is the operation strong or weak?
[ ] Are graph invariants tested?
[ ] Are coordinate invariants tested?
[ ] Are cache invalidation rules tested?
[ ] Are stereo references tested if relevant?
[ ] Are molecule-level props preserved or documented otherwise?
[ ] Are atom/bond-level props remapped, dropped, recomputed, or rejected?
[ ] Does the source molecule remain unchanged for value-style APIs?
[ ] Does this require RDKit parity tests?
[ ] Does this require golden regeneration?
[ ] Does this add or update a known-failure case?
[ ] If known-failure is used, is it executable xfail rather than skip?
[ ] Does batch behavior preserve order and error indices?
[ ] Does this introduce or modify unsafe code?
[ ] Does this introduce unsupported behavior that should fail explicitly?
[ ] Does this require documentation updates?
[ ] Does this require updating policy_invariants.md or testing_contract.md?
```

---

## 23. Current Priority Roadmap

Recommended near-term priorities:

```text
1. Implement unified topology invariant runner.
2. Add strong/weak topology operation classification.
3. Add executable known-failure mechanism.
4. Add cache invalidation tests using force-cache-then-transform pattern.
5. Add stereo index validity checks after topology operations.
6. Add atom/bond property remap/drop/error tests once such props are public.
7. Add topology-operation I/O roundtrip matrix.
8. Add unsafe soundness comments for existing FFI wrappers.
```

---

## 24. Summary

The testing contract can be summarized as:

```text
Every policy should map to a test.
Every topology operation should map to an invariant matrix.
Every supported RDKit-compatible behavior should map to parity tests.
Every known mismatch should remain executable and visible.
Every unsupported behavior should fail explicitly.
Every public behavior should be documented and testable.
```

Policy without tests becomes a wish. Tests without policy become scattered checks. COSMolKit should keep both aligned.
