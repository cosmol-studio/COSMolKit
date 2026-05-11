# COSMolKit Policy Invariants

This document defines the **policy invariants** of COSMolKit.

A policy invariant is a behavior rule that COSMolKit intentionally promises to users and developers. Unlike low-level structural invariants, which are required for a molecule object to be valid at all, policy invariants describe the design contract of this library.

These rules are not optional preferences. Any new feature, refactor, optimization, or API change should preserve these invariants unless the invariant is explicitly revised in this document.

This document is a **project policy**. It does not mean every policy is already fully covered by tests today. The current test coverage and required test mapping are tracked separately in [`testing_contract.md`](./testing_contract.md).

---

## 1. Core Principle

COSMolKit is not a blind clone of RDKit.

COSMolKit aims to provide a Rust-native, Python-friendly, batch-friendly, and AI-ready molecular toolkit, while using RDKit parity as the correctness floor for supported cheminformatics behavior.

The core development principle is:

> Correctness first, explicit behavior second, performance third.

Performance optimizations must not weaken correctness, parity guarantees, error visibility, object consistency, or value semantics.

---

## 2. Support Status Policy

Every feature should have a clear support status.

Public-facing feature status should use one of the following categories:

| Status | Meaning |
|---|---|
| `supported` | Fully implemented and covered by COSMolKit tests |
| `supported_with_rdkit_parity` | Implemented and tested against the pinned RDKit version |
| `preserved_only` | Parsed or stored, but not semantically interpreted |
| `experimental` | Available, but behavior may change |
| `unsupported` | Not implemented and should fail explicitly |

The status `unknown` may be used only for internal triage. It should not appear as a public user-facing feature guarantee.

In code, public feature status should be represented by a structured
`SupportStatus` and `FeatureSpec`, not by scattered comments. A feature
descriptor should include:

```text
feature name
feature category
support status
whether the feature is RDKit-parity-sensitive
short public documentation
```

Unsupported features must return a structured `UnsupportedFeature` error. They
must not panic through `unimplemented!`, return empty results, or produce
chemically meaningful-looking placeholder output.

### Policy invariant

Unsupported or partially supported behavior must not silently produce chemically meaningful-looking results.

If a feature is not supported, COSMolKit should prefer:

```text
UnsupportedFeature
```

over:

```text
silent fallback
best-effort guessing
partial behavior without warning
```

Operation support and feature support are related but not identical:

```text
FeatureSpec describes the public feature capability.
MoleculeOpSpec describes a registered molecule operation using that feature.
Support matrices connect features to operations.
```

Every public operation must be traceable to either a `FeatureSpec` or a clear
non-feature internal role.

---

## 3. Value-Style API Policy

COSMolKit public APIs should prefer value-style transformations.

For Python-facing molecule operations, the default behavior should be:

```python
mol2 = mol.operation()
```

rather than:

```python
mol.operation_in_place()
```

The original molecule should remain unchanged unless the API name clearly indicates mutation.

Examples of value-style operations:

```python
mol.with_hydrogens()
mol.without_hydrogens()
mol.with_kekulized_bonds()
mol.with_2d_coords()
mol.sanitize()
```

### Policy invariant

Public transformation APIs must not mutate the original molecule unless explicitly documented as mutable.

Example:

```python
mol0 = Molecule.from_smiles("CCO")
mol1 = mol0.with_hydrogens()

assert mol0.num_atoms() == 3
assert mol1.num_atoms() > 3
```

---

## 4. Explicit Mutation Policy

Mutation of existing molecules is allowed only through registered operations.

Examples are registered value-style operations such as
`mol.without_hydrogens(...)`, `mol.with_hydrogens(...)`, or future
operation-registry methods that return a new molecule.

### Policy invariant

Mutation-capable APIs must be explicit in naming, type, or usage pattern.

Avoid ambiguous public Python APIs such as:

```python
mol.add_hydrogens()
mol.remove_atom(3)
mol.set_bond_type(...)
```

unless the method name or documentation clearly states that it mutates the object in place.

Preferred names or patterns for explicit mutation include:

```python
mol.add_hydrogens_in_place()
mol.remove_atom_in_place(...)
```

---

## 5. Copy-on-Write Ownership Policy

COSMolKit may use shared internal storage such as `Arc<T>` and copy-on-write behavior.

Copy-on-write is an implementation strategy, not a public API promise.

The public contract is value semantics:

```text
clone + transform must not visibly mutate the source molecule
```

Internal storage sharing should be preserved when practical, but exact `Arc` block detachment behavior is an optimization target, not a stable public contract.

### Policy invariant

Copy-on-write must preserve value semantics.

This means:

```text
cloning and transforming a molecule must not alter the source molecule
```

Implementation-level sharing is allowed. User-visible aliasing bugs are not allowed.

---

## 6. Strong and Weak Topology Operation Policy

Topology-related operations should be classified as either **strong topology-changing operations** or **weak topology-state operations**.

### Strong topology-changing operations

Strong topology-changing operations change atom count, bond count, atom ordering, bond-table ordering, or atom/bond index mapping.

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

These operations may require atom and bond remapping.

### Weak topology-state operations

Weak topology-state operations keep atom and bond indices stable, but modify chemical graph state.

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

These operations may require cache invalidation, stereo updates, aromaticity updates, valence updates, or RDKit parity checks, but they do not necessarily require atom remapping.

### Policy invariant

Every topology-related operation must define whether it is strong or weak.

Strong operations must handle dependent data through remapping, recomputation, invalidation, explicit dropping, or explicit error.

Weak operations must invalidate or recompute all derived state affected by the changed chemical graph state.

The operation registry must also record the operation domain. Strong/weak is a
topology-index classification; some operations, such as coordinate generation,
do not change topology but still require operation-contract checks.

Recommended operation domains:

```text
Topology
Coordinate
```

Coordinate-domain operations should still declare affected derived state,
support status, parity policy, and invariant requirements.

---

## 7. Topology / Coordinate Alignment Policy

A molecule may contain topology, conformers, and properties as separate internal components.

Whenever atom indices change, coordinate rows must be updated, remapped, or removed consistently.

### Policy invariant

After any strong topology-changing operation:

```text
number of coordinate rows == number of atoms
```

for every conformer.

If atoms are removed, remaining atom coordinates must follow the correct atom mapping.

If atoms are reordered, conformer rows must be reordered using the same mapping.

No stale coordinate rows are allowed.

For weak topology-state operations that do not change atom indices, coordinates should remain unchanged unless the operation explicitly documents otherwise.

---

## 8. Atom and Bond Index Remapping Policy

Any operation that removes, inserts, or reorders atoms or bonds must handle dependent data carefully.

Dependent data may include:

```text
coordinates
stereochemistry
atom properties
bond properties
ring information
cached adjacency
cached perception results
drawing annotations
atom maps
```

### Policy invariant

Strong topology-changing operations must either:

1. correctly remap dependent data,
2. recompute dependent data,
3. explicitly invalidate dependent data,
4. explicitly drop unsupported dependent data,
5. or return an explicit error.

They must not leave stale indices behind.

For example, after removing an atom:

```text
No bond may reference the removed atom.
No stereo center may reference the removed atom.
No conformer row may remain for the removed atom.
No cache may continue using the old atom index space.
```

---

## 9. Cache Invalidation Policy

COSMolKit may cache derived molecular data for performance.

Examples of cacheable data:

```text
adjacency lists
ring information
aromaticity perception
valence perception
stereochemistry perception
distance matrices
drawing layouts
fingerprint intermediates
```

Caches are never allowed to become a hidden source of incorrect behavior.

### Policy invariant

Any operation that changes data underlying a cache must invalidate or refresh that cache.

Tests for cache-sensitive operations should follow this pattern:

```text
1. Build or access cache.
2. Perform topology-related operation.
3. Access cache-derived API again.
4. Compare with fresh recomputation or a known-valid result.
```

A stale cache is considered a correctness bug, not a performance bug.

---

## 10. Molecule-Level Property Policy

Molecule-level properties represent metadata about the molecule as a whole.

Examples:

```text
name
source_id
supplier_id
SDF data fields
dataset split
custom user properties
```

### Policy invariant

Single-source value-style transformations should preserve molecule-level properties by default.

Example:

```python
mol0 = Molecule.from_smiles("CCO")
mol0 = mol0.with_prop("ID", "mol_001")

mol1 = mol0.with_hydrogens()
mol2 = mol1.without_hydrogens()

assert mol2.prop("ID") == "mol_001"
```

Multi-source or splitting operations must define operation-specific property behavior.

Examples:

```text
combine(A, B)
fragment(mol)
merge_fragments(...)
split_components(...)
```

For such operations, the implementation must explicitly choose one of the following:

```text
preserve left-side props
preserve right-side props
merge props
namespace props
drop conflicting props
return an error on conflict
```

The chosen behavior must be documented and tested.

---

## 11. Atom-Level and Bond-Level Property Policy

Atom-level and bond-level properties are tied to specific atom or bond indices.

When topology changes, these properties require remapping.

SDF property lists are not ordinary molecule-level metadata once expanded.
The raw SDF data field should remain available as a molecule-level SDF field,
but any atom-list or bond-list interpretation must be represented as typed
property-list state or applied to the corresponding atom/bond properties with a
documented remapping policy.

### Policy invariant

Atom-level and bond-level properties must not silently attach to the wrong atom or bond after topology changes.

After atom or bond remapping, COSMolKit must choose one of the following:

```text
remap property
drop property
recompute property
return UnsupportedFeature or another explicit error
```

Silent misalignment is forbidden.

SDF property-list expansion must not silently ignore row-count mismatches. A
property list targeting atoms must match the atom table after RDKit-compatible
empty-value handling; a property list targeting bonds must match the bond table.
If the target count cannot be matched, strict parsing must return a structured
parse error and non-strict parsing must either preserve the raw field without
expansion or report an explicit unsupported path.

---

## 12. Stereochemistry Policy

Stereochemistry is high-risk and must be treated as first-class molecular state.

Stereochemical data may depend on:

```text
atom order
bond order
explicit hydrogens
implicit hydrogens
wedge/dash bonds
3D coordinates
SMILES ligand ordering
RDKit-compatible chiral tags
V3000 enhanced stereo COLLECTION groups
```

### Policy invariant

Any operation that may affect stereochemistry must either:

1. preserve stereochemistry correctly,
2. remap stereochemistry correctly,
3. recompute stereochemistry,
4. explicitly remove invalid stereochemistry,
5. or return an explicit error.

It must never leave stereochemistry pointing to invalid atom or bond indices.

For example, after `without_hydrogens()`:

```text
A tetrahedral center must not reference a deleted hydrogen atom.
```

If stereo preservation is not supported for a case, the operation should fail loudly or mark the behavior as unsupported.

Enhanced stereo groups are typed topology state. They must not be stored only as
string properties. RDKit collection ids may be preserved as metadata, but they
must not become COSMolKit atom, bond, or stereo row ids.

### Substance Group and Molfile Extended State Policy

Molfile/SDF SGroups, attach points, SGroup data fields, SGroup brackets,
SGroup cstates, and V3000 collection data are chemical or drawing semantics,
not generic molecule props.

Policy invariant:

```text
RDKit row, bookmark, sequence, and external ids must never be used as
COSMolKit row ids.
```

Parsers may keep temporary maps such as:

```text
rdkit_atom_bookmark -> AtomId
rdkit_bond_bookmark -> BondId
rdkit_sgroup_sequence_id -> SubstanceGroupId
```

Persistent molecule state must use COSMolKit ids for references and may preserve
RDKit ids only as explicit metadata fields such as `rdkit_sequence_id` or
`external_id`.

The following must be typed state when supported:

```text
SGroup atom membership
SGroup bond membership
SGroup parent atoms
SGroup parent group
SGroup attach points
SGroup cstate vectors
SGroup brackets and data display location
DAT SGroup field name/type/units/query metadata/data values
V3000 enhanced stereo groups
```

The following may remain in raw props only when they do not currently affect
behavior and the parser preserves enough information for a later typed upgrade:

```text
original RDKit labels
source sequence ids
vendor-specific display hints
unknown non-semantic metadata
```

Agents must not flatten SGroups, enhanced stereo, attach points, or property
list expansion into string props to make parsing easier. If the typed model is
missing or ambiguous, stop and add the typed state first, or ask the human
author to confirm a design exception.

---

## 13. RDKit Parity Policy

RDKit is the correctness oracle for supported RDKit-compatible behavior.

The RDKit version used for parity must be pinned and documented.

### Policy invariant

For features marked as `supported_with_rdkit_parity`, COSMolKit behavior should match the pinned RDKit version.

The operation registry must not use an ambiguous boolean such as
`rdkit_parity`. It should use a structured parity policy:

```text
NotApplicable
RequiredWhenSupported
RequiredNow
```

`RequiredWhenSupported` means the feature is intended to be RDKit-compatible,
but may currently be `unsupported`. Once its support status becomes
`experimental`, `supported`, or `supported_with_rdkit_parity`, it must have a
parity matrix entry and a parity profile.

`RequiredNow` means the operation currently claims RDKit-compatible behavior.
Missing parity corpus, golden output, or comparison schema is then a policy
violation.

`NotApplicable` is for COSMolKit-native features that do not have an RDKit
behavior oracle.

Parity-sensitive features include:

```text
SMILES parsing
SMILES writing
SDF/MolBlock parsing
SDF/MolBlock writing
hydrogen addition/removal
kekulization
sanitization
aromaticity behavior
formal charge handling
stereochemistry
Morgan fingerprints
distance geometry bounds
2D depiction behavior
```

When mismatch occurs, COSMolKit must not hide it by weakening tests.

Allowed responses to mismatch:

```text
fix COSMolKit
update expected behavior with clear RDKit source evidence
mark feature as unsupported
add known-failure case with explanation
```

Disallowed responses:

```text
relax assertion without explanation
skip only the failing field
add molecule-specific hacks
silently accept both outputs
guess chemistry behavior without source or tests
```

---

## 14. Known-Failure Policy

Known failures are allowed only as executable expected failures.

A known failure is not a skipped case and not permission to hide a mismatch.

### Policy invariant

A known-failure case must:

```text
still run
record a stable case_id
record the expected failure kind
record the reason
record the affected feature
fail if it unexpectedly passes
fail if the failure mode changes unexpectedly
```

Known-failure records should use common fields where possible and add domain-specific fields where needed.

Common fields:

```text
case_id
input or input_reference
feature
expected_failure_kind
reason
date_added or created_at
```

Domain-specific examples:

```text
operation
invariant
rdkit_version
branch
expected_error_kind
```

Known failures should be reviewed periodically. A fixed known failure should be removed from the known-failure list and converted into a normal regression test.

---

## 15. Golden Test Policy

Golden files are allowed, but they must be reproducible.

### Policy invariant

Golden files must be generated from deterministic scripts and a pinned reference environment.

A golden output should be traceable to:

```text
RDKit version
generation script
input corpus
output format
feature under test
```

This traceability may be provided by one or more of:

```text
generator script version assertions
golden README metadata
generated metadata file
CI environment pinning
file-level metadata
```

Each individual JSONL row does not need to repeat the RDKit version if the version is recorded elsewhere.

Golden outputs should not be manually edited unless explicitly documented.

If a golden output changes, the reason must be clear:

```text
RDKit version changed
COSMolKit bug fixed
test corpus changed
output format intentionally changed
```

---

## 16. Batch API Policy

Batch APIs are a core COSMolKit feature, not an afterthought.

Batch operations should be:

```text
ordered
traceable
parallel-safe
error-aware
deterministic when possible
```

### Policy invariant

Batch operations must preserve input order unless explicitly documented otherwise.

For example:

```python
batch = MoleculeBatch.from_smiles_list(smiles)
out = batch.with_hydrogens()

assert len(out) == len(batch)
assert out[i] corresponds to batch[i]
```

If an item fails, the failure should be traceable to the original input index.

Silent dropping of failed molecules is forbidden unless explicitly requested by the user through a documented error policy such as:

```python
errors="skip"
```

---

## 17. Error Visibility Policy

COSMolKit should prefer visible, structured errors over silent failure.

### Policy invariant

Failure modes should expose a meaningful error kind and enough context for debugging.

Preferred error context includes:

```text
error kind
operation name
input molecule identifier if available
batch index if available
atom index or bond index if relevant
underlying reason
```

Plain string wrappers are acceptable temporarily, but they must not erase the underlying error category when that category is known.

Preferred error categories include:

```text
UnsupportedFeature
ParseError
SanitizeError
ValenceError
StereoError
CoordinateError
ParityMismatch
InvalidMoleculeState
```

Avoid vague errors that provide no actionable context.

---

## 18. No Silent Chemistry Guessing Policy

Chemical behavior should not be guessed casually.

If the correct behavior is unclear, COSMolKit should prefer:

```text
unsupported
known mismatch
explicit TODO
parity investigation
```

over an unverified heuristic.

### Policy invariant

Heuristic chemistry behavior must be marked as heuristic and covered by tests.

A heuristic must not be presented as RDKit-parity behavior unless verified.

---

## 19. Performance Policy

Performance is important, but it must not compromise correctness.

Optimizations are encouraged when they preserve:

```text
molecule validity
RDKit parity
value semantics
error visibility
cache correctness
thread safety
```

### Policy invariant

An optimization must not change public behavior unless the behavior change is intentional, documented, and tested.

Unsafe code, cache reuse, in-place buffers, or parallelism must be justified by tests that protect correctness.

Implementation-level sharing and allocation behavior may be tested as performance targets, but should not be treated as public user-facing contracts unless explicitly documented.

---

## 20. Unsafe Code Policy

Rust `unsafe` may be used only when necessary.

### Policy invariant

Unsafe code must be localized, documented, and covered by tests.

Every unsafe block or unsafe wrapper should explain:

```text
why unsafe is needed
what assumptions must hold
which invariants protect those assumptions
why the operation is sound
```

For simple centralized FFI wrappers, a module-level soundness comment is acceptable if it clearly covers all unsafe calls in that wrapper.

Unsafe code must not be used to bypass borrow-checking inconvenience when a safe design is reasonable.

---

## 21. Thread Safety and Parallelism Policy

COSMolKit is expected to support batch and parallel workflows.

### Policy invariant

Parallel operations must not introduce shared mutable state bugs.

Parallel batch operations should satisfy:

```text
same input -> same output
same order -> same output order
same errors -> same error positions
```

Global mutable state should be avoided unless carefully synchronized and tested.

---

## 22. I/O Roundtrip Policy

I/O operations are part of the chemical correctness boundary.

Supported formats may include:

```text
SMILES
SDF
MolBlock V2000
MolBlock V3000
PNG/SVG depiction outputs
future graph formats
```

### Policy invariant

For supported I/O paths, roundtrip behavior should be tested.

Examples:

```text
SMILES -> Molecule -> SMILES
SDF -> Molecule -> SDF
Molecule -> SDF -> Molecule
Molecule -> graph -> Molecule, if supported
```

Roundtrip tests may be strict or lossy, but the expected level must be documented.

Topology-related operations should be covered by roundtrip tests when the operation affects data represented by that format.

---

## 23. AI-Native Export Policy

COSMolKit may provide AI-native exports such as:

```text
graph tensors
node features
edge features
coordinate arrays
torsion features
distance bounds
chirality constraints
diffusion model inputs
```

This policy applies when such exports are introduced or stabilized.

### Policy invariant

AI-native exports must document their schema.

A schema should define:

```text
node ordering
edge ordering
feature dimensions
feature meanings
coordinate convention
chirality convention
missing-value convention
version
```

Any breaking schema change should change the schema version.

Example schema names:

```text
cosmol-graph-v1
cosmol-torsion-v1
cosmol-diffusion-v1
```

---

## 24. Test Matrix Requirement for Topology-Related Operations

Every topology-related operation should be tested against the invariant matrix appropriate to its operation class.

Minimum required checks for strong topology-changing operations:

```text
graph index validity
coordinate row alignment
atom/bond remapping or explicit dropping
stereo index validity
cache invalidation
property preservation/remapping
COW source object unchanged
I/O roundtrip if applicable
RDKit parity if supported
```

Minimum required checks for weak topology-state operations:

```text
graph index validity
coordinate preservation unless documented otherwise
stereo validity if affected
cache invalidation
property preservation unless documented otherwise
COW source object unchanged
I/O roundtrip if applicable
RDKit parity if supported
```

A new topology-related operation should not be merged without the minimum invariant tests or a documented reason why a check is not applicable.

---

## 25. Documentation Policy

Public behavior must be documented at the same level as it is supported.

### Policy invariant

If a behavior is user-visible, it should be documented.

Documentation should clarify:

```text
whether the operation mutates or returns a new molecule
whether RDKit parity is expected
which unsupported cases may fail
how errors are reported
whether props are preserved
whether stereochemistry is preserved
whether coordinates are preserved
```

---

## 26. Compatibility Policy

COSMolKit may intentionally differ from RDKit when the design goal is different.

However, intentional differences must be explicit.

### Policy invariant

Any intentional difference from RDKit must be documented.

The documentation should explain:

```text
what differs
why COSMolKit differs
whether the difference is stable
whether users can request RDKit-compatible behavior
```

Unintentional differences in supported parity features are bugs.

---

## 27. Contributor Checklist

Before merging a change, check:

```text
[ ] Does this change mutate molecule topology or topology state?
[ ] Is it a strong topology-changing operation or weak topology-state operation?
[ ] If atom/bond indices change, are coordinates remapped or invalidated correctly?
[ ] Are caches invalidated or recomputed?
[ ] Are stereo references still valid?
[ ] Are molecule-level props preserved unless documented otherwise?
[ ] Are atom/bond-level props remapped, dropped, recomputed, or rejected explicitly?
[ ] Does this preserve value-style API behavior?
[ ] Does this preserve source-object immutability?
[ ] Does this need RDKit parity tests?
[ ] Does this need a known-failure case?
[ ] If known-failure is added, is it an executable xfail rather than a skip?
[ ] Does this introduce unsupported behavior that should fail loudly?
[ ] Does this change an AI-native export schema?
[ ] Does this change public API documentation?
[ ] Does this introduce unsafe code that needs a soundness comment?
```

---

## 28. Summary

COSMolKit policy invariants can be summarized as:

```text
Do not silently corrupt molecule state.
Do not silently guess unsupported chemistry.
Do not mutate value-style objects unexpectedly.
Do not leave stale indices, stale coordinates, stale stereo, or stale caches.
Do not weaken RDKit parity tests to make failures disappear.
Do not convert known failures into hidden skips.
Do not optimize by sacrificing correctness.
Do not expose unstable behavior as stable API.
```

The preferred development pattern is:

```text
define behavior
define support boundary
classify operation as strong or weak if topology-related
write invariant tests
write RDKit parity tests if applicable
implement feature
preserve explicit errors for unsupported cases
document public behavior
update testing_contract.md
```

COSMolKit should grow by making supported behavior more reliable, not by silently expanding ambiguous behavior.
