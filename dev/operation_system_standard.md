# COSMolKit Operation System Standard

This document defines the operation-system standard for COSMolKit.

It is a binding project rule and target architecture standard. New core operation code must follow this standard unless the human author explicitly approves an exception.

Implementation gaps must be tracked in plans or checklists. This document describes the desired operation system, not temporary implementation status.

---

## 1. Purpose

The operation system exists to guarantee:

```text
value semantics
explicit mutation authority
copy-on-write efficiency
strong/weak topology discipline
derived-state correctness
source-auditable behavior
strict development-time checking
low-overhead release execution
```

It prevents topology edits, cache invalidation, stereo updates, property remapping, and topology mapping from becoming scattered handwritten conventions.

---

## 2. Policy Relationship

This standard implements the project-level invariants in `policy_invariants.md`.

It must preserve:

- public transforms do not visibly mutate the source molecule
- mutation is allowed only through registered operations
- copy-on-write is internal, not a public API promise
- topology operations are classified as strong or weak
- unsupported behavior must fail explicitly, not produce plausible placeholders

Correctness and explicit behavior take priority over convenience or performance.

---

## 3. Scope

This standard applies to molecule operations that mutate or derive:

```text
topology
coordinates
atom or bond state
molecule properties
derived caches
stereo state
topology mapping
```

Truly read-only accessors are outside this system.

BioStructure operations follow the same design through `bio_structure_ops!` and `BioOpParts`.

---

## 4. Operation Model

Every registered molecule operation has five layers:

```text
registry declaration
macro-generated public wrapper
OpParts capability object
domain implementation body
contract-checked finish step
```

Responsibilities are separated:

- the registry declares allowed and required behavior
- the wrapper constructs `OpParts` and calls the implementation
- `OpParts` controls copy-on-write mutation, remap, mapping, trace, and invalidation
- the operation body implements chemistry or domain logic
- `finish()` validates the contract in strict builds and returns the result

No layer may silently take over another layer’s responsibility.

---

## 5. Registry

All public topology-related or coordinate-related molecule operations must be registered through `molecule_ops!`.

The registry is the source of truth for:

```text
MOLECULE_OPS
SUPPORT_MATRIX
OPERATION_INVARIANT_MATRIX
PARITY_MATRIX
```

Do not maintain parallel handwritten operation lists, support matrices, invariant matrices, or parity tables.

Any public mutation-capable API must be traceable to a registered `MoleculeOpSpec`, unless it is explicitly documented as a non-operation internal helper.

---

## 6. Required Spec Fields

Each registered operation must define:

```text
method
impl_fn
domain
kind
access
may_mutate
derived_effects
requires_mapping
allows_noop
feature
parity
io_roundtrip
invariant_profile
parity_profile when parity is required
```

Strong topology operations must also define their migration surface, including fields such as:

```text
topology_edit
auto_remap
```

These fields are behavioral contract inputs, not documentation-only metadata.

`access` is the authoritative block-capability declaration. For every molecule
block that the operation may touch, the registry must declare exactly one
access mode:

```text
none
read
write
```

`may_mutate` is the write subset of `access` and exists for generated matrices,
strict checks, and compatibility with current operation metadata. It must not
grant additional authority beyond `access`.

A block must not be exposed through both an independent read capability and a
write capability. If an operation needs to inspect a write-owned block before
modifying it, it must begin the write-owned block and read from the same local
owned working value.

---

## 7. Strong And Weak Operations

### Strong topology operations

An operation is strong if it changes any of:

```text
atom count
bond count
atom ordering
bond ordering
atom identity mapping
bond identity mapping
```

Strong operations must record topology edits through operation-system APIs such as:

```text
record_topology_edit(...)
record_topology_mapping(...)
begin_topology_mut() / commit_topology(...)
begin_coordinates_mut() / commit_coordinates(...)
begin_properties_mut() / commit_properties(...)
```

These APIs are responsible for:

```text
registry-checked access
mapping
cache invalidation
value semantics
trace recording
```

Operation bodies may compute edit plans and mutate their owned local working
blocks, but must not directly mutate `Molecule` internals, bypass begin/commit,
or keep a separate read view for a write-owned block.

Appending atoms or bonds is a strong topology edit.

### Weak topology-state operations

A weak operation preserves atom and bond identity and ordering, but changes local graph state.

Examples include:

```text
kekulize
sanitize
set_aromaticity
local bond-order update
formal-charge update
stereo assignment from existing topology or coordinates
```

Weak operations may begin and commit `topology` only when atom and bond tables remain stable. They must still clear or recompute affected derived state.

---

## 8. Operation Lifecycle

The lifecycle is:

```text
public method entry
registry spec selection
OpParts::new(source, spec)
begin each write-owned block through OpParts
read and mutate the local owned working blocks
commit changed blocks back through OpParts
record topology edit kind when topology identity/state changed
return OpOutcome
parts.finish(outcome)
strict validation
return Molecule
```

Operation bodies must use registry-derived begin/commit block access rather
than free `read_parts()` and independent mutation accessors.

For a write-owned block, reads and writes come from the same local owned
working value:

```rust
let mut topology = parts.begin_topology_mut()?;
let mut coordinates = parts.begin_coordinates_mut()?;
let mut properties = parts.begin_properties_mut()?;

let assignment = topology.add_hs_assignment(&params)?;
apply_add_hs_assignment(
    &mut topology,
    &mut coordinates,
    &mut properties,
    &assignment,
)?;

parts.commit_topology(topology)?;
parts.commit_coordinates(coordinates)?;
parts.commit_properties(properties)?;
parts.record_topology_edit(TopologyEditKind::Appending)?;
parts.record_topology_mapping(mapping);
parts.prove_preserved(
    DerivedState::RINGS | DerivedState::RING_FAMILIES,
    PreservationProof::LeafAtomAppend,
)?;
parts.clear_cache(WITH_HYDROGENS_SPEC.needs_update());
```

The operation must not keep a separate read view of a write-owned block while
also owning the mutable working value for that block.

For a read-owned block, `OpParts` may expose a read-begin method. For a
write-owned block, `OpParts` exposes only `begin_*_mut()` and `commit_*()`.
For an inaccessible block, no begin method is legal.

Strong operations must record the declared topology edit kind after the edit.
Weak operations may record only local stable-index topology edits.

---

## 9. OpParts

`OpParts` is the only mutable capability object for molecule operation bodies.

It owns:

```text
cheap working clone creation
copy-on-write block detachment
registry-derived block access construction
mutation permission checks
registry-driven remap and topology mapping
cache invalidation
handled-state tracing
topology mapping storage
operation finalization
```

It must not contain chemistry rules, sanitization policy, stereo perception, ring perception, operation-specific branching, or source-library guesses.

Operation bodies may mutate molecule state only through `OpParts` methods or references obtained from them.

`OpParts` must not expose a whole-molecule read view that overlaps with an
actively begun write-owned block. Read access to a write-owned block is allowed
only through the local owned working value.

---

## 10. Copy-On-Write

`OpParts::new()` must be cheap.

It may clone the top-level molecule wrapper and share internal blocks through `Arc` or equivalent COW storage.

It must not eagerly deep-clone topology, coordinates, properties, conformers, or caches.

Within one operation, each write-owned block should be materialized at most once
unless a source-backed or memory-bound reason is documented.

Release optimization must not bypass `OpParts`, skip required remap, or weaken invalidation.

---

## 11. Begin/Commit APIs

The framework must provide operation-scoped begin/commit methods derived from
the registry declaration:

```rust
let mut topology = parts.begin_topology_mut()?;
parts.commit_topology(topology)?;
```

The begin/commit API is the only route to mutable operation state. It exposes
methods according to block access mode:

```text
none  -> no method
read  -> begin_*_read()
write -> begin_*_mut(), commit_*()
```

For write-owned blocks, the operation receives an owned working block. Reads
and writes both happen through that same local value.

The begin/commit API must:

- remain internal to the operation framework
- respect registry access and mutation permissions
- avoid chemistry logic
- preserve cheap `OpParts::new()`
- materialize only write-owned blocks
- keep mapping, remap, and invalidation under framework control

Strong topology operations must record the edit kind declared by the registry:

```text
record_topology_edit(TopologyEditKind::Appending)
record_topology_edit(TopologyEditKind::Compacting)
record_topology_edit(TopologyEditKind::Renumbering)
record_topology_edit(TopologyEditKind::Merge)
```

The operation body owns chemistry and row-level mutation of its write-owned
blocks, including applying a topology mapping to other write-owned local blocks.
The framework owns access validation, topology mapping artifact storage, cache
invalidation tracing, handled-state tracing, and operation finalization.

---

## 12. Mutation Surface

`access` defines the legal block capability surface. `may_mutate` is the write
subset of that surface.

In strict builds, mutating outside this surface is a developer error.

Rules:

- a block declared `read` must not expose a mutable accessor
- a block declared `write` must not also be exposed through an independent read
  view
- a block declared `none` must not be accessible
- direct `topology_mut()` is not an operation-body API
- strong topology edits must record their declared edit kind and required
  topology mapping artifact before `finish()`

If more mutation authority is needed, update the registry and framework API first.

---

## 13. Derived State

Every operation must explicitly classify affected derived state through
`derived_effects`.

`derived_effects` is the primary registry contract. `must_handle()` and
`needs_update()` are derived compatibility views and must not be declared
directly in new molecule-operation registry entries.

The contract dimensions are:

```text
invalidate
recompute
preserve
require_handle
unsupported
```

### `invalidate`

`invalidate` means old derived state becomes stale and must be cleared unless
the operation also recomputes a replacement.

Invalidated states contribute to the derived compatibility view
`needs_update()`.

### `recompute`

`recompute` means the operation must produce a fresh framework-visible value
for the state.

Recomputed states contribute to both derived compatibility views:

```text
needs_update()
must_handle()
```

### `preserve`

`preserve` means the old derived state remains valid after the operation.

Preservation is not a label-only shortcut. If a topology operation declares
preserved derived state, the operation body must call an approved
framework-checked preservation proof, for example:

```text
PreservationProof::LeafAtomAppend
```

The proof must validate objective structural conditions, such as old atom and
bond identity preservation plus appended degree-one leaf atoms. Silent
preservation without proof is invalid in strict builds.

### `require_handle`

`require_handle` means the operation cannot ignore that state even if the state
is not materialized as a cache update. It must be satisfied by real
framework-recognized handling such as recompute, validation, clear/drop policy,
or structured unsupported behavior. Trace metadata alone is not enough.

`require_handle` contributes to the derived compatibility view
`must_handle()`.

### `unsupported`

`unsupported` means the operation recognizes that a derived-state effect exists
but the modeled operation cannot currently provide the required behavior.

Unsupported behavior must be explicit and structured. It must not be hidden by
cache clearing, placeholder values, or label changes.

### Materialized vs invalidation-only state

Some states are stored caches. Others are invalidation-only downstream products, such as drawing or fingerprint output.

Invalidation-only states may be invalidated, but must not be marked recomputed unless the operation actually produces a new value.

---

## 14. Mapping And Remap

Strong topology operations must record topology mapping when required by the registry.

If `auto_remap` declares a block, the operation must apply the topology mapping
to that write-owned local block and record the mapping artifact through
`OpParts`; strict trace validation verifies that the registry-required remap
surface was covered.

Operation bodies must not mutate unbegun blocks or directly edit molecule
internals to satisfy remap requirements.

If a block cannot be remapped meaningfully, the operation must explicitly drop it, return structured unsupported behavior, or use an approved policy change.

---

## 15. Mapping Artifacts And Outcomes

Topology mapping is an internal operation artifact, not a separate report
module. If the registry requires mapping, `finish()` must verify that the
mapping artifact was recorded.

`OpOutcome` must accurately indicate whether the operation changed the molecule or returned an allowed no-op.

If `allows_noop` is false, a successful operation must not silently return unchanged output.

---

## 16. Unsupported Behavior

Unsupported behavior must fail through structured errors.

Operations must not:

- guess missing chemistry
- emit plausible placeholders
- silently skip required state handling
- downgrade strong edits into partial weak edits

If required source behavior depends on state not yet modeled, keep that path explicitly unsupported.

---

## 17. Source Porting

Operation ports from RDKit or other libraries must follow `source_reproduction_protocol.md`.

This means:

- source lines remain inline near the Rust body
- behavior and performance markers are reviewed line by line
- unsupported branches stay explicit
- marker status must match real implementation behavior

“Fully ported” means all relevant source branches for the modeled input space are represented, state handling is complete, and tests cover reproduced behavior.

---

## 18. Testing

Every registered operation must have tests derived from:

```text
kind
domain
access
may_mutate
derived_effects
requires_mapping
parity
io_roundtrip
```

Minimum requirements:

- strong operations need mapping, remap, invariant, source-unchanged, and cache-state tests
- weak operations need stable-index, derived-state, invariant, and source-unchanged tests
- parity-required operations need executable parity tests or executable known failures
- unsupported operations need explicit-error and source-unchanged tests

Invariant profiles must map to meaningful check sets. Different profiles must not silently collapse into one global default unless documented.

---

## 19. Strict And Release Builds

Core operation work must pass:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict
```

Strict mode catches:

```text
unauthorized mutation
missing handled state
missing invalidation
missing mapping
post-operation invariant failure
```

Release mode may remove checking overhead, but must use the same operation path, `OpParts` accessors, and COW model.

Only checks disappear. Semantics must not change.

---

## 20. Review

An operation change is review-complete only when reviewers confirm:

- registry spec matches behavior
- strong/weak classification is correct
- `may_mutate` is minimal and sufficient
- `derived_effects` matches invalidated, recomputed, preserved, handled, and
  unsupported state
- unsupported behavior is explicit
- source-port markers are accurate
- performance refactors preserve value semantics and framework ownership

Changing labels or metadata without behavior and test evidence is non-compliant.

---

## 21. Migration Direction

Operation-system evolution must move toward more explicit framework ownership.

Preferred directions:

```text
stronger registry semantics
complete topology edit helpers
operation-scoped working-set sessions
generation-based cache invalidation
metadata-derived invariant runners
```

Non-preferred directions:

```text
ad hoc public mutation helpers
parallel metadata systems
chemistry logic inside OpParts
performance shortcuts that weaken contracts
```

---

## 22. Canonical Decision Rule

When convenience conflicts with discipline, use this order:

```text
correctness first
explicit behavior second
operation-system contract third
performance inside that boundary
```

If a change appears to require bypassing the registered operation system, stop and ask the human author.
