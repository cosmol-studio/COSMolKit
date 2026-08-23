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

## 6. Registry Inputs And Generated Spec Fields

Every generated `MoleculeOpSpec` contains:

```text
method
impl_fn
domain
kind
topology_edit
access
may_mutate
auto_remap
derived_effects
semantic_preconditions
requires_mapping
support
parity
io_roundtrip
```

The `molecule_ops!` entry must explicitly provide `method`, `impl_fn`, `kind`,
`access`, `derived_effects`, `feature`, `parity`, and `invariant_profile`.
Every parity policy other than `not_applicable` also requires
`parity_profile`. The macro supplies defined defaults for optional inputs such
as `domain`, `topology_edit`, `may_mutate`, `auto_remap`,
`semantic_preconditions`, `requires_mapping`, and `io_roundtrip`; entries
should still spell out behaviorally meaningful values when omission would make
review ambiguous.

Strong topology operations must also define their migration surface, including fields such as:

```text
topology_edit
auto_remap
```

These fields are either runtime contract inputs or generated evidence
requirements. None may be treated as unowned documentation-only metadata.

`access` is the authoritative block-capability declaration. For every molecule
block that the operation may touch, the registry must declare exactly one
access mode:

```text
none
read
write
```

`may_mutate` must equal `access.write`. It exists for generated matrices,
strict checks, and compatibility with current operation metadata; it does not
grant additional authority beyond `access`.

A block must not be exposed through both an independent read capability and a
write capability. If an operation needs to inspect a write-owned block before
modifying it, it must begin the write-owned block and read from the same local
owned working value.

### Field execution ownership

Registry fields must have an identified enforcement or evidence owner. The
current ownership is:

| Field | Execution or evidence owner |
|---|---|
| `access` | `OpParts` capability construction; strict runtime access checks |
| `may_mutate` | strict runtime consistency and mutation-trace checks |
| `derived_effects` | strict cache APIs, preservation proofs, and `finish()` trace validation |
| `requires_mapping` | strict `finish()` validation of the mapping artifact |
| `semantic_preconditions` | `OpParts::new` / `new_in_place` entry validation in strict and default release builds |
| `support` | macro-generated public-wrapper rejection and `SUPPORT_MATRIX` |
| `parity` / `parity_profile` | macro-generated `PARITY_MATRIX`; separate parity tests and CI must provide behavioral evidence |
| `io_roundtrip` | registry metadata only today; operation-specific tests are required but are not selected by a universal field-driven runner |
| `invariant_profile` | macro-generated `OPERATION_INVARIANT_MATRIX` plus invariant tests and CI evidence |

This table must describe current enforcement honestly. `invariant_profile` is
present in the generated matrix, but the current `OperationInvariantEntry`
construction maps all profiles to the same required check set;
profile-specific execution must not be claimed until such a runner exists.

Removed fields such as `allows_noop`, `must_handle`, `require_handle`, and
`derived_effects.unsupported` are not current contract inputs.

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
return current blocks through the scoped mutation capability on success or error
record topology edit kind when topology identity/state changed
return `Result<(), OperationError>`
parts.finish()
strict validation
return Molecule
```

Fallible operation bodies must use registry-derived scoped mutation
capabilities rather than free `read_parts()` and independent mutation
accessors. The scoped capability owns the begin/commit lifecycle and returns
the current block to the working molecule on both `Ok` and `Err`.

For a write-owned block, reads and writes come from the same local owned
working value:

```rust
let changed = parts.with_topology_coordinates_properties_mut(
    |parts, topology, coordinates, properties| {
        let assignment = parts.with_block_read_parts(
            topology.clone(),
            coordinates.clone(),
            properties.clone(),
            |read| read.add_hs_assignment(&params).map_err(map_add_hs_error),
        )?;
        apply_add_hs_assignment(
            parts,
            topology,
            coordinates,
            properties,
            &assignment,
        )
    },
)?;

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
write-owned block, operation bodies use scoped `with_*_mut()` capabilities.
Low-level `begin_*_mut()` and `commit_*()` methods implement and test that
lifecycle, but a block obtained from `begin_*_mut()` must not remain checked
out across a fallible return. For an inaccessible block, no begin or scoped
mutation method is legal.

Strong operations must record the declared topology edit kind after the edit.
Weak operations may record only local stable-index topology edits.

---

## 9. OpParts

`OpParts` is the only mutable capability object for molecule operation bodies.

`OpParts` is defined in the private molecule-operation runtime module. That
module is the only module allowed to own or touch the internal working
`Molecule`. `#[mol_op_body]` implementation functions live in sibling modules
and receive only `&mut OpParts`; they must not share a module with the runtime
state.

It owns:

```text
cheap working clone creation
copy-on-write block detachment
registry-derived block access construction
mutation permission checks
registry-driven remap and topology mapping
cache invalidation
derived-effect tracing
topology mapping storage
operation finalization
```

It must not contain chemistry rules, sanitization policy, stereo perception, ring perception, operation-specific branching, or source-library guesses.

Operation bodies may mutate molecule state only through `OpParts` methods or references obtained from them.

`OpParts` must not expose a whole-molecule read view that overlaps with an
actively begun write-owned block. Read access to a write-owned block is allowed
only through the local owned working value.

Operation-body helpers must not accept raw `Molecule` or `&Molecule` as a
convenience path. Use `MoleculeReadParts`, slices, coordinate blocks, or typed
assignment/update plans instead. If an algorithm historically consumed a whole
molecule, operations must use a two-stage shape: read-only calculation from
narrowed inputs, followed by block-level writeback through `OpParts`.

Canonical operation-body shape:

```rust
#[mol_op_body(with_example_state, parts)]
fn with_example_state_impl(args: ExampleArgs) -> Result<(), OperationError> {
    parts.with_topology_mut(|parts, topology| {
        let plan = parts.with_topology_read_parts(topology.clone(), |read| {
            crate::example_domain::compute_plan(read, &args).map_err(|source| {
                OperationError::Example {
                    operation: &WITH_EXAMPLE_STATE_SPEC,
                    source,
                }
            })
        })?;
        crate::example_domain::apply_plan(topology, &plan);
        Ok(())
    })?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(DerivedState::DRAWING);
    Ok(())
}
```

Forbidden operation-body shape:

```rust
fn with_example_state_impl(...) -> Result<(), OperationError> {
    let molecule = parts.working.clone();
    helper_that_accepts_whole_molecule(&molecule);
    Ok(())
}
```

The operation source tree carries guard tests for this boundary. If a new
operation appears to need direct `parts.working`, a raw molecule escape, or an
operation body inside the runtime module, treat that as a design exception and
stop for human-author approval instead of weakening the guards.

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
invalidation tracing, derived-effect tracing, and operation finalization.

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

`derived_effects` is the primary registry contract. `needs_update()` is a
derived compatibility view and must not be declared directly in new
molecule-operation registry entries.

The contract dimensions are:

```text
recompute
preserve
invalidate
operation_defined
```

### `invalidate`

`invalidate` means old derived state becomes stale and must be cleared.

Invalidated states contribute to the derived compatibility view
`needs_update()`.

### `recompute`

`recompute` means the operation must produce a fresh framework-visible value
for the state, or explicitly clear it when reproduced source behavior leaves
no materialized replacement.

Recomputed states contribute to the derived compatibility view:

```text
needs_update()
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

### `operation_defined`

`operation_defined` is an explicit escape hatch for a source-required state
transition that cannot be classified truthfully as preserve, recompute, or
invalidate. It delegates the transition mechanism, not the correctness
obligation. The operation may use the normal cache set, validity-update, and
clear APIs, while strict finalization still requires the declared state to be
updated or cleared. Source alignment, focused regression coverage, and the
declared parity profile remain responsible for proving the value semantics.

This exception is currently allow-listed only for `valence` in the
hydrogen-removal operation family. The registry macro, strict runtime, and a
registry-wide test reject every other use. Expanding the allow-list requires an
explicit design decision and coordinated guardrail changes; it is not a
general-purpose alternative to selecting one of the standard categories.

### Read authority is separate

The effect categories do not grant cache-read authority. Reading derived cache
state requires block-level authority from `access.read: [derived_cache]` (or a
write-owned derived-cache block). `preserve` is a proof obligation, not read
permission.

The current capability model controls the derived-cache block as a whole. It
does not enforce separate read permissions for rings, valence, aromaticity,
stereo, or other individual cache entries.

`needs_update()` remains the only molecule-operation compatibility view and is
derived as `recompute | invalidate | operation_defined`. Molecule operations have no
`must_handle()` view or `must_handle`/`require_handle` registry input.

Unsupported behavior is not a derived effect. It must be returned as a
structured operation error and must not be encoded as cache metadata.

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

## 15. Mapping Artifacts And Completion

Topology mapping is an internal operation artifact, not a separate report
module. If the registry requires mapping, `finish()` must verify that the
mapping artifact was recorded.

Operation bodies return `Result<(), OperationError>` and do not classify a
successful source execution as changed or unchanged. Contract obligations are
derived from structural trace facts such as claimed and committed writable
blocks, recorded topology edits, mappings, remaps, and derived-state effects.
Whether the source algorithm happened to produce an equal value is not a
generic operation-contract dimension and belongs in source-parity tests when it
is externally observable.

---

## 16. Unsupported Behavior

Unsupported behavior must fail through structured errors.

This is a capability-boundary design rule, not a parity-test disposition. It
applies when an independently identified API, option family, or state model is
outside the declared supported surface. A mismatching input inside a supported
surface is a bug and must not be relabeled unsupported.

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
semantic_preconditions
requires_mapping
support
parity
io_roundtrip
invariant_profile
parity_profile when parity is required
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
cargo test -p cosmolkit-core --release --features op-contracts-strict
```

Small focused test filters may use the default debug profile during iteration.
Large local runs, parity suites, and CI test runs should use release mode with
the same strict feature set. Release-mode testing must not relax the checks
listed below; they are controlled by `op-contracts-strict`, not by Cargo's
optimization profile.

Strict mode catches:

```text
unauthorized mutation
missing derived-effect handling
missing mapping
post-operation invariant failure
```

The Cargo optimization profile and the contract feature set are independent.
`--release --features op-contracts-strict` retains strict checks. A default
published release build without the strict feature may omit checking overhead,
but it must use the same operation path, `OpParts` accessors, and COW model.

Only feature-gated checks disappear. Semantics must not change.

---

## 20. Review

An operation change is review-complete only when reviewers confirm:

- registry spec matches behavior
- strong/weak classification is correct
- `may_mutate` is minimal and sufficient
- `derived_effects` matches invalidated, recomputed, preserved, and explicitly
  allow-listed operation-defined state
- derived-cache reads are covered by block-level `access`
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
