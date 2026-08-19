# Topology Operation Design

This document describes the current topology-operation design in COSMolKit.
The binding cross-cutting rules are defined in
[`operation_system_standard.md`](./operation_system_standard.md); the
derived-state model is defined in
[`derived_effects_permission_model.md`](./derived_effects_permission_model.md).

## Purpose

Topology and graph-state operations can affect atom and bond identity,
coordinates, properties, stereo references, mappings, and derived caches. The
operation system keeps those effects inside one registry-controlled path while
preserving value-style and in-place APIs.

Public mutation-capable topology and coordinate APIs must be traceable to a
registered `MoleculeOpSpec`. Internal read-only algorithms and operation-local
assignment plans are not separate operations.

## Source Of Truth

`molecule_ops!` is the source of truth for registered molecule operations. It
generates:

```text
Molecule methods
in-place Molecule methods when requested
MoleculeOpSpec values
MOLECULE_OPS
SUPPORT_MATRIX
OPERATION_INVARIANT_MATRIX
PARITY_MATRIX
```

Do not maintain a second `TopologyOpSpec`, `DependentData` list, handwritten
operation table, or parallel wrapper path.

## State Blocks

The operation capability model controls these molecule blocks:

```text
topology
coordinates
properties
derived_cache
```

Registry `access` assigns each block one mode:

```text
none
read
write
```

`may_mutate` must equal `access.write`. It is retained for generated metadata,
strict trace checking, and compatibility with the current public spec shape; it
does not grant authority beyond `access`.

A write-owned block is read and mutated through the same local owned value. An
operation must not keep an independent read view that aliases a write-owned
block.

## Strong And Weak Operations

### Strong topology operations

A strong operation changes at least one of:

```text
atom count
bond count
atom ordering
bond ordering
atom identity mapping
bond identity mapping
```

Its registry entry declares an exact `topology_edit` kind such as appending,
compacting, renumbering, or merge. When required, it must record a
`TopologyMapping` and remap every block listed in `auto_remap`.

Appending atoms or bonds is a strong edit even when every old atom index is
preserved.

### Weak topology-state operations

A weak operation preserves atom and bond identity and ordering but changes
local graph state, for example bond orders, aromatic flags, formal charges,
radicals, or atom stereo tags.

Weak operations use `TopologyEditKind::Local` when they touch topology. They do
not manufacture a strong-operation mapping and must not reorder or resize atom
or bond tables.

## Registry Shape

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

Registry-only inputs such as `feature`, `invariant_profile`, and
`parity_profile` generate the support, invariant, and parity matrices rather
than becoming additional runtime spec fields.

The macro supplies defined empty/default values for optional declarations.
Current registry entries should still spell out behaviorally meaningful fields
when omission would make review ambiguous.

## Operation Lifecycle

The macro-generated value-style wrapper performs:

```text
support check
OpParts::new(source, spec)
operation body
parts.finish()
return new Molecule
```

The in-place wrapper performs the same source-port body through
`OpParts::new_in_place`. It does not keep a full old molecule solely for
rollback. Once the in-place working state has been taken, a later
operation-body or finalization error may leave already completed mutations in
place, but the molecule's internal block storage must be complete and
structurally usable. This error policy is defined in
[`inplace_operation_api_design.md`](./inplace_operation_api_design.md).

Operation bodies return:

```rust
Result<(), OperationError>
```

They do not classify successful execution as changed or unchanged.

## Mutation Capabilities

Fallible operation bodies use scoped capabilities such as:

```text
with_topology_mut(...)
with_coordinates_mut(...)
with_topology_and_properties_mut(...)
with_topology_coordinates_properties_mut(...)
```

These methods return the current owned block to the working molecule on both
success and error. Low-level `begin_*_mut()` and `commit_*()` methods implement
and test the lifecycle, but an operation body must not leave a block checked
out across a fallible return.

Read-only calculation uses the narrow `MoleculeReadParts` view supplied by the
corresponding capability method. Helpers called by operation bodies accept
narrow inputs or explicit assignment/update plans, not a raw `Molecule` as an
escape hatch.

## Mapping And Remap

Strong compacting, renumbering, and merging edits must make old/new identity
relationships explicit through `TopologyMapping` when the registry requires a
mapping.

The operation body applies the mapping to its write-owned local coordinate or
property blocks and records the mapping through `OpParts`. Strict finalization
checks:

```text
required mapping exists
every auto_remap block was recorded as remapped
the recorded topology edit matches the registry
every begun writable block was committed
```

If source behavior cannot preserve or remap a dependent state, the operation
must drop it according to declared policy or return a structured error. It must
not leave stale indices.

## Derived State

Topology operations declare exactly three derived-effect categories:

```text
recompute
preserve
invalidate
```

They are pairwise disjoint. Cache write/clear methods and preservation proofs
record the actual handling trace. `needs_update()` is derived as
`recompute | invalidate`; Molecule operations have no `must_handle` field.

Derived-cache read authority comes from block-level `access`, not from an
effect category. Unsupported source behavior is a structured operation error,
not a derived effect.

## Value And In-Place APIs

The default public Rust and Python operation returns a new molecule and leaves
the source unchanged. Rust obtains value semantics through cheap Arc-backed
block cloning followed by copy-on-write detachment of write-owned blocks.

An explicitly named trailing-underscore method is the in-place variant. It
moves the target into the same operation machinery so unique blocks can be
reused without a defensive full clone. The two public forms execute the same
operation implementation body.

## Testing Requirements

Every registered operation needs evidence appropriate to its metadata:

- strong operations: mapping, remap, source-unchanged, block alignment, cache,
  and structural invariant tests;
- weak operations: stable identity/order, derived-state, source-unchanged, and
  invariant tests;
- parity-required operations: pinned-source corpus and option-branch parity;
- unsupported paths: structured errors without placeholder output;
- in-place operations: success, body-error, and post-error structural-validity
  tests.

Strict development and CI validation uses:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --release --features op-contracts-strict
```

The Cargo release profile does not disable strict features. Published default
release builds omit development checks but retain the same wrappers,
capability methods, operation bodies, storage semantics, and chemical behavior.
