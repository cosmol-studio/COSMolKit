# In-Place Operation API Design

## Goal

Expose explicit in-place operation methods without adding a third public API
family. Normal methods keep value semantics; trailing-underscore methods mutate
the receiver through the same registered operation system.

```rust
let next = mol.with_hydrogens()?;
mol.add_hydrogens_()?;
mol.sanitize_()?;
mol.kekulize_()?;
```

## Naming

- A trailing underscore means "mutates `self`".
- Every public `Molecule` in-place operation must end with `_`.
- No public `Molecule` method may use a trailing `_` for any meaning other than
  in-place mutation.
- Default generated in-place name: `{method}_`.
- Builder-style names may override the in-place name:
  - `with_hydrogens` -> `add_hydrogens_`
  - `without_hydrogens` -> `remove_hydrogens_`

## Execution Model

The ops macro should generate two entry points for eligible operations:

```rust
pub fn with_hydrogens(&self) -> Result<Molecule, OperationError>;
pub fn add_hydrogens_(&mut self) -> Result<(), OperationError>;
```

Both entry points call the same operation implementation. The difference is in
`OpParts` construction:

- Borrowed mode keeps current behavior: source is preserved and writable blocks
  are cloned into a working molecule.
- In-place mode works on the receiver and avoids the working-copy clone when a
  block is uniquely owned.

## Clone Boundary

In-place methods only promise to avoid the operation-system copy needed for
non-mutating value semantics.

- If a block has no other `Arc` owners, the operation may take and mutate the
  original block directly.
- If a block is shared with another `Molecule`, the operation must clone before
  mutation to preserve value semantics.
- Clones required by the operation implementation itself remain the
  implementation's responsibility and are not hidden by the in-place API.

## Error Semantics

In-place operations provide a basic failure guarantee, not transactional
rollback:

- If an in-place operation returns an error, the receiver is not guaranteed to
  equal its pre-call value. It may retain the partial changes made before the
  source algorithm reported the error.
- The receiver's internal storage must remain complete after a returned error.
  Every block checked out from the working molecule must be returned, and an
  operation-system placeholder or default block must never escape into the
  receiver as an error-path artifact.
- Internal completeness does not mean that sanitization, kekulization,
  stereochemistry assignment, or another requested chemistry operation
  completed successfully. Callers must still handle the returned error.
- Derived state affected by partial changes must be invalidated or updated
  according to the registered operation contract.
- The operation system must not preserve a full old `Molecule`, clone writable
  blocks solely for rollback, or otherwise turn the in-place path into the
  non-mutating COW path.

Fallible operation bodies must use scoped mutable block capabilities. These
capabilities return the current owned blocks to the working molecule on both
`Ok` and `Err` before the public method returns. A raw `begin_*_mut()` value
must not cross a fallible `?` path before its matching `commit_*()`.

Users that require failure-preserving value semantics should call the
non-mutating method. An individual operation may document and test a stronger
transactional guarantee when all of its fallible work happens before mutation,
but that is not the general in-place contract.

## Contract State

In-place execution must not keep a full old `Molecule` as the source, because
that would reintroduce shared block ownership. Contract checks that need old
state should use a lightweight source snapshot containing only the fields
required by the operation contract. This snapshot is contract-only state and
must be compiled only when `op-contracts` is enabled.

## Registry Shape

Eligible operations should declare whether an in-place method is generated and,
when needed, its public name.

```rust
inplace: true,
inplace_method: add_hydrogens_,
```

The generated in-place method must still go through operation registry and state
migration machinery. Development-only permission, invariant, and contract checks
run under `op-contracts-strict`; default release builds must not keep
contract-only source snapshots, preservation-proof checks, registry access
checks, topology-edit declaration checks, or begin/commit lifecycle checks.
Public APIs must not expose mutable block storage directly.
