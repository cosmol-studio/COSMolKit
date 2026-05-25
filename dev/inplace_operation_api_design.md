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
- If an in-place operation returns an error, the receiver is not guaranteed to
  equal its pre-call value. Users that require failure-preserving value
  semantics should call the non-mutating method.

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
