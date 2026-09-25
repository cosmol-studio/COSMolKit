# COSMolKit Operation System Standard

This document defines the operation-system standard for COSMolKit.

It is a binding project rule and target architecture standard. New operation code must follow this standard unless the human author explicitly approves an exception.

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
Crate ownership and the final capability/transaction shape are defined by
[crate_architecture.md](./crate_architecture.md); public names are defined by
[public_api_design.md](./public_api_design.md). This document owns operational
contracts, not a second runtime architecture.

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
- the operation body delegates chemistry to its unique detached algorithm owner
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
output
result_type
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

`output` declares wrapper cardinality. Its default is `single`; `multiple`
generates an ordered `Vec<Molecule>` wrapper and cannot be combined with an
in-place wrapper. Cardinality changes wrapper shape only and grants no block or
derived-state authority.

Multiple-output operations that need status, provenance, or other domain
metadata may declare `result_type` together with `assemble_fn`. The operation
body returns metadata separately, all molecule values are finalized by
`MultiOutputOpParts`, and the assembler receives only those finalized values
plus the metadata. Supplying only one of these fields on a multiple-output
operation is a compile-time error. Single-output operations must not declare
`assemble_fn`.

Single-output value operations may register `result_type` without an assembler
to return a `MoleculeResult`-derived result with one marked molecule field
(`M` or `Option<M>`). Its public default is `Molecule`; the body returns the
same result parameterized by private `PendingMolecule<Access>`. Only that
registered operation receives the generated `pending_molecule()` capability.
It seals the already staged detached blocks after the body's mapping/effect
bookkeeping; it does not construct a live molecule. Further capability access,
duplicate sealing, foreign-transaction pending values, and dropped pending
values are rejected. Pending containers have no public constructor, extraction,
clone, dereference, or molecule methods.

The generated wrapper alone invokes result finalization. It uses the existing
transaction's validation and commit path, then the derived field conversion
replaces the pending value with the finished molecule. No domain assembler or
registry closure receives a live molecule. An absent optional field still
requires the transaction to complete successfully. This mechanism grants no
new block permissions and is not available to in-place or multiple-output
operations; their existing lifecycles are unchanged. Validation feature gates
do not disable pending ownership, single-use, or module-privacy enforcement.

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
| `output` | macro-selected single or multiple capability lifecycle; compile-time rejection of multiple-output in-place wrappers |
| `result_type` / `assemble_fn` | registered single-output pending-field finalization, or multiple-output assembly after runtime candidate validation |
| `access` | marker-specific generated capabilities and module privacy; additional strict runtime access checks |
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

The canonical transaction shapes are defined in
[Crate Architecture, Parent Operation Flow](./crate_architecture.md#5-parent-operation-flow).
The generated wrapper owns transaction creation, finalization, and abort;
the operation body uses only generated marker-specific capabilities to
extract or read authorized values, call the unique domain owner, stage results,
and record mappings and effects. It cannot construct or commit a live molecule.

Read and write access to a write-owned block use the same local working value,
not an independent read view. Every required block must be returned through
the authorized lifecycle before a fallible return. Scoped mutation restores
runtime-owned blocks on both success and error; raw checkout/install helpers
are not an escape hatch for operation bodies.

### Multiple-output lifecycle

The domain owner produces an ordered collection of detached candidates.
`MultiOutputOpParts<'_, Access>` validates each candidate and applies the
declared effects before constructing public outputs. There is no parallel
branch-ID graph, parent-candidate derivation framework, or domain-side live
molecule constructor. Ordering, duplication, empty-result behavior, and domain
metadata follow the selected source behavior and its tests. Multiple-output
operations do not have an in-place form.

### In-place execution and failure semantics

Eligible value and in-place entry points share one registered implementation.
The in-place form may reuse uniquely owned storage; shared blocks require COW
to preserve other values. Algorithm-internal copies remain the algorithm's
responsibility, not a promise eliminated by the operation wrapper.

The general in-place guarantee is basic failure safety, not rollback:

- A returned error may leave completed partial changes in the receiver.
- Internal block storage must be complete and structurally usable; checked-out,
  placeholder, or default replacement blocks must not escape.
- Storage completeness does not imply successful sanitization, kekulization,
  or stereochemistry assignment. The error remains visible.
- Affected derived state must be invalidated or updated under the contract.
- Do not keep a full old `Molecule` or clone writable blocks solely for rollback.

Fallible operation bodies use generated scoped mutation capabilities. A raw
begin/check-out value cannot cross a fallible return without its corresponding
return of ownership. Runtime cleanup and unwind behavior remain framework
responsibilities; domain algorithms do not gain abort authority.

Callers needing source-preserving failure semantics use the value form.
An operation may guarantee stronger atomicity when its fallible work precedes
mutation, but that does not redefine the general in-place contract.

Contract validation in in-place mode uses only the lightweight old-state fields
it needs, not a full old molecule that forces shared ownership. Contract-only
snapshots are compiled only with `op-contracts`. Default release omits these
snapshots and diagnostic preservation/access/edit/lifecycle assertions, not
required correctness or compile-time capability boundaries.

Eligible declarations use `inplace: true` and, when needed, an explicit
`inplace_method`. Naming is governed by the public API standard; neither form
exposes mutable storage.

---

## 9. OpParts

The runtime boundary is `ops::runtime::{context,multiple,registry}` inside
`cosmolkit`. Operation bodies live outside that semantic module subtree and
receive marker-specialized capabilities generated from their declaration.
Physical file placement alone is not a Rust privacy boundary.

The runtime owns COW, access enforcement, mapping artifacts, derived-state
tracking, validation, and finalization. It contains no chemistry or
operation-specific algorithm. Domain crates receive explicit detached values,
slices, or assignments, never `Molecule`, `OpParts`, runtime views, or commit
authority.

Bodies cannot reach runtime fields, unrestricted read/write primitives, or
constructor/finish/abort methods. Helpers cannot recover a raw molecule or
broader capability than the declaration grants. Real-module compile-pass and
compile-fail tests must enforce this boundary in default and strict builds.

---

## 10. Copy-On-Write

`OpParts::new()` must be cheap.

It may clone the top-level molecule wrapper and share internal blocks through `Arc` or equivalent COW storage.

It must not eagerly deep-clone topology, coordinates, properties, conformers, or caches.

Within one operation, each write-owned block should be materialized at most once
unless a source-backed or memory-bound reason is documented.

Release optimization must not bypass `OpParts`, skip required remap, or weaken invalidation.

---

## 11. Authorized Block Lifecycle

The operation declaration generates the accessible read/write surface.
A block declared `none` has no capability; a read-owned block cannot be
mutated; a write-owned block is inspected and changed through the same
authorized working value. Concrete generated extraction/staging APIs follow
the architecture's capability projection, not an independent handwritten list.

Block lifecycle machinery remains internal, materializes only write-owned
blocks, and keeps access, remapping, and finalization under runtime control.
A strong operation records its declared appending, compacting, renumbering,
or merge edit and required mapping; a weak topology edit records only local
stable-identity changes. The domain owner computes row changes and detached
mappings, while the thin body stages them and records the required artifacts.

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

[derived_effects_permission_model.md](./derived_effects_permission_model.md)
defines the four pairwise-disjoint effect categories, permitted cache actions,
proof requirements, materialized versus invalidation-only state, and the
narrow `operation_defined` allow-list. Keep that contract in one place.

Every operation declares its affected state through `derived_effects`.
Cache read authority comes from block access, not from an effect label.
The runtime records and validates the operation's actual handling; metadata
alone does not prove preservation or successful recomputation. Unsupported
capabilities remain structured errors, never effect categories.

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

## 19. Strict And Release Builds

Core algorithm work must pass:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --release --features op-contracts-strict
```

Runtime, operation integration, and macro work also require affected-crate
checks with the affected capability features enabled:

```bash
cargo check -p cosmolkit --features op-contracts-strict
cargo test -p cosmolkit --release --features op-contracts-strict
cargo test -p cosmolkit --release --test migration_run_privacy
```

Final cross-crate validation uses:

```bash
cargo test --workspace --release --features cosmolkit/op-contracts-strict,cosmolkit-core/op-contracts-strict
```

Core strict alone is not runtime validation. Preserve exact commands, exits,
nonzero counts, and failures; only plan-prescribed stage exclusions are allowed.
A stage pass is not a workspace pass.

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
