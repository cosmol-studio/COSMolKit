# COSMolKit Operation Registry and Macro-Controlled State Migration Design

This document defines the canonical design for COSMolKit molecule operations.

This is not a transitional sketch. New core operation code must be built around
this registry/macro/parts/trace model. If an implementation appears to require a
manual shortcut, an unregistered operation, direct mutable access to `Molecule`,
or a non-macro operation wrapper, the agent must stop and confirm the design
exception with the human author before changing code.

The goal is to make every molecule operation:

```text
declared
permission-controlled
state-aware
debug-checkable
release-zero-check-cost
```

COSMolKit must not rely on every developer or agent remembering which internal `Arc` block to clone, which cache to invalidate, which dependent state to remap, or which operation is allowed to touch which part of the molecule.

Instead, all topology-related operations must be registered in a central operation registry, and their implementation bodies must be wrapped by macro-generated capability control.

The old parity-first implementation may be used as algorithm reference only. It
must not define the new operation architecture.

---

## 1. Core Idea

COSMolKit uses a value-style `Molecule`.

Public operations must normally return a new molecule:

```rust
let mol2 = mol.with_kekulized_bonds(true)?;
let mol3 = mol2.without_hydrogens()?;
```

The original molecule must not be silently mutated.

Internally, each operation is controlled by:

```text
1. A global operation registry
2. A macro-generated wrapper
3. Capability-limited operation parts
4. Operation outcome reporting
5. Debug/test contract checking
6. Release-mode check removal
```

The central philosophy is:

```text
Registry declares what the operation is allowed and required to do.
Macro generates the wrapper and capability type.
Operation body only receives allowed parts.
Real action functions automatically record handling.
Finish checks the contract in dev/test.
Release removes the checks.
```

This architecture is mandatory for topology-related operations. A hand-written
public operation method is allowed only for methods that are explicitly
classified as non-topology, read-only accessors. Any operation that can change
topology, coordinates, stereo, atom/bond properties, molecule properties, or
derived caches must be registered.

---

## 2. Required Macro Style

COSMolKit operation bodies must use this form:

```rust
#[mol_op_body(kekulize, parts)]
fn kekulize_impl() -> Result<OpOutcome, OperationError> {
    let topology = parts.topology_mut();

    kekulize_topology(topology)?;

    parts.recompute_aromaticity()?;
    parts.recompute_valence()?;

    parts.invalidate(DerivedState::Drawing | DerivedState::Fingerprint);

    Ok(OpOutcome::Changed)
}
```

The attribute macro:

```rust
#[mol_op_body(kekulize, parts)]
```

must rewrite the function into an internal implementation equivalent to:

```rust
fn kekulize_impl(parts: &mut TopologyOnlyParts) -> Result<OpOutcome, OperationError> {
    ...
}
```

The `parts` variable is injected by the macro and has the capability type determined by the operation registry.

For example:

```text
kekulize may_mutate = [Topology]
```

therefore `parts` is a `TopologyOnlyParts`.

The body can call:

```rust
parts.topology_mut()
```

but cannot call:

```rust
parts.conformers_mut()
parts.props_mut()
```

because those methods do not exist on `TopologyOnlyParts`.

This is intentional.

---

## 3. Operation Registry

All molecule operations must be registered in a central registry.

Example:

```rust
molecule_ops! {
    op kekulize {
        method: kekulize,
        impl_fn: kekulize_impl,
        kind: Weak,

        may_mutate: [Topology],

        must_handle: [
            Aromaticity,
            Valence,
        ],

        invalidates: [
            Drawing,
            Fingerprint,
        ],

        allows_noop: true,

        support: Unsupported,
        parity: RequiredWhenSupported,
        io_roundtrip: true,
    }

    op without_hydrogens {
        method: without_hydrogens,
        impl_fn: without_hydrogens_impl,
        kind: Strong,

        may_mutate: [
            Topology,
            Conformers,
        ],

        must_handle: [
            Coordinates,
            Stereo,
            AtomProps,
            BondProps,
            Caches,
        ],

        invalidates: [
            Adjacency,
            Rings,
            Valence,
            Stereo,
            Drawing,
            Fingerprint,
        ],

        requires_mapping: [
            Atom,
            Bond,
        ],

        allows_noop: true,

        support: Unsupported,
        parity: RequiredWhenSupported,
        io_roundtrip: true,
    }

    op remove_atom(atom: AtomId) {
        method: remove_atom,
        impl_fn: remove_atom_impl,
        kind: Strong,

        may_mutate: [
            Topology,
            Conformers,
            Props,
        ],

        must_handle: [
            Coordinates,
            Stereo,
            AtomProps,
            BondProps,
            Caches,
        ],

        invalidates: [
            Adjacency,
            Rings,
            Valence,
            Stereo,
            Drawing,
            Fingerprint,
        ],

        requires_mapping: [
            Atom,
            Bond,
        ],

        allows_noop: false,

        support: Unsupported,
        parity: NotApplicable,
        io_roundtrip: true,
    }
}
```

This registry is the single source of truth.

Do not duplicate operation metadata in scattered constants, comments, tests, or documentation.

The registry must generate the public operation specs and the testing matrices.
Hand-written matrices are not allowed because they drift from operation
metadata. The generated surfaces include:

```rust
MOLECULE_OPS
SUPPORT_MATRIX
OPERATION_INVARIANT_MATRIX
PARITY_MATRIX
```

`MOLECULE_OPS` lists every registered operation.

`SUPPORT_MATRIX` connects public `FeatureSpec` entries to operations.

`OPERATION_INVARIANT_MATRIX` lists every registered operation and the invariant
profile/check set it must satisfy, including operations that currently return
`UnsupportedFeature`.

`PARITY_MATRIX` lists operations whose parity policy requires a parity profile.
For `RequiredWhenSupported` operations, the matrix may exist before the feature
is supported; tests should assert explicit unsupported behavior until the
feature moves to a supported status.

The registry must also describe operation parameters. Parameter types are part of
the operation contract because generated public methods, operation body binding,
test registries, and trace diagnostics all depend on the same signature. An
operation that accepts an atom id, bond id, atom order, selection, options
object, or policy object must declare those parameters in the registry entry.

---

## 4. Registry Field Meaning

### `method`

The public method generated on `Molecule`.

Example:

```rust
method: kekulize
```

generates:

```rust
impl Molecule {
    pub fn kekulize(&self) -> Result<Molecule> {
        ...
    }
}
```

---

### `impl_fn`

The operation body function.

Example:

```rust
impl_fn: kekulize_impl
```

This function must be annotated:

```rust
#[mol_op_body(kekulize, parts)]
fn kekulize_impl() -> Result<OpOutcome, OperationError> {
    ...
}
```

---

### operation parameters

Operations may declare parameters directly in the registry:

```rust
op remove_atom(atom: AtomId) {
    method: remove_atom,
    impl_fn: remove_atom_impl,
    kind: Strong,
    ...
}

op renumber_atoms(order: AtomOrder) {
    method: renumber_atoms,
    impl_fn: renumber_atoms_impl,
    kind: Strong,
    ...
}

op set_formal_charge(atom: AtomId, charge: i8) {
    method: set_formal_charge,
    impl_fn: set_formal_charge_impl,
    kind: Weak,
    ...
}
```

The registry macro must generate the public method with exactly these
parameters:

```rust
impl Molecule {
    pub fn remove_atom(&self, atom: AtomId) -> Result<Molecule> {
        ...
    }
}
```

The attribute macro must pass the same parameters to the operation body through
the generated internal signature. The operation body must not recover operation
parameters by reading global state or by receiving `self`.

Parameter policy:

```text
atom and bond identities use typed ids, not raw usize in public APIs
multi-atom inputs use typed selections or order objects
option sets use explicit option structs
RDKit compatibility option structs live in compat modules
```

If an operation requires a parameter shape that the registry cannot express,
the registry must be extended first. Do not hand-write a special public wrapper.

---

### `kind`

Operation kind:

```rust
kind: Strong
kind: Weak
```

A `Strong` operation changes atom/bond index space.

Examples:

```text
add_atom
remove_atom
with_hydrogens
without_hydrogens
renumber_atoms
fragment
combine
```

A `Weak` operation keeps atom/bond index space stable but changes chemical graph state.

Examples:

```text
kekulize
sanitize
set_bond_type
set_formal_charge
set_aromaticity
update_valence
```

---

### `may_mutate`

The internal storage blocks this operation is allowed to mutate.

Examples:

```rust
may_mutate: [Topology]
may_mutate: [Topology, Conformers]
may_mutate: [Conformers]
may_mutate: [Props]
```

This controls the generated `parts` type.

For example:

```text
may_mutate = [Topology]
```

means the body receives a topology-only capability object.

The body cannot access conformers or props mutably.

---

### `must_handle`

The dependent states this operation must explicitly handle.

Examples:

```rust
must_handle: [
    Aromaticity,
    Valence,
]
```

Handling does not necessarily mean mutation.

Valid handling actions include:

```text
Remap
Recompute
Invalidate
Drop
Preserve
NotPresent
NoOp
Error
```

Operation bodies must not manually mark these actions by calling a generic `handle()` method.

Instead, they must call real action methods such as:

```rust
parts.recompute_valence()?;
parts.recompute_aromaticity()?;
parts.remap_conformers(&mapping)?;
parts.remap_stereo_or_recompute(&mapping)?;
parts.drop_stereo();
parts.preserve_coordinates_verified();
```

These real action methods automatically record handling in dev/test builds.

---

### `invalidates`

The derived states or caches that must be invalidated by the operation.

Example:

```rust
invalidates: [
    Drawing,
    Fingerprint,
]
```

The operation body must call:

```rust
parts.invalidate(DerivedState::Drawing | DerivedState::Fingerprint);
```

`invalidate()` is a real action.

It must actually update cache invalidation state and cache generation counters.

---

### `requires_mapping`

Only relevant for strong operations.

Example:

```rust
requires_mapping: [
    Atom,
    Bond,
]
```

If the operation outcome is `Changed`, the operation must record the required mappings.

Example:

```rust
let mapping = parts.remove_atoms_with_mapping(&removable)?;
parts.record_atom_mapping(&mapping);
```

If the operation returns `NoOp`, mapping is not required.

---

### `allows_noop`

Whether the operation may return `OpOutcome::NoOp`.

Example:

```rust
allows_noop: true
```

This is important because some operations may be called on molecules that already satisfy the target state.

For example:

```text
without_hydrogens on a molecule with no removable hydrogens
kekulize on a molecule already kekulized
sanitize on an already sanitized molecule
```

No-op is not an error if the operation declares `allows_noop: true`.

---

### `support`

The current public support state of the operation's feature.

```rust
support: Unsupported
```

The operation support state must be derived from a public `FeatureSpec` whenever
the operation is user-facing. Unsupported operations must return structured
`UnsupportedFeature` errors. They must not panic, silently no-op, or produce
placeholder chemistry.

### `parity`

The operation's RDKit parity obligation.

```rust
parity: RequiredWhenSupported
```

Allowed values:

```text
NotApplicable
RequiredWhenSupported
RequiredNow
```

`RequiredWhenSupported` means the operation is intended to match RDKit once
implemented, but the current support status may still be `Unsupported`.

`RequiredNow` means the operation currently claims supported RDKit-compatible
behavior and must appear in the parity matrix with a pinned RDKit profile,
corpus, golden output, and comparison schema.

`NotApplicable` means RDKit is not the behavior oracle for this operation.

---

### `io_roundtrip`

Whether this operation participates in I/O roundtrip tests.

```rust
io_roundtrip: true
```

---

## 5. Operation Outcome

Every operation body must return:

```rust
Result<OpOutcome, OperationError>
```

Definition:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OpOutcome {
    Changed,
    NoOp {
        reason: &'static str,
    },
}
```

Unsupported or invalid paths are not outcomes. They are structured errors:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OperationError {
    Unsupported {
        operation: &'static MoleculeOpSpec,
        reason: &'static str,
    },
    InvalidInput {
        operation: &'static MoleculeOpSpec,
        message: String,
    },
    Chemistry {
        operation: &'static MoleculeOpSpec,
        message: String,
    },
    InvariantViolation {
        operation: &'static MoleculeOpSpec,
        failure: InvariantFailure,
    },
}
```

An unsupported path must never be represented as `NoOp`. `NoOp` means "the
operation was valid and the molecule already satisfied the requested state".
Unsupported means "the operation cannot correctly execute for this input or
parameter set".

### `Changed`

The operation actually changed molecular state.

For `Changed`, dev/test checks are strict:

```text
must_handle states must be handled
invalidates states must be invalidated
strong operations must record required mappings
weak operations must preserve atom/bond index space
result invariants must hold
```

---

### `NoOp`

The operation was valid but did not change molecular state.

Example:

```rust
Ok(OpOutcome::NoOp {
    reason: "no removable hydrogens",
})
```

For `NoOp`, dev/test checks are different:

```text
allows_noop must be true
source molecule must remain unchanged
output must be equivalent to input
must_handle does not need remap/recompute if no state was affected
no stale cache may be introduced
```

No-op must not be inferred only from hash equality.

The operation body must explicitly return `OpOutcome::NoOp`.

---

## 6. No Manual Generic `handle()`

Operation bodies must not manually call a generic method like:

```rust
parts.handle(DependentState::Valence, HandleAction::Recompute);
```

This is discouraged because it allows an operation to claim that a state was handled without actually performing the handling.

Instead, operation bodies must call real action methods:

```rust
parts.recompute_valence()?;
parts.recompute_aromaticity()?;
parts.remap_conformers(&mapping)?;
parts.invalidate(...);
```

Each real action method must internally record the corresponding handling action in dev/test builds.

Example:

```rust
impl TopologyOnlyParts<'_> {
    pub fn recompute_valence(&mut self) -> Result<()> {
        recompute_valence_impl(self.topology_mut())?;

        #[cfg(any(debug_assertions, feature = "op-contracts"))]
        self.trace.mark_handled(
            DependentState::Valence,
            HandleAction::Recompute,
        );

        Ok(())
    }
}
```

This ensures that handling is coupled to real work.

---

## 7. Cache Generation

Derived caches must use generation counters.

Example:

```rust
#[derive(Debug, Clone, Default)]
pub struct CacheGenerations {
    adjacency: u64,
    rings: u64,
    valence: u64,
    aromaticity: u64,
    stereo: u64,
    drawing: u64,
    fingerprint: u64,
}
```

A derived cache must conceptually look like:

```rust
pub struct DerivedCache<T> {
    value: Option<T>,
    generation: u64,
}
```

Invalidation must:

```text
clear the cached value
increase the generation counter
record invalidation in trace under op-contracts
```

Example:

```rust
impl TopologyData {
    pub fn invalidate(&mut self, states: DerivedState) {
        if states.contains(DerivedState::Drawing) {
            self.cache.drawing = None;
            self.cache_gen.drawing += 1;
        }

        if states.contains(DerivedState::Fingerprint) {
            self.cache.fingerprint = None;
            self.cache_gen.fingerprint += 1;
        }
    }
}
```

`parts.invalidate(...)` must call the real invalidation method:

```rust
impl TopologyOnlyParts<'_> {
    pub fn invalidate(&mut self, states: DerivedState) {
        self.topology_mut().invalidate(states);

        #[cfg(any(debug_assertions, feature = "op-contracts"))]
        self.trace.mark_invalidated(states);
    }
}
```

Dev/test contract checks must verify that states declared in `invalidates` were actually invalidated and that their generation counters changed when applicable.

---

## 8. Capability Parts

The registry field `may_mutate` controls which `Parts` type the operation body receives.

Examples:

```text
may_mutate = [Topology]
  -> TopologyOnlyParts

may_mutate = [Conformers]
  -> ConformersOnlyParts

may_mutate = [Props]
  -> PropsOnlyParts

may_mutate = [Topology, Conformers]
  -> TopologyConformersParts

may_mutate = [Topology, Props]
  -> TopologyPropsParts

may_mutate = [Conformers, Props]
  -> ConformersPropsParts

may_mutate = [Topology, Conformers, Props]
  -> AllParts
```

Each parts type only exposes mutable accessors for the blocks it owns.

For example, `TopologyOnlyParts` exposes:

```rust
parts.topology()
parts.topology_mut()
parts.recompute_valence()
parts.recompute_aromaticity()
parts.invalidate(...)
```

It does not expose:

```rust
parts.conformers_mut()
parts.props_mut()
```

This means an operation cannot accidentally mutate a block it did not declare.

---

## 9. Operation-Level COW

COSMolKit must use operation-level COW.

At operation start, the generated wrapper detaches only the blocks listed in `may_mutate`.

For example:

```text
kekulize may_mutate = [Topology]
```

The generated wrapper clones/detaches only topology:

```rust
let mut parts = TopologyOnlyParts::new(self, &KEKULIZE_SPEC);
```

Internally:

```rust
TopologyOnlyParts {
    source: self,
    topology: self.topology.as_ref().clone(),
}
```

Conformers and props remain shared.

For:

```text
without_hydrogens may_mutate = [Topology, Conformers]
```

the wrapper detaches topology and conformers:

```rust
TopologyConformersParts {
    source: self,
    topology: self.topology.as_ref().clone(),
    conformers: self.conformers.as_ref().clone(),
}
```

Props remain shared.

This avoids repeated manual calls such as:

```rust
let mut topology = (*self.topology).clone();
let mut conformers = (*self.conformers).clone();
```

inside every operation.

The operation body works on owned mutable parts, and `finish()` assembles the new `Molecule`.

---

## 10. Macro-Generated Public Method

Given:

```rust
molecule_ops! {
    op kekulize {
        method: kekulize,
        impl_fn: kekulize_impl,
        kind: Weak,
        may_mutate: [Topology],
        must_handle: [Aromaticity, Valence],
        invalidates: [Drawing, Fingerprint],
        allows_noop: true,
        support: Unsupported,
        parity: RequiredWhenSupported,
        io_roundtrip: true,
    }
}
```

or a parameterized operation:

```rust
molecule_ops! {
    op remove_atom(atom: AtomId) {
        method: remove_atom,
        impl_fn: remove_atom_impl,
        kind: Strong,
        may_mutate: [Topology, Conformers, Props],
        must_handle: [Coordinates, Stereo, AtomProps, BondProps, Caches],
        invalidates: [Adjacency, Rings, Valence, Stereo, Drawing, Fingerprint],
        requires_mapping: [Atom, Bond],
        allows_noop: false,
        support: Unsupported,
        parity: NotApplicable,
        io_roundtrip: true,
    }
}
```

the registry macro generates something equivalent to:

```rust
impl Molecule {
    pub fn kekulize(&self) -> Result<Molecule> {
        let mut parts = TopologyOnlyParts::new(self, &KEKULIZE_SPEC);

        let outcome = kekulize_impl(&mut parts)?;

        parts.finish(outcome)
    }
}
```

and for the parameterized operation:

```rust
impl Molecule {
    pub fn remove_atom(&self, atom: AtomId) -> Result<Molecule> {
        let mut parts = TopologyConformersPropsParts::new(self, &REMOVE_ATOM_SPEC);

        let outcome = remove_atom_impl(&mut parts, atom)?;

        parts.finish(outcome)
    }
}
```

The actual user-facing API remains simple:

```rust
let mol2 = mol.with_kekulized_bonds(true)?;
```

---

## 11. Macro-Generated Spec

The registry macro also generates:

```rust
pub(crate) const KEKULIZE_SPEC: MoleculeOpSpec = MoleculeOpSpec {
    name: "kekulize",
    domain: OperationDomain::Topology,
    kind: OpKind::Weak,
    may_mutate: MutBlock::Topology,
    must_handle: DependentState::Aromaticity | DependentState::Valence,
    invalidates: DerivedState::Drawing | DerivedState::Fingerprint,
    requires_mapping: MappingRequirement::empty(),
    allows_noop: true,
    support: SupportStatus::Unsupported {
        reason: "real kekulization chemistry has not been ported",
    },
    parity: ParityPolicy::RequiredWhenSupported,
    io_roundtrip: true,
};
```

This generated spec is used by:

```text
wrapper finish checks
invariant test registry
support matrix
parity matrix
I/O roundtrip registry
documentation generation if needed
```

---

## 12. Macro-Generated Test Matrices

The registry macro must generate operation test matrices.

Example:

```rust
pub(crate) static OPERATION_INVARIANT_MATRIX: &[OperationInvariantEntry] = &[
    OperationInvariantEntry {
        operation: &KEKULIZE_SPEC,
        profile: "weak_topology_state",
        required_checks: InvariantCheckSet::GRAPH_INDEX_VALIDITY
            | InvariantCheckSet::CACHE_INVALIDATION
            | InvariantCheckSet::SOURCE_UNCHANGED
            | InvariantCheckSet::UNSUPPORTED_IS_EXPLICIT,
    }
];
```

The macro must also generate the support, invariant, and parity matrices:

```rust
pub(crate) static SUPPORT_MATRIX: &[SupportMatrixEntry] = &[
    SupportMatrixEntry {
        feature: &KEKULIZE_FEATURE,
        operation: Some(&KEKULIZE_SPEC),
    },
];

pub(crate) static OPERATION_INVARIANT_MATRIX: &[OperationInvariantEntry] = &[
    OperationInvariantEntry {
        operation: &KEKULIZE_SPEC,
        profile: "weak_topology_state",
        required_checks: InvariantCheckSet::GRAPH_INDEX_VALIDITY
            | InvariantCheckSet::CACHE_INVALIDATION
            | InvariantCheckSet::SOURCE_UNCHANGED
            | InvariantCheckSet::UNSUPPORTED_IS_EXPLICIT,
    },
];

pub(crate) static PARITY_MATRIX: &[ParityMatrixEntry] = &[
    ParityMatrixEntry {
        operation: &KEKULIZE_SPEC,
        feature: &KEKULIZE_FEATURE,
        profile: "kekulize_clear_aromatic_flags",
        rdkit_version: None,
    },
];
```

If `io_roundtrip: true`, it must register:

```rust
pub(crate) static IO_ROUNDTRIP_OPS: &[&MoleculeOpSpec] = &[
    &KEKULIZE_SPEC,
];
```

This allows tests to be driven by the same registry as the implementation.

The generated test registries are not optional. Every registered topology
or coordinate operation must appear in `OPERATION_INVARIANT_MATRIX`. Fields such
as `parity` and `io_roundtrip` only control additional parity/roundtrip
obligations; they do not opt the operation out of structural operation contract
tests. Unsupported operations still require tests proving that they fail
explicitly and do not mutate the source molecule.

The invariant runner has two layers:

```text
structural molecule invariants
  atom ids match atom table rows
  bond ids match bond table rows
  bond endpoints are valid
  no self-loop bonds unless explicitly supported by a future policy
  coordinate row counts match atom count
  stereo references are valid
  cache generations are internally consistent

operation contract invariants
  registry may_mutate was obeyed
  operation body touched only allowed blocks
  must_handle states were handled by real action methods
  invalidates states were actually invalidated
  generation counters changed for invalidated caches
  Strong Changed operations recorded required mappings
  Weak Changed operations preserved atom/bond index space
  NoOp is allowed, explicit, and semantically equivalent to input
  source molecule remains unchanged
```

The registry macro must provide enough metadata for tests to run the operation
in a uniform way. Parameterized operations must use test input generators
registered alongside the operation metadata, not ad hoc hand-written loops.

---

## 13. Attribute Macro Responsibility

The required attribute macro form is:

```rust
#[mol_op_body(kekulize, parts)]
fn kekulize_impl() -> Result<OpOutcome, OperationError> {
    ...
}
```

The attribute macro must:

```text
1. Find the operation name: kekulize
2. Find the injected parts variable name: parts
3. Determine the expected parts type from the registry
4. Rewrite the function signature to accept the correct parts parameter
5. Make the parts variable available inside the body
```

Conceptual expansion:

```rust
fn kekulize_impl(parts: &mut TopologyOnlyParts) -> Result<OpOutcome, OperationError> {
    ...
}
```

For a parameterized operation:

```rust
#[mol_op_body(remove_atom, parts)]
fn remove_atom_impl(atom: AtomId) -> Result<OpOutcome, OperationError> {
    let mapping = parts.remove_atom_with_mapping(atom)?;
    parts.record_atom_mapping(&mapping.atom_mapping);
    parts.record_bond_mapping(&mapping.bond_mapping);
    parts.remap_conformers(&mapping.atom_mapping)?;
    parts.remap_stereo_or_recompute(&mapping)?;
    parts.remap_or_drop_atom_props(&mapping.atom_mapping)?;
    parts.remap_or_drop_bond_props(&mapping.bond_mapping)?;
    parts.invalidate(
        DerivedState::Adjacency
            | DerivedState::Rings
            | DerivedState::Valence
            | DerivedState::Stereo
            | DerivedState::Drawing
            | DerivedState::Fingerprint,
    );
    Ok(OpOutcome::Changed)
}
```

the expansion is conceptually:

```rust
fn remove_atom_impl(
    parts: &mut TopologyConformersPropsParts,
    atom: AtomId,
) -> Result<OpOutcome, OperationError> {
    ...
}
```

For:

```rust
#[mol_op_body(without_hydrogens, parts)]
fn without_hydrogens_impl() -> Result<OpOutcome, OperationError> {
    ...
}
```

if the registry says:

```rust
may_mutate: [Topology, Conformers]
```

then the expansion is conceptually:

```rust
fn without_hydrogens_impl(
    parts: &mut TopologyConformersParts,
) -> Result<OpOutcome, OperationError> {
    ...
}
```

The operation body must not receive `self` directly.

This prevents bypassing the capability system.

---

## 14. Registry and Attribute Macro Coordination

The operation registry remains the source of truth.

The attribute macro must not define operation metadata itself.

The attribute macro only binds the implementation body to a registered operation.

The registry defines:

```text
operation name
method name
operation kind
allowed mutable blocks
required handled states
invalidated caches
mapping requirements
test participation
```

The attribute macro defines:

```text
this function is the body for that operation
inject the correct parts variable
enforce that the body uses the registry-defined parts type
```

If an operation body refers to an operation not present in the registry, compilation must fail.

If a registry entry names an `impl_fn` that has no matching `#[mol_op_body(...)]`, compilation must fail.

---

## 15. Contract Trace

In dev/test builds, parts must record an operation trace.

Trace must include:

```rust
pub struct OperationTrace {
    touched_blocks: MutBlock,
    handled: HandledStates,
    invalidated: DerivedState,
    recorded_mappings: MappingRequirement,
    outcome: Option<OpOutcome>,
}
```

Trace data is the source for operation contract checks. It must be written only
by capability parts and real action methods, not directly by operation bodies.
An operation body may request work:

```rust
parts.remap_conformers(&mapping)?;
parts.recompute_valence()?;
parts.invalidate(DerivedState::Fingerprint);
```

but it must not directly mutate trace fields. If a new handling action is
needed, add a real action method to the appropriate parts type and make that
method record the trace.

This trace is only compiled when:

```rust
#[cfg(any(debug_assertions, feature = "op-contracts"))]
```

Release builds without `op-contracts` must not include trace storage or trace checks.

---

## 16. Feature Flags and Check Strength

COSMolKit must provide feature tags controlling check strength.

Recommended features:

```toml
[features]
default = []
op-contracts = []
op-contracts-strict = ["op-contracts"]
op-contracts-hash = ["op-contracts"]
op-contracts-paranoid = ["op-contracts-strict", "op-contracts-hash"]
```

### No feature, release mode

```text
No trace
No contract checks
No snapshot hashing
No invariant runner inside operation wrapper
Only normal operation code remains
```

This is the performance mode.

---

### `op-contracts`

Basic dev/test checking.

Checks:

```text
may_mutate capability obeyed
must_handle states recorded by real action methods
invalidates states recorded
OpOutcome is valid
NoOp only if allows_noop = true
Weak operation preserves atom/bond counts
Strong Changed operation records required mappings
source molecule remains unchanged
```

---

### `op-contracts-strict`

Stronger structural checking.

Additional checks:

```text
graph invariants
coordinate row alignment
stereo index validity
cache generation changed for invalidated states
may_mutate-excluded blocks unchanged by cheap structural snapshot
```

---

### `op-contracts-hash`

Hash-based checking.

Additional checks:

```text
topology hash
conformer hash
props hash
cache generation snapshot
source hash
output equivalence for NoOp
```

Hash checks must not assume that `Changed` always changes a hash.

Hash checks are used to verify consistency with `OpOutcome`, not to infer outcome.

---

### `op-contracts-paranoid`

Maximum checking for CI or debugging.

Additional checks:

```text
strict checks
hash checks
I/O roundtrip if enabled
optional RDKit parity hook if test environment supports it
full invariant runner
```

This feature is not intended for release wheels.

---

## 17. Release Cost Requirement

Release builds without contract features must have no operation-contract checking cost.

Specifically:

```text
OperationTrace must not be compiled.
Contract check code must not be compiled.
Debug snapshots must not be compiled.
Hash checks must not be compiled.
Invariant runner inside finish must not be compiled.
```

The generated release method must be equivalent to:

```rust
pub fn kekulize(&self) -> Result<Molecule> {
    let mut parts = TopologyOnlyParts::new_release(self);
    let outcome = kekulize_impl(&mut parts)?;
    parts.finish_release(outcome)
}
```

The only unavoidable runtime work is the real operation work:

```text
detach/clone blocks declared in may_mutate
perform chemistry operation
remap if needed
invalidate real caches if needed
assemble new Molecule
```

No debug-only abstraction must remain in release.

---

## 18. Strong Operation Contract

For `Strong` operations:

```text
atom/bond index space may change
```

If outcome is `Changed`, the operation must:

```text
record required atom/bond mapping
remap or drop coordinates
remap, recompute, drop, or error on stereo
remap, drop, or error on atom/bond props
invalidate affected derived states
preserve molecule-level props unless documented otherwise
produce a valid molecule
```

If outcome is `NoOp`, the operation must:

```text
declare allows_noop = true
leave molecule equivalent to input
avoid creating stale cache state
```

---

## 19. Weak Operation Contract

For `Weak` operations:

```text
atom/bond index space must remain stable
```

If outcome is `Changed`, the operation must:

```text
not change atom count
not change bond count
not reorder atom table
not reorder bond table
preserve coordinates unless documented otherwise
handle all must_handle states
invalidate all declared derived states
produce a valid molecule
```

If outcome is `NoOp`, the operation must:

```text
declare allows_noop = true
leave molecule equivalent to input
avoid creating stale cache state
```

---

## 20. Example: Kekulize

Registry:

```rust
molecule_ops! {
    op kekulize {
        method: kekulize,
        impl_fn: kekulize_impl,
        kind: Weak,

        may_mutate: [Topology],

        must_handle: [
            Aromaticity,
            Valence,
        ],

        invalidates: [
            Drawing,
            Fingerprint,
        ],

        allows_noop: true,

        support: Unsupported,
        parity: RequiredWhenSupported,
        io_roundtrip: true,
    }
}
```

Body:

```rust
#[mol_op_body(kekulize, parts)]
fn kekulize_impl() -> Result<OpOutcome, OperationError> {
    if parts.topology().is_kekulized() {
        return Ok(OpOutcome::NoOp {
            reason: "molecule is already kekulized",
        });
    }

    let topology = parts.topology_mut();

    kekulize_topology(topology)?;

    parts.recompute_aromaticity()?;
    parts.recompute_valence()?;

    parts.invalidate(DerivedState::Drawing | DerivedState::Fingerprint);

    Ok(OpOutcome::Changed)
}
```

Expected generated public method:

```rust
impl Molecule {
    pub fn kekulize(&self) -> Result<Molecule> {
        let mut parts = TopologyOnlyParts::new(self, &KEKULIZE_SPEC);

        let outcome = kekulize_impl(&mut parts)?;

        parts.finish(outcome)
    }
}
```

---

## 21. Example: Without Hydrogens

Registry:

```rust
molecule_ops! {
    op without_hydrogens {
        method: without_hydrogens,
        impl_fn: without_hydrogens_impl,
        kind: Strong,

        may_mutate: [
            Topology,
            Conformers,
        ],

        must_handle: [
            Coordinates,
            Stereo,
            AtomProps,
            BondProps,
            Caches,
        ],

        invalidates: [
            Adjacency,
            Rings,
            Valence,
            Stereo,
            Drawing,
            Fingerprint,
        ],

        requires_mapping: [
            Atom,
            Bond,
        ],

        allows_noop: true,

        support: Unsupported,
        parity: RequiredWhenSupported,
        io_roundtrip: true,
    }
}
```

Body:

```rust
#[mol_op_body(without_hydrogens, parts)]
fn without_hydrogens_impl() -> Result<OpOutcome, OperationError> {
    let removable = find_removable_hydrogens(parts.topology());

    if removable.is_empty() {
        return Ok(OpOutcome::NoOp {
            reason: "no removable hydrogens",
        });
    }

    let mapping = parts.remove_atoms_with_mapping(&removable)?;

    parts.record_atom_mapping(&mapping.atom_mapping);
    parts.record_bond_mapping(&mapping.bond_mapping);

    parts.remap_conformers(&mapping.atom_mapping)?;
    parts.remap_stereo_or_recompute(&mapping)?;
    parts.remap_or_drop_atom_props(&mapping.atom_mapping)?;
    parts.remap_or_drop_bond_props(&mapping.bond_mapping)?;

    parts.invalidate(
        DerivedState::Adjacency
            | DerivedState::Rings
            | DerivedState::Valence
            | DerivedState::Stereo
            | DerivedState::Drawing
            | DerivedState::Fingerprint,
    );

    Ok(OpOutcome::Changed)
}
```

---

## 22. Hash and Snapshot Rules

Hash checks are optional and controlled by:

```text
op-contracts-hash
op-contracts-paranoid
```

Hash checks must obey these rules:

```text
NoOp must preserve semantic equivalence.
Changed does not require every hash to change.
Some real recomputations may produce identical hashes.
Hash changed does not prove chemical correctness.
Hash unchanged does not prove no work was done.
```

Hash checks are auxiliary.

Correctness still depends on:

```text
operation trace
cache generation
structural invariants
RDKit parity
I/O roundtrip
curated corpus
```

---

## 23. Required Implementation Shape

The final implementation must include both macros:

```text
molecule_ops!
  declares operation metadata, parameter lists, generated specs, public methods,
  and test registries

#[mol_op_body(op_name, parts)]
  binds an operation body to a registry entry and injects the capability-limited
  parts variable
```

Registered operations must not be implemented through manually written public
methods that bypass generated wrappers. A manually written operation wrapper can
silently drift from the registry and is therefore not allowed.

The only permitted hand-written public methods on `Molecule` are:

```text
read-only accessors
builder/editor entry points
methods explicitly classified as non-topology and non-derived-state-changing
```

If a future operation does not fit the macro model, the macro model must be
extended. Do not create a one-off exception without human-author approval.

---

## 24. No Direct Self Access in Operation Body

Operation bodies must not receive `self`.

They must only receive the injected `parts` variable.

This is forbidden:

```rust
#[mol_op_body(kekulize, parts)]
fn kekulize_impl(self: &Molecule) -> Result<OpOutcome, OperationError> {
    ...
}
```

This is also forbidden:

```rust
fn kekulize_impl(parts: &mut TopologyOnlyParts, mol: &Molecule) -> Result<OpOutcome, OperationError> {
    ...
}
```

The body must be capability controlled:

```rust
#[mol_op_body(kekulize, parts)]
fn kekulize_impl() -> Result<OpOutcome, OperationError> {
    ...
}
```

The macro provides `parts`.

---

## 25. Design Summary

The final architecture is:

```text
molecule_ops! registry
  = single source of truth

#[mol_op_body(op_name, parts)]
  = operation body binding

may_mutate
  = compile-time capability control

must_handle
  = dev/test responsibility contract

real action methods
  = perform work and automatically record handling

invalidate
  = real cache invalidation + generation bump

OpOutcome
  = explicit Changed / NoOp reporting

finish(outcome)
  = assemble new Molecule and check contract

feature tags
  = control check strength

release without op-contracts
  = no trace, no check, no snapshot, no hash cost
```

The shortest rule:

```text
Registry declares.
Macro restricts.
Parts mutate.
Real actions record.
Outcome explains.
Finish checks.
Release erases checks.
```

---

## 26. Final Principle

COSMolKit molecule operations must not rely on developer memory.

Every operation must be declared once in the registry.

Every operation body must be capability-limited.

Every dependent state must be handled by real action methods.

Every cache invalidation must update generation.

Every no-op must be explicit.

Every debug/test build must catch contract violations.

Every release build must pay only for real molecular work.
