# Crate Architecture

This document defines the current target architecture for the COSMolKit
workspace. It is a normative design document for new crate boundaries and for
future migrations of existing implementation code. Historical plans in
`dev/plans/` and `dev/archive/` are not rewritten by this document.

The concrete final crate decomposition and the complete layer diagram are
defined in [`final_target_architecture.md`](./final_target_architecture.md).

## 1. Ownership Model

`cosmolkit` is the public API and molecule-runtime crate. It is intentionally
more than a thin facade. It owns the stateful object and the rules required to
keep that object valid:

- public `Molecule` and `MoleculeBuilder` APIs;
- the private operation runtime, including `OpParts`;
- operation specifications, generated registries, and contract matrices;
- capability checks, mapping checks, derived-state effects, and invariants;
- extraction and commit of detached working blocks;
- public value-style and trailing-underscore in-place wrappers;
- facade-level re-exports and feature selection.

`OpParts` belongs in `cosmolkit` for ownership purposes, but it is not a
general public capability object. It must remain private or `pub(crate)` and
must only be reachable by generated operation wrappers and runtime code.
Users must use `Molecule` APIs, never construct or commit `OpParts` directly.

`cosmolkit-core` is the source-backed algorithm crate. It provides chemistry
and structural algorithms that operate on explicit model values and return
typed values, assignments, reports, or transformed blocks. It must not depend
on `cosmolkit`, `Molecule`, `OpParts`, the operation registry, or the runtime
cache authority.

The core crate is a capability layer, not a one-feature-per-crate rule. Tightly
coupled foundational chemistry implementations may remain together there,
including valence calculation, sanitization, and hydrogen transforms, when
they share model-level machinery. A feature such as hydrogen handling therefore
does not need a standalone crate. Separately reusable families such as
fingerprints, descriptors, search, or file IO may be sibling domain crates;
their dependency on core is determined by actual algorithm reuse, not by a
public feature alias.

The name `cosmolkit-core` is retained for package continuity; its role is no
longer "the crate that owns Molecule". It is a reusable implementation crate,
similar in dependency direction to a fingerprint or search crate, and may be
depended on by other domain crates when the dependency is acyclic and
explicit.

## 2. Dependency Direction

The target dependency graph is:

```text
Python bindings and other language adapters
                    |
                    v
               cosmolkit
                /     \
               v       v
   cosmolkit-core   cosmolkit-model
          |
          v
   cosmolkit-model
```

The diagram is illustrative rather than an exhaustive tree. The normative
rule is that all model and implementation crates form an acyclic dependency
graph below `cosmolkit`; dependencies point toward value types and lower-level
algorithms. `cosmolkit-core` must never point upward to `cosmolkit` to obtain a
`Molecule` or operation capability. This prevents a dependency cycle and makes
the runtime boundary enforceable by the type graph.

`cosmolkit` may depend on both `cosmolkit-model` and `cosmolkit-core`. The
facade crate is therefore the place where a `Molecule` operation is assembled:
it extracts authorized values, calls the algorithm crate, validates the
result, records contract effects, and commits the result.

`cosmolkit-core` must not depend on `cosmolkit-cx`. CX is shared notation syntax
owned by the SMILES and SMARTS/search layers; only those target-representation
crates parse `ParsedCxExtensions` and lower them into a molecule or query. Any
current core-to-CX edge is migration residue from colocated notation code and
must be removed when those parsers are extracted.

## 3. Model Crate

`cosmolkit-model` contains stable value types shared by algorithm crates. These
types enforce local structural invariants, but they do not own molecule-runtime
lifecycle authority. The initial scope is:

- atom, bond, element, bond-order, and identifier value types;
- query AST and primitive query predicates needed by atom and bond values;
- adjacency and graph value structures;
- 2D and 3D conformer value types;
- topology, coordinate, and ordinary molecule-property blocks;
- validated value-level mappings and assignments where they are independent
  of the runtime lifecycle.

The model crate must not own:

- `Molecule` or its `Arc`-backed state ownership;
- `OpParts`, operation specs, registries, capability trackers, or traces;
- derived-cache validity and invalidation authority;
- computed-property cache locks;
- operation commit or rollback/invalid-state policy;
- source-library dispatch or facade-specific compatibility policy.

Detached model values may be edited by algorithms. Such editing must not
provide a path to mutate or replace the state of a live `Molecule`, and the
model API must preserve local structural invariants before a block is returned
or passed to another layer. Cross-crate algorithms receive immutable views or
owned detached blocks and may use source-shaped batch editing where the model
validates the resulting block. Moving a type to `cosmolkit-model` does not
authorize `Molecule::topology_mut()`, unrestricted live-state replacement, or
runtime commit access.

`DerivedCacheBlock` is deliberately excluded from the initial model scope.
It is contract-managed state, not an ordinary algorithm input. It remains
owned by the molecule runtime until a separate cache-state design explicitly
proves that moving it cannot widen mutation or read authority.

## 4. Algorithm Boundary

Domain algorithms in `cosmolkit-core` and future sibling crates must expose
the narrowest input that matches the source behavior and the operation access
declaration:

```rust
pub fn calculate(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
) -> Result<Descriptor, Error>
```

For a topology-authorized transformation, an owned detached block is preferred:

```rust
pub fn transform(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    params: &Params,
) -> Result<TopologyTransform, Error>
```

Read-only algorithms receive immutable views or slices. Assignment-producing
algorithms return assignments. Algorithms that are authorized to transform a
whole block may return an owned transformed block, together with the mapping
or report required by the parent contract. They must not return a complete
`Molecule` and must not perform runtime cache or capability bookkeeping.

The signature is a compile-time data-visibility boundary, not a replacement
for semantic correctness. The parent runtime still enforces registry access,
mapping obligations, derived effects, preconditions, source-preservation rules,
unsupported errors, and final invariants.

Search follows the same boundary while its implementation remains in the
current core package. Predicate evaluation is generic over a read-only
`SearchTargetAccess` view whose data comes from model blocks and explicit ring
or valence assignments. The current `Molecule` implementation of that trait is
only a migration adapter; it is not an algorithm-owned dependency and must be
removed when search moves to its own crate. New search code must use the
detached view or the block-level constructor, never add another `Molecule`
parameter.

### Value and behaviour facades

When a domain value has several operations that interpret it, keep one
canonical value type and place the domain behaviour behind an operator facade.
The value type may retain methods for local construction, accessors, and local
validation, but it must not acquire dependencies on parsers, matchers,
serializers, planners, or runtime state merely to expose those behaviours.

The facade borrows the canonical value and groups related operations:

```rust
// model/value crate
pub struct QueryGraph { /* data and local validation */ }

// domain implementation crate
pub struct QueryGraphOperator<'a> {
    inner: &'a QueryGraph,
}
```

This rule applies to query graphs and to future data/behaviour pairs with the
same dependency shape. It does not authorize a second wrapper type with the
same public domain name, a duplicate value model, or long-lived conversion
adapters. Construction that returns an owned value remains an explicit domain
entrypoint, while interpretation operations are methods on the operator
facade. For example, SMARTS construction is `search::parse_smarts(...) ->
QueryGraph`; it must not be an associated function on `QueryGraphOperator`
because no existing `inner` value participates in parsing. Existing top-level
facade APIs may remain as thin forwarders to preserve source compatibility.
For example, `QueryGraph::to_smarts()` may delegate to
`self.operator().to_smarts()` without making the model crate depend on the
SMARTS writer. Such forwarders belong to the top-level or domain facade, never
to the lower-level canonical value implementation.

### Parser construction boundary

SMARTS parsing uses a private construction state, not a second public query
model and not a query-bearing `Molecule` as its canonical intermediate:

```text
SMARTS
  -> parser state + QueryGraphBuilder
  -> one controlled lowering (`finish`)
  -> QueryGraph
```

`ParserState` owns syntax bookkeeping such as branch stacks, pending bonds, and
ring-closure tables. `QueryGraphBuilder` owns only query graph construction
state and source metadata needed to construct the final value. Ring-closure
records are consumed when their bond is emitted and are not retained as a
second graph representation. `finish` validates index alignment, assigns
stable IDs, builds adjacency, and returns the sole `QueryGraph` value.

RDKit's parser creates query atoms and bonds directly in an `RWMol`; the builder
is therefore a COSMolKit construction aid for preserving the same semantic
action order, not a claim that RDKit has a `ParsedSmartsGraph` type. Source
post-processing must preserve RDKit's ordering and behavior. Any remaining
adapter that consumes a concrete molecule is an explicitly bounded integration
dependency and must not become a public parser/model conversion path.

## 5. Parent Operation Flow

Every public molecule operation follows this shape:

```text
Molecule
  -> cosmolkit runtime extracts authorized owned blocks
  -> cosmolkit-core (or a domain crate) runs the source-backed algorithm
  -> cosmolkit runtime validates mappings, effects, and invariants
  -> cosmolkit runtime commits the complete result
  -> Molecule
```

An algorithm error must not install a partially transformed real molecule.
Passing detached owned values preserves this property for the algorithm call;
the existing in-place error and panic policy remains the responsibility of
the `cosmolkit` runtime and must not be delegated to a child crate. Contract
obligations and final invariants must be validated before the transformed state
becomes the authoritative, externally observable `Molecule` state.

`cosmolkit` is the sole crate allowed to own or accept a live `Molecule`.
`cosmolkit-core` and every other domain crate operate on detached
`cosmolkit-model` values and return detached values, assignments, matches, or
structured errors. They must not accept `Molecule`, `MoleculeBuilder`,
`OpParts`, derived-cache blocks, or operation runtime types. Any crate-private
adapter left during migration is transitional and must be removed before the
architecture is complete; it cannot contain a second algorithm
implementation.

### Operation capability projection

The generated operation declaration is the single source for both runtime
metadata and compile-time block visibility. It generates a private zero-sized
access marker and specialises the one transaction type directly:

```rust
#[mol_op_body(with_hydrogens, parts)]
fn with_hydrogens_impl(
    parts: &mut OpParts<'_, WithHydrogensAccess>,
    params: AddHsParams,
) -> Result<(), OperationError> {
    let (topology, coordinates, properties) = parts.extract_all_writable()?;
    let (topology, coordinates, properties) =
        cosmolkit_core::hydrogens::transform(topology, coordinates, properties, &params)?;
    parts.commit_all_writable(topology, coordinates, properties)
}
```

There is no per-operation context wrapper around `OpParts`. The marker is the
capability projection: generated inherent methods expose only the declared
read and write blocks. Runtime internals may retain unrestricted helpers, but
they remain private and cannot be passed to an algorithm crate.

Multiple-output operations use the same rule with one distinct runtime object:

```rust
#[mol_multi_op_body(tautomers, parts)]
fn tautomers_impl(
    parts: &mut MultiOutputOpParts<'_, TautomersAccess>,
    params: TautomerOptions,
) -> Result<(), OperationError> {
    let candidates = cosmolkit_core::tautomer::enumerate(
        parts.topology()?,
        parts.coordinates()?,
        parts.properties()?,
        &params,
    )?;
    parts.emit_all(candidates)
}
```

The algorithm returns an ordered
`Vec<(TopologyBlock, CoordinateBlock, MoleculeProperties)>`, not a
`Molecule`, draft object, branch handle, or `OpParts`. `MultiOutputOpParts`
validates every candidate, applies the operation’s contract effects, and only
then constructs each public output molecule. It intentionally has no generic
branch-id or source-derivation abstraction: candidate enumeration belongs to
the algorithm, while authoritative installation belongs to the runtime.

## 6. Migration Constraints

This is a staged architecture migration, not a directory move.

1. Extract model values and the query AST without changing chemistry behavior.
2. Move one low-coupling algorithm family behind a value-based API.
3. Implement the corresponding `cosmolkit::Molecule` wrapper and preserve the
   existing public behavior and operation contract.
4. Add source markers and parity tests at the new algorithm boundary.
5. Only after the vertical slice is stable, migrate additional domains.

The query AST must be lowered before moving `Atom` and `Bond`, because those
types currently reference query types. Search parser/matcher code may then
depend on the model crate without creating a model-to-search cycle.

Conformer generation, distance geometry, force fields, alignment, and tautomer
enumeration have coupled dependencies and should be migrated after the model
and one simpler algorithm family establish the boundary. Feature flags remain
useful for selecting optional algorithms and dependencies, but they are not a
substitute for the crate-level ownership boundary.

Feature selection follows two distinct layers:

- Fine-grained features gate one public capability and its optional
  implementation dependency, such as `smiles`, `tautomer`, or `conformer`.
- User-facing bundles such as `common_api`, `chemistry_api`, `3d_api`, and
  `full` are explicit compositions of fine-grained features.

A fine-grained feature must not depend on another domain's public feature just
  because its implementation happens to reuse code. Shared implementation
  belongs in a lower-level internal crate or an explicit algorithm dependency;
  it must not make unrelated facade methods appear. The operation registry may
  retain complete contract metadata for review, while generated public
  wrappers, operation bodies, and optional domain dependencies are gated by
  their owning fine-grained feature.

While implementations remain colocated in `cosmolkit-core`, Cargo's lack of
private features is handled with double-underscore implementation features
(`__io_impl`, `__fingerprint_impl`, `__depict_impl`, `__hashing_impl`, and
`__batch_impl`). Public features select only their own implementation feature;
an internal edge records shared code reuse without enabling the other domain's
public feature or its facade re-exports. These aliases are migration
scaffolding and should disappear when the corresponding implementation crates
are separated.

## 7. Public Compatibility

The supported user entry point is `cosmolkit`. New documentation, examples,
Python bindings, and integration tests must import `Molecule` and public
operations from that crate. `use cosmolkit_core::Molecule` is not a supported
target architecture and must not be preserved by making `cosmolkit-core`
depend on the facade.

During migration, compatibility shims may exist only as deliberate, temporary
transitions with an explicit removal decision. They must not duplicate a
second `Molecule` implementation or a second operation runtime.

## 8. Non-Negotiable Invariants

- There is exactly one authoritative `Molecule` representation.
- There is exactly one operation runtime and registry, in `cosmolkit`.
- Child algorithm crates cannot obtain `Molecule`, `OpParts`, capabilities, or
  derived-cache authority.
- Public in-place molecule methods retain the trailing-underscore convention.
- Contract validation and invariant checks happen in the parent runtime before
  commit.
- Source ports remain source-backed and follow
  `dev/source_reproduction_protocol.md`.
