# Crate Architecture

This document defines the current target architecture for the COSMolKit
workspace. It is a normative design document for new crate boundaries and for
future migrations of existing implementation code. Historical plans in
`dev/plans/` and `dev/archive/` are not rewritten by this document.

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
