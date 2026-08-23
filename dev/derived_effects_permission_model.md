# Derived Effects Permission Model

This document defines the current derived-state contract for registered
`Molecule` operations. It is a companion to
[`operation_system_standard.md`](./operation_system_standard.md).

The BioStructure operation framework has its own state model. In particular,
its current `must_handle` field is not part of the molecule derived-effects
model described here.

## Two Independent Axes

The operation system separates two questions:

1. Which molecule blocks may the operation read or write?
2. What must the operation do with affected derived state?

`access` answers the first question. `derived_effects` answers the second.
Neither declaration implicitly grants authority belonging to the other.

```text
access: {
    read:  [topology, derived_cache],
    write: [properties],
}

derived_effects: {
    recompute:         [aromaticity],
    preserve:          [],
    invalidate:        [drawing, fingerprint],
    operation_defined: [valence],
}
```

## Block-Level Read Access

An operation may read the derived-cache block only when its registry access
includes `derived_cache` as a readable or write-owned block. `OpParts` converts
the registry declaration into `MoleculeReadParts` block capabilities.

The current model has only block-level derived-cache read authority. It does
not independently authorize reads of individual cache entries such as rings,
valence, aromaticity, or stereo. Once a strict read view has
`derived_cache` access, every exposed entry in that block is readable.

Consequently:

- `preserve` is not read permission;
- `recompute` is not read permission and does not itself forbid reading;
- `invalidate` is not read permission;
- `operation_defined` is not read permission;
- per-cache-item read permissions do not currently exist.

Strict builds check the block capability. Default release builds use the same
operation and capability API path but omit the development-time access
assertions.

## Derived-Effect Categories

The only valid `derived_effects` categories are:

| Category | Required meaning | Permitted framework action |
|---|---|---|
| `recompute` | Produce a fresh framework-visible value, or explicitly clear it when the operation's source behavior leaves it unavailable | `set_*_cache`, validity-marker update, or `clear_cache` |
| `preserve` | Prove that the old value remains valid across the edit | An approved `prove_preserved` proof; no cache mutation authority |
| `invalidate` | Declare that the old value is stale | `clear_cache`; no cache write authority |
| `operation_defined` | Reproduce a source-required state transition that cannot be expressed truthfully as preserve, recompute, or invalidate | `set_*_cache`, validity-marker update, or `clear_cache`, as required by the source-port implementation |

The four categories must be pairwise disjoint. Strict `finish()` rejects an
overlapping declaration.

### `recompute`

`recompute` grants strict cache-write permission for the named state. It also
permits clearing that state when the reproduced source behavior does not leave
a materialized replacement.

After an operation has touched a block, strict finalization requires every
declared recompute state to have been updated or cleared through `OpParts`.
Writing a cache state not declared in `recompute` is a programming error.

### `preserve`

`preserve` states that an existing derived value remains valid after the
operation. It does not authorize reading or writing the cache.

After an operation has touched a block, strict finalization requires an
approved structural proof for every declared preserved state. For example,
`PreservationProof::LeafAtomAppend` checks the old graph identity and the shape
of appended atoms before rings or ring families may be retained.

### `invalidate`

`invalidate` states that an old derived value is stale and must be cleared.
After an operation has touched a block, strict finalization requires every
declared invalidated state to have been cleared through the framework.

Invalidation-only downstream products, including drawing and fingerprint
state, must not be marked recomputed unless the operation actually produces a
fresh framework-visible value.

### `operation_defined`

`operation_defined` is a narrow escape hatch for a source-required state
transition that does not fit the three standard mechanisms. It delegates the
transition mechanism, not the correctness obligation: the operation may use
the ordinary cache set, validity-update, and clear APIs, and strict
finalization still requires the declared state to be updated or cleared.
Correctness of the chosen transition remains established by the line-for-line
source port, focused regression tests, and the declared parity profile.

This category grants no derived-cache read authority. Read access continues to
come only from the operation's block-level `access` declaration.

The current contract permits exactly one use: `valence` in the hydrogen-removal
operation family. The macro rejects every other registry declaration, strict
runtime validation repeats the allow-list check, and a registry test locks the
two generated hydrogen-removal specifications to that single state. Adding a
second use requires an explicit design decision and corresponding updates to
all three guardrails.

## Compatibility View

`MoleculeOpSpec::needs_update()` is a derived compatibility view:

```text
needs_update() = recompute | invalidate | operation_defined
```

It is not an independent registry input.

The molecule operation system has no `must_handle()` compatibility view and no
`must_handle` or `require_handle` registry field. Those former concepts must
not be described as current molecule-operation contract dimensions.

## Unsupported Behavior

Unsupported behavior is not a derived-cache effect. An operation that cannot
reproduce a required source branch returns a structured operation error, such
as `OperationError::Unsupported` or `OperationError::UnsupportedFeature`.

This rule governs explicitly out-of-scope capabilities. It cannot be used to
classify a parity difference inside a supported operation boundary.

It must not encode unsupported behavior by clearing a cache, producing a
placeholder value, or adding an `unsupported` entry to `derived_effects`.

## Enforcement

The operation body records real actions through `OpParts`:

```text
set_*_cache(...) or validity-marker update -> updated cache trace
clear_cache(...)                           -> cleared cache trace
prove_preserved(...)                       -> preservation proof trace
```

With `op-contracts` enabled, cache mutation APIs check that the corresponding
effect permits the action, and `finish()` checks category disjointness and
completion of the declared obligations after a block was touched. For
`operation_defined`, the runtime verifies trace completion but does not infer
the chemically correct value; that obligation is unchanged from ordinary
source-port code that writes a recomputed cache.

The default release build executes the same wrapper, source-port body,
`OpParts` mutation methods, COW/in-place storage path, and cache updates. It
omits permission assertions and contract-trace validation. Optimization mode
must not select a different chemical or state-migration algorithm.

## Non-Operation Read APIs

Truly read-only APIs outside a registered operation may build or reuse derived
data according to their own documented contract. They do not receive mutation
authority merely because they can read a `Molecule`.

Public mutation-capable topology and coordinate APIs remain subject to the
operation registry and capability rules in the operation-system standard.
