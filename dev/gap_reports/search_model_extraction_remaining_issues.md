# Search and Model Extraction: Remaining Issues

Status: open review baseline, 2026-09-04

This report records the implementation gaps found after introducing the
`cosmolkit-types` and `cosmolkit-model` crates and making `QueryGraph` a
first-class search value. It is a problem inventory, not an implementation
plan. Each item is intended to be analyzed and resolved separately.

The current worktree already contains unrelated and earlier migration changes.
This report does not classify those changes as part of the fixes below.

## Target boundary

The intended dependency and ownership boundary is:

```text
cosmolkit-types  -> public chemical vocabulary
cosmolkit-model  -> shared concrete values and local structural validity
cosmolkit-core   -> source-backed chemistry algorithms
cosmolkit        -> Molecule owner, operation runtime, and public assembly
```

For query functionality, the intended conceptual boundary is:

```text
QueryGraph / QueryAtom / QueryBond / predicates -> query data model
SMARTS parsing / query writing / compilation / matching -> search behavior
```

There must be one authoritative implementation for each concrete value and one
operation runtime. Temporary migration adapters are acceptable only when they
have an explicit removal decision and do not create a second public model.

## Priority summary

| ID | Area | Current state | Priority |
|---|---|---|---|
| M-01 | Concrete model ownership | Shared Atom/Bond, coordinate, property, SGroup, stereo-group, topology, and mapping values are owned by `cosmolkit-model`; runtime cache and Molecule remain in core | Resolved for this slice |
| M-02 | Query ownership | Query values and predicate vocabulary are now canonically defined by `cosmolkit-model`; core retains only search behavior and compatibility re-exports | Resolved for this slice |
| S-01 | Query-to-molecule projection | no production path projects `QueryGraph` to `Molecule`; query fixtures are query-native and concrete fixtures no longer manufacture query-bearing molecules | Resolved |
| S-02 | SMARTS parser post-processing | ordinary parsing, CXSMILES lowering for modeled records, and directional-bond stereo operate directly on `QueryGraph`; concrete-only CX records return structured unsupported errors | Resolved for current model scope |
| F-01 | Feature isolation | QueryGraph no longer exposes a projection-backed fingerprint facade; fingerprint implementation remains feature-gated and no-default core compilation passes | Resolved for this slice |
| Q-01 | Query wrappers | `QueryAtom`/`QueryBond` compose concrete values with mandatory predicates; legacy `Deref`, optional `query()` access, and query-bearing `Molecule` test fixtures are removed | Resolved |
| Q-02 | Compiled query | `CompiledQuery` now retains a pre-built VF2 query graph and RDKit node order; focused parity coverage is in place | Resolved for this slice |
| M-03 | Model invariants | topology and coordinate value validation now covers IDs, endpoints, adjacency, stereochemical references, group references, coordinate rows, and conformer IDs; runtime commit boundaries invoke it | Resolved for this slice |
| F-02 | Capability graph | public features no longer select other public domain features; colocated reuse is expressed through explicit `__*` implementation features | Resolved for this slice |
| R-01 | Runtime ownership | `Molecule`, `OpParts`, and the operation registry/runtime still live in `cosmolkit-core`, although the target owner is `cosmolkit` | High |

## M-01: Shared concrete model ownership

### Evidence

The complete `Atom`, `AtomSpec`, `Bond`, and `BondSpec` value definitions now
live in `crates/cosmolkit-model/src/atom.rs` and `bond.rs`, including the
property, stereo, and PDB metadata carried by the former core model.
`AtomId` and `BondId` are owned by the same crate. Core's `model::atom` and
`model::bond` modules contain only compatibility re-exports of those concrete
types. There is no core-owned duplicate concrete atom or bond implementation
and no conversion layer between the two identities.

The M-01 concrete value ownership slice is complete. `cosmolkit-model` owns
the shared topology, coordinate, property, SGroup, stereo-group, and
topology-mapping values. Core still retains `Molecule`, the operation
lifecycle, the derived cache, and source-specific remapping orchestration;
that runtime relocation is tracked separately by R-01. Query data ownership
is tracked by M-02 and wrapper cleanup by Q-01.

### Why this was a real architectural problem

This was not merely an unfinished directory move. It created two incompatible
type universes with the same conceptual names. The completed model extraction
now gives downstream code and core algorithms one concrete atom/bond identity,
without a lossy adapter or a second implementation.

### Remaining boundary

The concrete value extraction is complete. The remaining ownership question is
the relocation of `Molecule`, builders, and operation runtime, which is tracked
as R-01 rather than as another concrete-model implementation.

### Constraints for the fix

- Do not add bidirectional conversion traits as a permanent solution.
- Do not leave duplicate `Atom`/`Bond` implementations under different paths.
- Do not remove fields merely to make the types line up; preserve represented
  source behavior or explicitly record unsupported state.

## M-02: Query data ownership is resolved for this slice

### Resolution

Resolved for the current extraction slice. The canonical definitions now live
in `crates/cosmolkit-model/src/query.rs`: `QueryNode`, the atom and bond
predicate vocabulary, range queries, recursive query data, `QueryAtom`,
`QueryBond`, `QueryGraph`, and `QueryGraphError`. `cosmolkit-core/src/model/query.rs`
and `core::search` only re-export those values for source compatibility.

Search-side code retains `QueryGraphOperator`, parser/matcher/writer behavior,
and no concrete-molecule lowering helper. Remaining behavior work is tracked
by Q-02; none of these implementations makes the model depend on search.

### Why this mattered

Before this slice, the apparent model module depended on search behavior
instead of the other way around. A crate that wanted only query vocabulary
still pulled in the search implementation. The model now owns the query data
without importing parser, matcher, writer, valence, or molecule-runtime code.

### Remaining boundary

The remaining work is behavior-side: complete the search operator boundary and
compiled-query planning. Those issues are tracked by S-02 and Q-02; no second
query data model should be introduced.

### Constraints for the fix

- Query data must not re-export its canonical definitions from search.
- Search behavior may wrap or operate on the model graph, but the model must
  not depend on search.
- Do not add a second query graph solely to bridge the old and new modules.

## S-01: Remove legacy `QueryGraph -> Molecule` projection

### Evidence

The former `query_graph_to_molecule()` helper has been removed. SMARTS and
CXSMARTS writers, parser post-processing, matcher entry points, and production
fingerprint APIs consume `QueryGraph` or concrete `Molecule` values according
to their own input contract. Query fixtures construct the canonical query
graph directly. Concrete IO and chemistry fixtures no longer manufacture
query-bearing molecules; tests for unsupported concrete query records were
removed or reduced to their concrete-state assertions.

### Why this matters

The query model was introduced specifically to prevent query semantics from
being represented as a concrete molecule. The old test adapter restored that
path for query-target fingerprint fixtures, so it was a migration adapter and
not part of the target design. It is now gone. A future query-target Pattern
fingerprint parity surface would still require an explicit QueryGraph input
contract; it must not be implemented by recreating a projection helper.

### Resolution

The former `QueryGraphOperator::pattern_fingerprint()` facade and the
test-only projection helper were removed. Pattern and Layered fingerprints
retain their concrete-molecule contracts; query parsing, writing, and matching
remain query-native. Any future query-target fingerprint implementation must
accept `&QueryGraph` explicitly and own its query semantics without converting
to `Molecule`.

## S-02: Query-native SMARTS post-processing

### Evidence

The ordinary `parse_smarts()` path performs:

```text
preprocess -> parser state -> QueryGraphBuilder::finish() -> QueryGraph
```

Query-H merging, parser-state cleanup, and directional-bond stereo inference
are performed directly on the query graph. The directional path reproduces
RDKit's `setBondStereoFromDirections` traversal over query adjacency and
concrete bond state without constructing a `Molecule`.

CXSMILES extensions use the independent `cosmolkit-cx` syntax parser and a
QueryGraph-native lowerer. Query predicates, atom metadata, coordinates,
stereo, and radicals are applied directly to the query graph. Concrete-only
LN/SGroup/variable-attachment records return a structured `UnsupportedFeature`
error rather than being projected through a molecule. Plain suffix names retain
the existing `parse_name` and strict-mode behavior.

### Why this matters

The directional path is no longer coupled to the concrete molecule lifecycle.
The remaining work is parity expansion for less common CX records and a future
decision about query-native representations for concrete-only metadata; the
production path no longer has a concrete projection escape hatch.

### Resolution questions

1. Which CX extension records are valid query data and which remain concrete
   molecule metadata?
2. How should CX coordinates, enhanced stereo, labels, and query annotations
   be represented in the model without widening the `Molecule` boundary?
3. Which source cleanup operations are query semantics and which are concrete
   molecule normalization?

The completed directional migration preserves the RDKit call order and source
markers. Future CX work must use the same source-backed lowering discipline;
replacing the unsupported boundary with ad hoc heuristic parsing is not
acceptable.

## F-01: Feature-disabled core isolation

### Evidence

The former `QueryGraphOperator::pattern_fingerprint()` facade has been removed.
The remaining Pattern fingerprint implementation is independently gated behind
the `fingerprints` feature, so the core crate remains checkable without that
domain capability.

Verification command:

```text
cargo check -p cosmolkit-core --no-default-features
```

The command now passes; warnings do not indicate a feature-isolation failure.

### Remaining boundary

The gate preserves feature isolation. Query-specific fingerprint fixtures no
longer use a concrete projection; a direct query-target fingerprint contract
would be a separate future capability.

## Q-01: Query atom and bond wrappers use explicit projections

### Evidence

`QueryAtom` stores a concrete `Atom` plus a mandatory predicate and exposes the
two parts through explicit `atom()` and `predicate()` accessors. `QueryBond`
does the same for `Bond` and its mandatory predicate. The compatibility
`Deref` implementations and optional `query()` accessors have been removed.
Production search, parser, template, and CrystalFF consumers now use explicit
projections, and all test fixtures use the same boundary.

### Why this mattered

The former implicit deref mirrored the old query-bearing atom/bond model and
encouraged consumers to treat a query value as a concrete atom or bond. The
explicit projections make the data boundary visible at every call site while
retaining the concrete state required by RDKit-compatible matching.

### Resolution

The wrapper remains a deliberate data composition: query matching needs the
concrete atom/bond attributes as well as the query predicate tree, but neither
part is implicit through deref. `predicate()` is mandatory and consumers must
choose `atom()` or `bond()` when they need concrete state. This keeps the query
model explicit without introducing a second atom/bond representation.

Any change must preserve isotope, charge, radical, stereo, property, and dummy
atom matching semantics already covered by parity tests.

## Q-02: Compiled query planning

### Evidence

`CompiledQuery::compile()` now builds and retains the immutable VF2 query
adjacency view and the RDKit `SortNodesByFrequency` atom order. Its `matches()`
path passes both values into the canonical matcher, so repeated matches do not
rebuild the query-side graph plan.

### Why this matters

The plan is deliberately limited to source-backed preprocessing already used
by the VF2 implementation. It does not invent a heuristic predicate compiler
or alter result ordering. A focused test compares compiled and ordinary match
results on the same query and target.

### Resolution

Resolved for this slice. The compiled object owns one query graph plus its
detached VF2 topology view and stable source-defined node order. Recursive
query preparation and target-derived context remain per-match because they
depend on the target molecule and match parameters.

## M-03: Model-level local invariants

### Evidence

`TopologyBlock::validate()` now checks atom/bond table IDs, bond endpoints and
self-loops, stereo atom references, adjacency reconstruction, SGroup and
StereoGroup references, and parent/group IDs. `CoordinateBlock::validate_for_atom_count()`
checks coordinate row counts and duplicate conformer IDs independently for 2D
and 3D values. Constructors remain usable for detached source-shaped values;
validation is explicit because a conformer does not carry an atom count by
itself.

### Why this matters

These checks make model-level structural validity observable before an
algorithm hands a block to another layer. Runtime topology and coordinate
commit paths invoke the same checks, while cache validity and operation
contracts remain runtime-owned.

### Resolution

Resolved for this slice. Local checks live in `cosmolkit-model`; cross-block
atom-count checks are invoked by the runtime with the authoritative topology
size. No runtime cache, operation contract, or mapping authority was moved
into the model crate.

## F-02: Feature capability graph

### Evidence

The previous core feature graph included relationships such as:

```toml
batch = ["dep:indicatif", "dep:rayon", "depict", "fingerprints", "io"]
depict = ["dep:glam", "dep:resvg", "dep:tiny-skia", "dep:usvg", "io"]
fingerprints = ["dep:rayon", "io"]
hashing = ["fingerprints"]
```

### Why this matters

Implementation reuse was expressed as public capability implication. Enabling
one feature could expose unrelated APIs and pull in unrelated domain
dependencies, contrary to the intended separation between fine-grained
capabilities and user-facing bundles.

### Resolution

Resolved for this colocated implementation slice. `cosmolkit` public features
now select only their matching `cosmolkit-core` feature. Core records necessary
same-crate reuse with explicit `__*` implementation features, so hashing can
reuse fingerprint internals and batch can reuse IO/depict/fingerprint code
without enabling the corresponding public facade capabilities. These aliases
are intentionally temporary and will be removed when the implementation
crates are split.

## R-01: Molecule and operation runtime remain in core

### Evidence

The target architecture assigns the public `Molecule`, `MoleculeBuilder`,
`OpParts`, operation specifications, generated registries, and contract runtime
to `cosmolkit`. In the current workspace, `Molecule` and `MoleculeBuilder`
remain in `cosmolkit-core/src/model`, while `OpParts` and the operation runtime
remain in `cosmolkit-core/src/operations/ops`.

### Why this matters

The shared value model has been extracted, but the authoritative live-state
owner and operation lifecycle have not yet moved to the public runtime crate.
Until that migration is complete, `cosmolkit-core` is still more than a pure
algorithm crate and the dependency boundary described by the architecture
document is not yet enforceable by the crate graph.

### Constraints for the fix

- Keep exactly one `Molecule` representation and one operation runtime.
- Move ownership in vertical slices; do not add a second facade-side runtime.
- Preserve the existing public `cosmolkit` API while changing the internal
  owner and dependency direction.

## Recommended resolution order

The issues should be handled in this order:

1. Keep the completed M-01 model extraction as the base and do not reintroduce
   duplicate concrete values.
2. ~~Move query data definitions out of search and make M-02 real rather than a
   re-export facade.~~ Completed in the current extraction slice; the model is
   the canonical owner and core re-exports are compatibility-only.
3. ~~Keep the S-01 production boundary intact. Replace the remaining
   test-only Pattern/Layered fingerprint projection when the matcher gains a
   direct QueryGraph target boundary.~~ S-01 is complete: no QueryGraph to
   Molecule projection remains. S-02's modeled CX/stereo post-processing is
   complete.
4. ~~Fix F-01~~ Completed in the current extraction slice; decide the final
   ownership of pattern fingerprinting together with S-01.
5. ~~Resolve Q-01's compatibility wrappers.~~ Production paths now use the
   explicit `atom()`/`bond()` and mandatory `predicate()` APIs; migrate legacy
   test fixtures separately, then continue with Q-02.
6. Add model-local validation for M-03.
7. Rebuild the feature graph under F-02 once actual crate dependencies are
   known.
8. Migrate `Molecule` and the operation runtime under R-01 after the value and
   query boundaries are stable.

## Acceptance boundary

This report is closed only when:

- one concrete `Atom`/`Bond`/block model is used by core and exposed model
  consumers;
- query data has a canonical model definition independent of search behavior;
- SMARTS and matcher code operate on the canonical query model without a
  general legacy molecule projection;
- the target live-state owner and operation runtime are located in `cosmolkit`;
- `cargo check -p cosmolkit-core --no-default-features` succeeds;
- no duplicate operation or query implementation remains;
- feature bundles do not accidentally expose unrelated domain capabilities;
- source-backed behavior and existing parity contracts remain intact.
