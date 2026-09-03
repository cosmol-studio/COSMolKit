# Search and Model Extraction: Remaining Issues

Status: open review baseline, 2026-09-01

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
| M-02 | Query ownership | `QueryGraph` and predicates remain implemented in `core::search`, while `core::model::query` only re-exports them | Blocking |
| S-01 | Query-to-molecule projection | `QueryGraph::to_molecule()` remains a legacy internal projection used by pattern fingerprinting | High |
| S-02 | SMARTS parser path | parser still builds a temporary `Molecule`, performs cleanup, then converts it to `QueryGraph` | High |
| F-01 | Feature isolation | `QueryGraph::pattern_fingerprint()` is unconditional and breaks `--no-default-features` | High |
| Q-01 | Query wrappers | `QueryAtom`/`QueryBond` contain concrete `Atom`/`Bond`, implement `Deref`, and expose an always-present `query()` | Medium |
| Q-02 | Compiled query | `CompiledQuery` is currently a transparent wrapper and performs no compilation or planning | Medium |
| M-03 | Model invariants | public detached blocks and conformers have weak or absent local validation | Medium |
| F-02 | Capability graph | public feature definitions still leak implementation dependencies between domain capabilities | Medium |

## M-01: Two concrete molecule models coexist

### Evidence

The complete `Atom`, `AtomSpec`, `Bond`, and `BondSpec` value definitions now
live in `crates/cosmolkit-model/src/atom.rs` and `bond.rs`, including the
property, stereo, and PDB metadata carried by the former core model.
`AtomId` and `BondId` are owned by the same crate. Core's `model::atom` and
`model::bond` modules contain only aliases that specialize the model's opaque
query payload with the existing search-layer query node. There is no
core-owned duplicate concrete atom or bond implementation and no conversion
layer between the two identities.

The M-01 concrete value ownership slice is now complete. `cosmolkit-model`
owns the shared topology, coordinate, property, SGroup, stereo-group, and
topology-mapping values. Core retains `Molecule`, operation lifecycle, derived
cache, and source-specific remapping orchestration. Query payload decoupling is
tracked separately by M-02/Q-01; the current generic payload specialization
keeps this slice independent of a second concrete model.

### Why this was a real architectural problem

This was not merely an unfinished directory move. It created two incompatible
type universes with the same conceptual names. The alias-based migration now
gives downstream code and core algorithms one concrete atom/bond identity,
without a lossy adapter or a second implementation.

### Resolution questions

1. Which fields belong in the authoritative shared model, including property,
   stereo, PDB, and source-specific state?
2. Which fields are concrete value state and which remain runtime-owned?
3. Can core's current model be moved into `cosmolkit-model` in one vertical
   slice, with `Molecule` remaining in the facade, without introducing a
   second public representation?
4. Which builders and internal constructors must move with the values, and
   which must remain facade/runtime-only?

### Constraints for the fix

- Do not add bidirectional conversion traits as a permanent solution.
- Do not leave duplicate `Atom`/`Bond` implementations under different paths.
- Do not remove fields merely to make the types line up; preserve represented
  source behavior or explicitly record unsupported state.

## M-02: Query data still belongs to search implementation

### Evidence

`crates/cosmolkit-core/src/model/query.rs` only re-exports
`AtomQueryPredicate`, `BondQueryPredicate`, `QueryNode`, and query graph types
from `crate::search`.

The actual definitions remain in `search/query.rs` and
`search/query_graph.rs`. `QueryGraph` directly imports `Molecule`, `Atom`,
`Bond`, parser parameters, writer functions, matcher functions, and fingerprint
types.

### Why this matters

The current layering makes the apparent model module depend on search behavior
instead of the other way around. A crate that wants only query vocabulary still
pulls in the search implementation, and the model cannot become an independent
dependency without moving the actual definitions.

### Resolution questions

1. Which query fields and predicate variants are pure data and can move to the
   shared model without parser imports?
2. Should `QueryGraph` contain only graph data and local endpoint validation,
   with parser/matcher/writer behavior implemented by a search-side operator?
3. Which concrete atom fields must remain available in `QueryAtom` without
   reintroducing query state into `Atom`?

### Constraints for the fix

- Query data must not re-export its canonical definitions from search.
- Search behavior may wrap or operate on the model graph, but the model must
  not depend on search.
- Do not add a second query graph solely to bridge the old and new modules.

## S-01: Legacy `QueryGraph -> Molecule` projection remains

### Evidence

`QueryGraph::to_molecule()` in `search/query_graph.rs` rebuilds a core
`Molecule` from query-bearing atom and bond records. The method is crate-private,
but `QueryGraph::pattern_fingerprint()` calls it before invoking the existing
pattern fingerprint implementation.

### Why this matters

The query model was introduced specifically to prevent query semantics from
being represented as a concrete molecule. This projection restores that path
inside the implementation and requires the query graph to retain enough legacy
concrete state to rebuild a molecule. It is therefore a migration adapter, not
the target design.

### Resolution options to analyze

- Move pattern fingerprinting to consume query graph data directly, if the
  source algorithm requires query semantics.
- Move only a narrowly defined, private pattern-fingerprint input projection
  into the fingerprint implementation, with an explicit removal issue if the
  source algorithm genuinely requires concrete graph storage.
- Remove `QueryGraph::pattern_fingerprint()` from the query data type and expose
  it as a search/fingerprint operation.

The preferred direction is to keep data types free of domain operations and
avoid a general-purpose `to_molecule()` escape hatch.

## S-02: SMARTS parsing still uses a temporary concrete molecule

### Evidence

`mol_from_smarts()` currently performs:

```text
preprocess -> parse to temporary Molecule -> CX/name handling
-> merge query H -> stereo cleanup -> QueryGraph::from_query_molecule
```

The implementation and its local complexity note explicitly acknowledge the
builder-based reconstruction and additional cloning.

### Why this matters

This keeps the parser coupled to the concrete molecule lifecycle and makes the
first-class query graph a post-processing result rather than the parser's
canonical output. It also makes parser cleanup depend on concrete-molecule
helpers that were designed for a different state model.

### Resolution questions

1. Which parser reductions can write directly to `QueryGraph` data?
2. Which cleanup operations are query semantics and which are concrete-molecule
   normalization?
3. Can CXSMILES metadata and query-H merging be represented as model-level
   graph edits without constructing a `Molecule`?

The fix must preserve the RDKit call order and source markers; replacing the
temporary molecule with ad hoc heuristic parsing is not acceptable.

## F-01: Feature-disabled core does not compile

### Evidence

`QueryGraph::pattern_fingerprint()` is compiled unconditionally but references
`PatternFingerprintParams`, `Fingerprint`, `FingerprintError`, and the
`fingerprint` module, all gated behind the `fingerprints` feature.

Reproduced command:

```text
cargo check -p cosmolkit-core --no-default-features
```

It fails because those fingerprint items are unavailable.

### Resolution options to analyze

1. Gate the method and its imports with `#[cfg(feature = "fingerprints")]`.
2. Move the operation to a fingerprint/search operator module that is already
   feature-gated.
3. Make fingerprinting a mandatory dependency, which is contrary to the
   current fine-grained capability policy and should not be chosen without an
   explicit architecture decision.

The first or second option preserves feature isolation. The final API shape
should be decided together with S-01 rather than adding another wrapper method.

## Q-01: Query atom and bond wrappers retain compatibility baggage

### Evidence

`QueryAtom` stores a concrete `Atom` plus a predicate and implements `Deref<Target
= Atom>`. `QueryBond` does the same for `Bond`. Both expose `query()` methods
that always return `Some(&self.predicate)`.

### Why this matters

This mirrors the old query-bearing atom/bond model and makes query values look
like concrete values with an extra field. It also encourages consumers to rely
on concrete methods through deref instead of using an explicit query model.

### Resolution questions

- Should `QueryAtom` contain only query attributes and an explicitly modeled
  concrete projection, or should it contain a shared atom payload with query
  predicate data as separate fields?
- Which concrete metadata is semantically required by SMARTS/MCS and which is
  only present because of the legacy molecule representation?
- Should `query()` be removed in favor of the non-optional `predicate()` API?

Any change must preserve isotope, charge, radical, stereo, property, and dummy
atom matching semantics already covered by parity tests.

## Q-02: `CompiledQuery` is currently a naming wrapper

### Evidence

`CompiledQuery::compile()` only stores the original `QueryGraph`. Its `matches()`
method calls the same ordinary substructure matcher used by an uncompiled graph.

### Why this matters

The type promises a compiled or reusable match plan, but currently has no
precomputed candidate order, predicate lowering, adjacency view, or other
compiled state. This adds an abstraction name without changing behavior or
execution cost.

### Resolution options to analyze

- Remove `CompiledQuery` until a real source-backed compilation boundary is
  implemented.
- Keep it as an explicit immutable ownership wrapper, but rename/document it
  as such and avoid claiming compilation.
- Implement the actual RDKit/search-side preprocessing needed for repeated
  matching, with focused parity and performance validation.

Do not invent a heuristic optimizer merely to justify the type.

## M-03: Model-level local invariants are not enforced consistently

### Evidence

`TopologyBlock` exposes public `atoms`, `bonds`, and `adjacency` fields and
derives `Default`. `Conformer2D::new()` and `Conformer3D::new()` accept arbitrary
coordinate lengths without checking the associated atom count.

### Why this matters

The architecture document says the model owns local structural validity, but
these constructors permit detached values that are internally inconsistent.
Runtime validation cannot be the only protection if algorithms exchange model
values directly.

### Resolution questions

- Which invariants can be checked without knowing the owning molecule's
  operation lifecycle?
- Should blocks use private fields plus validated constructors, or public
  fields plus an explicit `validate()` boundary?
- Where should atom-count/conformer-count relationships be checked when the
  block does not itself carry atom count?

The eventual design must distinguish local structural invariants from
runtime cache and operation-contract invariants.

## F-02: Feature graph still leaks implementation dependencies

### Evidence

Current core features include relationships such as:

```toml
batch = ["dep:indicatif", "dep:rayon", "depict", "fingerprints", "io"]
depict = ["dep:glam", "dep:resvg", "dep:tiny-skia", "dep:usvg", "io"]
fingerprints = ["dep:rayon", "io"]
hashing = ["fingerprints"]
```

### Why this matters

Implementation reuse is currently expressed as public capability implication.
Enabling one feature can expose unrelated APIs and pull in unrelated domain
dependencies. This contradicts the intended separation between fine-grained
capabilities and user-facing bundles.

### Resolution questions

1. Which dependencies are genuinely shared lower-level implementations?
2. Which current feature edges only exist because code remains colocated in
   core?
3. Which bundles should be defined in `cosmolkit`, while core features remain
   independently selectable?

This should be addressed after the concrete model boundary is stable, because
moving implementation can change the true dependency graph.

## Recommended resolution order

The issues should be handled in this order:

1. Define the complete authoritative concrete model and eliminate M-01's
   duplicate `Atom`/`Bond` types.
2. Move query data definitions out of search and make M-02 real rather than a
   re-export facade.
3. Resolve S-01 and S-02 together so query operations no longer depend on a
   general `QueryGraph -> Molecule` projection.
4. Fix F-01 and decide the final ownership of pattern fingerprinting.
5. Reassess Q-01 and Q-02 after the data/behavior split is real.
6. Add model-local validation for M-03.
7. Rebuild the feature graph under F-02 once actual crate dependencies are
   known.

## Acceptance boundary

This report is closed only when:

- one concrete `Atom`/`Bond`/block model is used by core and exposed model
  consumers;
- query data has a canonical model definition independent of search behavior;
- SMARTS and matcher code operate on the canonical query model without a
  general legacy molecule projection;
- `cargo check -p cosmolkit-core --no-default-features` succeeds;
- no duplicate operation or query implementation remains;
- feature bundles do not accidentally expose unrelated domain capabilities;
- source-backed behavior and existing parity contracts remain intact.
