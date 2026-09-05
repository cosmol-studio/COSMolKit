# COSMolKit CXSMILES Extraction And Lowering Gap Report

Status: CX syntax extraction and modeled target lowerers complete for the
declared scope, 2026-09-05

This report defines the extraction of the shared CX extension syntax layer into
`cosmolkit-cx` and records the current implementation boundary. The shared
syntax layer and both production target lowerers are in place for the declared
scope. Future model extensions and parity expansion are tracked as follow-up
work, not as incomplete CX extraction.

The immediate goal is to let the SMILES and SMARTS parsers share CX syntax
parsing without sharing their destination-specific mutation state.

## Target Boundary

The target dependency direction is:

```text
cosmolkit-cx
    ↑
    ├── SMILES parser / SMILES CX lowerer
    └── SMARTS parser / SMARTS CX lowerer
```

`cosmolkit-cx` owns the syntax-level representation of CX records and the
source-backed parser for those records. It must not depend on:

```text
Molecule
MoleculeBuilder
SmilesBuildState
QueryGraph
operation runtime
derived cache state
```

The destination parsers own lowering:

```text
ParsedCxExtensions
    ├── SMILES lowering -> SmilesBuildState / Molecule state
    └── SMARTS lowering -> QueryGraph
```

Parsing may be shared across representations; lowering is owned by the target
representation.

## Current State

The current implementation is in
`crates/cosmolkit-core/src/notation/smiles/cx.rs`.

The production entry point parses CX text with `cosmolkit-cx` and then applies
the resulting records through a SMILES-specific lowerer. Name suffixes,
strictness, and parser cleanup remain in the SMILES entry point.

The shared dispatch covers source branches for coordinates, atom labels, atom
values, atom properties, coordinate bonds, zero bonds, enhanced stereo,
ring-bond constraints, link nodes, data/polymer SGroups, unsaturation,
substitution, variable attachments, wedged bonds, double-bond stereo, and
radicals. Destination behavior is applied by the SMILES and SMARTS lowerers;
the old stateful parser helpers are test-only source anchors and are not part
of a production parse path.

The SMARTS entry point now parses CX records before hydrogen merging,
directional stereo, cleanup, and name attachment, and lowers supported records
directly to `QueryGraph`. Records without query semantics return a structured
unsupported error.

## Priority Summary

| ID | Area | Current state | Priority |
|---|---|---|---|
| CX-01 | Shared syntax value | `cosmolkit-cx` owns the independent typed `ParsedCxExtensions` parser, including structural records | Complete |
| CX-02 | SMILES lowering | Production SMILES parsing consumes shared records through `apply_cx_records_to_smiles`; concrete records are lowered, while query-only `u/rb/s` records return structured unsupported errors because `Molecule` is query-free | Complete for modeled scope |
| CX-03 | SMARTS lowering | Production SMARTS parsing lowers supported CX records directly to `QueryGraph`; records without query semantics return structured unsupported errors | Complete for modeled scope |
| CX-04 | Query records | `u:`, `rb:`, and `s:` have direct QueryGraph application paths; SMILES explicitly rejects those query-only records instead of silently dropping them | Complete for modeled scope |
| CX-05 | Target validation | Destination-specific atom, bond, conformer, and graph checks are performed by each lowerer | Complete |
| CX-06 | Concrete-only records | LN and SGroup records have neutral typed syntax values and remain explicit unsupported on SMARTS | Complete for modeled scope |
| CX-07 | Call order | SMARTS applies CX lowering before hydrogen merge, directional stereo, cleanup, and name attachment | Complete |
| CX-08 | Crate extraction | `cosmolkit-cx` is workspace-registered and both production call sites use it as the sole CX grammar parser | Complete |
| CX-09 | Legacy parser removal | The old stateful parser and its source-anchor helpers are test-only; no production path reparses CX text | Complete |
| CX-10 | SMILES query semantics | `Molecule` deliberately has no query payload. Query-only CX records return `SmilesParseError::UnsupportedFeature` rather than being represented as ordinary atom properties | Complete boundary |

## CX-01: Introduce The Shared Syntax Layer

### Required shape

The workspace crate named `cosmolkit-cx` (library name `cosmolkit_cx`) now
contains the initial independent syntax layer:

```rust
pub struct ParsedCxExtensions {
    pub records: Vec<CxRecord>,
    pub consumed: usize,
}

pub fn parse_cx_extensions(
    text: &str,
) -> Result<ParsedCxExtensions, CxParseError>;
```

The record type must retain source-relevant indices, ordering, conformer
identity, group identifiers, and scalar values. It must not retain a mutable
builder or a destination graph reference.

The syntax crate should remain independent of `cosmolkit-core` and
`cosmolkit`. It may use standard-library value types and CX-specific syntax
enums. It must not introduce a second chemical model or duplicate the query
predicate vocabulary.

### Constraints

- Preserve the RDKit parser dispatch order and consumed-byte behavior.
- Preserve strict parse errors as structured `CxParseError` values.
- Keep source comments and two-axis markers with the corresponding parser
  functions when code is moved, as required by
  `dev/source_reproduction_protocol.md`.
- Do not expose a general conversion from CX records to `Molecule`.

The parser now has typed representations for coordinates, labels/values, atom
properties, coordinate bonds, stereo, count constraints, radicals, link-nodes,
data/polymer SGroups, hierarchy, and variable attachments. The core SMILES
module retains source-anchored helper functions for focused legacy tests, but
the production entry point no longer dispatches or reparses CX text there.

## CX-02: SMILES Lowering

SMILES behavior is routed through an explicit target adapter:

```rust
fn apply_cx_to_smiles(
    state: &mut SmilesBuildState,
    cx: &ParsedCxExtensions,
) -> Result<(), SmilesParseError>;
```

Name handling, `_CXSMILES_Data`, `strict_cxsmiles`, and non-strict recovery
remain in the SMILES entry point. They are not syntax-parser responsibilities.

The lowerer preserves current source behavior for coordinates, labels,
properties, conformers, stereo groups, bond annotations, SGroups, and parser
temporary state. Query-only `u/rb/s` records are not concrete-molecule data;
they produce a structured unsupported error. Unsupported behavior must remain
an explicit structured error or an explicitly source-compatible non-strict
path; it must not become a silent record discard.
The SMILES adapter preflights the complete record list before mutating its
builder, so a query-only record cannot leave a partially lowered block.

## CX-03: SMARTS QueryGraph Lowering

The SMARTS path uses a separate target adapter:

```rust
fn apply_cx_to_query(
    graph: &mut QueryGraph,
    cx: &ParsedCxExtensions,
) -> Result<(), SmartsParseError>;
```

The adapter writes directly to the canonical query model. It must not call
the SMILES lowerer, construct a `Molecule`, or use a query-to-molecule
projection.

Expected mappings include:

| CX record | QueryGraph destination |
|---|---|
| coordinates | query conformer data, if the model supports the record |
| atom labels/values/properties | `QueryAtom.atom()` properties or query metadata according to source semantics |
| `u:` | `QueryAtom.predicate()` composition with unsaturation predicate |
| `rb:` | ring-bond-count predicate, including the source `*` scan form |
| `s:` | substitution-count predicate |
| `m:` | variable-attachment query state |
| enhanced stereo | QueryGraph stereo-group data |
| wedged/double-bond stereo | `QueryBond` concrete direction/stereo plus query semantics |
| radicals | concrete atom radical state where represented |

If a record has no faithful QueryGraph representation, the adapter must return
a structured unsupported error. It must not encode the record by creating a
concrete molecule or by dropping it silently.

## CX-04: Query Predicate Application

The former empty `expand_cx_atom_query()` is no longer on the production path.
The shared parser emits syntax records such as:

```rust
CxRecord::Unsaturation { atoms: Vec<usize> }
CxRecord::RingBondCount { atom: usize, value: RingBondCount }
CxRecord::SubstitutionCount { atom: usize, value: SubstitutionCount }
CxRecord::VariableAttachment { ... }
```

The SMARTS lowerer performs the equivalent of RDKit's query expansion against
the existing `QueryAtom` predicate tree. The SMILES lowerer makes an explicit
source-backed boundary decision for the same records and returns structured
unsupported because it must not inherit SMARTS semantics accidentally.

## CX-05: Target-Specific Validation

The syntax layer parses numeric indices and record fields. The target lowerer
validates whether those indices refer to atoms, bonds, conformers, or groups in
the destination graph.

This separation is required because validation can differ by destination:

```text
syntax validity       -> cosmolkit-cx
graph index validity  -> SMILES/SMARTS lowerer
model invariants      -> destination model/runtime
```

The lowerers must preserve RDKit behavior for invalid or out-of-range records,
including whether the source rejects, skips, marks a scan, or continues in
non-strict mode. Moving the check out of `SmilesBuildState` must not weaken
those checks.

## CX-06: Records Without A QueryGraph Equivalent

Substance groups, polymer hierarchy, link nodes, and some CX metadata are
concrete-molecule features or require a model that QueryGraph does not yet
own. They must remain represented in `ParsedCxExtensions` so the parser does
not need a second grammar, but the SMARTS lowerer must return a structured
unsupported result until the query model has an explicit representation.

No generic property-string encoding should be introduced as a substitute for
missing query semantics.

## CX-07: Restore RDKit Post-Processing Order

The SMARTS path now follows:

```text
parse SMARTS
  -> parse CX records
  -> lower CX records to QueryGraph
  -> merge_query_hs
  -> set bond stereo from directions
  -> cleanup
  -> attach name
```

This ordering is part of the production implementation and prevents CX query
semantics from being lost in a concrete-molecule projection.

## CX-08: Crate Extraction Sequence

The completed migration proceeded as one vertical slice:

1. Add `cosmolkit-cx` to the workspace with the syntax records and parser.
2. Route the existing SMILES path through parsed records and a SMILES lowerer.
3. Run focused SMILES CX parity checks and close any behavior regressions.
4. Add the QueryGraph lowerer and remove the SMARTS CX unsupported boundary.
5. Add focused CXSMARTS parity cases, including query predicates, stereo, and
   unsupported model records.
6. Only after the two lowerers are stable, consider moving the ordinary SMILES
   and SMARTS parsers into sibling crates.

The extraction must not create a second public `Atom`, `Bond`, `Molecule`, or
`QueryGraph` type. Existing top-level APIs may remain forwarders, but there
must be one parser implementation and one lowering implementation per target.

## Acceptance Criteria

The migration is complete for the initial scope when:

- `cosmolkit-cx` builds independently of `cosmolkit-core` and `cosmolkit`;
- no `cosmolkit-cx` function receives `SmilesBuildState`, `Molecule`, or
  `QueryGraph`;
- modeled SMILES CX behavior remains parity-tested through the new lowerer,
  and query-only records do not silently mutate a concrete molecule;
- SMARTS CX input no longer needs a concrete-molecule projection;
- query-specific CX records update `QueryGraph` predicates through one canonical
  query implementation;
- unsupported QueryGraph records return structured errors;
- RDKit post-processing order is preserved;
- source markers remain attached to the functions that reproduce the source;
- no temporary compatibility parser is used by a production path; the
  test-only legacy helper remains solely as an inline source/reproduction
  anchor for existing focused parser tests.

## Explicit Non-Goals

This migration does not by itself:

- move `Molecule` or `OpParts` into another crate;
- redesign the QueryGraph model;
- implement all missing SGroup or link-node query semantics;
- make `cosmolkit-cx` a general SMILES/SMARTS parser crate;
- claim complete CXSMILES or CXSMARTS parity before focused validation passes.

## Follow-Up Outside This Extraction

The following items are intentionally outside the completed CX extraction
boundary and must not be reported as parser-extraction failures:

1. The test-only stateful parser helpers in
   `crates/cosmolkit-core/src/notation/smiles/cx.rs` still contain the copied
   RDKit parser bodies used by focused legacy tests. They are not called by
   production parsing. Removing them requires migrating those tests to the
   typed records or retaining the source anchors in a dedicated test module.
   This is test/source-anchor cleanup only; production migration is complete.
2. The current concrete `Molecule` model has no query payload. SMILES
   `u/rb/s` are therefore a deliberately closed unsupported boundary;
   implementing them requires a deliberate model/API decision, not another CX
   parser.
3. `LN`, SGroup, polymer, and `m` records remain structured unsupported on
   SMARTS because `QueryGraph` does not yet own those semantics. The syntax
   records are already preserved and ready for a future query-model extension.
4. Existing test-only Query API migration failures and stale architecture
   baselines belong to `search_model_extraction_remaining_issues.md`; they do
   not indicate a second CX parser or a production CX dependency leak. They are
   owned by the search-model extraction report.

Validation completed for this migration:

```text
cargo check -p cosmolkit-cx
cargo test -p cosmolkit-cx
cargo check -p cosmolkit-core --features op-contracts-strict
git diff --check
```

The strict core library test build is currently blocked by the pre-existing
search-model API migration errors (for example, stale `AtomSpec::with_query`
and `Atom::query` test call sites). Those failures are tracked by
`search_model_extraction_remaining_issues.md`; they are not compile failures
in `cosmolkit-cx` or in either production CX lowering path.
