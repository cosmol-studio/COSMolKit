# RDKit SDF Query-Graph Contract Migration

## Outcome

`cosmolkit-io` now owns an explicit detached result boundary for query-bearing
V2000 and V3000 MolBlocks:

- `read_mol_block_detached()` returns `MolBlockRecord::{Concrete, Query}`;
- `read_sdf_graph_record_detached()` adds SDF data fields without lowering the
  graph to concrete topology;
- `QueryMolBlockRecord` carries the canonical `cosmolkit-model::QueryGraph`,
  typed substance groups, molecule properties, and the source coordinate
  dimension; and
- the legacy concrete readers remain fail-closed and direct callers to the
  query-aware APIs.

No query predicate is coerced into an `AtomSpec`, `BondSpec`, or concrete
`Molecule`.

## Source-Backed Coverage

The detached owner now preserves the modeled RDKit query state from:

- V2000 query atom symbols (`*`, `A/AH`, `Q/QH`, `X/XH`, `M/MH`, `R`, `R#`,
  numbered R groups, `L`, and `LP`);
- V2000 old-style and `M  ALS` atom lists, including new-list replacement of
  old atomic-number predicates while preserving other query components;
- V2000 R-group labels (`M  RGP`), ring-bond counts (`M  RBC`), substitution
  counts (`M  SUB`), and unsaturation (`M  UNS`);
- V2000 bond types 5--8, unknown bond types using RDKit's any-query fallback,
  and ring/non-ring topology constraints;
- V3000 atom lists, negated lists, complex atom symbols, wildcard/R-group
  symbols, `RGROUPS` label/isotope/null-query state, `HCOUNT`, `RBCNT`,
  `UNSAT`, and query expansion by `CHG`/`MASS`;
- V3000 bond query types and `TOPO`; and
- coordinates, enhanced-stereo collections, typed SGroups, MolBlock
  properties, and SDF data fields attached to query records.

The copied RDKit parser bodies remain beside the Rust implementations with
two-axis status markers. The new query-property helpers remain behavior-axis
`❗` until broader fixture-oracle parity is run; local inspection supports
equivalent asymptotic behavior and allocation shape for the modeled records.

## Explicit Residual Boundaries

- V2000 `M  MRV SMA` remains structured unsupported in `cosmolkit-io` because
  compiling its recursive SMARTS requires a deliberate `cosmolkit-search`
  integration boundary; `cosmolkit-io` does not depend on the search crate.
- V3000 template attachment-order lists and `OBJ3D` remain
  structured unsupported where the detached model lacks the required state.
- Query MolBlock writing is not claimed by this slice. Existing V2000/V3000
  writers accept concrete `TopologyBlock` values only.
- Runtime chemistry finalization and public `Molecule` construction are a
  narrow adapter concern. A query record is preserved by the detached API and
  rejected structurally by the concrete runtime API instead of being lowered.

## Validation

```text
cargo check -p cosmolkit-io
cargo test -p cosmolkit-io --release
75 passed, 0 failed
```

Focused cases cover V2000 and V3000 atom/bond query trees, query-scan marker
preservation, as-drawn substitution degree, invalid unsaturation values,
unknown-bond any fallback, concrete-reader rejection, coordinate dimension,
and query SDF data-field propagation.
