# SMILES RDKit Remaining Source Scan

Date: 2026-05-14

## Result

Frozen-scope coverage is **not** at `100%`.

`dev/smiles_rdkit_full_port_checklist.md` Step 167 requires a zero-gap rewrite
only if every frozen-baseline function is fully covered. That condition is not
met in the current tree, so this report records the remaining frozen-scope gaps
instead of claiming closure.

## Frozen Baseline Used For Audit

This audit uses the frozen baseline defined by the current checklist pass:

- shared parser/core helpers from SMILES parse and postprocessing
- parser functions in `crates/cosmolkit-core/src/notation/smiles.rs`
- writer functions in `crates/cosmolkit-core/src/notation/smiles_write.rs`
- canon-ranking functions in `crates/cosmolkit-core/src/notation/canon_rank.rs`
- ring functions in `crates/cosmolkit-core/src/chemistry/rings.rs`
- kekulization functions in `crates/cosmolkit-core/src/chemistry/kekulize.rs`
- valence/sanitize functions in `crates/cosmolkit-core/src/chemistry/valence.rs`
  and `crates/cosmolkit-core/src/operations/ops.rs`

The source files compared remain:

- `third_party/rdkit/Code/GraphMol/SmilesParse/SmilesParse.cpp`
- `third_party/rdkit/Code/GraphMol/SmilesParse/SmilesParseOps.cpp`
- `third_party/rdkit/Code/GraphMol/SmilesParse/CXSmilesOps.cpp`
- `third_party/rdkit/Code/GraphMol/SmilesParse/SmilesWrite.cpp`
- `third_party/rdkit/Code/GraphMol/Chirality.cpp`
- `third_party/rdkit/Code/GraphMol/Canon.cpp`
- `third_party/rdkit/Code/GraphMol/new_canon.cpp`
- `third_party/rdkit/Code/GraphMol/Kekulize.cpp`
- `third_party/rdkit/Code/GraphMol/FindRings.cpp`
- `third_party/rdkit/Code/GraphMol/MolOps.cpp`
- `third_party/rdkit/Code/GraphMol/MolOps.h`

## Audit Standard

For this frozen scope, `100% coverage` requires all of the following:

- every frozen-baseline function is behaviorally closed for the currently
  modeled input state space,
- copied-source markers are no longer left at unresolved or known-gap status
  inside those frozen functions,
- targeted tests exist where the checklist required them,
- the final source rescan finds no remaining direct function/helper gap within
  the frozen scope.

Any remaining `RDKit❌❌`, `RDKit❗❗`, `RDKit✔️❌`, `RDKit❗✔️`, or `RDKit✔️❗`
marker in the frozen-scope implementations is evidence that zero-gap closure
has not been reached.

## Current Remaining Gaps

### 1. Parser path is not frozen-scope complete

`crates/cosmolkit-core/src/notation/smiles.rs` still contains unresolved marker
classes inside frozen-scope parser and helper implementations:

- `RDKit❌❌`: 1
- `RDKit❗❗`: 2
- `RDKit✔️❌`: 14
- `RDKit❗✔️`: 2037

This means the parser side cannot be described as fully covered.

### 2. Writer path is not frozen-scope complete

`crates/cosmolkit-core/src/notation/smiles_write.rs` still contains unresolved
marker classes inside frozen-scope writer implementations:

- `RDKit❗❗`: 2
- `RDKit❗✔️`: 582

The checklist-targeted writer helpers added in recent steps are covered, but
broader writer orchestration and copied-source closure are still incomplete.

### 3. Canon ranking is not frozen-scope complete

`crates/cosmolkit-core/src/notation/canon_rank.rs` still contains unresolved
marker classes inside the frozen canon-ranking set:

- `RDKit✔️❌`: 119
- `RDKit❗✔️`: 83

So the canon-ranking block is not in a zero-gap state.

### 4. Ring perception is not frozen-scope complete

`crates/cosmolkit-core/src/chemistry/rings.rs` still contains unresolved marker
classes inside the frozen ring-perception functions:

- `RDKit❌❌`: 7
- `RDKit❗✔️`: 78

The explicit non-URF branch remains outside the frozen-scope closure standard,
and aggregate `findSSSR` / `symmetrizeSSSR` copied-source closure is still open.

### 5. Kekulization is not frozen-scope complete

`crates/cosmolkit-core/src/chemistry/kekulize.rs` still contains unresolved
marker classes inside the frozen kekulization functions:

- `RDKit❗✔️`: 365

Recent checklist work closed several targeted helpers and regressions, but the
full frozen kekulization block is not yet marker-closed.

### 6. Valence/sanitize orchestration is not frozen-scope complete

`crates/cosmolkit-core/src/chemistry/valence.rs` is close but not closed:

- `RDKit❗✔️`: 13

`crates/cosmolkit-core/src/operations/ops.rs` still carries unresolved frozen
sanitize/cleanup orchestration markers:

- `RDKit❗❗`: 3
- `RDKit✔️❌`: 180
- `RDKit❗✔️`: 33

The frozen valence/sanitize scope therefore remains incomplete.

## Representative Blocking Examples

Examples observed during this audit that prevent a zero-gap claim:

- `notation/smiles.rs`: unresolved copied-source blocks remain throughout the
  frozen parser path, including lexer/source-port sections and broader parse/CX
  orchestration.
- `notation/smiles_write.rs`: unresolved copied-source blocks remain in top-level
  writer orchestration such as `MolToRandomSmilesVect` / `MolToSmiles`-aligned
  logic.
- `notation/canon_rank.rs`: frozen ranking helpers still carry `RDKit✔️❌` and
  `RDKit❗✔️` markers, so the block is not performance/behavior closed.
- `chemistry/rings.rs`: `find_ring_families` still includes an explicit
  `RDKit❌❌` non-URF branch, and the aggregate `findSSSR` /
  `symmetrizeSSSR` surfaces still carry `RDKit❗✔️` markers.
- `operations/ops.rs`: frozen sanitize orchestration still contains unresolved
  cleanup/property/kekulize-related blocks.

## Conclusion

The frozen baseline remains useful as the Step 177 final-audit scope, but the
current implementation state is **not** zero-gap and must **not** be described
as `100% coverage`.
