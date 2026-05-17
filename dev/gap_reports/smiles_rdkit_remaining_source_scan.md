# SMILES RDKit Remaining Source Scan

Date: 2026-05-16

## Result

Frozen-scope coverage is **not** at `100%`.

This Step 219 audit rescans the current tree against the frozen baseline named
by `dev/smiles_rdkit_full_port_checklist.md`. A zero-gap rewrite is allowed
only when every frozen-baseline function is marker-closed and no direct frozen
helper/function gap remains. That condition is still false.

## Frozen Baseline Used For Audit

The frozen baseline for this checklist pass remains:

- parser/core helpers in `crates/cosmolkit-core/src/notation/smiles.rs`
- writer helpers in `crates/cosmolkit-core/src/notation/smiles_write.rs`
- canon-ranking helpers in `crates/cosmolkit-core/src/notation/canon_rank.rs`
- ring helpers in `crates/cosmolkit-core/src/chemistry/rings.rs`
- kekulization helpers in `crates/cosmolkit-core/src/chemistry/kekulize.rs`
- valence/sanitize helpers in `crates/cosmolkit-core/src/chemistry/valence.rs`
  and `crates/cosmolkit-core/src/operations/ops.rs`

The RDKit comparison set remains:

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

## Zero-Gap Standard

For this frozen scope, `100% coverage` requires all of the following:

- every frozen-baseline function is behaviorally closed for the currently
  modeled input state space,
- no copied-source line inside the frozen scope remains at
  `RDKit❌❌`, `RDKit❗❗`, `RDKit❗✔️`, `RDKit✔️❗`, or `RDKit✔️❌`,
- the checklist-required targeted tests have passed,
- a direct source rescan finds no remaining frozen helper/function gap.

Any remaining marker in those classes blocks a zero-gap claim.

## Current Frozen-Scope Status

Full frozen-scope marker scan command:

```text
rg -n "RDKit(❌❌|❗❗|❗✔️|✔️❗|✔️❌):" \
  crates/cosmolkit-core/src/notation/smiles.rs \
  crates/cosmolkit-core/src/notation/smiles_write.rs \
  crates/cosmolkit-core/src/notation/canon_rank.rs \
  crates/cosmolkit-core/src/chemistry/rings.rs \
  crates/cosmolkit-core/src/chemistry/kekulize.rs \
  crates/cosmolkit-core/src/chemistry/valence.rs \
  crates/cosmolkit-core/src/operations/ops.rs
```

Current unresolved marker counts by frozen-scope file:

- `crates/cosmolkit-core/src/notation/smiles.rs`: `0`
- `crates/cosmolkit-core/src/notation/smiles_write.rs`: `0`
- `crates/cosmolkit-core/src/notation/canon_rank.rs`: `0`
- `crates/cosmolkit-core/src/chemistry/rings.rs`: `0`
- `crates/cosmolkit-core/src/chemistry/kekulize.rs`: `292`
- `crates/cosmolkit-core/src/chemistry/valence.rs`: `4`
- `crates/cosmolkit-core/src/operations/ops.rs`: `214`

These frozen-scope files are now marker-closed and are no longer Step 219
blockers:

- `crates/cosmolkit-core/src/notation/smiles.rs`
- `crates/cosmolkit-core/src/notation/smiles_write.rs`
- `crates/cosmolkit-core/src/notation/canon_rank.rs`
- `crates/cosmolkit-core/src/chemistry/rings.rs`

The frozen scope is still blocked by unresolved markers in the kekulize,
valence, and sanitize/orchestration paths.

## Exact Remaining Blockers

### 1. Parser path has zero remaining frozen-scope markers

File: `crates/cosmolkit-core/src/notation/smiles.rs`

Parser-only marker scan command:

```text
rg -n "RDKit(❌❌|❗❗|❗✔️|✔️❗|✔️❌):" crates/cosmolkit-core/src/notation/smiles.rs
```

Result:

```text
no matches
```

Conclusion:

- there are zero remaining parser-scope non-closed RDKit markers in
  `crates/cosmolkit-core/src/notation/smiles.rs`
- no parser continuation steps are required at this time

### 2. Kekulize path is still marker-gap blocked

File: `crates/cosmolkit-core/src/chemistry/kekulize.rs`

Kekulize-only marker scan command:

```text
rg -n "RDKit(❌❌|❗❗|❗✔️|✔️❗|✔️❌):" crates/cosmolkit-core/src/chemistry/kekulize.rs
```

Aggregate unresolved marker counts in the current tree:

- `RDKit✔️❌`: 292

Remaining blocking copied-source blocks in the current tree:

- `kekulize_fragment_assignment`
  - ring-info/wedged-ring ordering and filtered ring-copy block
  - aromatic-flag clearing block
  - pyrrolic N/P explicit-H reset block
- `kekulize_worker_matching`
  - full worker loop block remains `RDKit✔️❌`
- `permute_dummies_and_match`
  - `QuestionEnumerator` copied block remains `RDKit✔️❌`
  - `permuteDummiesAndKekulize` copied block remains `RDKit✔️❌`

Closed since the previous scan:

- `kekulize_assignment`
- `kekulize_fused_assignment`
- `mark_double_bond_candidates`
- `kekulize_if_possible_assignment`

Conclusion:

- `crates/cosmolkit-core/src/chemistry/kekulize.rs` is still a frozen-scope
  blocker
- the remaining blockers are now concentrated in fragment orchestration,
  worker matching, and dummy-permutation handling

### 3. Valence path still has an explicit performance gap

File: `crates/cosmolkit-core/src/chemistry/valence.rs`

Valence-only marker scan command:

```text
rg -n "RDKit(❌❌|❗❗|❗✔️|✔️❗|✔️❌):" crates/cosmolkit-core/src/chemistry/valence.rs
```

Aggregate unresolved marker counts in the current tree:

- `RDKit✔️❌`: 4

Remaining blocking function:

- `ValenceContext::new`

The current marker state is deliberate: when the adjacency cache is absent,
COSMolKit builds an ephemeral adjacency snapshot instead of reading graph
adjacency directly from persistent molecule storage as RDKit does.

### 4. Sanitize/property/cleanup orchestration still has explicit performance gaps

File: `crates/cosmolkit-core/src/operations/ops.rs`

Sanitize/orchestration marker scan command:

```text
rg -n "RDKit(❌❌|❗❗|❗✔️|✔️❗|✔️❌):" crates/cosmolkit-core/src/operations/ops.rs
```

Aggregate unresolved marker counts in the current tree:

- `RDKit✔️❌`: 214

Remaining blocking functions:

- `add_hs_terminal_coords`
- `sanitize_after_remove_hs_removal`
- `run_sanitize_pipeline`
- `sanitize_nitrogens_cleanup_assignment`
- `sanitize_phosphorus_cleanup_assignment`
- `sanitize_halogen_cleanup_assignment`
- `sanitize_cleanup_incident_bonds`
- `sanitize_organometallic_cleanup_assignment`
- `sanitize_metal_bond_cleanup_assignment`
- `sanitize_is_hypervalent_nonmetal`
- `sanitize_organometallic_single_bonded_metals`

These are no longer behaviorally unresolved in the audited subset, but the
copied-source comments still document known performance/representation gaps
relative to RDKit, so Step 219 cannot claim frozen-scope closure.

## Conclusion

The current tree is materially ahead of the previous report: the frozen parser,
writer, canon-ranking, and ring-perception files are now marker-closed. Even
so, the overall frozen baseline is still **not** zero-gap because unresolved
markers remain in `chemistry/kekulize.rs`, `chemistry/valence.rs`, and
`operations/ops.rs`.

Do not describe the frozen scope as `100% coverage` until those remaining
function-level blockers are closed and the final strict validation steps pass.

## Strict Validation Snapshot

Formatting command:

```text
cargo fmt --all
```

Strict check command:

```text
cargo check -p cosmolkit-core --features op-contracts-strict
```

Result:

- passed
- warnings remain elsewhere in the crate, but no strict-check failure blocks
  this checklist pass

Checklist Step 375 used split equivalent targeted commands because Cargo accepts
only one test-name filter per invocation:

```text
cargo test -p cosmolkit-core --features op-contracts-strict notation::smiles::tests
cargo test -p cosmolkit-core --features op-contracts-strict notation::smiles_write::tests
cargo test -p cosmolkit-core --features op-contracts-strict canon_rank
cargo test -p cosmolkit-core --features op-contracts-strict rings
cargo test -p cosmolkit-core --features op-contracts-strict kekulize
cargo test -p cosmolkit-core --features op-contracts-strict valence
cargo test -p cosmolkit-core --features op-contracts-strict sanitize
```

Result:

- all targeted strict test groups passed

Full strict test command:

```text
cargo test -p cosmolkit-core --features op-contracts-strict
```

Result:

- passed: `823` unit tests
- passed: `1` doctest
- failed: `0`

Final Step 380 status for this checkpoint:

- full strict validation passes
- frozen-scope coverage is still **not** `100%`
- exact remaining blockers are the unresolved marker groups already listed in:
  - `crates/cosmolkit-core/src/chemistry/kekulize.rs`
  - `crates/cosmolkit-core/src/chemistry/valence.rs`
  - `crates/cosmolkit-core/src/operations/ops.rs`
- the checklist stays open until those blocker groups are closed and the full
  frozen-scope scan returns zero unresolved markers

## Step 456 Final Frozen-Scope Rescan Update

Date: 2026-05-16

Full frozen-scope marker scan command:

```text
rg -n "RDKit(❌❌|❗❗|❗✔️|✔️❗|✔️❌):" \
  crates/cosmolkit-core/src/notation/smiles.rs \
  crates/cosmolkit-core/src/notation/smiles_write.rs \
  crates/cosmolkit-core/src/notation/canon_rank.rs \
  crates/cosmolkit-core/src/chemistry/rings.rs \
  crates/cosmolkit-core/src/chemistry/kekulize.rs \
  crates/cosmolkit-core/src/chemistry/valence.rs \
  crates/cosmolkit-core/src/operations/ops.rs
```

Result:

- no matches

Final frozen-scope unresolved marker counts by file:

- `crates/cosmolkit-core/src/notation/smiles.rs`: `0`
- `crates/cosmolkit-core/src/notation/smiles_write.rs`: `0`
- `crates/cosmolkit-core/src/notation/canon_rank.rs`: `0`
- `crates/cosmolkit-core/src/chemistry/rings.rs`: `0`
- `crates/cosmolkit-core/src/chemistry/kekulize.rs`: `0`
- `crates/cosmolkit-core/src/chemistry/valence.rs`: `0`
- `crates/cosmolkit-core/src/operations/ops.rs`: `0`

The final two blockers closed in this continuation were:

- `ValenceContext::new`
- `add_hs_terminal_coords`

Strict validation commands:

```text
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict
```

Strict validation result:

- `cargo check` passed
- `cargo test` passed: `839` unit tests, `1` doctest, `0` failures

Final conclusion:

- frozen-scope coverage is exactly `100%`
- the frozen-scope checklist has no remaining unresolved RDKit markers
- no further continuation steps are required from the Step 448 rescan
