# SDF Query Finalization: Pending Optimizations

Status: **Partially addressed; remaining optimizations pending**

Review target: commit `b2e4766`, associated with version `0.5.0-rc.4`.

This document records concerns from the supplied review. It is not an
independent code audit, an acceptance claim, or an alternative execution plan.
Proposed changes require source-backed verification and incorporation into
`dev/plans/crate_architecture_completion_plan.md` before implementation.
The architecture, public API rules, and source-reproduction protocol remain
authoritative. No code, equality semantics, or ownership changes are authorized
by this record alone.

## Design properties to preserve

- Keep live concrete `Molecule` and detached `QueryGraph` separate.
- Preserve typed `QueryPredicateOrigin::{Explicit, CarrierDerived}` rather than
  inferring query identity from `_MolFileAtomQuery` or `_MolFileBondQuery`.
- Preserve query state across chemistry transforms through validated topology
  mappings; do not silently lose explicit predicates.
- Keep SGroup bond references and three-dimensional geometry in typed model
  fields rather than magic string properties.
- Keep chemistry algorithms in their domain owners and live-object creation
  and commit authority in the runtime.

## 1. Commit size and source-port auditability

Status: **Pending optimization**

The reviewed commit combines query provenance and transport, chemistry
consumers, legacy stereochemistry, MolBlock finalization, SGroups, CX handling,
binding contracts, regressions, adapter changes, and release version updates.
This breadth makes source-level review and regression attribution difficult.

The actual commit summary is 68 files, 15,633 insertions and 2,361 deletions;
the supplied review's approximate line counts must not be used as verified
repository statistics.

For future work, prefer independently reviewable, dependency-ordered changes
for provenance, transport, query-aware chemistry, legacy stereo, finalization,
typed SGroups, regression coverage, and release metadata. Preserve source and
validation evidence for each logical change. This recommendation does not
authorize rewriting the existing commit or fragmenting complete behaviors
into falsely completed partial implementations.

## 2. Query state is an overlay, not a current carrier snapshot

Status: **Access boundary corrected; representation redesign deferred**

The review identifies `QueryStateRef` as borrowing `QueryAtom` and `QueryBond`
rows while consumers primarily need predicates and their origins. Following
a topology transform, their embedded carriers may describe an earlier state.
Validation of counts, IDs, and bond endpoints does not establish equality of
all carrier chemistry fields with the current topology.

Audit consumers for accidental reads of stale embedded carriers. Document
which state is authoritative throughout the transform lifecycle and verify
the exact guarantees of `try_for_topology()`.

Potential designs include a more explicit overlay name or predicate/origin-only
borrowed access. `QueryOverlayRef`, `QueryAtomState`, and `QueryBondState` are
review suggestions, not approved names or requirements to introduce additional
storage types. Prefer narrowing the existing abstraction over duplicating it.

Resolution evidence should include a transformed-topology case that preserves
explicit predicates while proving consumers read current carrier chemistry
from the authoritative topology.

## 3. Provenance and public equality

Status: **Representation-equality contract documented and tested**

The review notes that derived equality includes `predicate_origin`, which can
make otherwise identical carrier/predicate rows unequal. Define whether
public equality means exact representation equality, operational-state
equality, or some explicitly bounded query-semantic equivalence.

Do not simply remove provenance from equality: `Explicit` and `CarrierDerived`
can produce different `hasQuery()` behavior and different chemistry pipeline
results. Identical predicate syntax alone does not prove semantic equivalence.
Conversely, representation equality should not be documented as general
chemical or matching equivalence.

Audit observable equality uses and document the chosen contract. Add examples
and regressions for equal carriers/predicates with different origins. Separate
comparison methods, if justified, require normal API review and registration;
the suggested names `structural_eq()` and `semantic_eq()` are not approved APIs.

## 4. Preserve typed finalization errors

Status: **Pending optimization**

The reported `MolPostError::Processing(String)` conversions flatten sanitization,
hydrogen-removal, stereochemistry, and other errors into display strings.
This may lose typed causes and make precise Rust or binding error mapping
depend on fragile message parsing.

Audit each conversion and consider typed variants with retained source chains,
while preserving existing source-defined failure behavior. Keep dependency
direction acyclic and avoid exposing runtime internals. Validate that callers
can distinguish failure categories without parsing strings and that public
error contracts remain consistent with binding registration.

## 5. Public SDF concrete/query boundary

Status: **Boundary fixed in detached code and public registration; full public reader integration pending**

Registration of `Molecule::from_sdf` and its parameterized variant does not
establish implementation or support. A query-bearing `MolBlockRecord` cannot
be returned as a concrete `Molecule` without losing semantics.

Preserve the previously discussed direction: concrete-only molecule readers
return a structured error for query-bearing results; a complete reader uses
an IO record container with an explicit concrete/query payload. Reuse existing
models instead of introducing a third chemistry model or query-bearing live
`Molecule`. Exact public names and signatures must follow the sole plan and
binding contract.

Classification must account for finalization, not just parsing: an operation
may introduce query semantics after a concrete record has been parsed. Tests
must cover ordinary concrete input, encoded query input, and query introduced
during finalization, preserving record metadata in every applicable path.

## 6. Attachment expansion and appended query rows

Status: **Pending implementation under the sole plan**

The supplied review reports an explicit attachment-expansion error and
`remap_query_rows()` rejection of appended rows. Verify the precise current
support boundary; do not infer that enabling an option always changes graph
category when there are no applicable attachment points.

Full source-backed expansion needs coherent atom/bond additions, coordinates,
query predicates and origins, topology mappings, and associated metadata.
Existing-row remapping alone cannot define the semantics of appended rows.
Preserve failure atomicity and explicit errors until the complete behavior is
implemented; neither silent dropping nor a heuristic partial expansion is
acceptable completion. An unsupported boundary is not evidence that the
required capability is finished.

Validate concrete-to-query promotion, already-query records, additions and
index alignment, coordinate handling, option branches, and no-op cases against
the pinned source/reference before claiming complete support.

## 7. Release version scope

Status: **Intent confirmed; documentation clarity pending**

The rc.4 scope intentionally covered the 21 previously published current Rust
crates. Python, WASM, historical macros, test support, and tooling were excluded
from package-version upgrades and publication. Internal dependency requirements
were synchronized where needed.

Keeping `cosmolkit-py` at rc.2 and WASM at rc.3, with WASM explicitly pinned
instead of inheriting the workspace version, was deliberate for this release.
Do not automatically synchronize unpublished adapter versions as a cleanup.
Future release documentation should distinguish package version, dependency
version, publication scope, and adapter implementation status.

## Acceptance boundary

Before claiming complete public SDF support, resolve the provenance/equality
contract, preserve typed error semantics, complete the public concrete/query
boundary, and implement and validate required attachment expansion. Audit the
overlay's carrier-access risks alongside those changes.

This record makes no percentage-complete claim and does not certify the
review's statement that the internal architecture is already fully validated.
Track executable work and completion evidence only in the sole migration plan.

## Implemented boundary correction (2026-09-21)

The subsequent request authorized direct correction of the public SDF boundary,
overlay carrier access, and equality semantics. It did not authorize completion
of attachment expansion or a new overlay representation.

### Overlay consumer audit and correction

| Consumer | Actual use and correction |
|---|---|
| core aromaticity | Predicate inspection and alignment validation; use `validate_for_topology` instead of exposing rows. |
| core kekulize | Query identity/predicate inspection and alignment validation; same narrow validation change. |
| core CIP ranks | Query identity and alignment validation; same narrow validation change. |
| core legacy stereo | Alignment validation; current topology remains the carrier authority. |
| core sanitize | Alignment validation; no carrier chemistry obtained from query rows. |
| core hydrogen transforms | Identity/predicate tests and validation; RemoveHs previously cloned whole overlay rows for transport. It now materializes transport rows through the existing identity mapping and `remap_query_rows`, using current topology carriers. |
| IO mol-post | Constructs aligned overlays and remaps them with current topology; no raw overlay accessor was required. |

The two raw `QueryStateRef::atoms/bonds` accessors were removed. The existing
type and storage remain; there is no new overlay class or parallel graph.
Compile-fail doctests prevent restoring the raw accessor call paths unnoticed.
Model regression coverage uses changed atom charge and bond order to verify
current-carrier remapping without altering explicit predicates or provenance.
A core RemoveHs regression contrasts stale ordinary-H carrier state with
current isotope-labelled-H topology and checks default versus remove-isotope
behavior. Existing consumers still use the source-backed chemistry owners.

### Equality decision

Derived equality is deliberately retained, including provenance. Public rustdoc
on QueryAtom, QueryBond and QueryGraph now distinguishes stored representation
from matching equivalence. Tests compare equal carriers and predicate trees
with different origins, verify inequality and clone equality, and check the
graph-level consequence. No semantic comparison API was introduced.

### SDF result boundary

The existing detached `SdfGraphRecord` now owns one `into_concrete` conversion
that rejects a Query payload with `SdfReadError::QueryRecord`. Single-record
and multi-record concrete SDF readers reuse it; the latter preserves record
index and offsets around the typed error. Query-preserving readers still
return the original query payload and data fields. This is a CK type boundary,
not an upstream unsupported chemistry claim.

The binding registry now reserves `SdfRecord`, `SdfGraph` and the two record
constructors alongside the concrete Molecule constructors. All remain
registered/unsupported until full public integration is accepted. Comments
freeze classification after finalization, metadata preservation, and exclusive
runtime authority to construct Molecule. Registry tests protect the distinct
return types and honest implementation status.

The sole plan's downstream SDF implementation and public test steps were
updated accordingly without renumbering or changing existing completion marks.
This correction does not expose a new public parser, implement attachment
expansion, or certify full SDF support. MolPostError source-chain cleanup and
all unrelated review items remain deferred.

### Validation evidence

All exits below are actual process return codes, not inferred from log text.
Counts summarize executed test results; zero-test targets are not evidence of
behavioral coverage. Temporary logs are diagnostic conveniences, not required
inputs for reproducing these commands.

| Command | Exit | Result |
|---|---|---|
| `cargo test -p cosmolkit-model --release --test migration_mod_query` | 0 | 16 passed, 0 failed/ignored. |
| `cargo test -p cosmolkit-model --release --doc` | 0 | 2 compile-fail tests passed. |
| `cargo test -p cosmolkit-core --release --features op-contracts-strict --test migration_h_remove_candidates query_overlay` | 0 | 1 passed, 15 filtered out. |
| `cargo test -p cosmolkit-io --release --features cosmolkit-core/op-contracts-strict --test migration_io_sdf_read sdf_concrete_boundary` | 0 | 1 passed, 8 filtered out. |
| `cargo test -p cosmolkit --release --features io,op-contracts-strict --test migration_run_schema` | 0 | 7 passed, 0 failed/ignored. |
| `cargo check -p cosmolkit-core --features op-contracts-strict` | 0 | Compile check passed. |
| `cargo check -p cosmolkit --features io,op-contracts-strict` | 0 | Compile check passed. |
| `cargo test -p cosmolkit-model --release` | 0 | 104 passed including 2 compile-fail doctests, 0 failed/ignored. |
| `cargo test -p cosmolkit-core --release --features op-contracts-strict` | 0 | Initial run: 455 passed; final rerun including the new consumer regression: 456 passed, 0 failed/ignored. |
| `cargo test -p cosmolkit-io --release --features cosmolkit-core/op-contracts-strict` | 0 | 402 passed, 0 failed, 1 existing opt-in oracle ignored; no claim that the ignored oracle ran. |
| `cargo test -p cosmolkit --release --features io,op-contracts-strict` | 0 | 260 passed, 0 failed/ignored; includes both real-module-layout privacy tests with their default/strict probes. |
| `cargo test --workspace --release --features cosmolkit/op-contracts-strict,cosmolkit-core/op-contracts-strict` | 101 | Compilation blocked by existing Python (664 errors), WASM (98 errors), and tautomer-oracle (2 errors) API gaps; workspace tests did not pass. |
| `cargo fmt --all -- --check` | 0 | Formatting passed after Rust edits. |
| `.venv/bin/python dev/tools/architecture/check_plan.py --plan dev/plans/crate_architecture_completion_plan.md` | 0 | Plan validation passed without checker changes. |

The correction has targeted and affected-crate evidence, not full workspace
acceptance, full SDF support, or renewed upstream parity claims. No version bump,
publication, commit, or language-adapter repair was part of this correction.
