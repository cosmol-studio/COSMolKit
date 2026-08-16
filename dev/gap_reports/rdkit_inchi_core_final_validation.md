# RDKit InChI Core Final Validation

## Scope

This report closes the narrowed source-port plan for these four scalar COSMolKit APIs:

- `molecule.to_inchi()`;
- `molecule.to_inchi_key()`;
- `inchi_to_key(inchi)`;
- `Molecule.from_inchi(inchi)`.

The official target is GCC/Linux `libinchi` with `COMPILE_ANSI_ONLY`,
`TARGET_API_LIB`, and `fPIC`. The RDKit adapter oracle is pinned to RDKit
`2026.03.1`; the official engine source is InChI `1.07.5`. Exact-result claims
below apply only to source-defined behavior for which the complete recorded
fields and active branches compare exactly. They do not include the authorized
official-C undefined path.

## Source Frames

- The configured official target contains `1364` active definitions.
- The five selected official roots are `GetINCHI`, `FreeINCHI`,
  `GetINCHIKeyFromINCHI`, `GetStructFromINCHI`, and `FreeStructFromINCHI`.
- Compiler direct-edge and address-taken callback evidence produces a selected
  callback-complete closure of `1050` active definitions: `1014` direct-closure
  definitions plus `36` callback targets.
- All `1050` selected official source locations have complete in-function
  source frames and complete Rust implementations. `1049` locations have
  closed behavioral first axes. `NormalizeAndCompare` is the one selected
  location with an open first axis.
- All `30` functions in the selected RDKit adapter closure have complete
  source frames, complete Rust implementations, and closed behavioral first
  axes.
- There are no remaining selected source implementation Port units.

The compiler derivation, root mapping, per-function source locations, and
frame checks are recorded in
`dev/gap_reports/rdkit_inchi_core_scope_audit.md`. Historical completed-prefix
evidence remains in
`dev/gap_reports/rdkit_inchi_completed_prefix_parity_audit.md`.

## Focused Rust Tests

Every selected official and RDKit adapter Port unit has a function-specific
focused Rust test that executed in the main test harness; a provenance harness
with zero selected tests was never counted as a focused pass. The per-function
commands and results are retained in the two audit reports above.

The remaining open-marker function was checked with:

```bash
cargo test -p cosmolkit-inchi \
  source_port__ichirvr4__normalizeandcompare__line_1274 -- --nocapture
```

The main harness executed one matching test and passed. The separately
authorized Rust behavior was checked with:

```bash
cargo test -p cosmolkit-inchi \
  source_policy__normalizeandcompare__undefined_initial_buffer_allocation_returns_structured_error \
  -- --nocapture
```

The main harness executed exactly one matching test and reported `1 passed; 0
failed; 1265 filtered out`. It verifies `SourceHeapError::AllocationFailed`,
caller-visible mutation order, allocation and ownership state, cleanup state,
and absence of fallback.

Integration-focused suites also passed for `inchi_core_scalar_api`, the strict
`inchi_core_bridge`, `inchi_rust_facade`, and `python/tests/test_inchi.py`.
These cover value semantics, graph/index and coordinate invariants,
atom/bond/stereo/isotope/charge mapping, options, diagnostics, allocation
errors, and structured unsupported-state rejection.

## Official C Oracle

The official-C behavior oracles are non-production test dependencies. The
production implementation does not call C through FFI or an external command.

The selected-root aggregate command was:

```bash
cargo test -p cosmolkit-inchi --release \
  'source::api::inchi_dll::tests::official_c_oracle__rdkit_core_roots__exact' \
  -- --exact --test-threads=1
```

One matching test executed and passed. Its `15` official-C records compare the
selected root return statuses, diagnostics, atom/bond/stereo/isotope-H/charge/
radical/coordinate fields, input preservation, and cleanup/reset behavior.

The dedicated `NormalizeAndCompare` command was:

```bash
cargo test -p cosmolkit-inchi \
  official_c_oracle__normalizeandcompare__exact -- --nocapture
```

One matching test executed and passed. It processed `241` provenance-checked
records. All `240` records with source-defined return, mutation, allocation,
cleanup, integer, comparison, and call-order behavior matched Rust
field-for-field with zero mismatches. The remaining record is undefined-path
evidence and has no asserted exact C result.

The function-specific official-C suites recorded in the audits passed their
complete comparison schemas. Therefore, exact official-C results are claimed
for the compared source-defined fields and active branches of the selected
closure, and nowhere beyond that boundary.

## Pinned RDKit Oracle

All `30` selected RDKit adapter functions passed their focused tests and pinned
RDKit `2026.03.1` C++ oracle tests. The four public-entry record sets contain:

| Entry | Complete records | Result |
|---|---:|---|
| `Molecule.to_inchi()` (`RDKit MolToInchi`) | `65` | every recorded source-defined field and branch matched exactly |
| `Molecule.to_inchi_key()` (`RDKit MolToInchiKey`) | `9` | every recorded source-defined field and branch matched exactly |
| `inchi_to_key()` (`RDKit InchiToInchiKey`) | `9` | every recorded source-defined field and branch matched exactly |
| `Molecule.from_inchi()` (`RDKit InchiToMol` / `MolFromInchi`) | `33` | every recorded source-defined field and branch matched exactly |

The schemas include return and exception paths, diagnostic ordering, graph
mutation and preservation, atom/bond/stereo/isotope/H/charge/radical fields,
coordinates, options, cleanup, toolkit call order, and source-molecule
preservation as applicable. The oracle runners are tests only and are absent
from production dependencies.

## Strict 5000 Public API Regression

The Schematic issue 11 reproduction initially found eight valid strict-corpus
rows that failed COSMolKit InChI generation with
`SourceHeapError::PointerOutOfBounds`: `440`, `904`, `2040`, `2556`, `3248`,
`3620`, `3811`, and `4944`. The source-level repair preserves the official
`ichi_bns.c:10951-10957` short-circuit order in `BalancedNetworkSearch`; the
`bIgnoreVertexNonTACN_group` call is no longer evaluated on states for which
the preceding official predicates are false.

The exact focused regression command executed one test:

```bash
cargo test -q -p cosmolkit-core \
  --features op-contracts-strict \
  --lib \
  inchi::tests::inchi_generation_matches_chematic_issue_11_pointer_regressions \
  -- --exact --nocapture
```

It reported `1 passed; 0 failed; 2573 filtered out` and checks the complete
InChI, warning return code, and InChIKey for every formerly failing row.

The installed Python extension was then rebuilt from the current Rust source
with `.venv/bin/maturin develop --manifest-path python/Cargo.toml`. A
fail-closed public-wrapper harness processed all `5000` rows of
`tests/smiles_5000.smi`, independently parsed each input with pinned RDKit
`2026.03.1` and COSMolKit, and compared `rdkit.Chem.MolToInchi` with
`cosmolkit.Molecule.to_inchi()` by exact string equality. Its result was:

```text
corpus_rows=5000
exact=5000
mismatch=0
rdkit_parse_error=0
cosmolkit_parse_error=0
rdkit_generation_error=0
cosmolkit_generation_error=0
failure_count=0
```

The harness exits nonzero unless every count has exactly the values above; no
row, warning, exception, or output field is skipped.

This corpus check is now also a persistent Rust integration test at
`crates/cosmolkit-core/tests/rdkit_inchi_parity.rs`. Its pinned-RDKit golden
generator is part of the unified generator as the `inchi` suite, and committed
goldens contain `150` small-profile records and `5000` strict-profile records.
The test checks corpus/golden row counts, contiguous row numbers, source SMILES
order, and every result. The small profile has `132` nonempty RDKit outputs,
which require complete InChI byte equality; `5` RDKit empty outputs, for which
Rust must return either an empty output or the policy-defined structured
`UnsupportedState`; and `13` RDKit parse failures, which must remain Rust parse
failures. No unrelated generation-error category is accepted. Every
strict-profile row has a nonempty RDKit InChI, so all `5000` strict comparisons
are exact byte equality.

The persisted commands both executed one matching test:

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_inchi_parity \
  inchi_matches_pinned_rdkit_for_every_active_profile_row -- --exact

COSMOLKIT_PARITY_PROFILE=smiles_5000 \
  cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_inchi_parity \
  inchi_matches_pinned_rdkit_for_every_active_profile_row -- --exact
```

The small profile reported `1 passed`, `0 failed`, `0 filtered out` after
processing all `150` rows. The strict profile reported `1 passed`, `0 failed`,
`0 filtered out` after processing all `5000` rows.

## Authorized Undefined-Behavior Mapping

`NormalizeAndCompare` at
`third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1274-1672` has one
active official-C undefined path. Initial `inchi_strbuf_init` allocation
failure leaves `strbuf->pStr == NULL`; the caller reaches
`strcpy(existing_formula, strbuf->pStr)`. C defines no return value, signal,
mutation state, or cleanup sequence for that execution.

The human-authorized Rust policy maps only this audited path to deterministic
`SourceHeapError::AllocationFailed`. This is structured Rust behavior, not an
official-C exact result. It does not authorize any inferred behavior for other
undefined, unverified, or unsupported paths. `NormalizeAndCompare` remains
`INCHI❗❌`.

## Remaining Selected Call Graph

No selected source implementation Port unit remains. The selected official
closure still contains one source location with an open behavioral first axis:

```text
GetStructFromINCHI
  -> GetStructFromINCHIEx
  -> ReadWriteInChI
  -> ConvertInChI2Struct
  -> AllInchiToStructure
  -> InChI2Atom
  -> OneInChI2Atom
  -> RestoreAtomMakeBNS
  -> RunBnsRestore1
  -> NormalizeAndCompare
```

All edges above come from the configured GCC IPA call graph. The open marker is
retained solely because the active C function contains the undefined path
described above. Its `240` source-defined oracle records pass exactly, and its
authorized Rust error path passes independently. No global all-path parity
claim is made.

## Frozen Outside-Scope Code

Outside the five-root closure, `108` active definitions are already ported and
remain frozen, while `206` active definitions are excluded and unscheduled.
Existing code is not deleted or promoted by this report. The complete
function-level disposition table is in the `Frozen And Excluded Active
Definitions` section of
`dev/gap_reports/rdkit_inchi_core_scope_audit.md`.

The excluded categories include MolBlock and SDF input, V3000, IXA, AuxInfo
conversion wrappers, incremental INCHIGEN, version query, extended polymer
input, and CLI/demo functionality. None is exposed through the narrowed four-
API public boundary.

## Unsupported States

The public boundary exposes only `mol_to_inchi`, `mol_to_inchi_key`,
`inchi_to_inchi_key`, and `mol_from_inchi`, plus their Rust-facade and Python
forms. Unmodeled molecule states, including unsupported atom/bond semantics,
query state, substance-group state, and mixed unsupported state, are rejected
with structured `UnsupportedState` errors. Invalid input, invalid source
output, toolkit failure, and allocation failure retain separate structured
error categories.

There is no heuristic result, placeholder chemistry, silent fallback,
SMILES/MolBlock transit, production FFI, or production external-command path.
Frozen entry points are not silently routed through the supported scalar APIs.

## Performance Markers

Behavior and performance remain independent. This validation did not batch-
promote second-axis markers. Existing `INCHI✔❌`, `RDKit✔❌`, and other
reviewed second-axis states remain unchanged wherever allocation, copying,
buffering, lookup, or other known performance differences remain.

In particular, `NormalizeAndCompare` remains `INCHI❗❌`: the first symbol
records the official undefined path, and the second records the known
performance gap. Passing source-defined behavior tests does not change either
fact.

## Full Validation

| Command or harness | Recorded result |
|---|---|
| `cargo fmt --all -- --check` | passed |
| `cargo test -p cosmolkit-inchi` | passed; main harness reported `1271 passed`, `0 failed`, `0 filtered out`; provenance harness reported `1 passed` |
| strict 5000 Python public-wrapper exact comparison | pinned RDKit `2026.03.1`; `5000` exact, `0` mismatches, `0` parse errors, `0` generation errors |
| persisted Rust InChI small-profile parity | `150` rows processed; `1 passed`, `0 failed`, `0 filtered out` |
| persisted Rust InChI strict-profile parity | `5000` exact nonempty InChI byte comparisons; `1 passed`, `0 failed`, `0 filtered out` |
| `cargo check -p cosmolkit-core --features op-contracts-strict` | passed |
| `cargo test -p cosmolkit-core --release --features op-contracts-strict` | passed; core harness reported `2529 passed`, `0 failed`, and `45` explicitly ignored tests |
| `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict` | exit `0`; InChI main harness reported `1271 passed`, `0 failed`, `0 filtered out`; strict workspace harnesses passed |
| `.venv/bin/pytest python/tests/test_inchi.py` | all selected InChI Python tests executed and passed |
| `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html` | completed successfully |
| `.venv/bin/basedpyright python/tests python/examples` | `0 errors`, `412 warnings`, `0 notes`; exit status was nonzero because warnings are configured to affect the status |

Before the final type check, `.venv/bin/maturin develop --manifest-path
python/Cargo.toml` refreshed the installed extension and stub. The repository
and installed stubs then had the same SHA-256,
`fdf851cc6442209092dc1c2b7ce0131ef61ed46658e864472d5bf0482503d960`.
The basedpyright result satisfies the plan's zero-error criterion, but it is
not represented as a warning-free or zero-exit command.

Final state: there are zero remaining selected source implementation Port
units; `1049` selected official source locations and all `30` selected RDKit
adapter functions retain their exact source-defined comparison results; one
selected source location remains `INCHI❗❌` because official C has an
undefined path and Rust has a known performance gap. No exact-result claim is
made for the authorized undefined path or for frozen and unsupported scope.
