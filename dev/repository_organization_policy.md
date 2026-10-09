# Repository Organization Policy

This policy defines where tests, test data, generated reference data, and
test-data tooling belong in COSMolKit. The terms **MUST**, **MUST NOT**,
**SHOULD**, and **SHOULD NOT** indicate requirement strength.

## 1. Test Placement

Tests MUST be placed according to the boundary they verify.

### 1.1 Inline tests

Tests MAY be defined inside `src/` when they verify private functions, local
algorithms, source-port functions, or module-local invariants. Inline tests
SHOULD use small inputs defined directly in test code.

A fixed source-backed regression MAY remain inline when moving it would
require widening production visibility. Its expected values must be fixed;
ordinary tests must not invoke reference implementations or prepare corpora.

Tests that exercise only public cross-module behavior, complete shared
corpora, or a public file-format workflow MUST NOT be implemented as inline
tests.

### 1.2 Rust crate tests

Tests that exercise a crate through its public interface MUST be placed under:

```text
crates/<crate>/tests/
```

This includes fixed regressions for public API behavior, cross-module
workflows and operation contracts. Fixed upstream-derived tables are fixtures,
not corpus parity merely because they contain many rows. Corpus preparation,
reference execution and corpus comparison belong to the top-level
current `parity-tests_fixed/` crate, not to owning crates.

Facade crates MUST NOT repeat complete suites already covered by an underlying
crate. They SHOULD test only facade-specific exports and behavior.

The explicitly designated 77-case structure-tag matrix is the special
source-regression exception described in `test_boundaries.md`. Its complete
reference-dependent comparison lives in `parity-tests_fixed/tests/`; ordinary local
validation/numeric tests stay in core. Test-only detached dependencies preserve
the source boundary without exposing a new production/public algorithm API.

### 1.3 Python tests

Python tests MUST be placed under `python/tests/`. They MUST focus on behavior
that can be changed by the Python boundary, including type conversion,
exception conversion, ownership, mutation, array shape and dtype, ordering,
serialization, and Python-visible exports.

Python tests MAY reuse representative Rust cases. They MUST NOT repeat a
complete Rust parity corpus unless the Python boundary can independently
change the complete tested result.

## 2. Test Data Ownership

For COSMolKit 0.5.0, corpus and designated special-regression inputs MUST be
stored under `parity-tests_fixed/testdata/`. Their generated references belong
under `parity-tests_fixed/expected/`, and reports under `parity-tests_fixed/reports/`.
The repository-level `testdata/` directory retains shared ordinary-regression
fixtures and historical provenance, not a second corpus preparation workflow.
Existing consumers are not moved by this documentation update.

Fixed regressions MAY use `include_str!`/`include_bytes!` or a path anchored
at `CARGO_MANIFEST_DIR` to access repository-level fixtures. They MUST NOT
depend on the process working directory or maintain convenience copies in
crate-local or Python-local fixture directories. A shared support crate is
not required for reading fixed files.

Repository tests MAY read fixtures directly from pinned `third_party/`
submodules. Do not copy an upstream fixture merely to change its location.
Record each selected file's submodule revision, relative path and checksum
in its fixture documentation. Missing inputs are errors, never skips.
Anchor runtime fixture paths at `CARGO_MANIFEST_DIR`, not the working directory.
Do not use compile-time `include_str!`/`include_bytes!` for these external
fixtures: read them at test runtime so packaged test sources remain compilable.

Production libraries, build scripts, examples and documentation examples MUST
NOT require these external test files. Release-profile repository tests MAY
use them; `--release` does not mean publication. Published package builds must
work without the source repository or initialized submodules. Coverage CI
first builds the production library without initialized submodules; only a
successful build permits fetching pinned third-party inputs and running tests.
The build and tests share coverage instrumentation, release profile and target
directory; do not clean compiled artifacts between these stages. This is a
build-time dependency check, not proof about runtime filesystem access or
package contents. Package contents MUST NOT include the external fixture tree.
Running repository-only fixture tests from a downloaded crate is not a
supported substitute for the full repository test suite; document their
required checkout instead of silently treating absent data as a pass.

This permission covers input/expected data, not compiling third-party code,
parsing source code into expectations or invoking oracles in regression tests.
Corpus/reference execution remains owned by `parity-tests_fixed`.

A test MAY create temporary files in a temporary directory. Test execution
MUST treat committed inputs and generated expected data as read-only and MUST
NOT write into `testdata/`.

### 2.1 Large externally distributed audit corpora

A complete externally distributed corpus MAY remain uncommitted when its size
or distribution terms make repository storage unsuitable and it is used only
by an explicit large-stress audit, not by ordinary tests. This exception does
not permit an undocumented local corpus. The repository MUST commit:

- the upstream release and stable source URL;
- the exact source checksum and expected record count;
- deterministic, atomic preparation and selection code;
- the shard assignment algorithm and output-manifest schema;
- the complete audit profile, reference-version pins, and acceptance rules;
  and
- a documented repository-owned command that validates every prepared shard
  before execution.

The source, prepared shards, and run outputs MUST remain outside tracked
repository data. A run manifest MUST bind the external corpus identity to the
Git state, installed implementation, reference environment, audit code, and
result checksums. The current ChEMBL 37 implementation of this exception is
[`tools/chembl_parity/`](tools/chembl_parity/README.md), relative to this
`dev/` directory.

## 3. Test Data Layout

Corpus inputs, generated references and reports MUST remain distinguishable:

```text
parity-tests_fixed/
  testdata/
  expected/
  reports/
  README.md

testdata/
  <format-or-domain>/fixtures/
```

Not every domain needs every directory. The current runner owns corpus
selection, expected-data naming and preparation. Domain fixture READMEs keep
provenance and purpose only; they MUST NOT define separate preparation flows.
Preserve historical fixture paths while their consumers still need them.

Examples:

```text
parity-tests_fixed/testdata/
parity-tests_fixed/expected/corpus/smiles_5000/num_heavy_atoms_smiles/reference.jsonl
testdata/mol2/fixtures/
testdata/molblock/fixtures/
```

Directory names MUST use lower snake case unless an upstream filename must be
preserved. Permanent paths MUST NOT use migration-stage names such as `new`,
`old`, `phase_1`, `final_fix`, or `agent_output`.

## 4. Test Data Terms

### 4.1 Fixture

A fixture is a committed input file used by one or more tests. Fixtures MUST
NOT be generated during normal test execution. Externally derived fixtures
MUST record their source project, source path, source version or commit,
selection method, license context, and file checksum in a nearby README or
manifest.

### 4.2 Corpus

A corpus is a named collection of inputs processed as a suite. Its provenance,
selection, filtering, ordering, and checksum MUST be documented. Corpus cases
MUST NOT be silently removed to make parity pass.

### 4.3 Expected data

Expected data are outputs used for comparison. They MAY be generated locally
or in CI and MAY remain uncommitted and ignored by Git.

The repository MUST commit every fixture, corpus, generator, schema,
generation option, and reference-version pin needed to reproduce expected
data. Generated expected data MUST include a machine-readable identity
manifest and output checksums.

Each expected domain/profile directory MUST contain `manifest.json`. The
manifest MAY list multiple outputs from that domain, but every output entry
MUST carry the generator, input, option, schema, record-count, and checksum
identity needed to validate that output. Preparation of a narrower suite MAY
publish a manifest containing only the outputs selected by that suite; tests
for other outputs must then fail with the preparation command for their
domain.

A cached expected-data family MAY be reused only when every identity field
matches exactly. A missing, incomplete, corrupt, or stale family MUST be
regenerated before its tests run. Required tests MUST fail explicitly when
valid expected data cannot be prepared; they MUST NOT skip, use stale data, or
fall back to weaker assertions.

### 4.4 Cache identity

The identity manifest MUST include at least:

```text
reference implementation name and version or commit
generator source checksum
input corpus checksums
fixture checksums when used
output schema version
generation profile and options
expected record count
generated output checksums
platform identity when output is platform-dependent
```

Cache keys are an optimization. A cache hit MUST still be validated against
the identity manifest. Directory existence alone is not validation.

CI cache keys MUST include an explicit cache-schema version in addition to
the generated-data identity. Because GitHub Actions caches are immutable, a
workflow MUST NOT depend on overwriting an invalid exact-key cache. The
supported repair pattern is a stable identity restore prefix plus a unique
run/attempt suffix for saves. Preparation validates a restored candidate; if
it regenerates any domain, the workflow saves the validated replacement under
the current unique key. The next restore selects the newest matching valid
candidate. Increment the cache-schema version when the cache layout or restore
protocol changes, not merely when generator inputs change.

### 4.5 External reference implementation

RDKit, Gemmi, official InChI, and similar implementations MAY generate or
verify expected data. They MUST NOT become production runtime dependencies.
Reference dependencies MAY be required by the top-level parity pipeline,
but ordinary regression test binaries MUST NOT invoke generators or oracles.

### 4.6 Known failures

Known failures MUST be stored separately from test logic and MUST remain
executable. Tests MUST NOT hide failures through filtering, broad exception
handling, loop-local skipping, reduced comparison schemas, or silent case
removal.

## 5. Expected Data Preparation

Preparation and Cargo comparison are separate stages:

1. Select all registered tasks by default, or the explicitly requested tasks.
2. Validate the complete selected input/variant set.
3. Reuse only reference generations whose full identity and checksums match.
4. Call the selected Python generators sequentially with corpus, task parameters
   and concurrency; each function parallelizes cases internally. Rust labels,
   validates and atomically publishes each generation. Preserve invalid evidence.
5. Save the prepared function/corpus selection and preflight ALL its references.
6. Run independent Cargo tests for function/corpus-type pairs; Cargo parallelizes
   tests, which use checked snapshots and report every selected result.

A failure in preparation or global preflight prevents all selected operation
calls. A comparison mismatch must never regenerate references to fit CK.
Preparation must finish before corpus tests. Preflight is read-only and never
starts Python; Cargo tests also fail rather than generate missing data. Ordinary
regressions remain independent of references. Tests never modify expectations.
Corpus type is explicit in registration and source selection, not inferred from
file suffixes. The same function on SMILES and SDF is two separate tests.

Reference adapters and existing preparation scripts are implementation tools,
not independent task registries. Keep task/variant selection and result schemas
in Rust. Generated corpora, caches and reports stay in ignored output paths;
committed fixture snapshots remain read-only. The current runner's actual
coverage and limits are documented in [parity-tests_fixed/README.md](../parity-tests_fixed/README.md).

Validate generated output schemas, input identity, reference version,
generator identity, exact case counts and checksums before comparison.
Load validated snapshots so subsequent file changes cannot alter a run's
expected results. Missing data must never first be discovered halfway through
chemistry execution. Do not duplicate these mechanisms in domain crates.

## 6. Rust Parity Pipeline Ownership

`parity-tests_fixed/` (`cosmolkit-parity-tests-fixed`, unpublished) is the
current 0.5.0 home for corpus and special-regression orchestration.
The current runner's chemistry dependency is public
`cosmolkit` with `full`, not individual domain crates. Its Rust registry
declares operations, variants, typed inputs/results and reference identity.
Do not create a parallel test-support crate, registry or owner-local corpus
runner. Fixed regressions remain independent of the parity runner.

Special source regressions explicitly designated in `test_boundaries.md` reuse
the same prepare CLI and the `special_regression` Cargo test target.
They are not corpus-task registrations or corpus passes. This exception does
not move ordinary fixed owner regressions or permit a second test framework.

Python/JS binding verification uses representative cases for FFI/WASM and
language-boundary behavior. It does not repeat complete chemistry corpora
without a specific boundary need. Complete 0.5.0 validation is pending;
historical results do not establish current coverage.

## 7. Cross-Layer Coverage

The same fixture, corpus, or expected output MAY be consumed by multiple
layers, but the same behavior SHOULD NOT be exhaustively retested at every
layer.

- The top-level Rust parity runner owns complete chemistry parity corpora.
- Python tests SHOULD use representative cases for binding-specific behavior.
- Facade crates SHOULD verify exports without repeating core suites.
- Private source-port tests SHOULD verify private branch and field behavior at
  the smallest stable boundary.

## 8. Naming And Changes

Test files and test names MUST describe verified behavior. Permanent names
MUST NOT be based only on migration history.

Repository-wide data moves SHOULD be separated from behavior changes. A change
that adds externally derived inputs MUST also add provenance, the reference
version, selection rules, checksums, and the supported preparation command.

An agent MUST NOT invent another local fixture or expected-data convention. It
must use the closest existing domain or update this policy first.
