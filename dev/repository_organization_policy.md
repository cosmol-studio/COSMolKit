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
`parity-tests/` crate, not to owning crates.

Facade crates MUST NOT repeat complete suites already covered by an underlying
crate. They SHOULD test only facade-specific exports and behavior.

### 1.3 Python tests

Python tests MUST be placed under `python/tests/`. They MUST focus on behavior
that can be changed by the Python boundary, including type conversion,
exception conversion, ownership, mutation, array shape and dtype, ordering,
serialization, and Python-visible exports.

Python tests MAY reuse representative Rust cases. They MUST NOT repeat a
complete Rust parity corpus unless the Python boundary can independently
change the complete tested result.

## 2. Test Data Ownership

Committed reusable test inputs MUST be stored under the repository-level
`testdata/` directory. This includes fixtures, corpora, known-failure
declarations, schemas, source manifests, and generation configuration.

Fixed regressions MAY use `include_str!`/`include_bytes!` or a path anchored
at `CARGO_MANIFEST_DIR` to access repository-level fixtures. They MUST NOT
depend on the process working directory or maintain convenience copies in
crate-local or Python-local fixture directories. A shared support crate is
not required for reading fixed files.

`third_party/` and submodule working trees MUST NOT be used as test-data
locations. Tests MUST remain runnable when optional submodules and external
source checkouts are absent. An upstream fixture needed by tests MUST have a
committed test copy under `testdata/` with provenance and checksum metadata.

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

`testdata/` MUST group data by stable format or domain. Inputs, corpora,
expected outputs, and metadata MUST remain distinguishable.

```text
testdata/
  <format-or-domain>/
    fixtures/
    corpus/
    expected/
      <reference-implementation>/
    README.md
```

Not every domain needs every directory. Shared SMILES corpora belong under
`testdata/smiles/corpus/`; expected results derived from those corpora belong
under the behavior domain being verified, not beside the corpus.

Examples:

```text
testdata/smiles/corpus/smiles_small.smi
testdata/smiles/corpus/smiles_5000.smi
testdata/inchi/expected/rdkit/smiles_small/inchi.jsonl
testdata/mol2/fixtures/
testdata/molblock/expected/rdkit/
testdata/bio/expected/gemmi/
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

The top-level Rust parity `run` command owns preparation before comparison:

1. Select all registered tasks by default, or the explicitly requested tasks.
2. Validate the complete selected input/variant set.
3. Reuse only reference generations whose full identity and checksums match.
4. Generate missing or invalid references with the pinned oracle, publishing
   each validated generation atomically. Preserve invalid evidence.
5. Preflight ALL selected references before the first COSMolKit operation.
6. Execute Rust comparisons and report every selected result.

A failure in preparation or global preflight prevents all selected operation
calls. A comparison mismatch must never regenerate references to fit CK.
Standalone preparation/preflight commands are optional diagnostics, not
prerequisites users must manually sequence. Ordinary `cargo test` regressions
never generate or modify committed fixtures or reference results.

Reference adapters and existing preparation scripts are implementation tools,
not independent task registries. Keep task/variant selection and result schemas
in Rust. Generated corpora, caches and reports stay in ignored output paths;
committed fixture snapshots remain read-only. The current runner's actual
coverage and limits are documented in [parity-tests/README.md](../parity-tests/README.md).

Validate generated output schemas, input identity, reference version,
generator identity, exact case counts and checksums before comparison.
Load validated snapshots so subsequent file changes cannot alter a run's
expected results. Missing data must never first be discovered halfway through
chemistry execution. Do not duplicate these mechanisms in domain crates.

## 6. Rust Parity Pipeline Ownership

`parity-tests/` (`cosmolkit-parity-tests`, unpublished) is the sole home for
new corpus parity orchestration. Its chemistry dependency is public
`cosmolkit` with `full`, not individual domain crates. Its Rust registry
declares operations, variants, typed inputs/results and reference identity.
Do not create a parallel test-support crate, registry or owner-local corpus
runner. Fixed regressions remain independent of the parity runner.

Python/JS binding verification may later consume the complete Rust-validated
case set or an explicitly selected subset, with extra coverage for FFI/WASM
and language-boundary differences. This is a design requirement, not a claim
that these projections or large-corpus execution are already implemented.

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
