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

A focused source-level parity test MAY remain inline when moving it would
require widening production visibility. This includes official C/C++ oracle
tests for private ported functions. Shared file inputs and generated expected
data used by such tests MUST still be resolved through the shared test-support
API.

Tests that exercise only public cross-module behavior, complete shared
corpora, or a public file-format workflow MUST NOT be implemented as inline
tests.

### 1.2 Rust crate tests

Tests that exercise a crate through its public interface MUST be placed under:

```text
crates/<crate>/tests/
```

This includes public API behavior, cross-module workflows, operation
contracts, complete shared corpora, and public comparison with external
reference outputs. A test belongs to the crate whose public behavior it
verifies.

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

Crates and bindings MUST reference these files through shared test support and
MUST NOT maintain convenience copies in crate-local or Python-local fixture
directories.

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
Reference dependencies MAY be required by an explicit preparation command or
oracle suite, but ordinary test binaries MUST NOT silently invoke generators.

### 4.6 Known failures

Known failures MUST be stored separately from test logic and MUST remain
executable. Tests MUST NOT hide failures through filtering, broad exception
handling, loop-local skipping, reduced comparison schemas, or silent case
removal.

## 5. Expected Data Preparation

Normal test execution MUST NOT generate, replace, or repair expected data.
Preparation is an explicit step that runs before tests locally or in CI.

Each reference implementation for which the repository defines generated
expected data MUST expose one documented public preparation entrypoint. The
currently defined generated family is:

```text
tools/testdata/rdkit/generate_all.py
```

If generated Gemmi or official InChI expected-data families are added, their
public entrypoints MUST be `tools/testdata/gemmi/generate_all.py` and
`tools/testdata/official_inchi/generate_all.py`, respectively. A live oracle
alone does not require a generated-data entrypoint.

Feature-specific helpers MAY exist but MUST be private implementation details,
for example `_generate_inchi.py`. Documentation and agent instructions MUST
point to the family entrypoint.

Preparation MUST:

1. compute the complete expected cache identity;
2. validate an existing manifest and every output checksum;
3. reuse the family only after exact validation;
4. otherwise generate into a temporary sibling directory;
5. validate schema, record counts, and output checksums;
6. atomically replace the stale family only after complete success.

Complete-family checksum validation belongs to preparation, not to ordinary
test lookup. A preparation process MUST read and validate every output in the
family once before publishing or reusing it. A test process MUST validate only
the manifest entry and output file that the test actually consumes. It SHOULD
combine checksum and record validation with its normal sequential read when
the loader API permits that. It MUST NOT rescan unrelated outputs in the same
family. In particular, separate test binaries MUST NOT each hash a complete
multi-gigabyte expected-data family.

Input checksums determine whether preparation may reuse a cached family.
Output checksums detect incomplete or corrupt generated files. Ordinary tests
rely on the successfully prepared family identity and validate their requested
output; they do not repeat preparation's full input-and-output audit.

A generator MUST NOT modify fixtures, corpora, known-failure declarations,
production source, or unrelated expected-data families.

Live official-source oracle runners are not expected-data generators. They
belong under `tools/oracles/<reference>/` and MAY be called by focused tests
when the required compiler and vendored oracle source are available.

## 6. Shared Test Support

Repository-root discovery, `testdata/` path resolution, profile selection,
manifest validation, checksum validation, and shared loading MUST have one Rust
implementation in a `publish = false` test-support crate.

The support crate MUST be used only as a dev-dependency, MUST NOT depend on a
production chemistry crate, and MUST NOT become a production runtime
dependency. Inline tests MAY use it under `cfg(test)`.

Python MAY use a separate thin loader when required by the language boundary,
but it MUST resolve and validate the same data and manifests.

## 7. Cross-Layer Coverage

The same fixture, corpus, or expected output MAY be consumed by multiple
layers, but the same behavior SHOULD NOT be exhaustively retested at every
layer.

- Rust core tests MAY run complete chemistry parity corpora.
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
