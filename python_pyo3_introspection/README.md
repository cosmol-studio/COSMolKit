# Python introspection experiment

An isolated PyO3 experiment using real `cosmolkit` calls. Production bindings,
chemistry and release workflows are unchanged.

Python-specific rules are defined in [BINDING_CONTRACT.md](BINDING_CONTRACT.md).

## Conventions

- Rust owns behavior. Binding methods forward to its public facade.
- Official PyO3 introspection supplies normal input/output types.
- Native bindings decide whether a callable accepts Params only, direct keywords,
  or both. The generator never flattens configuration classes or changes calls.
- Shared generation groups explicitly marked overload branches, refines
  integer/slice overloads, applies explicit return-type declarations and adjusts
  property setters for BasedPyright.
- Keep these transformations general; do not add per-business stub patches.
- Follow these conventions directly. Representative runtime and static usage
  checks validate the binding sample; registry reconciliation is independent.

## Run

Fetch dependencies once when setting up a fresh checkout:

```bash
cargo fetch --manifest-path python_pyo3_introspection/Cargo.toml
```

Build and check the complete representative binding sample:

```bash
.venv/bin/python python_pyo3_introspection/check.py \
    --extractor target/pyo3_introspection_extractor/debug/pyo3-introspection
```

The extractor can instead be on PATH as `pyo3-introspection`. It must come from
the same PyO3 Git revision pinned in `Cargo.toml` as the binding library.

The command does five things, once each:

1. Build the experimental extension and read its actual filename from Cargo.
2. Extract the official stub.
3. Generate the Python package and its public `__init__.pyi`.
4. Run the ordinary calling examples in `tests/runtime.py`.
5. Type-check `tests/positive.py` against the public stub.

There is one explicit interpreter setting for the Rust build:
`PYO3_PYTHON=sys.executable`. Runtime examples receive the generated package
directory as an argument, and the type checker uses that same interpreter.
The runner does not replace PYTHONPATH or assume a Linux library filename.

Results live in ignored `artifacts/smoke/`. The original stub stays under
`official/` as generation input; the generated package does not contain a
second `_native.pyi`. Static usage checks do not recheck the original stub
or the generated Python re-export module. Rust build artifacts are reused
under `target/typed_facade_projection/`.

This command does **not** build, install or validate a wheel. No publication
workflow is changed by this experiment.

The command passes all 55 runtime checks across 15 public exports and reports
zero public-usage type errors/warnings. Precise NumPy return annotations come
from explicit binding metadata, not automatic upstream NumPy inference.

## How the shared generation works

Params classes use a small Rust declaration:

```rust
parameters!(MaccsFingerprintParams, ck::MaccsFingerprintParams, {
    n_bits: usize,
});
```

It provides the constructor and writable fields only. It does not imply that a
callable accepts those fields as separate keywords: some configurations require
special initialization or methods and must remain parameter objects.

The native method explicitly declares its supported call forms. For example:

```rust
#[pyo3(signature = (params=None, *, n_bits=None))]
fn fingerprint_maccs(
    &self,
    params: Option<&MaccsFingerprintParams>,
    n_bits: Option<usize>,
) -> PyResult<Fingerprint>
```

The main entrypoint rejects mixed non-None Params/keywords and calls one of two
real, registered Rust branch functions. Params is forwarded by reference; the
keyword branch builds the configuration and calls the Params branch, which owns
the single facade call. None means omitted. No branch has a dummy body.

Official introspection emits one combined signature, not two overloads:

```python
def fingerprint_maccs(
    self, /, params: MaccsFingerprintParams | None = None,
    *, n_bits: SupportsIndex | None = None,
) -> Fingerprint: ...
```

The same official extraction includes both registered branches. Each carries
its target as machine-readable documentation, for example:

```rust
#[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
#[pyo3(signature = (params=None))]
fn _fingerprint_maccs_params(
    &self, params: Option<&MaccsFingerprintParams>,
) -> PyResult<Fingerprint>
```

`#[doc]` is supported Rust/PyO3 documentation metadata, not a new PyO3 attribute.
The official stub retains it. Keep internal branch names underscore-prefixed;
module targets use the function name without a class prefix. The generator
reads the marked declarations, groups them by target, renames their signatures,
adds `@typing.overload`, and replaces the broad entrypoint declaration. It hides
the helper names and metadata from the public stub, without reading Params
fields or reconstructing types/defaults. Targets and branches share this
grouping rule without changing the generator.

The public result is:

```python
@typing.overload
def fingerprint_maccs(
    self, /, *, n_bits: SupportsIndex | None = None,
) -> Fingerprint: ...

@typing.overload
def fingerprint_maccs(
    self, /, params: MaccsFingerprintParams | None = None,
) -> Fingerprint: ...
```

Both valid forms type-check. A one-off check of mixed Params/non-None keywords
reports two expected `reportCallIssue` diagnostics for the method and module
sample; runtime also rejects those mixtures. The diagnostic is retained in
ignored `artifacts/mixed_calls.py` and is separate from the ordinary check command.

Generated `__init__.py` only re-exports native objects, without Python argument
wrappers or method replacement. Internal branch methods remain callable through
their native underscored names; they are implementation details, not additional
public APIs.

The existing collection bindings delegate to shared Rust integer/slice indexing.
Their stub union is refined into two overloads. Runtime class objects are
re-exported, not duplicated; the public stub is a derived declaration.

NumPy bindings directly return `Bound<PyArray1<u8>>` / `Bound<PyArray2<u8>>` and
use Rust NumPy's `IntoPyArray`. Return annotations use the shared metadata rule.
Official Rust
NumPy 0.30.0 does not implement PyO3's introspection TYPE_HINT for PyArray, so
the original official stub still says `-> Incomplete`. To supply the public
annotation explicitly, the actual native method declares:

```rust
#[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int], numpy.dtype[numpy.uint8]]"}"#]
fn to_numpy<'py>(
    &self, py: Python<'py>,
) -> PyResult<Bound<'py, numpy::PyArray1<u8>>>
```

The shared rule imports the named modules, parses the declared Python type
expression and replaces only that callable's public return annotation. It
removes the metadata docstring from the public stub. The matrix method declares
`tuple[int, int]` in the same way. Other external return types use the same role;
adding an interface, dtype or dimensionality does not add generator branches.

This is a handwritten annotation owned by the binding, not Rust type inference:
Rust does not prove that its string agrees with the actual return type. Runtime
dtype/shape/copy checks and the existing precise NumPy type examples both pass.
The JSON remains visible in native runtime `__doc__`; generation does not alter
runtime documentation, calls or the official input stub.

NumPy comes directly from crates.io; its source is not adapted or patched.
The `[patch.crates-io]` entries in `Cargo.toml` route its PyO3 dependencies to
the binding library's pinned Git source; they do not modify NumPy itself.

Structured error conversion remains a shared Rust helper. Error variant
projections are explicitly declared, and unprojected source causes become
RuntimeError messages; this is not universal automatic error projection.

## Separate registry check

```bash
.venv/bin/python python_pyo3_introspection/registry_check.py
```

This checks all canonical records enabled by facade `core` and `fingerprints`.
It reads the already built registry reader and generated package; it never
compiles, generates or selects a hand-maintained interface scope.

The sample implements only a small part of the enabled registry, so this command
must report missing bindings rather than pass. Complete diagnostics are retained
in ignored `artifacts/registry.json`.

Current result: 795 canonical records yield 807 expected paths; runtime and stub
each match 18 and miss 789. Four sample exports are not registered. This is an
incomplete binding sample.

Registry reconciliation checks presence, field names and assignability; it does
not fully verify signatures, defaults or setter effects. Assigning a field its
existing value does not prove that its setter updates the underlying value.
