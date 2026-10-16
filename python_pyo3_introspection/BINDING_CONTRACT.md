# Python Binding Contract

Bindings may adapt representation and calling syntax, not Rust business behavior.
This contract supplements [public API design](../dev/public_api_design.md).
It defines requirements, not a claim that the experiment already satisfies all of them.

## Forwarding and configuration

- Call the canonical public Rust implementation; do not duplicate algorithms,
  introduce heuristics, retry failures, or substitute empty/default results.
- Allow only language adaptation: typed extraction, explicit enum-string
  conversion, configuration dispatch, Python indexing and result conversion.
- Bindings explicitly choose supported Params/keyword forms. Preserve supplied
  Params state; do not flatten and rebuild it. Omitted options use Rust defaults
  or Rust initialization methods, not a second Python policy.
- Reject mixed Params/keywords with TypeError. Define whether None, omission and
  empty values differ. Constructors/setters share conversion rules; failed
  assignments preserve the previous value.

## Errors

Raise exceptions for Rust failures; never return an error-shaped success value.
Each Rust error type has one consistent binding projection:

| Rust error shape | Python projection |
|---|---|
| Variant without fields | Typed exception with readonly domain and enum kind |
| Variant with fields | Same, plus readonly fields preserving Rust names, types and meanings |
| Error containing a cause | Same, with the concrete cause recursively projected through __cause__ |

- str(error) uses Rust Display. Expose message only when it is an actual Rust
  payload field; never infer kind or fields by parsing a message.
- Preserve structured cause types and fields, not just their text. System and
  third-party mappings must be explicit; do not blanket-convert causes to RuntimeError.
- Python call/type errors remain language exceptions, separate from domain errors:
  TypeError for incompatible inputs, ValueError for unknown enum strings, and
  IndexError for Python container index bounds.

### Example 1: variant without fields

Rust definition and returned value (excerpt):

```rust
enum FingerprintError {
    Unsupported,
}
// Err(FingerprintError::Unsupported)
```

The exception caught as `except ck.FingerprintError as error` exposes:

```python
error.domain  # "fingerprint"
error.kind    # ck.FingerprintErrorKind.Unsupported
str(error)    # Rust Display text; no invented payload fields
```

### Example 2: variant with fields

```rust
enum FingerprintError {
    SparseIndexOutOfRange { index: u64, size: u64 },
}
// Err(FingerprintError::SparseIndexOutOfRange { index: 9, size: 8 })
```

```python
fp = ck.Fingerprint.from_on_bits(8, [1, 4])
try:
    fp.get_bit(9)
except ck.FingerprintError as error:
    assert error.kind == ck.FingerprintErrorKind.SparseIndexOutOfRange
    assert error.index == 9
    assert error.size == 8
```

Fields are read from the owned Rust error, not extracted from its message.
For a shared exception class, fields absent from the active variant return None.

### Example 3: nested cause

```rust
enum SmilesError {
    Parse(SmilesParseError),
}
enum SmilesParseError {
    Syntax { offset: usize, message: String },
}
// Err(SmilesError::Parse(SmilesParseError::Syntax {
//     offset: 3, message: "unclosed branch".into(),
// }))
```

Required Python projection:

```python
try:
    ck.mol_from_smiles("C(C")
except ck.SmilesError as error:
    assert error.kind == ck.SmilesErrorKind.Parse
    cause = error.__cause__
    assert isinstance(cause, ck.SmilesParseError)
    assert cause.kind == ck.SmilesParseErrorKind.Syntax
    assert cause.offset == 3
    assert cause.message == "unclosed branch"
```

This last example defines the required contract; the current experiment's
RuntimeError cause fallback does not yet satisfy it.

## Values and ownership

- Public types must be explicit; do not hide known types behind PyAny/Any.
  Preserve optional values, integer ranges, precision, order and batch error positions.
- NumPy conversions declare dtype, shape and copy/sharing behavior. Convert
  representation without recomputing or altering domain results.
- Preserve Rust value/in-place semantics and failure atomicity. Do not mutate
  inputs from value calls or clone whole state merely to implement a wrapper.

## Declarations and checks

- Official introspection supplies ordinary declarations. Explicit overload
  branches are registered, executed functions; the native main entrypoint dispatches.
- return_type metadata is an explicit binding-owned annotation, not Rust inference.
  Its type expression must match runtime behavior. Doc metadata remains visible
  in native __doc__ unless a separate documentation policy changes that.
- Generation only applies shared declaration rules; no business-specific patches,
  Python dispatch wrappers or runtime method replacement.
- Runtime examples check actual calls; typing examples check public declarations.
  Reuse representative cases for shared mechanisms rather than duplicate tests
  for every new API. Registry reconciliation is separate and is not a second API definition.

### Example: native NumPy return with an explicit annotation

The actual Rust binding still returns a native NumPy array:

```rust
#[doc = r#"{
    "role": "return_type",
    "imports": ["numpy"],
    "type": "numpy.ndarray[tuple[int], numpy.dtype[numpy.uint8]]"
}"#]
fn to_numpy<'py>(
    &self,
    py: Python<'py>,
) -> PyResult<Bound<'py, numpy::PyArray1<u8>>> {
    let values = (0..self.inner.n_bits())
        .map(|bit| self.inner.get_bit(bit).map(u8::from))
        .collect::<Result<Vec<_>, _>>()
        .map_err(errors::FingerprintError::from_rust)?;
    Ok(values.into_pyarray(py))
}
```

Generation adds the declared import and replaces only the public return type:

```python
import numpy

def to_numpy(self) -> numpy.ndarray[tuple[int], numpy.dtype[numpy.uint8]]: ...
```

For a matrix, the binding declares `tuple[int, int]` and returns PyArray2<u8>.
The same generation rule handles either declaration; it does not recognize
NumPy APIs by name. Runtime conversion remains in Rust, and the original official
stub remains unchanged. The annotation is explicit, not automatically inferred.

### Example: two executed overload branches

Inside Molecule's PyO3 methods, register a main entrypoint and two real branches:

```rust
#[pyo3(signature = (params=None, *, n_bits=None))]
fn fingerprint_maccs(
    &self,
    params: Option<&MaccsFingerprintParams>,
    n_bits: Option<usize>,
) -> PyResult<Fingerprint> {
    if let Some(params) = params {
        if n_bits.is_some() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "parameter object and configuration keywords are mutually exclusive",
            ));
        }
        self._fingerprint_maccs_params(Some(params))
    } else {
        self._fingerprint_maccs_keywords(n_bits)
    }
}

#[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
#[pyo3(signature = (params=None))]
fn _fingerprint_maccs_params(
    &self,
    params: Option<&MaccsFingerprintParams>,
) -> PyResult<Fingerprint> {
    let defaults;
    let effective = if let Some(params) = params {
        &params.inner
    } else {
        defaults = ck::MaccsFingerprintParams::default();
        &defaults
    };
    self.inner.fingerprint_maccs_with_params(effective)
        .map(|inner| Fingerprint { inner })
        .map_err(errors::MaccsFingerprintError::from_rust)
}

#[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
#[pyo3(signature = (*, n_bits=None))]
fn _fingerprint_maccs_keywords(&self, n_bits: Option<usize>) -> PyResult<Fingerprint> {
    let params = MaccsFingerprintParams {
        inner: ck::MaccsFingerprintParams {
            n_bits: n_bits.unwrap_or_else(|| ck::MaccsFingerprintParams::default().n_bits),
        },
    };
    self._fingerprint_maccs_params(Some(&params))
}
```

The main entrypoint calls these functions; neither is a signature-only placeholder.
Generation groups their official signatures by target, hides the internal names
and emits (SupportsIndex is the official input type for usize):

```python
class Molecule:
    @typing.overload
    def fingerprint_maccs(
        self, /, params: MaccsFingerprintParams | None = None,
    ) -> Fingerprint: ...

    @typing.overload
    def fingerprint_maccs(
        self, /, *, n_bits: SupportsIndex | None = None,
    ) -> Fingerprint: ...
```

```python
mol.fingerprint_maccs()                  # Rust defaults
mol.fingerprint_maccs(params)            # Preserve the supplied Params
mol.fingerprint_maccs(n_bits=166)        # Explicit keyword branch
mol.fingerprint_maccs(params, n_bits=166) # Type-check error and runtime TypeError
```

Here None means omitted. Generation does not infer fields from Params or add
runtime dispatch. A new target or branch uses the same role/target rule.
