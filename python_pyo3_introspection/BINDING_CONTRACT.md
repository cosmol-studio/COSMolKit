# Python Binding Contract

Bindings may adapt representation and calling syntax, not Rust business behavior.
This contract supplements [public API design](../dev/public_api_design.md).
It governs binding implementation, generated types and package behavior.

## Forwarding and configuration

- Call the canonical public Rust implementation; do not duplicate algorithms,
  introduce heuristics, retry failures, or substitute empty/default results.
- Allow only language adaptation: typed extraction, explicit enum-string
  conversion, configuration dispatch, Python indexing and result conversion.
- Bindings explicitly choose supported Params/keyword forms. Preserve supplied
  Params state; do not flatten and rebuild it. Omitted options use Rust defaults
  or Rust initialization methods, not a second Python policy.
- Reject mixed Params/keywords with TypeError based on argument presence, not
  non-None values. Explicit None still selects a call form: supplying params
  together with any configuration keyword is an error even if either is None.
  Preserve omission separately using a typed native argument-presence extractor;
  do not use Option::is_some to decide whether an argument was supplied.
  Use `Argument<Option<T>>` at the entry point and `configuration_arguments!`
  before unwrapping values; keep shared arguments outside its configuration list.
  None may retain its declared default/value meaning within one selected form.
  Independent arguments shared by the declared forms are not conflicting fields.
  Constructors/setters share conversion rules; failed
  assignments preserve the previous value.

An optional Python constructor input over concrete Rust storage uses
that Rust field's actual default for omission/None. The getter returns the
stored concrete value; an explicit empty collection remains empty.
Constructor and setter apply the same conversion. This rule does not apply to
genuinely optional Rust fields, where None and a present empty value stay distinct.

Configuration follows the common public API contract. Rust defines the input
values, composition, defaults and execution lifecycle. Bindings only adapt their
representation. The wrapper holds the actual Rust configuration as its sole
input storage. Constructors/setters use the same typed extraction and replace
native fields only after conversion succeeds. Sequence getters return tuples;
map getters return readonly mappings, recursively projecting their contents.
Do not expose list.append/dict item assignment that silently edits a copy.
Nested configuration getters return independent values: edit and assign back
explicitly. Native callbacks, random sources and executors retain their actual
state and ownership rather than being rebuilt as input records. Capture native
invocation inputs before execution or callbacks. Keyword paths construct inputs
directly, without allocating an intermediate editable Python configuration.

Configuration bindings use native `Configuration` metadata and
`ck::project_configuration!(ck::TypeName)`. Do not repeat their fields or types
in binding declarations. Rust-resolved field aliases select shared typed
projections; generated native struct literals enforce completeness. Required
inputs, None/default interpretation, composition and registered factories are
explicit native annotations. Fixed arrays return tuples; setters retain native
length and element-type checks.

Stateful/accessor records require an explicit ownership adapter, for example
`ck::project_configuration!(ck::EmbedParams, using embedding_parameters)`.
The native declaration marks fields handled as state or by that adapter;
ordinary fields still come from native metadata. Keep native initialization,
callbacks, diagnostics and random-source ownership in the existing adapter;
never rebuild them from a Default record. Unsupported shapes fail compilation.

```python
p = ck.SmartsParseParams()
source = {"{A}": "O"}
p.replacements = source        # Receives an independent value.
source["{A}"] = "C"            # Does not change p.
rules = dict(p.replacements)    # Explicit editable copy of the readonly mapping.
rules["{A}"] = "N"
p.replacements = rules         # Replaces the actual Rust field.
# p.replacements["{A}"] = "N" raises TypeError.
```

## Callbacks

Callable configuration uses concrete native Rust wrappers whose official PyO3
input/output type is a concrete `Callable[[...], bool]` or a void notification
`Callable[[], None]`, never public Any/PyAny. The
wrapper checks callable extraction and forwards actual Rust callback fields.
Opaque CPython handles are private storage, not a second algorithm.

Callback errors belong to one invocation. Rust's bool callback receives false
after a Python failure; the binding then restores the original Python exception
object before returning the operation result. Do not stringify or replace that
exception. Void progress callbacks retain the Rust owner's one-notification-per-row
behavior even after an earlier callback failure, retaining the first error.
Release the calling Python interpreter while native parallel workers may call
Python; otherwise they cannot acquire it to deliver notifications.
A final-match callback receives an independent owned COW Molecule value for
the corresponding matching input. Every callback gets its own value; Python
object identity with the source is not preserved. It may retain or mutate that
value without changing the input, later callbacks or generated products.
Only its boolean return accepts/rejects the candidate match. Original input
receivers remain subject to normal shared/exclusive borrowing, including during
in-place reaction execution. Never manufacture a reentrant loan or substitute
value execution/writeback for an actual in-place call. Configuration
construction/callability alone does not establish actual callback execution.
Configuration classes retaining Python callbacks expose those references to
Python's cyclic collector and clear them through the native GC protocol;
retaining a callback must not create an uncollectable configuration cycle.
Compound reaction values retain matching configuration callbacks through value
and removal results, with real collector-visible ownership. A getter returns a
configuration value copy, not a mutable alias of the reaction's stored options.
Reaction matching operations must install those callbacks per invocation using
the corresponding native COW input values and fresh error scope as ordinary matching;
retaining a configuration alone is not proof of executing a reaction callback.

Random-bits providers use the actual shared native source. Integer seeds and
unsigned results use PyO3's native arbitrary-integer conversion, never decimal
text or fixed-width truncation. Each synchronous native iteration call has a
thread-local stacked error scope: nested calls and threads cannot steal another
invocation's original Python exception. Scope creation does not rebuild or
reset the retained provider. A negative provider result raises the actual
unsigned-conversion OverflowError; the binding invents no numerical fallback.
The Python random source is the single collector-visible callback owner;
native provider clones refer weakly back to it. Options retain that actual
Python source and return it from the getter. Every public iterator/operation
retaining the native provider must also retain its Python owner for the entire
use. An improperly retained bare native clone raises ReferenceError after its
owner expires; this is not a normal enumeration outcome or a numerical fallback.

Tautomer callbacks receive live, callback-scoped source/progress objects.
Candidate views and iterators derived from live progress share its expiry;
access after return or exception raises ReferenceError. Explicit snapshot() or
to_owned() retains a native owned value: candidate snapshots share their actual
record Arc and coordinate owner, while a borrowed source snapshot explicitly
copies its source blocks. Ordinary access does not request a snapshot or force
candidate COW. Progress builds a reference-only index on its first entry lookup;
subsequent lookups/iterator steps are O(1), without a persistent mutable loan.
Scoring must act on the actual native short loan, not a copied graph plus cache
copy-back. Atom/Bond/Properties query results follow the existing detached-value
projection; they are independent copies, not expiring views. In particular,
properties() copies the property collection and is not a zero-copy operation.
Only live loans are callback-thread-affine: attempting to read one from another
thread raises RuntimeError before touching native borrowed storage. Explicit
owned snapshots retain the native value's cross-thread usability, including
their candidate views and iterators. Holding a Python handle does not extend a
live callback loan's lifetime.

```python
retained = []
def inspect(source: ck.TautomerMoleculeView, progress: ck.TautomerProgress) -> bool:
    key, candidate = next(progress.entries())
    retained.append(candidate.snapshot())  # Shared native candidate ownership.
    retained.append(progress.snapshot())   # Explicit, not performed on callback entry.
    return True
# Saving candidate itself instead would retain only an expired callback loan.
```

```python
params = ck.SubstructMatchParams(final_match=lambda target, indices: len(indices) > 1)
# final_match: Callable[[Molecule, Sequence[int]], bool] | None
params.final_match = None
```

## Errors

Raise exceptions for Rust failures; never return an error-shaped success value.
Each Rust error type has one consistent binding projection:

| Rust error shape | Python projection |
|---|---|
| Variant without fields | Typed exception with readonly domain and enum variant |
| Variant with fields | Same, plus readonly fields preserving Rust names, types and meanings |
| Error containing a cause | Same, with the concrete cause recursively projected through __cause__ |

- str(error) uses Rust Display. Expose message only when it is an actual Rust
  payload field; never infer variant or fields by parsing a message.
- Preserve structured cause types and fields, not just their text. System and
  third-party mappings must be explicit; do not blanket-convert causes to RuntimeError.
- Python call/type errors remain language exceptions, separate from domain errors:
  TypeError for incompatible inputs, ValueError for unknown enum strings, and
  IndexError for Python container index bounds.

Python enum errors use readonly `variant` as their discriminant. Preserve actual
Rust `kind` payload fields and canonical `kind()` methods unchanged; the
synthetic discriminant must not overwrite them. Companion enum type names may
retain their `ErrorKind` spelling. Record errors invent no `variant`.
Actual Rust fields that are Python keywords receive only the normal syntax
escape, for example RingDecomposerError's `from` is exposed as readonly `from_`.
This is declared by the native binding, never repaired in generated stubs.

```python
try:
    ck.expand_one_letter_sequence("!", "aa")
except ck.ResidueSequenceError as error:
    assert error.variant == ck.ResidueSequenceErrorKind.UnexpectedLetter
    assert error.kind == "peptide"  # Actual Rust field, not the discriminant.
    assert error.letter == "!" and error.source_code == 33
```

The registry identifies Rust error types using its error role. A Python discriminator
enum is that registered error's discriminant projection, identified by the
official typed `variant` getter declaration. Registry reconciliation requires the
referenced companion class in both runtime and stub; it does not accept names
merely because they end in `Kind`, nor maintain a second enum inventory.

Canonical callable `error_type` declarations also identify direct error classes.
The reader resolves existing registered aliases, or the unambiguous declared
Rust type name, and rejects naming conflicts. No redundant type inventory is
maintained. For example, the existing `EmbedParams.with_json -> Result<EmbedParams,
ConformerError>` declaration requires a real typed `ConformerError` exception
and its discriminant/payload projection even without a separate type record.

A Rust record error without a discriminant stays a record error. The native
exception declares `#[doc = r#"{"role":"error_record"}"#]`; the official class
metadata identifies this shape to the checker, not to the generator. Its actual
readonly payload accessors and real public Rust fields must remain typed, with
Rust field names and native Display preserved. Public field getters do not
require a duplicate manually maintained accessor registry. Exhaustive Rust
destructuring guards the public record shape; never expose or invent private
fields. A source-only record instead requires a concrete typed readonly
cause getter. No empty record
exception is accepted merely because it declares this role. For example:

```rust
struct PropertyValueError { expected: PropertyValueKind, actual: PropertyValueKind }
// Canonical public methods: expected() and actual(), no kind().
```

```python
try:
    ck.PropertyValue.Int(7).as_bool()
except ck.PropertyValueError as error:
    assert error.expected() == ck.PropertyValueKind.Bool
    assert error.actual() == ck.PropertyValueKind.Int
    # No fabricated error.variant or PropertyValueErrorKind.
```

Starting only from canonical registered/direct error classes, reconciliation
also follows official native readonly getters explicitly marked
`#[doc = r#"{"role":"error_source"}"#]`, using their concrete exception
annotation (including precise unions), recursively. Use `source` when available,
or `cause` if a real Rust payload already uses `source`. The checker does not
guess from the field name. For example, ProteinProjectionError.source remains
its integer entity-ID payload, while its marked cause getter returns the actual
BioStructureError cause or None. Both spellings refer to the actual __cause__
object, never a reconstructed copy; never override Python's built-in __cause__
descriptor. These companions represent
the real Rust source graph, not another API inventory. Every referenced class
must exist as a real exception at runtime and in the stub and satisfy the same
typed error shape rules. Broad Exception/BaseException/Any annotations, naming
suffixes and manually maintained source lists are not proof. Native conversion
must preserve the actual Rust Error::source graph and available payloads.

```rust
struct BioSelectionMatchError { cause: BioRowTraverseError }
// Rust Error::source returns the stored BioRowTraverseError; no kind().
```

```python
try:
    structure.selected_atom_ids(selection)
except ck.BioSelectionMatchError as error:
    assert error.source is error.__cause__  # same actual exception, not a copy
    assert isinstance(error.source, ck.BioRowTraverseError)
    assert isinstance(error.source.source, (ck.BioRowModelError, ck.BioRowChainError))
    # No invented BioSelectionMatchErrorKind; str(error) remains Rust Display.
```

This exception remains a real native exception owning the Rust error. Record
metadata does not exempt arbitrary classes from payload validation, and enum
errors still require their typed discriminant companion.

A concrete standard-library Display-only, source-free leaf may declare native
`role=error_leaf`, but only as a companion reached through a real marked
source/payload edge. It owns the actual Rust error, exposes readonly domain and
Rust Display, and invents no kind or fields. It cannot replace an error with
public payloads or a further source; `std::io::Error` is not such a leaf.

```rust
// XyzReadError::Coordinate { source: std::num::ParseFloatError, ... }
// ParseFloatError has no public kind/payload or Error::source.
```

```python
cause = xyz_error.source
assert isinstance(cause, ck.ParseFloatError)
assert cause is xyz_error.__cause__ and cause.__cause__ is None
# cause.domain is readonly; str(cause) is the original Rust Display.
# No fabricated ParseFloatErrorKind.
```

An actual readonly native error payload may declare
`#[doc = r#"{"role":"error_payload"}"#]`. Reconciliation follows its official
concrete value/error annotation from an accepted error root, including precise
optional/union/container types, and requires matching native readonly descriptor
metadata. This allows real enum discriminants such as ValencePhase and concrete
nested error payloads; it does not register arbitrary business classes, functions
or constructors. Unregistered simple enums must have the official native
enum-constant shape and matching runtime variants, without business methods.
Rich non-error records/tagged enums declare native class metadata
`#[doc = r#"{"role":"error_value"}"#]` on a frozen PyO3 class. They expose actual
typed readonly native getters, no constructors, setters or business methods.
Reconciliation recursively follows those getter types, including recursive
values inside concrete built-in containers, starting only from a real accepted
error-payload edge. An unrelated marked class is still an unregistered export.
This role neither makes the value an exception nor creates an Error::source edge.
Any/object annotations and missing/mismatched payload classes fail. Payload
metadata never changes Rust Error::source or invents a cause:

```rust
enum ValenceError {
    InvalidValence { phase: ValencePhase, /* actual remaining Rust fields */ },
}
// The readonly native phase getter returns Option<ValencePhase> and carries
// role=error_payload. A non-source enum field does not become __cause__.
```

```python
try:
    ck.mol_from_smiles("[CH5]", sanitize=False).atom_metadata()
except ck.ValenceError as error:
    assert error.phase == ck.ValencePhase.Explicit
    assert error.__cause__ is None
```

For example, an XML failure's real `position` payload may be a readonly
`TextPos` value exposing integer `row` and `col`. Its class carries `error_value`,
and the exception's native `position` getter carries `error_payload`. The value
remains data, not `error.__cause__`. A recursive query predicate similarly keeps
its actual typed discriminant and child predicates rather than returning Debug
text or a dictionary with unspecified values.

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
error.variant    # ck.FingerprintErrorKind.Unsupported
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
    assert error.variant == ck.FingerprintErrorKind.SparseIndexOutOfRange
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
    assert error.variant == ck.SmilesErrorKind.Parse
    cause = error.__cause__
    assert isinstance(cause, ck.SmilesParseError)
    assert cause.variant == ck.SmilesParseErrorKind.Syntax
    assert cause.offset == 3
    assert cause.message == "unclosed branch"
```

The outer SmilesError preserves the concrete Parse, Hydrogen, Sanitize, Stereo
and Construction causes. Construction exposes the complete OperationError;
its real source edges remain distinct from readonly error payload fields.

## Values and ownership

- Molecule pickle uses the unchanged reconstruction path
  `cosmolkit._molecule_from_binary(bytes)`. `__reduce__` returns that native
  callable and one binary argument; `__reduce_ex__` validates Python's index
  protocol and uses the same reduction. Binary encoding/decoding belongs only
  to the canonical Rust API. An explicit `private_export` declaration re-exports
  the actual native function at package level, attributes its module identity
  to that public package, and keeps it outside public `__all__`. No callable is
  wrapped/replaced and native module initialization remains statically described.
- Consuming Rust resource methods transfer their resource exactly once. After
  `stream.batches(...)` enters the native consuming call, every original stream
  reference is unusable and raises RuntimeError, including if native parameter
  validation fails. The returned iterator owns the current cursor and continues
  from it. Python argument conversion failures before transfer retain the stream.
  Never reopen files, duplicate cursors or copy native validators in the binding.
- Public types must be explicit; do not hide known types behind PyAny/Any.
  Preserve optional values, integer ranges, precision, order and batch error positions.
- NumPy conversions declare dtype, shape and copy/sharing behavior. Convert
  representation without recomputing or altering domain results.
- Preserve Rust value/in-place semantics and failure atomicity. Do not mutate
  inputs from value calls or clone whole state merely to implement a wrapper.
- Python `Molecule.copy()` and `copy.copy(molecule)` forward Rust Clone directly:
  independent value semantics, initially shared molecular storage, writes
  isolated through normal operations. No eager deep-copy or universal O(1)
  promise: native Clone also copies its descriptor memo.

## Declarations and checks

- Official introspection supplies ordinary declarations. Explicit overload
  branches are registered, executed functions; the native main entrypoint dispatches.
- Indexing overloads require explicit native documentation
  `{"role":"indexing","protocol":"integer_slice_self"}` on `__getitem__`.
  This declares integer indexing returns the item and slicing returns the same
  collection type. The extracted signature must accept exactly an integer-like
  index or slice and return concrete item types or that collection. Generation
  rejects incompatible marked signatures; unmarked indexing stays unchanged.
  The shared Rust indexing implementation owns actual runtime behavior. Other
  indexing shapes use their ordinary signature or existing explicit overload
  branches, not new per-collection generator rules.
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
#[pyo3(signature = (params=Argument::Missing, *, n_bits=Argument::Missing))]
fn fingerprint_maccs(
    &self,
    py: Python<'_>,
    params: Argument<Option<PyRef<'_, MaccsFingerprintParams>>>,
    n_bits: Argument<Option<usize>>,
) -> PyResult<Fingerprint> {
    let keyword_form = n_bits.is_supplied();
    crate::arguments::configuration_arguments!(params; n_bits);
    if keyword_form {
        self._fingerprint_maccs_keywords(n_bits)
    } else {
        self._fingerprint_maccs_params(py, params.as_deref())
    }
}

#[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
#[pyo3(signature = (params=None))]
fn _fingerprint_maccs_params(
    &self,
    py: Python<'_>,
    params: Option<&MaccsFingerprintParams>,
) -> PyResult<Fingerprint> {
    let effective = params.map(|params| params.to_rust(py)).transpose()?.unwrap_or_default();
    self.inner.fingerprint_maccs_with_params(&effective)
        .map(|inner| Fingerprint { inner })
        .map_err(errors::MaccsFingerprintError::from_rust)
}

#[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
#[pyo3(signature = (*, n_bits=None))]
fn _fingerprint_maccs_keywords(&self, n_bits: Option<usize>) -> PyResult<Fingerprint> {
    let params = ck::MaccsFingerprintParams {
        n_bits: n_bits.unwrap_or_else(|| ck::MaccsFingerprintParams::default().n_bits),
    };
    self.inner.fingerprint_maccs_with_params(&params)
        .map(|inner| Fingerprint { inner })
        .map_err(errors::MaccsFingerprintError::from_rust)
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
mol.fingerprint_maccs(None)              # Params form, Rust defaults
mol.fingerprint_maccs(n_bits=None)        # Keyword form, Rust defaults
mol.fingerprint_maccs(params, n_bits=166) # Type-check error and runtime TypeError
mol.fingerprint_maccs(params, n_bits=None) # Same error: both forms were supplied
mol.fingerprint_maccs(None, n_bits=166)    # Same error: explicit None is not omission
```

Call-form presence and optional values are separate: None does not erase a
supplied argument. The native dispatcher forwards to one actual branch only.
Generation does not infer fields from Params, add missing None combinations
or implement runtime dispatch. A new target or branch uses the same role/target
rule without extending the generator.
