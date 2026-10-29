# JavaScript Binding Contract

Bindings project the canonical Rust facade. Language adaptation must not change
chemistry, defaults, value/in-place behavior or native ownership. The shared
[public API design](../dev/public_api_design.md) remains authoritative.

## Registry projection

Canonical `javascript_projection` metadata distinguishes ordinary runtime
exports (the default), `type_only` declarations, `native("TypeScript type")`
values, `error_details` accessors and `host_filesystem` APIs. These are language
representation/platform boundaries, never exemptions for unfinished bindings.
Filesystem-only types carry that boundary for their members; text IO remains
required. No checker-local API exclusion list is allowed. Registered setters
remain callable functions even when a writable property delegates to the same
native accessor. Factory names use camelCase, without case-variant aliases.

Numeric enums remain numeric values. Approved native enum methods can project as
domain-qualified module functions, recorded by `module_function("exportName")`
on the canonical callable. `residueInfoKindName(value)` and
`coordinateZPolicyFromName(text)` call their Rust methods directly; checked inputs,
native errors and automatically generated declarations use the existing adapters.
Do not replace these enums with classes or patch their generated declarations.

## Checked value inputs

The `macros/` crate generates binding declarations from canonical Rust
configuration metadata, and feature-selected TypeScript unions. Defaults come from
Rust; generated accessors delegate to shared input/storage adapters. It contains
no chemistry, business-type mapping list or independent ownership/error policy.
Unsupported field projections must fail compilation, not fall back to `any`.

Conditional union members use the native branch's `#[cfg(...)]` in
`type_union!`; the compiler selects the declaration before wasm-bindgen emits
the package. Do not remove unavailable types by editing the generated `.d.ts`.

- Function arguments, constructors, setters and callback return values use the
  same fallible Tsify/Serde input path (`input::read` or `input::decode`). Preserve
  the original JS value until decoding; do not first accept a coercing primitive
  WASM ABI parameter. Output primitives do not require input wrappers.
- Reuse Serde's type and integer-range checks. Fixed-size arrays decode the
  complete sequence and use Rust's checked array conversion: the upstream
  WASM tuple decoder may otherwise ignore extra elements. Do not implement a
  parallel numeric validator or repeat checks in individual business bindings.
  Reject fractional/out-of-range integers and implicit string/boolean-to-number
  conversions with TypeError. Floating-point inputs retain the native type's
  semantics; do not invent a blanket ban on NaN or infinity.
- Explicit representation checks may narrow Serde's accepted forms to the
  declared JS type: bigint inputs must remain bigint, and Uint8Array inputs must
  be Uint8Array. Never convert through Number or construct a TypedArray before
  validation to hide an invalid original input.
- Complete conversion before changing native state. A failed setter leaves the
  previous value unchanged. Native Params, Molecule, executors and callback
  handles retain their native ownership adapters; do not serialize/reconstruct
  them to obtain input validation.
- Dynamic native-class extraction uses `input::with_instance`; concrete
  wasm-bindgen bridges perform class/pointer validation. Bindings supply only
  their native read/write access. Do not repeat visitor/result/error bookkeeping
  or introduce a business-type lookup table. Direct typed native arguments keep
  their normal ABI; ordinary serializable inputs continue to use Tsify.
- Shared-only native configuration accessors return owner-retaining readonly
  views, typed `Readonly<Params>`, not editable detached copies. Writes fail
  before mutation. Explicit `copy()` produces independent editable configuration;
  replacement uses the owner's canonical setter.
- Generate declarations from actual binding/Tsify types. Transparent input
  wrappers remain ordinary JS values, not user-created wrapper objects. No
  hand-edited `.d.ts`, per-function validation patches or framework business-type
  mapping lists. Validate this mechanism with a few framework samples, not a
  growing per-interface consumer inventory.

```rust
#[wasm_bindgen(js_name = fromAtomicNumber)]
pub fn from_atomic_number(
    number: Ts<input::U8Input>,
) -> Result<Option<Element>, JsValue> {
    let number = input::read::<u8>(number)?;
    Ok(ck::Element::from_atomic_number(number).map(Element))
}
```

```javascript
Element.fromAtomicNumber(6);     // Carbon.
Element.fromAtomicNumber(262);   // TypeError, never wraps to 6.
Element.fromAtomicNumber(6.9);   // TypeError, never truncates to 6.
Element.fromAtomicNumber("6");   // TypeError, no string coercion.
```

The generated declaration uses `U8Input = number`; runtime decoding enforces
the range that TypeScript's `number` cannot express.

## Errors

Domain failures are concrete `Error` subclasses owning the native Rust error.
`name` identifies the error class; `message` and `toString()` preserve Rust
Display; `variant` is the exhaustive enum discriminant when one exists.
All native payload fields live in the strongly typed readonly `details` value,
including fields named `name`, `message`, `kind` or `source`. No per-error aliases.
`cause` retains the actual recursively projected Rust source and stable object
identity, separate from payloads which merely contain another error value.

Native synchronous iterators expose `[Symbol.iterator]` and typed `next()`
results. Iterator handles share the native cursor; they do not eagerly collect
or reset it. Completion is `{done: true, value: undefined}`.

```javascript
error.message;                  // Native Display.
error.details.message;          // Actual native message field, when present.
error.cause.details.offset;     // Typed native cause payload.
```

Class-extraction failures are TypeError at the shared input boundary. Native
operation failures and user callback exceptions keep their actual type and
identity; do not classify errors by matching their messages.
