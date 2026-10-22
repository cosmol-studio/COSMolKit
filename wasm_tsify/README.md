# WASM binding protocols

Small executable samples for wasm-bindgen + Tsify. Existing `wasm/` and release
entrypoints are unchanged. Chemistry calls use only the public `cosmolkit` facade.

## Run

Requires the `wasm32-unknown-unknown` Rust target, wasm-bindgen-cli 0.2.128,
Node.js and npm.

```bash
.venv/bin/python wasm_tsify/build.py
node wasm_tsify/tests/runtime.cjs
npx --yes --package typescript@5.8.3 tsc -p wasm_tsify/tests/tsconfig.json
```

The build compiles one extension and generates `pkg/` including its `.d.ts`.
There is no declaration rewriting. Runtime and typing checks are separate commands.

## Working boundaries

| Boundary | Representation and responsibility |
|---|---|
| Native objects | wasm-bindgen classes retain Rust state; borrowed arguments do not consume JS handles. `copy()` calls native COW copy. |
| Scalar Params | `parameters.rs` consumes native Configuration metadata. Field names/types/defaults are not redeclared in the binding. |
| Nested Params | Construction/assignment copies the value; the child getter shares its actual owned state. The sample child survives parent disposal. |
| Container fields | Frozen readonly snapshots on get; explicit validated whole-field assignment on set. Failed conversion leaves state unchanged. |
| Value objects | Tsify derives schema and Serde converts inside the function through `Ts<T>`; invalid input returns an exception. |
| Arrays | Fixed Rust arrays become TS tuples, numeric buffers become TypedArrays, readonly containers have explicit readonly declarations plus runtime freezing. |
| Enums | String enums and tagged payload unions derive from the same Rust declaration. |
| Callbacks | Explicit function signature on the binding; native COW snapshot passed to JS; callback exceptions propagate unchanged. A non-boolean result is rejected. |
| Errors | The transport sample preserves structured details and the actual JS cause on an `Error`. `failure_details()` validates/narrows an unknown error payload. |
| Large integers | `u64` uses JS `bigint`; it is not converted through `number`. |

```typescript
const params = new RingSearchParams();
params.include_dative_bonds = true;
const molecule = new Molecule('C1CC1');
const count = molecule.num_rings(params);
const nested = new NestedParams(params);
nested.child.include_hydrogen_bonds = true;
const selection = new Selection([0, 1]);
selection.indices = [2]; // selection.indices.push(2) is rejected
```

## Limits before full binding adoption

- Tsify handles serializable values, not native class ownership, methods or
  arbitrary Rust callbacks. Those remain explicit bindings. The callback sample
  is synchronous; it does not establish worker/thread or async callback support.
- The generated Params macro currently covers native scalar fields. Nested and
  container samples establish behavior, but are not yet automatic projections
  of arbitrary Configuration field shapes. Associated native-to-JS type adapters
  must live beside their binding declarations, never in a framework business map.
- No Params/options overload protocol is claimed here. These samples use a
  native Params argument or a value object, not an inferred union of both.
- Structured error **transport** is exercised, not complete CK error projection:
  the two chemistry demonstration calls currently expose their error messages.
  Their domain error types/fields must be projected before production adoption.
- wasm-bindgen scalar ABI conversions follow JS coercion rules; TypeScript is
  not runtime validation. Tsify value inputs use fallible Serde validation.
- Owned TypedArray output copies data. No zero-copy view across memory growth
  or release is promised. Finalization is nondeterministic; explicit `free()`
  and generated disposal remain available.

Tests contain six protocol groups and six negative typing witnesses, not a
business-interface inventory. New bindings reusing a proven protocol must not
expand these consumers. This prototype proves the listed representations work;
it does not yet freeze every protocol needed by the full Python inventory.

Upstream: [Tsify](https://github.com/madonoharu/tsify),
[fallible Ts conversion](https://docs.rs/tsify/0.5.8/tsify/struct.Ts.html).
