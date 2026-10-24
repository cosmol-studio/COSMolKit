# WASM binding protocols

Bindings and small protocol samples using wasm-bindgen + Tsify. Existing `wasm/`
and release entrypoints are unchanged. Chemistry calls use the public Rust facade.

Implementation rules: [Binding contract](BINDING_CONTRACT.md).

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

## Registry check

```bash
.venv/bin/python wasm_tsify/build.py --registry
.venv/bin/python wasm_tsify/registry_check.py
```

Build with `--registry` once to include inspection metadata from the canonical
Rust registry compiled into that exact WASM. Normal builds omit it. The checker
does not compile or run chemistry; it checks registered JS names, static/instance
placement and runtime/declaration shapes using TypeScript's parser. Missing
entries fail the command; the complete report is `target/registry.json`.
There is no second interface list. Canonical JS projection metadata identifies
type-only/native values, error payloads in `details` and host-filesystem APIs.
Filesystem exclusions and native projections remain visible in the report;
they are not counted as matched runtime exports. Python-only metadata is not
inferred as JS metadata. This is not behavioral or full Python-coverage
acceptance, and does not establish defaults or signature equivalence.

## Working boundaries

| Boundary | Representation and responsibility |
|---|---|
| Native objects | wasm-bindgen classes retain Rust state; borrowed arguments do not consume JS handles. `copy()` calls native COW copy. |
| Params/options | Canonical Configuration metadata generates classes and options. Calls borrow native Params or decode explicit options; names/types/defaults are not redeclared. |
| Checked inputs | `input.rs` wraps primitive values in transparent Tsify types. Serde validates original JS values before native conversion; generated setters use the same mechanism. |
| Nested Params | Construction/assignment copies the value; the child getter shares its actual owned state. A shared-only native accessor returns a readonly view; `copy()` explicitly creates editable configuration. Views survive parent disposal. |
| Container fields | Frozen readonly snapshots on get; explicit validated whole-field assignment on set. Failed conversion leaves state unchanged. |
| Value objects | Tsify derives schema and Serde converts inside the function through `Ts<T>`; invalid input returns an exception. |
| Arrays | Fixed Rust arrays become TS tuples, numeric buffers become TypedArrays, readonly containers have explicit readonly declarations plus runtime freezing. |
| Enums | Closed native enums retain their numeric values; serializable tagged unions use Tsify. |
| Callbacks | Explicit function signature on the binding; native COW snapshot passed to JS; callback exceptions propagate unchanged. A non-boolean result is rejected. |
| Errors | Concrete Error classes retain native Display, readonly typed `details` and actual recursively projected `cause`. |
| Large integers | `u64` uses JS `bigint`; it is not converted through `number`. |

```typescript
const params = new RingSearchParams();
params.includeDativeBonds = true;
const molecule = new Molecule('C1CC1');
const count = molecule.symmetrizedSssr(params).numRings;
const nested = new NestedParams({child: params});
nested.child.includeHydrogenBonds = true;
const selection = new Selection({indices: [0, 1]});
selection.indices = [2]; // selection.indices.push(2) is rejected
```

## Limits before full binding adoption

- Tsify handles serializable values, not native class ownership, methods or
  arbitrary Rust callbacks. Those remain explicit bindings. The callback sample
  is synchronous; it does not establish worker/thread or async callback support.
- Associated native-to-JS type adapters live beside their binding declarations,
  never in a framework business map. Remaining native callback/executor shapes
  are implemented as complete behavior groups, not declared supported by counts.
- Host filesystem datasets need an explicit host boundary; no virtual filesystem
  is implied. Binary archives are excluded from WASM by the shared API contract.
- Owned TypedArray output copies data. No zero-copy view across memory growth
  or release is promised. Finalization is nondeterministic; explicit `free()`
  and generated disposal remain available.

Tests contain seven protocol groups and ten negative typing witnesses, not a
business-interface inventory. New bindings reusing a proven protocol must not
expand these consumers. This prototype proves the listed representations work;
it does not yet freeze every protocol needed by the full Python inventory.

Upstream: [Tsify](https://github.com/madonoharu/tsify),
[fallible Ts conversion](https://docs.rs/tsify/0.5.8/tsify/struct.Ts.html).
