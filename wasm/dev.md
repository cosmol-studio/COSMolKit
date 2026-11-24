# WASM development

Run every command below from the repository root. Requires Rust,
Python 3.11+ and Node.js/npm.

We build one ES Module package with `--target web`, usable in both browsers
and Node. Importing the npm entry initializes WASM automatically. wasm-pack
generates the bindings unchanged; the small `js/index.js` entry loads them.
Release builds also run wasm-opt.

## Install build tools

```bash
cargo install wasm-pack --version 0.15.1 --locked
rustup target add wasm32-unknown-unknown
```

## Build the package

This builds the full package into `target/npm-package`. Use a fresh output
directory for a publication build. `--out-dir` is relative to `wasm/`.

```bash
wasm-pack build wasm --release --target web --out-dir ../target/npm-package \
    -- --locked --no-default-features --features full
```

Copy the automatic entry and license, then set the npm entry, name and label:

```bash
cp LICENSE target/npm-package/LICENSE
cp wasm/js/index.js wasm/js/index.d.ts target/npm-package/
npm pkg set --prefix target/npm-package \
    main=index.js module=index.js types=index.d.ts \
    name=@cosmol-studio/cosmolkit publishConfig.access=public \
    "publishConfig.tag=$(python3 wasm/presets.py tag full)"
npm pkg set --prefix target/npm-package --json \
    'files=["*.js","*.d.ts","*.wasm","snippets"]' \
    'sideEffects=["./index.js","./snippets/*"]'
```

## Create the npm archive

The first command lists the files; the second creates the `.tgz`. Neither
uploads anything.

```bash
npm pack --dry-run --json ./target/npm-package
npm pack ./target/npm-package
```

## Use the same package in Node

After installing the npm package, save this as `example.mjs` and run
`node example.mjs`. The same import works in browsers. No init call is needed:
the entry reads the WASM file in Node and loads its URL in browsers.

```javascript
import { molFromSmiles } from '@cosmol-studio/cosmolkit';

const mol = molFromSmiles('CCO');
console.log(mol.toSmiles());
mol.free();
```

## Run binding checks

Build a development copy in `wasm/pkg` with registry inspection and the few
framework test samples enabled. It is still the same `web` format, not a
Node-specific build. These check-only features must not be published.

```bash
python3 wasm/features.py
wasm-pack build wasm --dev --target web --out-dir pkg \
    -- --locked --no-default-features --features full,registry-check,protocol-tests
cp wasm/js/index.js wasm/js/index.d.ts wasm/pkg/
npm pkg set --prefix wasm/pkg main=index.js module=index.js types=index.d.ts
npm pkg set --prefix wasm/pkg --json private=true \
    'files=["*.js","*.d.ts","*.wasm","snippets"]' \
    'sideEffects=["./index.js","./snippets/*"]'
```

Run the checks separately; none recompiles Rust:

```bash
# Check exports and declarations against the compiled Rust registry.
python3 wasm/registry_check.py

# Run the framework samples in Node, loading the web package.
node wasm/tests/runtime.mjs

# Check TypeScript declarations and the typed samples.
npx --yes --package typescript@5.8.3 tsc -p wasm/tests/tsconfig.json
```

The registry report is `wasm/target/registry.json`.

## Build a smaller package

For example, `core-search` uses `core,search` instead of `full`:

```bash
wasm-pack build wasm --release --target web --out-dir ../target/npm-core-search \
    -- --locked --no-default-features --features core,search
cp LICENSE target/npm-core-search/LICENSE
cp wasm/js/index.js wasm/js/index.d.ts target/npm-core-search/
npm pkg set --prefix target/npm-core-search \
    main=index.js module=index.js types=index.d.ts \
    "name=$(python3 wasm/presets.py name core-search)" publishConfig.access=public \
    "publishConfig.tag=$(python3 wasm/presets.py tag core-search)"
npm pkg set --prefix target/npm-core-search --json \
    'files=["*.js","*.d.ts","*.wasm","snippets"]' \
    'sideEffects=["./index.js","./snippets/*"]'
npm pack --dry-run --json ./target/npm-core-search
```

`presets.py` only prints `features`, `name` or `tag` for a preset; it
never builds, checks or uploads. `--out-dir` is relative to `wasm/`;
absolute paths also work. `CARGO_TARGET_DIR` selects the compilation cache.
The standalone Cargo workspace retains its own lockfile.

After changing the Rust feature graph:

```bash
python3 wasm/features.py --write
```

## Test published prereleases

```bash
npm install @cosmol-studio/cosmolkit@rc
npm install @cosmol-studio/cosmolkit@dev
```

Choose one channel or an exact published version. Smaller packages use the
same labels. Follow [BINDING_CONTRACT.md](BINDING_CONTRACT.md) for binding
semantics; chemistry remains in the Rust facade.
