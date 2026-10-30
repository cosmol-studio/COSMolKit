# WASM development

## Test published prereleases

Dev and rc builds are for development testing:

```bash
npm install @cosmol-studio/cosmolkit@rc
npm install @cosmol-studio/cosmolkit@dev
```

Select one channel, or pin an exact published version for reproducibility.
The same tags apply to the smaller feature packages.

## Build from source

Requires Rust with `wasm32-unknown-unknown`, Python 3.11+, Node.js/npm and
`wasm-bindgen-cli 0.2.128`.

```bash
# Build a Node package, check the compiled registry and fixed framework samples.
.venv/bin/python wasm/build.py --check

# Build a publishable browser package, including generated declarations.
.venv/bin/python wasm/build.py --release --target web --out-dir target/npm-package
npm pack ./target/npm-package

# Smaller package; presets use the facade's feature combinations.
.venv/bin/python wasm/build.py --preset core-search --release --target web
```

The build compiles Rust, then wasm-bindgen generates JavaScript and `.d.ts`.
No declaration rewriting is performed. `--check` includes inspection metadata;
full check builds additionally include the fixed protocol samples. Such packages
are private and are not release artifacts. Smaller presets check their own
compiled registry. Production builds omit both inspection and protocol samples.

```bash
# Inspect an already-built check package without rebuilding.
.venv/bin/python wasm/registry_check.py

# After changing the Rust facade feature graph:
.venv/bin/python wasm/features.py --write
```

The registry report is `wasm/target/registry.json`. It checks registered names
and runtime/declaration shapes, not chemical correctness. Tests cover shared
binding protocols with a few representative values; new APIs using those
protocols do not require a second consumer inventory.

The standalone Cargo workspace keeps the WASM toolchain lock separate.
`CARGO_TARGET_DIR` can select its build cache. Publication uses the existing
npm package names and `dev`/`rc`/`latest` labels. The retained `wasm_old/` directory
is not used by builds, documentation or publication.

Follow [BINDING_CONTRACT.md](BINDING_CONTRACT.md) for input validation, ownership,
configuration, callbacks and errors. Only binding declarations and shared
adapters belong here; chemistry remains in the Rust facade.
