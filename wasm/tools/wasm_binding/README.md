# WASM Development And Build

`wasm/` is the binding-facing crate. It re-exports the Rust `cosmolkit` facade
and provides the ABI-safe `Molecule` projection used by Alef and
`wasm-bindgen`. Chemistry remains in the existing Rust crates; this directory
owns only the language-boundary shape and its checks.

## Prerequisites

Install the Rust target and the tools used by the binding check:

```bash
rustup target add wasm32-unknown-unknown
# Alef 0.83.3 and wasm-bindgen-cli 0.2.128 must be available on PATH.
# Bun is preferred for the runtime check; Node is accepted as a fallback.
```

The checked-in Alef input is [alef.toml](./alef.toml). It intentionally lives
under this tool directory. Do not run generation from the repository root:
Alef can create toolchain and provenance files beside the input config.

## Local Checks

These commands validate the source crate without generating JavaScript:

```bash
cargo test -p cosmolkit-wasm --release
cargo check -p cosmolkit-wasm --target wasm32-unknown-unknown
```

The complete boundary check generates an isolated crate, compiles it, runs
`wasm-bindgen`, executes the JavaScript API with Bun or Node, and type-checks
the generated declarations with TypeScript `strict` mode:

```bash
ALEF_BIN=/path/to/alef \
WASM_BINDGEN_BIN=/path/to/wasm-bindgen \
python3 wasm/tools/wasm_binding/run.py
```

When Alef and wasm-bindgen are already on `PATH`, the short form is:

```bash
python3 -B wasm/tools/wasm_binding/run.py
```

Set `BUN_BIN` or `NODE_BIN` to select a specific runtime. The runner uses
the pinned `typescript@5.8.3` package through `npx`; set `TSC_BIN` only
to explicitly test another compiler.

All generated crates, WASM binaries, JavaScript glue, declarations, and
TypeScript configuration are created below a temporary directory and removed
on exit. To retain the tested npm package, supply a new output directory:

```bash
python3 -B wasm/tools/wasm_binding/run.py --out-dir target/npm/cosmolkit
npm pack --dry-run --json ./target/npm/cosmolkit
```

The default is `--preset full`. Select a fixed smaller distribution with, for
example, `--preset core-bio --out-dir target/npm/core-bio`. Available presets
are `core`, `core-search`, `core-fingerprints`, `core-analysis`, `core-reaction`, `core-depict`, `core-3d`,
`core-bio`, `core-inchi`, and `full`. Each compiles without dependency defaults,
selects only its binding modules and applicable tests, and exports its own
version and npm tag. Domain-specific tests retain their assertions; mixed-domain
suites run when all their prerequisites are selected.

Export happens only after runtime and TypeScript checks pass. The package
includes its JavaScript snippets, WASM binary, declarations, README and license.

The Rust facade manifest is the only feature dependency tree. After changing
it, run `python3 wasm/tools/wasm_binding/features.py --write` to regenerate
`wasm/Cargo.toml`; the runner and unit tests reject a stale projection.
Binary serialization is excluded on WASM. Presets select feature names only;
they do not maintain a second set of prerequisite edges.
Generated files are build output, not committed source.

## Release Build Shape

The release package follows the same isolated sequence:

1. Alef extracts `wasm/src/lib.rs` into a temporary binding crate.
2. Cargo builds that crate for `wasm32-unknown-unknown` in release mode.
3. `wasm-bindgen --target web` emits the JavaScript module, declarations, and
   background `.wasm` file.
4. Runtime and TypeScript checks validate the generated package before export.
5. [publish.yml](../../../.github/workflows/publish.yml) publishes that artifact
   as `@cosmol-studio/cosmolkit`, using npm Trusted Publisher (OIDC). Each preset
   uses its preset name as dist-tag and a distinct artifact version; full RCs
   use `rc` and full stable releases use `latest`. Manual publishing is RC-only; stable releases
   require a matching version tag.

The committed Rust surface is the source of truth. Generated package files are
build output and must stay outside the repository tree or under `tmp/`.
