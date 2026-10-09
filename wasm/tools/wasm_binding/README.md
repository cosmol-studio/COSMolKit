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
package name. All packages share the Rust version and use `rc` or `latest` as
their release channel. Domain-specific tests retain their assertions; mixed-domain
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
2. Cargo builds that crate for `wasm32-unknown-unknown` with explicit release
   settings: `opt-level=3`, fat LTO, and `codegen-units=1`.
3. `wasm-bindgen --target web` emits the JavaScript module, declarations, and
   background `.wasm` file.
4. Runtime and TypeScript checks validate the package before export.
5. [publish.yml](../../../.github/workflows/publish.yml) publishes the full artifact
   as `@cosmol-studio/cosmolkit` and each smaller preset as
   `@cosmol-studio/cosmolkit-{preset}`, using npm Trusted Publisher (OIDC).
   Every package uses the same release version: RCs use `rc`, stable releases
   use `latest`. Manual CI publishing is RC-only; stable releases
   require a matching version tag.

The committed Rust surface is the source of truth. Generated package files are
build output and must stay outside the repository tree or under `target/`.

## First Publication

Build and validate each preset serially with a shared Cargo cache. Use a new
output directory for every build; the runner does not overwrite packages:

```bash
export CARGO_TARGET_DIR="$PWD/target/wasm"
mkdir -p target/npm-tarballs
for preset in core core-search core-fingerprints core-analysis core-reaction core-depict core-3d core-bio core-inchi full; do
    python3 -B wasm/tools/wasm_binding/run.py --preset "$preset" --out-dir "target/npm/$preset"
    npm pack "./target/npm/$preset" --pack-destination target/npm-tarballs
done
```

Check each tarball contains JavaScript, declarations, the `.wasm` binary and
snippets, not just `README.md` and `package.json`. The package manifest records
its preset, name, version and publication channel.

For the first RC publication, the package owner publishes the built tarballs
manually (authenticate with `npm login` first):

```bash
for package in target/npm-tarballs/*.tgz; do
    npm publish "$package" --access public --tag rc
done
```

Then configure **each package's** Settings → Trusted Publisher for GitHub
Actions: owner `cosmol-studio`, repository `COSMolKit`, workflow filename
`publish.yml`, no environment (unless one is also configured on the workflow
job). Allow direct `npm publish` if the form offers an action choice. Use npm
11.5.1 or newer and Node 22.14 or newer; CI uses Node 24 and already grants
`id-token: write`. No `NPM_TOKEN` is needed for these subsequent OIDC publishes.
See [npm's trusted-publisher documentation](https://docs.npmjs.com/trusted-publishers/).

A published name/version cannot be reused. After manually publishing this RC,
CI must publish a later version, not attempt the same version again.
