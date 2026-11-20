# COSMolKit

Development follows [dev/README.md](dev/README.md), [crate architecture](dev/crate_architecture.md) and [public API design](dev/public_api_design.md).

For upstream code ports, follow [dev/source_reproduction_protocol.md](dev/source_reproduction_protocol.md).

Build instructions: [Python](python/dev.md), [WASM](wasm/dev.md), [documentation](docs/README.md) and [parity tests](parity-tests/README.md).

## Project Layout

- `crates/cosmolkit/`: public Rust API, live Molecule and operation runtime.
- `crates/cosmolkit-types/`, `crates/cosmolkit-model/`: shared types and detached data.
- `crates/cosmolkit-core/` and domain crates: algorithms.
- `crates/cosmolkit-macros/`: declaration parsing and code generation.
- `python/`, `wasm/`: language bindings.
- `parity-tests/`: corpus and special-regression comparisons.
- `testdata/`: shared fixed fixtures and provenance.
- `docs/`, `docs-web/`: documentation and website.
- `dev/`: development standards, designs and plans.
