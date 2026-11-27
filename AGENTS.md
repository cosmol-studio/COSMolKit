# COSMolKit

Development follows [dev/README.md](dev/README.md), [crate architecture](dev/crate_architecture.md) and [public API design](dev/public_api_design.md).

For upstream code ports, follow [dev/source_reproduction_protocol.md](dev/source_reproduction_protocol.md).

Build instructions: [Python](python/dev.md), [WASM](wasm/dev.md), [documentation](docs/README.md) and [parity tests](parity-tests/README.md).

## Extending the Project

Agents tend to keep expanding a project to finish each local feature: add a CI
step, invent a test workflow, patch a generator, then repeat the same checks in
every language. Each addition may look useful alone, but together they make
every new API increase the project's maintenance burden. COSMolKit controls
this growth through the following division of responsibilities:

- Rust owns chemistry and behavior. Implement algorithms in their domain crate;
  expose them through the public facade, operation contracts and canonical
  registry. The Rust runtime enforces COW, state, cache and commit guarantees.
  Test those guarantees in runtime tests; APIs using that runtime inherit them.
  Individual operation tests cover behavior specific to that operation instead
  of repeating the runtime's COW, commit or cache tests.
- Python and JavaScript bind the same Rust API. Except for language-specific
  needs such as type conversion, indexing and callbacks, they directly forward
  Rust behavior. Shared projections generate declarations from actual binding
  types, so ordinary API additions reuse the framework and generation scripts.
- Tests follow the behavior's owner. Small Rust regressions cover local
  algorithms, boundaries and runtime behavior. Chemical consistency with the
  pinned upstream belongs in `parity-tests/`: ordinary corpus tasks exercise
  operations across a corpus; special regressions compare explicitly selected
  fixed reference cases. Both use the existing prepare/compare workflow.
- Python/JS tests use a few small samples to establish binding usability and
  language protocols. When adding a new type converter, indexing protocol or
  callback adapter, test it with representative samples once; additional APIs
  using the same binding code inherit its guarantees.
  For example, Tsify owns TypeScript type projection, shared input adapters own
  runtime type extraction, and shared error projections own exception conversion.
  Test these integration boundaries with representative samples; reuse their
  guarantees for each new API instead of repeatedly testing whether every
  projected field has a valid type or every argument rejects the same invalid
  input. Chemical consistency stays in Rust/parity.
- New features add implementations and registrations within this structure.
  CI discovers registered tasks through existing entrypoints, keeping the
  workflow independent of individual chemical features. This is why adding an
  API should leave mandatory commands and generic infrastructure stable.

When the existing runtime, binding framework or parity runner cannot support a
required capability, explain the
concrete gap and obtain maintainer approval for a shared solution first. Changes
to reference versions or numerical acceptance rules likewise require an explicit
decision. Follow [test boundaries](dev/test_boundaries.md) for test placement.

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
