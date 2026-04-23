# Repository Guidelines

## Project Structure & Module Organization
This repository is currently in bootstrap stage and contains only [README.md](/home/wangjingtong/code/COSMolKit/README.md). As implementation lands, keep a predictable Rust workspace layout:

- `crates/chem-core/` for molecular graph and chemistry perception logic
- `crates/bio-core/` for PDB/mmCIF parsing and structure hierarchy
- `crates/io/` for SDF/MOL/SMILES/PDB/mmCIF readers and writers
- `bindings/python/` for PyO3 bindings and packaging metadata
- `tests/` for integration/regression datasets and parity tests

Prefer small, focused crates over a monolith; keep public APIs in `lib.rs` and internal helpers in private modules.

## Build, Test, and Development Commands
Use standard Rust workflows from repo root once `Cargo.toml` is added:

- `cargo check` validates code quickly without producing release artifacts
- `cargo test` runs unit and integration tests
- `cargo fmt --all` applies formatting
- `cargo clippy --all-targets --all-features -D warnings` enforces lint quality

For Python tooling and bindings, use project-level `uv` env management:

- `uv sync --group dev` creates/updates root `.venv` with RDKit/testing/build tools
- `.venv/bin/maturin develop --manifest-path bindings/python/native/Cargo.toml` installs local extension into the shared env
- `.venv/bin/pytest` runs Python-facing tests

## Coding Style & Naming Conventions
Use Rust 2021 defaults: 4-space indentation, `snake_case` for functions/modules, `CamelCase` for types/traits, and `SCREAMING_SNAKE_CASE` for constants. Keep modules narrowly scoped and avoid large files with mixed responsibilities. Run `cargo fmt` and `cargo clippy` before opening a PR.

## Testing Guidelines
Place unit tests near code (`mod tests`) and integration tests under `tests/`. Name tests by behavior, e.g., `kekulize_handles_fused_aromatics`. Add regression fixtures for chemistry/biostructure edge cases and explicitly compare outputs against RDKit/Biopython where parity is a goal.

## Commit & Pull Request Guidelines
No commit history exists yet; adopt Conventional Commits from day one (e.g., `feat: add sdf streaming parser`, `fix: correct aromatic valence handling`). Keep commits atomic and buildable.

PRs should include: scope summary, validation steps (`cargo test`, fixtures used), linked issue(s), and notes on parity-impacting behavior changes.
