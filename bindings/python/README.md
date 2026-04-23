# Python Bindings

This directory hosts PyO3-based Python bindings for COSMolKit.

Use the repository-level Python environment:
- `uv sync --group dev`
- `.venv/bin/maturin develop --manifest-path bindings/python/native/Cargo.toml`

Planned layout:
- `pyproject.toml` for packaging
- Rust extension module crate (e.g., `bindings/python/native/`)
- Python shim package in `cosmolkit/`
