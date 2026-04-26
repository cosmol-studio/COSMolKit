# Python Bindings

This directory hosts PyO3-based Python bindings for COSMolKit.

Use the repository-level Python environment:
- `uv sync --group dev`
- `.venv/bin/python -m pip install -e python`

Current layout:
- `pyproject.toml` package metadata for the published Python distribution
- `native/` Rust extension module crate built with PyO3 + maturin
- no Python shim package; `cosmolkit` is provided directly by the compiled extension module

Current status:
- this is a placeholder package scaffold
- it currently only binds `cosmolkit-core` version metadata (no functional API binding yet)
- the extension currently exposes `placeholder()`, `rust_version()`, and `core_version()`

Local development:
- `uv sync --group dev`
- `.venv/bin/python -m pip install -e python`
- `.venv/bin/python -c "import cosmolkit; print(cosmolkit.placeholder())"`

Publishing:
- GitHub Actions workflow: `.github/workflows/python-publish.yml`
- expected PyPI secret: `PYPI_API_TOKEN`
- trigger: push a tag matching `python-v*`
