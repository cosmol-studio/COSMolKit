# Python Bindings

This directory hosts PyO3-based Python bindings for COSMolKit.

Use the repository-level Python environment:
- `uv sync --group dev`
- `.venv/bin/python -m pip install -e python`

Current layout:
- `pyproject.toml` package metadata for the published Python distribution
- `Cargo.toml` and `src/lib.rs` for the Rust extension module crate built with PyO3 + maturin
- no Python shim package; `cosmolkit` is provided directly by the compiled extension module

Current status:
- the package is still partial and not yet a full public API surface
- low-level bindings now expose:
  - `placeholder()`, `rust_version()`, and `core_version()`
  - `Molecule.from_smiles()`, `Molecule.read_sdf()`, `atoms()`, `bonds()`, and `find_chiral_centers()`
  - `tetrahedral_stereo_from_smiles()` for the ordered-ligand tetrahedral stereo representation
- `Molecule.ensure_conformer()` and the higher-level 3D pipeline remain intentionally unimplemented on the COSMolKit side

Examples:
- `python/examples/tetrahedral_stereo.py` shows the internal tetrahedral stereo representation exposed to Python
- `python/examples/io_and_properties.py` and the other example files still describe the intended longer-term API shape

Local development:
- `uv sync --group dev`
- `.venv/bin/python -m pip install -e python`
- `.venv/bin/python -c "import cosmolkit; print(cosmolkit.placeholder())"`

Publishing:
- GitHub Actions workflow: `.github/workflows/python-publish.yml`
- expected PyPI secret: `PYPI_API_TOKEN`
- trigger: see the workflow file for the current tag/release policy
