# Conda distribution

The `publish.yml` workflow publishes to the public Anaconda.org channel
`cosmol-studio`.

## One-time setup

1. Create or obtain upload access to the `cosmol-studio` organization on
   Anaconda.org, independently of the GitHub organization.
2. Create an Anaconda.org upload token with access to that organization.
3. Save it as the GitHub repository Actions secret `ANACONDA_API_TOKEN`.

Do not put tokens in workflow inputs or source control.

## Release flow

The `publish` job packages the wheel artifacts with `package_wheels.py` and
uploads the resulting `.conda` files. Stable versions go to the `main` label;
release candidates go to `rc`, and development builds go to `dev`.

Each of Linux x86_64, Linux ARM64, Windows x86_64, macOS Intel, and macOS ARM64
gets CPython 3.9–3.14 packages. CPython-specific wheels are used when present;
the ABI3 wheel supplies other supported versions. macOS universal2 wheels
are packaged separately for each conda architecture.

Packages retain the wheel's native libraries, `__init__.pyi`, and `py.typed`.
Conda supplies Python, NumPy, and platform runtime dependencies. The package
version follows the Python distribution's version.

## Install

```bash
conda install -c conda-forge cosmol-studio::cosmolkit
# Or use mamba with the same arguments.
```

For an explicitly selected release candidate:

```bash
conda install -c conda-forge -c cosmol-studio/label/rc cosmolkit
# Development builds:
conda install -c conda-forge -c cosmol-studio/label/dev cosmolkit
```
