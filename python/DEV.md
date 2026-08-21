# Python Dev Commands

This document lists the canonical commands for Python binding development in this repository.

## Environment

From repo root:

```bash
uv sync --group dev
```

## Generate `.pyi` Stubs (dev / abi3-py310)

From repo root:

```bash
cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen
```

Generated file:

```text
python/cosmolkit.pyi
```

## Build/Install Extension with maturin

### Dev install (editable)

From repo root:

```bash
.venv/bin/maturin develop --manifest-path python/Cargo.toml
```

### Release wheel build (abi3-py39)

From repo root:

```bash
.venv/bin/maturin build --release --manifest-path python/Cargo.toml --features release-abi3-py39 --out python/dist
```

## Build Python Documentation

From repo root:

```bash
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
rm -rf python/docs/build/html
.venv/bin/python -m sphinx -W --keep-going -E -b html python/docs/source python/docs/build/html
.venv/bin/python python/docs/check_seo.py python/docs/build/html
```

Generated HTML:

```text
python/docs/build/html/index.html
```

The strict build treats documentation warnings as failures. The final command
also verifies every public page's title, description, and canonical URL, plus
the generated sitemap, homepage structured data, shared project links on every
generated page, and copied ``robots.txt``.

## Type Checking

```bash
.venv/bin/basedpyright python/tests python/examples
```
