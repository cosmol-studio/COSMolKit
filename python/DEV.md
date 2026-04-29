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
