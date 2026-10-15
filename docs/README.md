# Documentation

The shared Sphinx project lives here, with Python pages in `source/python/`,
JavaScript / TypeScript pages in `source/javascript/`, and common styles,
templates and API extensions shared by both languages. Generated HTML belongs
in `build/html/`.

From the repository root, prepare the bindings and documentation tools:

```sh
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen
npm ci --prefix docs --no-audit --no-fund
.venv/bin/python wasm/tools/wasm_binding/run.py --preset full --out-dir target/docs-js-package
```

Reuse the built bindings and declarations for documentation-only edits. Build
both language guides and API references together:

```sh
.venv/bin/python -m sphinx -W --keep-going -E -b html docs/source docs/build/html
```

Python API signatures come from the native extension and checked stub;
JavaScript API signatures come from the public TypeScript declarations via
TypeDoc and sphinx-js. Do not maintain handwritten signature inventories.

The [Dioxus site](../docs-web/README.md) consumes this output and preserves
the public `/python`, `/python/api`, `/javascript` and `/javascript/api` routes.
Its route contract is `docs-web/routes.toml`; source paths are not public URLs.
