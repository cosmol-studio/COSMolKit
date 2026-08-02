# InChI Development Tools

These tools inspect the configured official `libinchi` source target or
regenerate owned Rust source types. They are development and oracle tooling;
they are not production dependencies.

Public entrypoints:

```bash
.venv/bin/python dev/tools/inchi/audit_inchi_active_call_graph.py
uv run --no-project --with tree-sitter --with tree-sitter-c \
  python dev/tools/inchi/generate_inchi_function_inventory.py
.venv/bin/python dev/tools/inchi/generate_inchi_source_types.py
```

The source-type generator crate is in `source_type_generator/`. Scripts that
rewrote the archived full-port plan are preserved under
[`../../archive/tools/inchi_plan_generation/`](../../archive/tools/inchi_plan_generation/).
