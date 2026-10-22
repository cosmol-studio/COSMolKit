# Python bindings

Bindings call the public Rust facade. PyO3 supplies types; the generator handles
declared overloads, indexing and NumPy annotations. See the
[binding contract](BINDING_CONTRACT.md) for implementation rules and examples.

Run all commands from the repository root.

## Typed atom properties

```python
tagged = mol.with_atom_property_int(0, "tracking_id", 42)
assert tagged.atom_property_int(0, "tracking_id") == 42
tagged.set_atom_property_string_(0, "label", "CA")
assert tagged.atom_property_string(0, "label") == "CA"
```

Use `string`, `int`, `uint`, `double`, `bool`, `int_vector` or `string_vector`
after `with_atom_property_`, `set_atom_property_` (with a final `_`) or
`atom_property_`. Reads return `None` only for missing keys; wrong types, invalid
atom indices and invalid UTF-8 raise `AtomPropertyReadError`. No numeric or
text coercion is performed. Generic enum-valued atom access is Rust-only.

## Embedding results

Explicit 3D coordinate inputs use Rust `Vec<[f64; 3]>`. Python lists, tuples
and NumPy rows are extracted directly by PyO3, which requires exactly three
components per row. Rust checks atom count and finite values. The official
stub describes a length-bearing sequence; it does not encode the length 3.

```python
result = mol.with_3d_conformers_result(8)
for conformer_id in result.conf_ids:  # successful IDs only, never -1
    print(conformer_id)
print(result.generated_count, result.termination)

single = mol.with_3d_conformer_result()
if single.conf_id is not None:
    print(single.conf_id)
```

`EmbedTermination` is `Completed`, `TimedOut`, or `Interrupted`, even when failure
tracking is disabled. `Completed` means normal completion, not that all requested
conformers survived embedding and pruning. Counts describe this call, not the
molecule's pre-existing conformers. Invalid inputs still raise their Rust error.

## Ring queries

```python
rings = mol.symmetrized_sssr()
rings = mol.symmetrized_sssr(include_dative_bonds=True)
print(rings.num_rings, rings.atom_rings, rings.bond_rings)
```

The result owns the symmetrized SSSR rows. The query does not modify the molecule
or its ring cache. It does not enumerate every possible cycle.

## Prepare

Install development tools and the extractor matching PyO3:

```bash
uv sync --group dev
cargo install pyo3-introspection --version 0.30.1 --locked --root target/pyo3-tools
```

## Build

```bash
.venv/bin/python python_pyo3_introspection/build.py \
    --extractor target/pyo3-tools/bin/pyo3-introspection
```

This command compiles the Rust extension, extracts types, generates the Python
package and `.pyi`, then runs the runtime samples and static type checks.

The package is written to
`python_pyo3_introspection/artifacts/smoke/package/cosmolkit/`.

## Check an existing build

Run runtime samples and static type checks without rebuilding:

```bash
.venv/bin/python python_pyo3_introspection/tests/runtime.py \
    python_pyo3_introspection/artifacts/smoke/package
.venv/bin/python -m basedpyright --pythonpath .venv/bin/python \
    --project python_pyo3_introspection/artifacts/smoke/typing.json
```

Use these commands when only consumer samples change. After changing bindings
or generation rules, run `build.py` again.

## Check registration

```bash
.venv/bin/python python_pyo3_introspection/registry_check.py
```

Compares the canonical registry with the built package and generated types,
including configuration fields and errors. Writes `artifacts/registry.json`.

## Build and install a wheel

After `build.py`, package the extension, `__init__.py`, `__init__.pyi` and `py.typed`:

```bash
.venv/bin/maturin build --manifest-path python_pyo3_introspection/Cargo.toml \
    --target-dir target/typed_facade_projection --features extension \
    --interpreter .venv/bin/python \
    --out python_pyo3_introspection/artifacts/wheels
uv pip install --python .venv/bin/python --reinstall \
    python_pyo3_introspection/artifacts/wheels/cosmolkit-*.whl
.venv/bin/python -c 'import cosmolkit as ck; print(ck.mol_from_smiles("CCO").num_atoms)'
```

The final command prints `3`. Add `--release` to `maturin build` for an optimized
wheel. Distribution name: `cosmolkit`; import name: `cosmolkit`.

Molecules support `pickle.dumps(mol)` and `pickle.loads(data)` through the stable
`cosmolkit._molecule_from_binary` reconstruction function and native binary format.

## Extending bindings

Configurations declare fields only in their Rust owner, using
`#[derive(cosmolkit_macros::Configuration)]`. Re-export its generated
hidden metadata alongside the native type; the facade also re-exports its
`__ck_configuration_TypeName` macro. The Python declaration names only that type:

```rust
ck::project_configuration!(ck::MaccsFingerprintParams);
```

The shared projection generates typed construction and writable properties over
the actual Rust configuration. Getters return tuples, readonly mappings or
nested configuration views. Child edits update the parent; containers are
replaced through whole-field assignment.
Conversion failures leave the Rust field unchanged, including fixed-array
length errors. No editable Python container mirror is retained.

```python
p = ck.MorganCallParams(from_atoms=[0])
p.from_atoms = [0, 2]          # Assigns the actual Rust field.
assert p.from_atoms == (0, 2)
# p.from_atoms.append(2)      # AttributeError: tuple has no append.
# p.from_atoms = ["wrong"]   # TypeError; the field stays (0, 2).
```

Native annotations specify construction and composition: `required` keeps
mandatory inputs, `none_default` accepts None as the native default, `flatten`
exposes explicitly composed options, and `factory` forwards the registered
`new()` method. Every ordinary parameter wrapper retains its native record.

Callbacks, random sources and private execution state use an explicit adapter:

```rust
ck::project_configuration!(ck::SubstructMatchParams, using match_parameters);
```

Their native declarations mark `state` and `adapter` fields; ordinary fields
still come from native metadata and assign the actual native fields.
Adapters retain actual native ownership, initialization and invocation errors.
Unmarked private fields, generic records and unsupported field types fail
compilation, rather than being inferred or downgraded.

Ordinary new APIs must not require generator, runner or consumer-sample changes.
`tests/runtime.py` and `tests/positive.py` contain only a few fixed framework
examples; they do not enumerate business APIs. New projection protocols require
an explicit design decision.
