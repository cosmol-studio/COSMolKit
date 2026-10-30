# COSMolKit for Python

COSMolKit provides Rust-backed molecular IO, fingerprints, ring queries,
conformer generation and structural biology APIs, with Python type annotations
and NumPy integration.

See the [documentation](https://kit.cosmol.org/) for the full API and
[examples/](examples/) for usage examples. Building or extending the bindings?
Read [dev.md](dev.md).

## Quick start

```bash
pip install cosmolkit
```

```python
import cosmolkit as ck

print(ck.__version__)
mol = ck.mol_from_smiles("CCO")
print(mol.num_atoms)  # 3
print(mol.to_smiles())
array = mol.fingerprint_morgan().to_numpy()
```

Non-mutating operations return new values. In-place molecule operations end
in `_`; `mol.copy()` shares storage until a modification requires a copy.

## Typed atom properties

```python
tagged = mol.with_atom_property_int(0, "tracking_id", 42)
assert tagged.atom_property_int(0, "tracking_id") == 42
tagged.set_atom_property_string_(0, "label", "CA")
assert tagged.atom_property_string(0, "label") == "CA"
```

Property APIs distinguish strings, signed and unsigned integers, doubles,
booleans, integer vectors and string vectors. Reads return `None` for a missing
key; invalid indices or incompatible value types raise an error.

## Conformers and rings

```python
result = mol.with_3d_conformers_result(8)
print(result.conf_ids)  # Successfully generated conformer IDs only.
print(result.generated_count, result.termination)

rings = mol.symmetrized_sssr()
print(rings.num_rings, rings.atom_rings, rings.bond_rings)
```

Conformer counts describe this call. `Completed` means the computation finished,
not that every requested conformer survived generation and pruning. Ring queries
do not modify the original molecule.

Molecules also support `pickle.dumps(mol)` and `pickle.loads(data)`.

## Filesystem paths

Paths are passed literally; `~` is not expanded automatically. Expand it
explicitly when needed:

```python
from pathlib import Path

batch = ck.mols_from_smiles_list(["CCO"])
batch.write_images(Path("~/images").expanduser())
```
