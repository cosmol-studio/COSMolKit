# COSMolKit Python

COSMolKit is a Python package for molecule graph workflows, SMILES/SDF IO,
coordinate access, molecule depiction, and high-throughput batch processing.

## API Model: Copy-On-Write (COW) Molecule Values

COSMolKit's Python API uses copy-on-write (COW) value semantics for molecule
transforms. Methods such as `with_hydrogens()`, `without_hydrogens()`,
`with_kekulized_bonds()`, and `with_2d_coords()` return a new `Molecule`; they
do not mutate the original object in place.

This is intentionally different from common RDKit Python workflows, where many
operations mutate an `RWMol` or update a molecule object directly. In COSMolKit,
keep the returned value:

```python
mol = Molecule.from_smiles("CCO")
mol_h = mol.with_hydrogens()

assert mol is not mol_h
print(mol.to_smiles())
print(mol_h.to_smiles())
```

## Installation

```bash
pip install cosmolkit
```

## Quick Start

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1O")
drawn = mol.with_2d_coords()

print(mol.to_smiles())
print(drawn.atoms()[0])

drawn.write_png("phenol.png", width=400, height=300)
```

Chiral tags are available directly on atoms, so code that works with the
SMILES/RDKit-style CW and CCW path does not need to switch to the ordered
tetrahedral view:

```python
chiral = Molecule.from_smiles("F[C@H](Cl)Br")

print(chiral.to_smiles())
print(chiral.to_smiles(isomeric_smiles=False))

for atom in chiral.atoms():
    if atom.chiral_tag() != "CHI_UNSPECIFIED":
        print(atom.idx(), atom.chiral_tag())
```

Use `Molecule.edit()` when you want an explicit editing workflow:

```python
editor = mol.edit()
cl = editor.add_atom("Cl")
editor.add_bond(0, cl, order="single")
mol2 = editor.commit()
```

## Batch Workflows

```python
from cosmolkit import MoleculeBatch

smiles = ["CCO", "c1ccccc1", "not-smiles"]
batch = MoleculeBatch.from_smiles_list(smiles, errors="keep")

prepared = batch.add_hydrogens(errors="keep").compute_2d_coords(errors="keep")
report = prepared.to_images("molecule_images", format="png", errors="skip")

print(prepared.valid_mask())
print(prepared.errors())
print(report)
```

The `errors` option controls invalid records:

- `errors="raise"` raises on the first batch validation failure.
- `errors="keep"` preserves failed records and exposes structured errors.
- `errors="skip"` omits failed records from the returned result or export.

Batch SMILES output accepts formatting options:

```python
canonical = prepared.to_smiles_list(canonical=True)
explicit = prepared.to_smiles_list(
    all_bonds_explicit=True,
    all_hs_explicit=True,
)
without_maps = prepared.to_smiles_list(ignore_atom_map_numbers=True)
```

## SDF and Arrays

```python
mol = Molecule.from_smiles("CCO").with_2d_coords()

sdf_text = mol.to_sdf_string(format="v2000")
restored = Molecule.read_sdf_record_from_str(sdf_text, coordinate_dim="2d")

coords = restored.coords_2d()
bounds = restored.dg_bounds_matrix()

print(coords.shape)
print(bounds.shape)
```

## Main Features

- SMILES parsing and writing with `Molecule.from_smiles()` and `to_smiles()`
- copy-on-write molecule value semantics for transforms
- SDF file and string IO
- atom and bond feature inspection
- hydrogen add/remove transforms
- Kekule bond representation
- CW/CCW chiral tags, chiral centers, and tetrahedral stereo inspection
- 2D coordinate generation and NumPy coordinate arrays
- distance-geometry bounds matrix export
- SVG and PNG molecule depictions
- ordered batch construction, transformation, filtering, and export
- explicit molecule editing with `Molecule.edit()`
