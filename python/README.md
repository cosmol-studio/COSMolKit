# COSMolKit

COSMolKit is a Rust-backed Python package for cheminformatics and molecular
graph workflows. The Python package currently exposes a strict, tested subset
of the Rust core while the chemistry implementation is expanded against RDKit
parity tests.

The public Python API is intentionally non-inplace by default: ordinary
operations return a new `Molecule` value instead of mutating the input object
implicitly. The intended optimization strategy for this API shape is Rust-side
copy-on-write storage, not ambiguous Python-layer mutation semantics.

## Installation

```bash
pip install cosmolkit
```

## Currently Available API

The current binding surface is useful for:

- SMILES parsing with `Molecule.from_smiles()`
- optional RDKit bridge with `Molecule.from_rdkit()` when `rdkit` is installed
- SDF loading with `Molecule.read_sdf()`
- atom and bond graph inspection with RDKit-like feature names
- chiral center inspection with `Molecule.find_chiral_centers()`
- hydrogen add/remove transforms
- Kekulization
- tetrahedral stereo inspection
- 2D coordinate generation and coordinate access
- molecule length / basic repr helpers
- SMILES export
- distance-geometry bounds matrix export
- SVG export with `Molecule.to_svg()`
- SVG file export with `Molecule.write_svg()`
- PNG file export with `Molecule.write_png()`
- SDF read/write and SDF string export
- SDF directory export with `Molecule.write_sdf_to_directory()`
- explicit edit context with `Molecule.edit()`
- module version helpers via `cosmolkit.rust_version()` and `cosmolkit.core_version()`

## Quick Start

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO")

print(mol)
print("atoms:", len(mol))

for atom in mol.atoms():
    print(atom.idx(), atom.atomic_num(), atom.total_valence())

for bond in mol.bonds():
    print(bond.begin_atom_idx(), bond.end_atom_idx(), bond.bond_type())
```

## Hydrogen Transforms and Kekulization

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1O")
mol_h = mol.add_hydrogens()
mol_no_h = mol_h.remove_hydrogens()
kek = mol.kekulize()

print(len(mol.atoms()), len(mol_h.atoms()), len(mol_no_h.atoms()))
print(kek.to_smiles())
```

## Tetrahedral Stereo and Chiral Tags

COSMolKit exposes an ordered tetrahedral stereo representation derived from the
internal chiral tags. Each record is `(center_atom_index, ordered_ligands)`.
Implicit hydrogen is represented as `None`.

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("F[C@H](Cl)Br")

print(mol.find_chiral_centers(include_unassigned=False))

for center, ligands in mol.tetrahedral_stereo():
    print("center:", center)
    print("ordered ligands:", ligands)
```

## 2D Coordinates and SDF I/O

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO").compute_2d_coords()

print("num conformers:", mol.num_conformers())
print("first position:", mol.conformer_positions()[0])

sdf_text = mol.to_sdf_string(format="v2000")
print("SDF length:", len(sdf_text))

mol.write_sdf("ethanol.sdf", format="v2000")
loaded = Molecule.read_sdf("ethanol.sdf")
print(loaded)
```

## SVG Rendering

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1O").compute_2d_coords()
svg = mol.to_svg(width=400, height=300)
mol.write_svg("phenol.svg", width=400, height=300)
mol.write_png("phenol.png", width=400, height=300)
```

`Molecule.to_svg()` returns the SVG string directly. `Molecule.write_svg()`
and `Molecule.write_png()` perform the file IO inside the Rust binding.
See `python/examples/image_exports.py` for a runnable file-export example.

## SMILES Export and DG Bounds

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("F[C@H](Cl)Br")

print("isomeric:", mol.to_smiles())
print("non-isomeric:", mol.to_smiles(False))

bounds = mol.dg_bounds_matrix()
print("matrix size:", len(bounds), "x", len(bounds[0]))
print("d(0,1):", bounds[0][1])
```

## Explicit Editing

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO")

editor = mol.edit()
cl = editor.add_atom("Cl")
editor.add_bond(0, cl, order="single")
editor.set_atom_charge(cl, -1)
mol2 = editor.commit()

print(mol2.to_smiles())
```

## RDKit Bridge

`Molecule.from_rdkit()` is available when `rdkit` is installed in the active
Python environment. The bridge copies basic graph information from an RDKit
molecule into COSMolKit and then uses COSMolKit logic for downstream derived
features.

```python
from rdkit import Chem
from cosmolkit import Molecule

rd_mol = Chem.MolFromSmiles("CCO")
mol = Molecule.from_rdkit(rd_mol)

print(mol.to_smiles())
```

## Planned APIs

The following API directions are planned but not currently available:

- sanitization pipeline methods such as `sanitize()`
- explicit ring/aromaticity perception methods
- substructure search
- fingerprint generation
- 3D embedding and force-field optimization
- alignment result APIs beyond the current minimal surface

## Development Status

COSMolKit is not a full RDKit replacement today. The Python package is being
built through strict parity testing against RDKit behavior where compatibility
is the goal. The currently usable binding surface is centered on deterministic
non-inplace molecule transforms, graph/stereo inspection, 2D coordinates,
SMILES/SDF IO, DG bounds export, SVG rendering, and explicit editing.
