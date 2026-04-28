# COSMolKit

COSMolKit is a Rust-backed Python package for cheminformatics and molecular
graph workflows. The project is currently in early development: the Python
package exposes a small, strict subset of the Rust core while the chemistry
implementation is expanded against RDKit parity tests.

The current package is useful for basic SMILES parsing, molecule graph
inspection, hydrogen expansion, Kekulization, and tetrahedral stereo
experiments. Higher-level workflows such as SDF IO, conformer generation,
substructure search, fingerprints, editing, and alignment are part of the
planned public API but are not implemented yet.

## Installation

```bash
pip install cosmolkit
```

## Quick Start

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO")

print("atoms")
for atom in mol.atoms():
    print(atom.idx(), atom.atomic_num())

print("bonds")
for bond in mol.bonds():
    print(bond.begin_atom_idx(), bond.end_atom_idx(), bond.bond_type())
```

## Add Explicit Hydrogens

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO")
mol_h = mol.add_hydrogens()

print(len(mol.atoms()))
print(len(mol_h.atoms()))
```

```python
# Not implemented yet
mol_no_h = mol_h.remove_hydrogens()
```

## Kekulization

`Molecule.kekulize()` is implemented for the currently supported aromatic
systems.

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1")
kekulized = mol.kekulize()

for bond in kekulized.bonds():
    print(bond.begin_atom_idx(), bond.end_atom_idx(), bond.bond_type())
```

## Tetrahedral Stereo

COSMolKit exposes an ordered tetrahedral stereo representation derived from the
internal chiral tags. Each record is `(center_atom_index, ordered_ligands)`.
Implicit hydrogen is represented as `None`.

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("[13CH3:7][C@H](F)Cl")

for center, ligands in mol.tetrahedral_stereo():
    print("center:", center)
    print("ordered ligands:", ligands)
```

You can also inspect RDKit-like chiral tags:

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("F[C@H](Cl)Br")
print(mol.find_chiral_centers(include_unassigned=False))
```

## Graph Data for ML Workflows

The current low-level API can be used to build simple graph tensors in user
code. This example intentionally ignores 3D coordinates because COSMolKit does
not implement conformer generation yet.

```python
import torch
from cosmolkit import Molecule


def graph_from_smiles(smiles: str):
    mol = Molecule.from_smiles(smiles)

    x = torch.tensor([[atom.atomic_num()] for atom in mol.atoms()], dtype=torch.long)

    edges = []
    edge_attr = []
    for bond in mol.bonds():
        i = bond.begin_atom_idx()
        j = bond.end_atom_idx()
        kind = bond.bond_type()
        edges.extend([(i, j), (j, i)])
        edge_attr.extend([kind, kind])

    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    return x, edge_index, edge_attr


x, edge_index, edge_attr = graph_from_smiles("CCO")
print(x)
print(edge_index)
print(edge_attr)
```

## Planned API Examples

The following examples show the intended direction of the public API. They are
included to document the target design, but they are not implemented yet and
will raise `NotImplementedError` in current releases.

### SDF IO

```python
from cosmolkit import Molecule

# Not implemented yet
mol = Molecule.read_sdf("ligand.sdf", sanitize=True)
mol.write_sdf("out.sdf")
```

### Sanitization and Valence Checks

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1", sanitize=False)

# Not implemented yet
report = mol.check_valence()
mol = mol.sanitize(strict=True)
mol = mol.perceive_rings().perceive_aromaticity()
```

### Substructure Search and Fingerprints

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1CCO")

# Not implemented yet
query = Molecule.query_from_smarts("c1ccccc1")
matches = mol.substructure_find(query)

# Not implemented yet
fp1 = mol.fingerprint_morgan(radius=2, n_bits=2048)
fp2 = Molecule.from_smiles("CCN").fingerprint_morgan(radius=2, n_bits=2048)
print(fp1.tanimoto(fp2))
```

### 3D Coordinates and Alignment

```python
from cosmolkit import Alignment, Molecule

mol = Molecule.from_smiles("CCO")

# Not implemented yet
mol3d = mol.add_hydrogens().embed_3d(seed=42, num_conformers=20)

# Not implemented yet
result = Alignment.find_most_similar_fragment(
    reference=mol3d,
    candidates=[mol3d],
    mutate_reference=False,
    mutate_candidates=False,
)
```

### Explicit Editing

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO")

# Not implemented yet
editor = mol.edit()
cl = editor.add_atom("Cl")
editor.add_bond(0, cl, order="single")
mol2 = editor.commit(sanitize=True)
```

## Development Status

COSMolKit is not a full RDKit replacement today. The implementation is being
built through strict parity testing against RDKit behavior where compatibility
is the goal.
