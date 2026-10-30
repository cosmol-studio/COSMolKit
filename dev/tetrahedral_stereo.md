# Tetrahedral Stereo

`tetrahedral_stereo` is a read-only query using the same tetrahedral-center
perception as `find_chiral_centers`. Python/JavaScript default
`include_unassigned` to true; Rust takes an explicit bool. False returns only
centers with specified configuration.

```rust
pub struct TetrahedralStereo {
    pub center: AtomId,
    pub ligands: Option<[LigandRef; 4]>,
}
```

Python returns `list[tuple[int, list[int | None] | None]]`:

- `(center, None)`: potential center with unspecified or unknown configuration.
- `(center, [a, b, c, None])`: specified configuration; inner `None` has no
  atom index (implicit hydrogen or lone pair). Rust distinguishes these as
  `LigandRef::ImplicitHydrogen` and `LigandRef::LonePair`.
- `[]`: no centers match the filter.

## Ligand Order

`center` and the integers in `ligands` are molecule atom indices. For
`F[C@H](Cl)Br`, these are `F=0, C=1, Cl=2, Br=3`, so
`(1, [0, 2, 3, None])` means center C with ordered positions F, Cl, Br, H.

The fourth position is the reference ligand. For matching nondegenerate 3D
geometry, when the first three positions have atom coordinates, the order
encodes positive orientation:

```text
det(pos(ligands[0]) - pos(center),
    pos(ligands[1]) - pos(center),
    pos(ligands[2]) - pos(center)) > 0
```

Thus `[0, 2, 3, None]` and `[0, 3, 2, None]` encode opposite handedness.
The query reads chemical stereo, not coordinates; it does not require or
repair a 3D conformer.

## Examples

```python
from cosmolkit import Molecule

Molecule.from_smiles("F[C@H](Cl)Br").tetrahedral_stereo()
# [(1, [0, 2, 3, None])]

Molecule.from_smiles("F[C@@H](Cl)Br").tetrahedral_stereo()
# [(1, [0, 3, 2, None])]

Molecule.from_smiles("CC(O)C(=O)O").tetrahedral_stereo()
# [(1, None)]

Molecule.from_smiles("CC(O)C(=O)O").tetrahedral_stereo(include_unassigned=False)
# []

Molecule.from_smiles("CCO").tetrahedral_stereo()
# []
```
