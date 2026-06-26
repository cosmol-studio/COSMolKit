# Tetrahedral Stereo

COSMolKit exposes tetrahedral stereo as:

```rust
pub struct TetrahedralStereo {
    pub center: usize,
    pub ligands: [LigandRef; 4],
}
```

`ligands` is not just the adjacency list. Its order is the stereochemical
value returned by `tetrahedral_stereo()`.

## Equality Rule

For one stereocenter, even permutations of the four ligands represent the same
handedness. COSMolKit canonicalizes them to the lexicographically smallest
ordered ligand list.

Odd permutations are not canonicalized together with even permutations because
they represent the opposite handedness.

For explicit four-ligand centers, all four positions participate in this rule.
For implicit-hydrogen centers, `None`/`ImplicitHydrogen` is kept in the fourth
slot so the first three real atoms remain usable for coordinate checks.

## `None`

`None` does not mean "no ligand". It means the fourth ligand is a hydrogen that
exists chemically but is implicit in the current molecule graph, so it has no
atom index. In Rust this is `LigandRef::ImplicitHydrogen`.

## Examples

```python
from cosmolkit import Molecule

Molecule.from_smiles("F[C@H](Cl)Br").tetrahedral_stereo()
# [(1, [0, 2, 3, None])]

Molecule.from_smiles("F[C@@H](Cl)Br").tetrahedral_stereo()
# [(1, [0, 3, 2, None])]

Molecule.from_smiles("F[C@](Cl)(Br)I").tetrahedral_stereo()
# [(1, [0, 2, 3, 4])]

Molecule.from_smiles("F[C@@](Cl)(Br)I").tetrahedral_stereo()
# [(1, [0, 2, 4, 3])]
```

RDKit/SMILES-style `CW/CCW` tags remain important I/O semantics. COSMolKit
maps them into the ordered ligand list instead of exposing a separate
`CW/CCW` field on `TetrahedralStereo`.
