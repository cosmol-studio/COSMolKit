"""COSMolKit usage: SMILES export + DG bounds matrix."""

import numpy as np

from cosmolkit import Molecule

mol = Molecule.from_smiles("F[C@H](Cl)Br", sanitize=True)

print("isomeric smiles:", mol.to_smiles())
print("non-isomeric smiles:", mol.to_smiles(False))
print(
    "chiral tags:",
    [
        (atom.idx(), atom.chiral_tag())
        for atom in mol.atoms()
        if atom.chiral_tag() != "CHI_UNSPECIFIED"
    ],
)

bounds = mol.dg_bounds_matrix()
print("bounds matrix shape:", bounds.shape)
print("d(0,1):", bounds[0, 1])
print("row sums:", bounds.sum(axis=1))
print("diagonal is zero:", np.allclose(np.diag(bounds), 0.0))
