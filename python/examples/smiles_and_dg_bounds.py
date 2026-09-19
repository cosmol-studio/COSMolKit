"""Public Python API example: SMILES export and DG bounds matrices."""

import numpy as np

import cosmolkit as ck

mol = ck.Molecule.from_smiles("F[C@H](Cl)Br", sanitize=True)

print("isomeric smiles:", mol.to_smiles())
print("non-isomeric smiles:", mol.to_smiles(False))
print(
    "chiral tags:",
    [
        (atom.idx(), atom.chiral_tag().name)
        for atom in mol.atoms()
        if atom.chiral_tag() != ck.ChiralTag.CHI_UNSPECIFIED
    ],
)

bounds = mol.dg_bounds_matrix()
print("bounds matrix shape:", bounds.shape)
print("d(0,1):", bounds[0, 1])
print("row sums:", bounds.sum(axis=1))
print("diagonal is zero:", np.allclose(np.diag(bounds), 0.0))
