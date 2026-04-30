"""COSMolKit usage: SMILES export + DG bounds matrix."""

from cosmolkit import Molecule

mol = Molecule.from_smiles("F[C@H](Cl)Br", sanitize=True)

print("isomeric smiles:", mol.to_smiles())
print("non-isomeric smiles:", mol.to_smiles(False))

bounds = mol.dg_bounds_matrix()
print("bounds matrix size:", len(bounds), "x", len(bounds[0]) if bounds else 0)
print("d(0,1):", bounds[0][1])
