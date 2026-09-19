"""Public Python API example: explicit sanitize and standardize steps."""

import cosmolkit as ck

raw = ck.Molecule.from_smiles("CN(=O)=O", sanitize=False)
sanitized = raw.sanitize(strict=True)

print("raw charges:", [atom.formal_charge() for atom in raw.atoms()])
print("sanitized charges:", [atom.formal_charge() for atom in sanitized.atoms()])
print("sanitized smiles:", sanitized.to_smiles())

kekule = ck.Molecule.from_smiles("c1ccccc1").with_kekulized_bonds(clear_aromatic_flags=True)

print("kekulized molecule smiles:", kekule.to_smiles(kekule=True))
