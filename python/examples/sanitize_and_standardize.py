"""Public Python API example: explicit sanitize and standardize steps."""

import cosmolkit as ck

raw = ck.mol_from_smiles("CN(=O)=O", sanitize=False)
sanitized = raw.sanitize()

# Inspect detached stored atom values without recalculating unsanitized valence.
print("raw charges:", [atom.formal_charge() for atom in raw.to_builder().atoms()])
print("sanitized charges:", [atom.formal_charge() for atom in sanitized.atoms()])
print("sanitized smiles:", sanitized.to_smiles())

kekule = ck.mol_from_smiles("c1ccccc1").with_kekulized_bonds(mark_atoms_bonds=True)

print("kekulized molecule smiles:", kekule.to_smiles(kekule=True))
