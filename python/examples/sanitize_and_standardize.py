"""COSMolKit usage: explicit sanitize and standardize steps."""

from cosmolkit import Molecule

raw = Molecule.from_smiles("CN(=O)=O", sanitize=False)
sanitized = raw.sanitize(strict=True)

print("raw charges:", [atom.formal_charge() for atom in raw.atoms()])
print("sanitized charges:", [atom.formal_charge() for atom in sanitized.atoms()])
print("sanitized smiles:", sanitized.to_smiles())

kekule = Molecule.from_smiles("c1ccccc1").with_kekulized_bonds(clear_aromatic_flags=True)

print("kekulized molecule smiles:", kekule.to_smiles(kekule=True))
