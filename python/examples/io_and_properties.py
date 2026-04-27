"""Ideal COSMolKit usage: IO + basic properties.

This script demonstrates the target public API shape only.
Current bindings are placeholders and will raise NotImplementedError.
"""

from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO", sanitize=True)
lig = Molecule.read_sdf("ligand.sdf", sanitize=True)

# Planned pythonic property access (not implemented yet):
# print(mol.formula, mol.exact_mass, mol.num_atoms, mol.num_bonds)
# for atom in mol.atoms:
#     print(atom.index, atom.element, atom.formal_charge, atom.total_valence)

_ = lig
