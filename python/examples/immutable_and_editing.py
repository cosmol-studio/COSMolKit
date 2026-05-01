"""COSMolKit usage: immutable-first transforms + explicit edit context."""

from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO", sanitize=True)

mol_h = mol.with_hydrogens()
mol_no_h = mol_h.without_hydrogens()
_ = mol_no_h

editor = mol.edit()
cl = editor.add_atom("Cl")
editor.add_bond(0, cl, order="single")
mol2 = editor.commit(sanitize=True)
_ = mol2
