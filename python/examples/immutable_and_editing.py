"""Ideal COSMolKit usage: immutable-first + explicit edit context.

This script demonstrates the target public API shape only.
Current bindings are placeholders and will raise NotImplementedError.
"""

from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO", sanitize=True)

mol_h = mol.add_hydrogens()
mol_no_h = mol_h.remove_hydrogens()
_ = mol_no_h

editor = mol.edit()
cl = editor.add_atom("Cl")
editor.add_bond(0, cl, order="single")
mol2 = editor.commit(sanitize=True)
_ = mol2
