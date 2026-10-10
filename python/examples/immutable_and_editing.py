"""Public Python API example: value-style transforms and explicit editing."""

import cosmolkit as ck

mol = ck.mol_from_smiles("CCO", sanitize=True)

mol_h = mol.with_hydrogens()
mol_no_h = mol_h.without_hydrogens()
_ = mol_no_h

builder = mol.to_builder()
cl = builder.add_atom(ck.AtomSpec(ck.Element.CL))
_ = builder.add_bond(ck.BondSpec(0, cl, ck.BondOrder.SINGLE))
mol2 = builder.build().sanitize()

assert mol.to_smiles() == "CCO"
assert mol2.num_atoms() == mol.num_atoms() + 1
print("source:", mol.to_smiles())
print("edited:", mol2.to_smiles())
