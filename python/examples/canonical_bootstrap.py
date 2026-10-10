"""Construct a molecule, write SMILES, fingerprint it, and generate 2D coordinates."""

import cosmolkit as ck

molecule = ck.mol_from_smiles("CCO")
print(ck.version(), molecule.num_atoms(), molecule.to_smiles())
print(molecule.fingerprint_morgan().on_bits())
params = ck.SmilesWriteParams(all_bonds_explicit=True)
print(molecule.to_smiles_with_params(params))
positioned = molecule.with_2d_coordinates()
print(positioned.coordinates_2d())
