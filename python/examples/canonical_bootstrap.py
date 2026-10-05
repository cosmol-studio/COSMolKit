"""Canonical facade-backed subset available in the default extension.

The complete 0.3.0 inventory is tracked separately; this example uses only
installed APIs and does not imply that all historical interfaces are ready.
"""

import cosmolkit

molecule = cosmolkit.Molecule.from_smiles("CCO")
print(cosmolkit.version(), molecule.num_atoms(), molecule.to_smiles())
print(molecule.morgan_fingerprint().on_bits())
params = cosmolkit.SmilesWriteParams(all_bonds_explicit=True)
print(molecule.to_smiles_with_params(params))
positioned = molecule.with_2d_coordinates()
print(positioned.coordinates_2d())
