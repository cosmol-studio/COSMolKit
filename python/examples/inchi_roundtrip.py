"""Public Python API example: object-oriented InChI conversion."""

from cosmolkit import Molecule, inchi_to_key

molecule = Molecule.from_smiles("CCO")
inchi = molecule.to_inchi()
molecule_key = molecule.to_inchi_key()

assert inchi_to_key(inchi) == molecule_key

restored = Molecule.from_inchi(inchi)
if restored is None:
    raise RuntimeError("the generated InChI did not produce a molecule")

print("InChI:", inchi)
print("InChIKey:", molecule_key)
print("restored SMILES:", restored.to_smiles())
