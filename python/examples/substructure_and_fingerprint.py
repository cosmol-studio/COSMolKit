"""Ideal COSMolKit usage: substructure search + fingerprint similarity.

This script demonstrates the target public API shape only.
Current bindings are placeholders and will raise NotImplementedError.
"""

from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1CCO", sanitize=True)
query = Molecule.query_from_smarts("c1ccccc1")
matches = mol.substructure_find(query)

fp1 = mol.fingerprint_morgan(radius=2, n_bits=2048)
fp2 = Molecule.from_smiles("CCN", sanitize=True).fingerprint_morgan(radius=2, n_bits=2048)
similarity = fp1.tanimoto(fp2)

_ = (matches, similarity)
