"""COSMolKit usage: Morgan fingerprints and Tanimoto similarity."""

from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1CCO", sanitize=True)
fp1 = mol.fingerprint_morgan(radius=2, n_bits=2048)
fp2 = Molecule.from_smiles("CCN", sanitize=True).fingerprint_morgan(
    radius=2, n_bits=2048
)
similarity = fp1.tanimoto(fp2)

result = mol.fingerprint_morgan_with_output(radius=2, n_bits=2048)
additional = result.additional_output()

print("fp1 bits:", fp1.on_bits()[:8])
print("similarity:", similarity)
print("atom counts:", additional.atom_counts())
print("bit info entries:", len(additional.bit_info_map()))
