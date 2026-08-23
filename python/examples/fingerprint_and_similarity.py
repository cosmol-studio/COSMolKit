"""Public Python API example: exact-parity fingerprint families."""

from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1CCO", sanitize=True)
fp1 = mol.fingerprint_morgan(radius=2, n_bits=2048)
fp2 = Molecule.from_smiles("CCN", sanitize=True).fingerprint_morgan(
    radius=2, n_bits=2048
)
similarity = fp1.tanimoto(fp2)

result = mol.fingerprint_morgan_with_output(radius=2, n_bits=2048)
additional = result.additional_output()

topological = mol.topological_fingerprint(fp_size=2048)
topological_output = mol.topological_fingerprint_with_output(
    fp_size=2048,
    atom_bits=True,
    bit_info=True,
)
avalon = mol.avalon_fingerprint(n_bits=512)
maccs = mol.maccs_fingerprint()
atom_pair = mol.fingerprint_atom_pair(n_bits=2048)
atom_pair_sparse_count = mol.fingerprint_atom_pair_sparse_count()
atom_pair_output = mol.fingerprint_atom_pair_with_output().additional_output()

print("fp1 bits:", fp1.on_bits()[:8])
print("similarity:", similarity)
print("atom counts:", additional.atom_counts())
print("bit info entries:", len(additional.bit_info_map()))
print("topological bits:", topological.on_bits()[:8])
print(
    "topological atom provenance counts:",
    [len(bits) for bits in topological_output.atom_bits()],
)
print("topological path entries:", len(topological_output.bit_info()))
print("Avalon bits:", avalon.on_bits()[:8])
print("MACCS bits:", maccs.on_bits()[:8])
print("AtomPair bits:", atom_pair.on_bits()[:8])
print(
    "AtomPair sparse counts:",
    list(atom_pair_sparse_count.nonzero_elements().items())[:8],
)
print("AtomPair atom counts:", atom_pair_output.atom_counts())
