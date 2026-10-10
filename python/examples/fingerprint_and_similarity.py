"""Public Python API example: exact-parity fingerprint families."""

import cosmolkit as ck

mol = ck.mol_from_smiles("c1ccccc1CCO", sanitize=True)
fp1 = mol.fingerprint_morgan(radius=2, fp_size=2048)
fp2 = ck.mol_from_smiles("CCN", sanitize=True).fingerprint_morgan(
    radius=2, fp_size=2048
)
similarity = fp1.tanimoto(fp2)

additional = ck.FingerprintAdditionalOutput()
additional.allocate_atom_counts()
additional.allocate_bit_info_map()
_ = mol.fingerprint_morgan(radius=2, fp_size=2048, additional_output=additional)

topological = mol.fingerprint_topological(fp_size=2048)
topological_output = mol.fingerprint_topological_with_output(
    fp_size=2048,
    atom_bits=True,
    bit_info=True,
)
avalon = mol.fingerprint_avalon(n_bits=512)
pattern = mol.fingerprint_pattern(n_bits=2048)
tautomeric_pattern = mol.fingerprint_pattern(n_bits=2048, tautomeric=True)
maccs = mol.fingerprint_maccs()
atom_pair = mol.fingerprint_atom_pair(fp_size=2048)
atom_pair_sparse_count = mol.fingerprint_atom_pair_sparse_count()
atom_pair_output = ck.FingerprintAdditionalOutput()
atom_pair_output.allocate_atom_counts()
_ = mol.fingerprint_atom_pair(additional_output=atom_pair_output)
layered_mask = mol.fingerprint_layered(fp_size=257)
layered = mol.fingerprint_layered_with_output(
    layers=0x3F,
    min_path=2,
    max_path=4,
    fp_size=257,
    atom_counts=[10] * mol.num_atoms(),
    set_only_bits=layered_mask,
    branched_paths=True,
    from_atoms=[0],
)

print("fp1 bits:", fp1.on_bits()[:8])
print("similarity:", similarity)
print("atom counts:", additional.atom_counts())
bit_info = additional.bit_info_map()
assert bit_info is not None
print("bit info entries:", len(bit_info))
print("topological bits:", topological.on_bits()[:8])
print(
    "topological atom provenance counts:",
    [len(bits) for bits in topological_output.atom_bits()],
)
print("topological path entries:", len(topological_output.bit_info()))
print("Avalon bits:", avalon.on_bits()[:8])
print("Pattern bits:", pattern.on_bits()[:8])
print("Tautomeric Pattern bits:", tautomeric_pattern.on_bits()[:8])
print("MACCS bits:", maccs.on_bits()[:8])
print("AtomPair bits:", atom_pair.on_bits()[:8])
print(
    "AtomPair sparse counts:",
    list(atom_pair_sparse_count.nonzero_elements().items())[:8],
)
print("AtomPair atom counts:", atom_pair_output.atom_counts())
print("Layered bits:", layered.fingerprint().on_bits()[:8])
print("Layered seeded atom counts:", layered.atom_counts())
