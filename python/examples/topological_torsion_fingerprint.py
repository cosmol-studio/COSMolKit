"""Modern, provenance, bulk, and legacy Topological Torsion fingerprints."""

import cosmolkit as ck


molecules = [
    ck.mol_from_smiles("CCCCO"),
    ck.mol_from_smiles("CCCCC"),
]
generator = ck.TopologicalTorsionFingerprintGenerator(
    params=ck.TopologicalTorsionParams(fp_size=2048)
)

sparse_count = molecules[0].fingerprint_topological_torsion_sparse_count_with_generator(generator)
sparse_bit = molecules[0].fingerprint_topological_torsion_sparse_with_generator(generator)
count = molecules[0].fingerprint_topological_torsion_count_with_generator(generator)
bit = molecules[0].fingerprint_topological_torsion_with_generator(generator)

print("sparse count:", sparse_count.nonzero_elements())
print("sparse bit:", sparse_bit.on_bits())
print("folded count:", count.nonzero_elements())
print("explicit bit:", bit.on_bits())

additional = ck.FingerprintAdditionalOutput()
additional.allocate_atom_to_bits()
additional.allocate_atom_counts()
additional.allocate_bit_paths()
additional.allocate_atoms_per_bit()
_ = molecules[0].fingerprint_topological_torsion_with_generator(generator, output=additional)
print("atom to bits:", additional.atom_to_bits())
print("atom counts:", additional.atom_counts())
print("bit paths:", additional.bit_paths())

bulk = generator.fingerprints(molecules, num_threads=2)
first = bulk[0]
assert first is not None
assert first.on_bits() == bit.on_bits()

legacy_unfolded = molecules[0].fingerprint_topological_torsion_sparse_count_legacy()
legacy_hashed_count = molecules[0].fingerprint_topological_torsion_count_legacy()
legacy_hashed_bit = molecules[0].fingerprint_topological_torsion_legacy()
print("legacy unfolded:", legacy_unfolded.nonzero_elements())
print("legacy hashed count:", legacy_hashed_count.nonzero_elements())
print("legacy hashed bit:", legacy_hashed_bit.on_bits())
