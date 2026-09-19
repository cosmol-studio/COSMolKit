"""Modern, provenance, bulk, and legacy Topological Torsion fingerprints."""

import cosmolkit as ck


molecules = [
    ck.Molecule.from_smiles("CCCCO"),
    ck.Molecule.from_smiles("CCCCC"),
]
generator = ck.get_topological_torsion_generator(fp_size=2048)

sparse_count = generator.get_sparse_count_fingerprint(molecules[0])
sparse_bit = generator.get_sparse_fingerprint(molecules[0])
count = generator.get_count_fingerprint(molecules[0])
bit = generator.get_fingerprint(molecules[0])

print("sparse count:", sparse_count.nonzero_elements())
print("sparse bit:", sparse_bit.on_bits())
print("folded count:", count.nonzero_elements())
print("explicit bit:", bit.on_bits())

additional = ck.AdditionalOutput()
additional.allocate_atom_to_bits()
additional.allocate_atom_counts()
additional.allocate_bit_paths()
additional.allocate_atoms_per_bit()
_ = generator.get_fingerprint(molecules[0], additional_output=additional)
print("atom to bits:", additional.atom_to_bits())
print("atom counts:", additional.atom_counts())
print("bit paths:", additional.bit_paths())

bulk = generator.get_fingerprints(molecules, num_threads=2)
assert bulk[0].on_bits() == bit.on_bits()

legacy_unfolded = ck.get_topological_torsion_fingerprint(molecules[0])
legacy_hashed_count = ck.get_hashed_topological_torsion_fingerprint(molecules[0])
legacy_hashed_bit = ck.get_hashed_topological_torsion_fingerprint_as_bit_vect(
    molecules[0]
)
print("legacy unfolded:", legacy_unfolded.nonzero_elements())
print("legacy hashed count:", legacy_hashed_count.nonzero_elements())
print("legacy hashed bit:", legacy_hashed_bit.on_bits())
