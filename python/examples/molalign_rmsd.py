"""Read-only RMSD measurement and explicit molecular alignment."""

import numpy as np

import cosmolkit as ck

reference_coordinates = np.array(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
)
probe_coordinates = reference_coordinates + np.array([3.0, -2.0, 1.0])

reference = ck.Molecule.from_smiles("CCC").with_only_3d_conformer(reference_coordinates)
probe = ck.Molecule.from_smiles("CCC").with_only_3d_conformer(probe_coordinates)
params = ck.AlignmentParameters(
    atom_map=[ck.AlignmentAtomMap(index, index) for index in range(3)]
)

transform_result = probe.alignment_transform_to(reference, params)
aligned_probe, applied_result = probe.with_alignment_to(reference, params)

print("read-only RMSD:", transform_result.rmsd())
print("applied RMSD:", applied_result.rmsd())
print("source unchanged:", np.array_equal(probe.coordinates_3d(), probe_coordinates))
print("aligned coordinates:\n", aligned_probe.coordinates_3d())
