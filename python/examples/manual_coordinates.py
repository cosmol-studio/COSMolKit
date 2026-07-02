import numpy as np

from cosmolkit import Molecule


mol = Molecule.from_smiles("CCO")

coords_2d = np.array(
    [
        [0.0, 0.0],
        [1.5, 0.0],
        [2.1, 1.2],
    ],
    dtype=np.float64,
)
mol_2d = mol.with_2d_coordinates(coords_2d)

coords_3d = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.1, 1.2, 0.4],
    ],
    dtype=np.float64,
)
mol_3d = mol.with_added_3d_conformer(coords_3d)
shifted = mol_3d.with_3d_coordinates(coords_3d + [0.0, 0.0, 1.0])
single = mol_3d.with_only_3d_conformer(coords_3d + [0.0, 0.0, 2.0])
cleared = mol_3d.with_cleared_3d_conformers()

editable = mol.with_2d_coordinates(coords_2d)
conf_id = editable.add_3d_conformer_(coords_3d)
editable.set_3d_coordinates_(coords_3d + [0.0, 0.0, 2.0], conformer_index=conf_id)
editable.clear_3d_conformers_()
only_conf_id = editable.set_only_3d_conformer_(coords_3d)

print(mol_2d.coordinates_2d())
print(mol_3d.num_conformers())
print(shifted.coordinates_3d())
print(single.num_conformers())
print(cleared.num_conformers())
print(editable.coordinates_3d(only_conf_id))
