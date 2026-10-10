"""Pinned original manual input conditions projected onto canonical parameter APIs.
Python returns owned NumPy rows with the runtime's dimension: XY or XYZ.
"""
import cosmolkit
import numpy as np
import pytest

def xyz(mol, index=0):
    return np.asarray(mol.coordinates_3d(index))

def test_empty_coordinate_arrays_keep_their_dimension():
    mol = cosmolkit.Molecule.from_smiles("")
    assert mol.coordinates_2d() is None
    two_d = mol.with_2d_coordinate_block(np.empty((0, 2)))
    three_d = two_d.with_only_3d_conformer(np.empty((0, 3)))
    for output, dimension in [
        (three_d.coordinates_2d(), 2),
        (three_d.coordinates_3d(), 3),
        (three_d.conformers_3d()[0].coordinates(), 3),
    ]:
        assert isinstance(output, np.ndarray)
        assert output.shape == (0, dimension)
        assert output.dtype == np.float64

def test_setting_2d_coordinates_is_value_style_and_validates_input():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    coords = np.array([[0.0, 0.0], [1.5, 0.0], [3.0, 0.0]], dtype=np.float32)

    with_coords = mol.with_2d_coordinate_block(coords)

    assert with_coords is not mol
    assert not mol.has_2d_coordinates()
    assert with_coords.has_2d_coordinates()
    output = with_coords.coordinates_2d()
    assert isinstance(output, np.ndarray)
    assert output.shape == (3, 2) and output.dtype == np.float64
    np.testing.assert_array_equal(output, coords)
    output[:] = 123
    np.testing.assert_array_equal(with_coords.coordinates_2d(), coords)

    with pytest.raises(ValueError, match="row count mismatch"):
        mol.with_2d_coordinate_block([[0.0, 0.0]])

    with pytest.raises(ValueError, match="non-finite"):
        mol.with_2d_coordinate_block([[0.0, 0.0], [1.0, np.nan], [2.0, 0.0]])

def test_setting_2d_coordinates_z_policy_and_in_place_update():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    coords3 = [[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [2.0, 0.2, 0.0]]

    strict = mol.with_2d_coordinate_block_with_params(coords3, cosmolkit.Coordinate2DInputParams(cosmolkit.CoordinateZPolicy.RequireZero))
    np.testing.assert_array_equal(strict.coordinates_2d(), np.asarray(coords3)[:, :2])

    with pytest.raises(ValueError, match="z_policy='error'"):
        mol.with_2d_coordinate_block_with_params(coords3, cosmolkit.Coordinate2DInputParams(cosmolkit.CoordinateZPolicy.Error))

    with pytest.raises(ValueError, match="require zero z"):
        mol.with_2d_coordinate_block_with_params(
            [[0.0, 0.0, 0.0], [1.0, 0.1, 0.25], [2.0, 0.2, 0.0]],
            cosmolkit.Coordinate2DInputParams(cosmolkit.CoordinateZPolicy.RequireZero),
        )

    assert mol.set_2d_coordinates_(coords3) is None
    assert mol.has_2d_coordinates()
    assert mol.coordinates_2d().shape == (3, 2)
    np.testing.assert_array_equal(mol.coordinates_2d(), np.asarray(coords3)[:, :2])

def test_adding_and_replacing_3d_coordinates_preserves_value_semantics():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    first = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.0, 1.0, 0.0]])
    second = first + np.array([0.0, 0.0, 1.0])
    replacement = first + np.array([0.25, 0.5, 0.75])

    one_conf = mol.with_added_3d_conformer(first)

    assert one_conf is not mol
    assert len(mol.conformers_3d()) == 0
    assert len(one_conf.conformers_3d()) == 1
    output = one_conf.coordinates_3d()
    assert isinstance(output, np.ndarray)
    assert output.shape == (3, 3) and output.dtype == np.float64
    np.testing.assert_array_equal(output.view(np.uint64), first.view(np.uint64))
    output[:] = 123
    np.testing.assert_array_equal(one_conf.coordinates_3d(), first)
    assert np.allclose(xyz(one_conf), first)

    two_confs = one_conf.with_added_3d_conformer(second.astype(np.float32))
    replaced = two_confs.with_3d_coordinates_with_params(replacement, cosmolkit.Replace3DCoordinatesParams(conformer_id=0))

    assert len(two_confs.conformers_3d()) == 2
    assert len(replaced.conformers_3d()) == 2
    assert np.allclose(xyz(two_confs, 0), first)
    assert np.allclose(xyz(replaced, 0), replacement)
    assert np.allclose(xyz(replaced, 1), second)

def test_3d_coordinate_in_place_api_returns_conformer_ids_and_validates_input():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    first = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.0, 1.0, 0.0]]
    replacement = [[0.1, 0.2, 0.3], [1.5, 0.2, 0.3], [2.1, 1.2, 0.3]]

    assert mol.add_3d_conformer_(first) == 0
    assert len(mol.conformers_3d()) == 1
    assert mol.set_3d_coordinates_(replacement) is None
    assert np.allclose(xyz(mol), replacement)

    with pytest.raises(ValueError, match="ConformerRowCount|row count mismatch"):
        mol.add_3d_conformer_([[0.0, 0.0, 0.0]])

    with pytest.raises(ValueError, match="no 3D conformer with id 7"):
        mol.set_3d_coordinates_with_params_(replacement, cosmolkit.Replace3DCoordinatesParams(conformer_id=7))

    with pytest.raises(ValueError, match="shape"):
        mol.with_added_3d_conformer([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])

    with pytest.raises(ValueError, match="non-finite"):
        mol.with_added_3d_conformer([[0.0, 0.0, 0.0], [1.0, float("inf"), 0.0], [2.0, 0.0, 0.0]])

def test_coordinate_ingress_rejects_oversized_broadcast_views_before_copying():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    oversized_3d = np.broadcast_to(
        np.zeros((1, 3), dtype=np.int8), (1_000_000, 3)
    )
    oversized_2d = np.broadcast_to(
        np.zeros((1, 2), dtype=np.int8), (1_000_000, 2)
    )

    with pytest.raises(
        ValueError,
        match=r"^3D coordinates row count mismatch: expected 3, got 1000000$",
    ):
        mol.with_only_3d_conformer(oversized_3d)

    with pytest.raises(
        ValueError,
        match=r"^2D coordinates row count mismatch: expected 3, got 1000000$",
    ):
        mol.with_2d_coordinate_block(oversized_2d)

def test_3d_conformer_clear_and_single_conformer_assignment_use_value_semantics():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    first = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.0, 1.0, 0.0]])
    second = first + np.array([0.0, 0.0, 1.0])
    replacement = first + np.array([0.25, 0.5, 0.75])

    multi = mol.with_added_3d_conformer(first).with_added_3d_conformer(second)
    cleared = multi.with_cleared_3d_conformers()
    single = multi.with_only_3d_conformer(replacement)

    assert len(multi.conformers_3d()) == 2
    assert len(cleared.conformers_3d()) == 0
    assert len(single.conformers_3d()) == 1
    assert np.allclose(xyz(single), replacement)

    with pytest.raises(ValueError, match="no 3D conformer"):
        cleared.coordinates_3d()

def test_3d_conformer_clear_and_single_conformer_assignment_in_place():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    first = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.0, 1.0, 0.0]]
    second = [[0.0, 0.0, 1.0], [1.4, 0.0, 1.0], [2.0, 1.0, 1.0]]

    mol.add_3d_conformer_(first)
    mol.add_3d_conformer_(second)
    assert len(mol.conformers_3d()) == 2

    assert mol.clear_3d_conformers_() is None
    assert len(mol.conformers_3d()) == 0

    assert mol.set_only_3d_conformer_(second) == 0
    assert len(mol.conformers_3d()) == 1
    assert np.allclose(xyz(mol), second)

    with pytest.raises(ValueError, match="row count mismatch"):
        mol.set_only_3d_conformer_([[0.0, 0.0, 0.0]])
