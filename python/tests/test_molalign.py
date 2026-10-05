import cosmolkit
import numpy as np
import pytest
from numpy.typing import NDArray


REFERENCE_COORDINATES = np.array(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
    dtype=np.float64,
)
PROBE_COORDINATES = REFERENCE_COORDINATES + np.array([3.0, -2.0, 1.0])


def molecule_with_coordinates(coordinates: NDArray[np.float64]) -> cosmolkit.Molecule:
    return add_conformer(cosmolkit.Molecule.from_smiles("CCC"), coordinates)


def identity_map() -> list[cosmolkit.AlignmentAtomMap]:
    return [cosmolkit.AlignmentAtomMap(index, index) for index in range(3)]


def test_alignment_parameters_and_results_keep_typed_source_state():
    atom_map = identity_map()
    params = cosmolkit.AlignmentParameters(atom_map=atom_map, weights=[1.0, 2.0, 3.0])

    assert params.probe_conformer_id == -1
    assert params.reference_conformer_id == -1
    assert params.max_iterations == 50
    assert params.atom_map is not None
    assert [(pair.probe_atom, pair.reference_atom) for pair in params.atom_map] == [
        (0, 0),
        (1, 1),
        (2, 2),
    ]

    probe = molecule_with_coordinates(PROBE_COORDINATES)
    reference = molecule_with_coordinates(REFERENCE_COORDINATES)
    result = probe.alignment_transform_to(reference, params)

    assert result.rmsd() == pytest.approx(0.0, abs=1.0e-10)
    assert np.asarray(result.transform().matrix()).shape == (4, 4)
    assert [(pair.probe_atom, pair.reference_atom) for pair in result.atom_map()] == [
        (0, 0),
        (1, 1),
        (2, 2),
    ]


def test_read_only_alignment_and_coordinate_rmsd_do_not_mutate_inputs():
    probe = molecule_with_coordinates(PROBE_COORDINATES)
    reference = molecule_with_coordinates(REFERENCE_COORDINATES)
    probe_before = coordinates_of(probe).copy()
    reference_before = coordinates_of(reference).copy()
    atom_map = identity_map()

    aligned = probe.best_alignment_to(
        reference,
        cosmolkit.BestAlignmentParameters(atom_maps=[atom_map]),
    )
    best_rmsd = probe.best_rmsd_to(
        reference,
        cosmolkit.BestAlignmentParameters(atom_maps=[atom_map]),
    )
    coordinate_rmsd = probe.coordinate_rmsd_to(
        reference,
        cosmolkit.CoordinateRmsdParameters(atom_maps=[atom_map]),
    )

    assert aligned.rmsd() == pytest.approx(0.0, abs=1.0e-10)
    assert best_rmsd == pytest.approx(0.0, abs=1.0e-10)
    assert coordinate_rmsd == pytest.approx(np.sqrt(14.0))
    assert np.array_equal(coordinates_of(probe), probe_before)
    assert np.array_equal(coordinates_of(reference), reference_before)


def test_explicit_value_and_in_place_alignment_have_unambiguous_mutation():
    probe = molecule_with_coordinates(PROBE_COORDINATES)
    reference = molecule_with_coordinates(REFERENCE_COORDINATES)
    params = cosmolkit.AlignmentParameters(atom_map=identity_map())
    source_before = coordinates_of(probe).copy()

    aligned, result = probe.with_alignment_to(reference, params)

    assert result.rmsd() == pytest.approx(0.0, abs=1.0e-10)
    assert np.array_equal(coordinates_of(probe), source_before)
    assert np.allclose(coordinates_of(aligned), coordinates_of(reference))

    in_place_result = probe.align_to_(reference, params)
    assert in_place_result.rmsd() == pytest.approx(0.0, abs=1.0e-10)
    assert np.allclose(coordinates_of(probe), coordinates_of(reference))


def test_conformer_alignment_returns_source_ordered_reports():
    molecule = add_conformer(molecule_with_coordinates(REFERENCE_COORDINATES), PROBE_COORDINATES)
    source_second = coordinates_of(molecule, 1).copy()

    pair_rmsds = molecule.all_conformer_best_rmsds(
        cosmolkit.AllConformerRmsdParameters(atom_maps=[identity_map()])
    )
    aligned, report = molecule.with_aligned_conformers(
        cosmolkit.ConformerAlignmentParameters(atom_indices=[0, 1, 2])
    )

    assert len(pair_rmsds) == 1
    assert pair_rmsds[0].probe_conformer_id() == 1
    assert pair_rmsds[0].reference_conformer_id() == 0
    assert pair_rmsds[0].rmsd() == pytest.approx(0.0, abs=1.0e-10)
    assert report.rmsds() == pytest.approx([0.0], abs=1.0e-10)
    assert np.array_equal(coordinates_of(molecule, 1), source_second)
    assert np.allclose(coordinates_of(aligned, 1), coordinates_of(aligned, 0))

    in_place_report = molecule.align_conformers_()
    assert in_place_report.rmsds() == pytest.approx([0.0], abs=1.0e-10)
    assert np.allclose(coordinates_of(molecule, 1), coordinates_of(molecule, 0))


def test_alignment_errors_are_visible_and_leave_state_unchanged():
    molecule = cosmolkit.Molecule.from_smiles("CCC")
    reference = molecule_with_coordinates(REFERENCE_COORDINATES)

    with pytest.raises(ValueError, match="no 3D conformers"):
        _ = molecule.best_rmsd_to(reference)

    before = coordinates_of(reference).copy()
    with pytest.raises(ValueError, match="probe atom index 7 is out of range"):
        _ = reference.align_to_(
            reference,
            cosmolkit.AlignmentParameters(
                atom_map=[cosmolkit.AlignmentAtomMap(7, 0)]
            ),
        )
    assert np.array_equal(coordinates_of(reference), before)


# Local test fixtures use the current canonical public coordinate constructor.
def add_conformer(molecule, coordinates):
    builder = molecule.to_builder()
    builder.add_3d_conformer(np.asarray(coordinates).tolist())
    return builder.build()


def coordinates_of(molecule, conformer_id=0):
    return np.asarray(next(c.coordinates() for c in molecule.conformers_3d() if c.id() == conformer_id))
