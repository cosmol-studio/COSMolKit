"""DRAW delivery proposals: canonical mutation/query/value and typed IO behavior."""
import errno
import numpy as np
import pytest
import cosmolkit


def test_default_and_configured_mutation_match_value_transform():
    source = cosmolkit.Molecule.from_smiles("CCO")
    assert not hasattr(source, "add_2d_coordinates")
    assert not hasattr(source, "compute_2d_coordinates_")
    assert not hasattr(source, "add_2d_coordinates_with_params")
    assert not hasattr(source, "compute_2d_coordinates_with_params_")
    value = source.with_2d_coordinates()
    assert source.coordinates_2d() is None and not source.has_2d_coordinates()
    assert value.has_2d_coordinates()
    assert source.add_2d_coordinates_() is None
    assert source.has_2d_coordinates()
    np.testing.assert_array_equal(source.coordinates_2d(), value.coordinates_2d())
    configured = cosmolkit.Molecule.from_smiles("CCO")
    assert configured.add_2d_coordinates_with_params_(cosmolkit.Coordinate2DParams()) is None
    np.testing.assert_array_equal(configured.coordinates_2d(), value.coordinates_2d())
    assert configured.to_smiles() == "CCO"


def test_with_2d_configuration_call_forms_share_one_operation():
    source = cosmolkit.Molecule.from_smiles("CCO")
    params = cosmolkit.Coordinate2DParams(canonical_orientation=True)
    configured = source.with_2d_coordinates_with_params(params)
    for value in (
        source.with_2d_coordinates(params),
        source.with_2d_coordinates(canonical_orientation=True),
    ):
        np.testing.assert_array_equal(value.coordinates_2d(), configured.coordinates_2d())
    assert source.coordinates_2d() is None
    assert source.add_2d_coordinates_(canonical_orientation=True) is None
    np.testing.assert_array_equal(source.coordinates_2d(), configured.coordinates_2d())


def test_empty_coordinate_presence_distinguishes_absence():
    molecule = cosmolkit.Molecule.new()
    assert not molecule.has_2d_coordinates() and molecule.coordinates_2d() is None
    assert molecule.add_2d_coordinates_() is None
    assert molecule.has_2d_coordinates()
    output = molecule.coordinates_2d()
    assert isinstance(output, np.ndarray) and output.shape == (0, 2)


@pytest.mark.parametrize("suffix", ["svg", "png"])
def test_write_error_has_registered_type_and_os_payload(tmp_path, suffix):
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    path = tmp_path / "missing" / ("drawing." + suffix)
    with pytest.raises(cosmolkit.DrawingWriteError) as caught:
        getattr(molecule, "write_" + suffix)(str(path), 120, 80)
    error = caught.value
    assert isinstance(error, OSError)
    assert error.domain == "drawing" and error.kind == "Io"
    assert error.errno == errno.ENOENT and error.filename == str(path)
    assert not molecule.has_2d_coordinates()
    assert not path.exists()
