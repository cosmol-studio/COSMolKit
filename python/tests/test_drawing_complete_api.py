"""DRAW delivery proposals: canonical mutation/query/value and typed IO behavior."""
import errno
import pytest
import cosmolkit


def test_default_and_configured_mutation_match_value_transform():
    source = cosmolkit.Molecule.from_smiles("CCO")
    value = source.with_2d_coordinates()
    assert source.coordinates_2d() is None and not source.has_2d_coordinates()
    assert value.has_2d_coordinates()
    assert source.compute_2d_coordinates_() is None
    assert source.has_2d_coordinates()
    assert source.coordinates_2d() == value.coordinates_2d()
    configured = cosmolkit.Molecule.from_smiles("CCO")
    assert configured.compute_2d_coordinates_with_params_(cosmolkit.Coordinate2DParams()) is None
    assert configured.coordinates_2d() == value.coordinates_2d()
    assert configured.to_smiles() == "CCO"


def test_empty_coordinate_presence_distinguishes_absence():
    molecule = cosmolkit.Molecule.new()
    assert not molecule.has_2d_coordinates() and molecule.coordinates_2d() is None
    assert molecule.compute_2d_coordinates_() is None
    assert molecule.has_2d_coordinates() and molecule.coordinates_2d() == []


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
