import cosmolkit as ck
import numpy as np
import pytest

def test_parameter_objects_keep_original_defaults_and_are_writable():
    for value, name, expected in [(ck.Coordinate2DInputParams(), 'z_policy', ck.CoordinateZPolicy.Ignore), (ck.Coordinate3DInputParams(), 'is_3d', True), (ck.Replace3DCoordinatesParams(), 'conformer_id', 0)]:
        assert getattr(value, name) == expected
        setattr(value, name, expected)
        assert getattr(value, name) == expected
        with pytest.raises(TypeError): setattr(value, name, None)
    assert ck.CoordinateZPolicy.from_name('ReQuIrE_ZeRo') == ck.CoordinateZPolicy.RequireZero
    with pytest.raises(ck.CoordinateInputError) as caught: ck.CoordinateZPolicy.from_name('unknown')
    assert (caught.value.domain, caught.value.kind, caught.value.value) == ('coordinate_input', 'UnknownZPolicy', 'unknown')
    with pytest.raises(OverflowError): ck.Replace3DCoordinatesParams(conformer_id=-1)

@pytest.mark.parametrize('z', [1e-12, -1e-12, 0.])
def test_original_inclusive_zero_z_threshold(z):
    source = ck.Molecule.from_smiles('CCO')
    result = source.with_2d_coordinate_block_with_params([[0., 0., z], [1., 2., z], [3., 4., z]], ck.Coordinate2DInputParams(ck.CoordinateZPolicy.RequireZero))
    np.testing.assert_array_equal(result.coordinates_2d(), [[0., 0.], [1., 2.], [3., 4.]])

@pytest.mark.parametrize('value', [float('inf'), float('-inf'), float('nan')])
def test_ignored_z_is_still_finite_checked_and_typed_failure_is_atomic(value):
    source = ck.Molecule.from_smiles('CCO')
    before = source.to_binary()
    with pytest.raises(ck.OperationError) as caught: source.set_2d_coordinates_([[0., 0., 0.], [1., 2., value], [3., 4., 0.]])
    assert caught.value.kind == 'CoordinateInput'
    cause = caught.value.__cause__
    assert isinstance(cause, ck.CoordinateInputError)
    assert (cause.domain, cause.kind, cause.dimension, cause.row, cause.column) == ('coordinate_input', 'NonFinite', '2D', 1, 2)
    assert source.to_binary() == before

def test_numeric_noncontiguous_protocol_does_not_alias_input_or_read_output():
    rows = np.arange(18, dtype=np.float32).reshape(3, 6)[:, ::2]
    expected = rows.copy()
    source = ck.Molecule.from_smiles('CCO')
    result = source.with_only_3d_conformer(rows)
    rows[:] = -99
    assert np.array_equal(result.conformers_3d()[0].coordinates(), expected)
    conformer = result.conformers_3d()[0]
    projected = conformer.coordinates()
    assert isinstance(projected, np.ndarray)
    assert projected.shape == (3, 3) and projected.dtype == np.float64
    np.testing.assert_array_equal(projected, result.coordinates_3d())
    projected[0][0] = -100
    assert np.array_equal(result.conformers_3d()[0].coordinates(), expected)
    np.testing.assert_array_equal(conformer.coordinates(), expected)
    assert source.conformers_3d() == []

def test_builder_existing_two_dimensional_projections_are_checked():
    builder = ck.Molecule.from_smiles('CCO').to_builder()
    assert builder.set_2d_coordinates([[0., 0.], [1., 2.], [3., 4.]]) is None
    np.testing.assert_array_equal(builder.build().coordinates_2d(), [[0., 0.], [1., 2.], [3., 4.]])
    with pytest.raises(ValueError, match='row count'): builder.set_2d_coordinates([[0., 0.]])
    with pytest.raises(ValueError, match='non-finite'): builder.set_2d_coordinates([[0., 0.], [1., float('nan')], [3., 4.]])
    with pytest.raises(ValueError, match='shape'): builder.set_2d_coordinates([[0., 0., 0.]] * 3)

def test_configured_inplace_non3d_append_only_and_all_typed_context():
    source = ck.Molecule.from_smiles('CCO')
    rows = [[0., 0., 0.], [1., 2., 3.], [4., 5., 6.]]
    params = ck.Coordinate3DInputParams(is_3d=False)
    assert source.add_3d_conformer_with_params_(rows, params) == 0
    assert source.add_3d_conformer_with_params_(rows, params) == 1
    assert [value.is_3d() for value in source.conformers_3d()] == [False, False]
    before = source.to_binary()
    with pytest.raises(ck.OperationError) as caught: source.set_3d_coordinates_with_params_(rows, ck.Replace3DCoordinatesParams(conformer_id=7))
    cause = caught.value.__cause__
    assert isinstance(cause, ck.CoordinateInputError)
    assert (cause.kind, cause.conformer_id, cause.count) == ('ConformerNotFound', 7, 2)
    assert source.to_binary() == before
    assert source.set_only_3d_conformer_with_params_(rows, params) == 0
    assert [value.id() for value in source.conformers_3d()] == [0]
    assert source.conformers_3d()[0].is_3d() is False

def test_xyz_missing_read_has_original_value_error_and_complete_context():
    source = ck.Molecule.from_smiles('CCO')
    before = source.to_binary()
    with pytest.raises(ck.Coordinate3DReadError, match='no 3D conformer') as caught: source.coordinates_3d()
    assert (caught.value.domain, caught.value.kind, caught.value.conformer_id, caught.value.count) == ('coordinate_read', 'ConformerNotFound', 0, 0)
    assert source.to_binary() == before
    rows = [[0., 0., 0.], [1., 2., 3.], [4., 5., 6.]]
    value = source.with_only_3d_conformer(rows)
    copied = value.coordinates_3d()
    copied[0][0] = 99
    np.testing.assert_array_equal(value.coordinates_3d(), rows)
    with pytest.raises(ck.Coordinate3DReadError) as caught: value.coordinates_3d(17)
    assert (caught.value.conformer_id, caught.value.count) == (17, 1)
