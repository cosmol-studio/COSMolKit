"""Local representation and binding regressions, not corpus parity."""

import numpy as np
import pytest

import cosmolkit as ck


def expected_array(fp):
    expected = np.zeros(fp.n_bits(), dtype=np.uint8)
    expected[fp.on_bits()] = 1
    return expected


@pytest.mark.parametrize("width,bits", [(0, []), (1, [0]), (8, [1, 4]), (65, [0, 63, 64])])
def test_scalar_numpy_preserves_bit_order_width_and_snapshot(width, bits):
    fp = ck.Fingerprint.from_on_bits(width, bits)
    actual = fp.to_numpy()
    assert actual.shape == (width,)
    assert actual.dtype == np.uint8
    assert actual.flags.c_contiguous and actual.flags.writeable
    np.testing.assert_array_equal(actual, expected_array(fp))
    actual[:] = 0
    np.testing.assert_array_equal(fp.to_numpy(), expected_array(fp))


@pytest.mark.parametrize("method", ["fingerprint_morgan", "fingerprint_avalon", "fingerprint_maccs",
                                   "fingerprint_pattern", "fingerprint_layered", "fingerprint_topological",
                                   "fingerprint_atom_pair", "fingerprint_topological_torsion"])
def test_scalar_numpy_is_shared_across_dense_fingerprint_algorithms(method):
    mol = ck.mol_from_smiles("CCO")
    fp = getattr(mol, method)()
    np.testing.assert_array_equal(fp.to_numpy(), expected_array(fp))
    assert fp.to_numpy().dtype == np.uint8
    assert mol.to_smiles() == "CCO"


@pytest.mark.parametrize("method", ["fingerprint_morgan_list", "fingerprint_pattern_list",
                                   "fingerprint_layered_list", "fingerprint_atom_pair_list",
                                   "fingerprint_topological_torsion_list"])
def test_batch_numpy_is_shared_across_dense_fingerprint_algorithms(method):
    batch = ck.mols_from_smiles_list(["CCO", "c1ccccc1"])
    fps = getattr(batch, method)()
    assert isinstance(fps, ck.FingerprintBatch)
    assert len(fps) == 2
    matrix = fps.to_numpy()
    assert matrix.shape == (2, fps[0].n_bits())
    assert matrix.dtype == np.uint8
    assert matrix.flags.c_contiguous and matrix.flags.writeable
    for index, fp in enumerate(fps):
        np.testing.assert_array_equal(matrix[index], fp.to_numpy())
    matrix[:] = 0
    np.testing.assert_array_equal(fps.to_numpy()[0], expected_array(fps[0]))
    assert batch.to_smiles_list() == ["CCO", "c1ccccc1"]


def test_batch_configured_morgan_and_generator_use_the_same_container():
    batch = ck.mols_from_smiles_list(["CCO", "CC"])
    fps = batch.fingerprint_morgan_list(radius=2, fp_size=65)
    assert isinstance(fps, ck.FingerprintBatch)
    assert fps.to_numpy().shape == (2, 65)
    params = ck.MorganFingerprintParams(generator=ck.MorganParams(radius=2, fp_size=65))
    configured = batch.fingerprint_morgan_list_with_params(params, ck.BatchQueryParams())
    np.testing.assert_array_equal(fps.to_numpy(), configured.to_numpy())
    generator = ck.MorganFingerprintGenerator(params=ck.MorganParams(radius=2, fp_size=65))
    generated = generator.fingerprints([ck.mol_from_smiles("CCO"), ck.mol_from_smiles("CC")])
    assert isinstance(generated, ck.FingerprintBatch)
    np.testing.assert_array_equal(generated.to_numpy(), fps.to_numpy())


def test_batch_index_slice_iteration_and_empty_width():
    fp = ck.Fingerprint.from_on_bits(8, [1, 4])
    fps = ck.FingerprintBatch([fp, None, fp])
    assert fps[0] is fp and fps[-1] is fp and fps[1] is None
    assert list(fps) == [fp, None, fp]
    assert list(fps[::-1]) == [fp, None, fp]
    assert isinstance(fps[::2], ck.FingerprintBatch)
    assert list(fps[1::2**63 - 1]) == [None]
    np.testing.assert_array_equal(fps[::2].to_numpy(), np.stack([fp.to_numpy()] * 2))
    assert fps[:0].to_numpy().shape == (0, 8)
    assert ck.FingerprintBatch([]).to_numpy().shape == (0, 0)
    assert ck.mols_from_smiles_list([]).fingerprint_morgan_list().to_numpy().shape == (0, 0)
    assert ck.FingerprintBatch([ck.Fingerprint.from_on_bits(0, [])]).to_numpy().shape == (1, 0)
    with pytest.raises(IndexError):
        _ = fps[3]
    with pytest.raises(IndexError):
        _ = fps[-4]
    with pytest.raises(TypeError):
        _ = fps[1.5]
    with pytest.raises(ValueError):
        _ = fps[::0]
    assert "FingerprintBatch" in repr(fps)


def test_batch_numpy_rejects_failed_rows_and_mismatched_widths():
    fp = ck.Fingerprint.from_on_bits(8, [1])
    with pytest.raises(ValueError, match="row 1 is None"):
        ck.FingerprintBatch([fp, None]).to_numpy()
    with pytest.raises(ValueError, match="row 1 has 9 bits; expected 8"):
        ck.FingerprintBatch([fp, ck.Fingerprint.from_on_bits(9, [])]).to_numpy()
    batch = ck.mols_from_smiles_list(["CCO", "C1CC"], errors="keep")
    fps = batch.fingerprint_morgan_list()
    assert len(fps) == 2 and fps[1] is None
    with pytest.raises(ValueError, match="row 1 is None"):
        fps.to_numpy()
