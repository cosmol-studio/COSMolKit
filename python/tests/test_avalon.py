"""Avalon configuration transport, ownership and typed source errors."""
import cosmolkit as ck
import pytest


def test_avalon_parameters_and_error_chain():
    mol = ck.Molecule.from_smiles("CCO")
    before = mol.to_binary()
    params = ck.AvalonFingerprintParams()
    params.n_bits = 64
    params.is_query = False
    params.bit_flags = 32767
    fp = mol.fingerprint_avalon(params)
    assert fp.on_bits() == [6, 14, 30, 31, 42]
    assert mol.fingerprint_avalon(n_bits=64).on_bits() == fp.on_bits()
    assert "n_bits=64" in repr(params)
    assert mol.to_binary() == before
    params.n_bits = 7
    with pytest.raises(ck.AvalonFingerprintError) as caught:
        mol.fingerprint_avalon(params)
    assert caught.value.kind == "InvalidArguments"
    assert isinstance(caught.value.__cause__, ck.AvalonEngineError)
    assert caught.value.__cause__.reason == "Avalon n_bits must be at least 8"
