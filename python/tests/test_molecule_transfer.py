"""Binary transfer boundary, independent of archive/pickle compatibility."""
import pickle
import struct

import cosmolkit as ck
import pytest


def test_transfer_owns_the_value_and_preserves_coordinates_and_original_archive():
    source = ck.mol_from_smiles("N[C@@H](C)C(=O)O").with_2d_coordinates()
    coordinates = [(i + 0.125, -0.0, -i * 0.25) for i in range(source.num_atoms)]
    source = source.with_only_3d_conformer(coordinates)
    archive = source.to_binary()
    payload = source.to_transfer_binary()
    assert isinstance(payload, bytes)
    restored = ck.Molecule.from_transfer_binary(payload)
    assert restored is not source
    assert restored.num_atoms == source.num_atoms
    assert restored.num_bonds == source.num_bonds
    assert restored.to_smiles() == source.to_smiles()
    assert restored.coordinates_2d.tobytes() == source.coordinates_2d.tobytes()
    actual = restored.coordinates_3d(source.conformers_3d[0].id)
    assert [struct.pack("ddd", *row) for row in actual] == [struct.pack("ddd", *row) for row in coordinates]
    assert restored.to_transfer_binary() == payload
    assert source.to_binary() == archive
    assert pickle.loads(pickle.dumps(source)).to_binary() == archive
    with pytest.raises(ck.MoleculeTransferError) as failure:
        ck.Molecule.from_transfer_binary(archive)
    assert isinstance(failure.value.__cause__, ck.PickleError)


def test_transfer_version_and_malformed_data_are_structured_errors():
    payload = ck.mol_from_smiles("CCO").to_transfer_binary()
    changed = bytearray(payload)
    changed[10] = ord("9")
    with pytest.raises(ck.MoleculeTransferError) as failure:
        ck.Molecule.from_transfer_binary(bytes(changed))
    error = failure.value
    assert error.variant == ck.MoleculeTransferErrorKind.VersionMismatch
    # The wire header uses Rust SemVer, not Python's PEP 440 package spelling.
    assert error.expected == ck.version()
    assert error.actual != error.expected
    assert error.__cause__ is None
    for malformed in (b"", payload[:-1], payload + b"\0"):
        with pytest.raises(ck.MoleculeTransferError) as failure:
            ck.Molecule.from_transfer_binary(malformed)
        assert failure.value.variant == ck.MoleculeTransferErrorKind.Codec
        assert isinstance(failure.value.__cause__, ck.PickleError)
