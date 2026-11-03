"""Internal viewer transfer: real native Protein ownership and typed wire errors."""
import gc
import struct

import cosmolkit as ck
import pytest

PDB = "ATOM      1  CA AALA A   7      12.210  13.912  10.555  0.50 21.00           C  \nEND\n"


def test_protein_viewer_transfer_and_typed_errors():
    source = ck.Protein.from_pdb(PDB)
    payload = source.to_transfer_binary()
    assert isinstance(payload, bytes)
    restored = ck.Protein.from_transfer_binary(payload)
    assert restored is not source
    assert (restored.num_models, restored.num_chains, restored.num_residues, restored.num_atoms) == (1, 1, 1, 1)
    atom = restored.atoms[0]
    assert atom.name.as_str().strip() == "CA"
    assert atom.element == ck.Element.C
    assert struct.pack("ddd", *atom.position) == struct.pack("ddd", *source.atoms[0].position)
    assert restored.to_transfer_binary() == payload
    del source
    gc.collect()
    assert restored.atoms[0].position == [12.21, 13.912, 10.555]
    changed = bytearray(payload)
    changed[10] = ord("9")
    with pytest.raises(ck.ProteinTransferError) as failure:
        ck.Protein.from_transfer_binary(bytes(changed))
    assert failure.value.variant == ck.ProteinTransferErrorKind.VersionMismatch
    # The wire header uses Rust SemVer, not Python's PEP 440 package spelling.
    assert failure.value.expected == ck.version()
    for malformed in (b"", payload[:-1], payload + b"\0"):
        with pytest.raises(ck.ProteinTransferError) as failure:
            ck.Protein.from_transfer_binary(malformed)
        assert failure.value.variant == ck.ProteinTransferErrorKind.InvalidData
    with pytest.raises(TypeError):
        ck.Protein.from_transfer_binary("not bytes")
