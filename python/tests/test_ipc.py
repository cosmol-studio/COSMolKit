"""Graph information content through the public Python bindings."""

import math

import cosmolkit as ck
import pytest


def test_distance_matrix_override_is_borrowed_and_checks_shape() -> None:
    mol = ck.Molecule.from_smiles("CC")
    separated = ck.Molecule.from_smiles("C.C").distance_matrix()
    assert (mol.ipc(), mol.avg_ipc()) == (2, 1)
    assert mol.ipc(separated) == mol.avg_ipc_with_distance_matrix(separated) == 0
    assert mol.avg_ipc([[0, 1], [0, 0]]) == 0
    assert mol.ipc_with_distance_matrix([[0, 0], [1, 0]]) == 2
    assert mol.avg_ipc(separated) == 0

    with pytest.raises(ck.IpcError) as error:
        mol.ipc([[0, 1], [0]])
    assert error.value.variant == ck.IpcErrorKind.InvalidDistanceMatrixShape
    assert (error.value.row, error.value.actual, error.value.expected) == (1, 1, 2)
    assert mol.ipc() == 2


def test_total_overflow_and_uncertified_reconstruction_are_numbers() -> None:
    raw = ck.SmilesParseParams(sanitize=False, remove_hs=False)
    path = ck.Molecule.from_smiles_with_params("C" * 1473, raw)
    assert path.ipc() == math.inf
    out_of_range = ck.Molecule.from_smiles_with_params("C" * 4000, raw)
    assert math.isnan(out_of_range.avg_ipc())
