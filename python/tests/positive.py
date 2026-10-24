"""Fixed consumer examples for the shared binding protocols.

HARD CONSTRAINT: Rust/PyO3 types and explicit projection metadata determine
public types. These examples are not a second business-API inventory.
Do not extend them for ordinary new APIs or move removed assertions elsewhere.
Changing this fixed protocol set requires an explicit framework design decision.
"""

from typing import assert_type
from collections.abc import Mapping

import numpy as np
import cosmolkit as ck


def native_enum_protocol() -> None:
    assert_type(ck.BondOrder.names, Mapping[str, ck.BondOrder])
    assert_type(ck.BondOrder.values, Mapping[int, ck.BondOrder])
    assert_type(ck.BondOrder.SINGLE.name, str)
    assert_type(ck.BondOrder(1), ck.BondOrder)


def facade_projection(molecule: ck.Molecule) -> None:
    assert_type(molecule.num_atoms, int)
    assert_type(molecule.with_hydrogens(), ck.Molecule)


def params_and_keywords(molecule: ck.Molecule) -> None:
    params = ck.MaccsFingerprintParams()
    params.n_bits = 166
    assert_type(params.n_bits, int)
    assert_type(molecule.fingerprint_maccs(), ck.Fingerprint)
    assert_type(molecule.fingerprint_maccs(params), ck.Fingerprint)
    assert_type(molecule.fingerprint_maccs(n_bits=None), ck.Fingerprint)


def integer_and_slice(batch: ck.FingerprintBatch, dataset: ck.SdfDataset) -> None:
    assert_type(batch[0], ck.Fingerprint | None)
    assert_type(batch[1:], ck.FingerprintBatch)
    assert_type(batch[[0, -1]], ck.FingerprintBatch)
    assert_type(batch[[True, False]], ck.FingerprintBatch)
    assert_type(dataset[0], ck.SdfRecord)
    assert_type(dataset[:1], ck.MoleculeBatch)


def owned_configuration(params: ck.SmartsParseParams) -> None:
    params.replacements = {"{A}": "O"}
    assert_type(params.replacements, Mapping[str, str])


def numpy_returns(fingerprint: ck.Fingerprint, batch: ck.FingerprintBatch) -> None:
    assert_type(
        fingerprint.to_numpy(), np.ndarray[tuple[int], np.dtype[np.uint8]]
    )
    assert_type(
        batch.to_numpy(), np.ndarray[tuple[int, int], np.dtype[np.uint8]]
    )


def precise_errors(text: str) -> None:
    try:
        ck.mol_from_smiles(text)
    except ck.SmilesError as error:
        assert_type(error.variant, ck.SmilesErrorKind)
        cause = error.source
        if isinstance(cause, ck.SmilesParseError):
            assert_type(cause.offset, int | None)
            assert_type(cause.message, str | None)
