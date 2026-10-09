"""Small regressions for public batch help and its executable examples."""

import inspect
from collections.abc import Callable
from typing import cast

import pytest
import cosmolkit as ck


@pytest.mark.parametrize(
    ("name", "details"),
    [
        ("from_smiles_list", ["input order", 'sanitize=True', 'errors="keep"', "SmilesParseParams"]),
        ("to_smiles_list", ["input order", "canonical, isomeric", "None", "SmilesWriteParams"]),
        ("with_hydrogens", ["source batch unchanged", 'errors="raise"', 'errors="keep"', "AddHsParams"]),
        ("fingerprint_morgan_list", ["input order", "radius 3, 2048 bits", "None", "MorganFingerprintParams"]),
    ],
)
def test_public_batch_help_preserves_native_documentation(name: str, details: list[str]):
    method = cast(Callable[..., object], getattr(ck.MoleculeBatch, name))
    doc = inspect.getdoc(method)
    assert doc is not None
    for detail in details:
        assert detail in doc
    assert "n_jobs" in doc
    native = cast(object, inspect.unwrap(method))
    assert inspect.getdoc(native) == doc
    configured = cast(object, getattr(ck.MoleculeBatch, name + "_with_params"))
    assert inspect.getdoc(configured)


def test_batch_help_examples_match_return_and_error_behavior():
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC"], errors="keep", n_jobs=1)
    before = batch.to_smiles_list()
    assert before == ["CCO", None]
    assert batch.to_smiles_list(canonical=False, n_jobs=1) == before
    assert [error.index() for error in batch.errors()] == [1]

    with pytest.raises(ck.BatchValidationError):
        _ = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC"])
    with pytest.raises(ck.BatchValidationError):
        _ = batch.with_hydrogens(errors="raise")
    assert batch.with_hydrogens().valid_mask() == [True, False]
    with_hs = batch.with_hydrogens(errors="keep", n_jobs=1)
    assert len(with_hs) == len(batch)
    assert with_hs[1] is None
    assert [error.as_dict() for error in with_hs.errors()] == [
        error.as_dict() for error in batch.errors()
    ]
    original, hydrogenated = batch[0], with_hs[0]
    assert original is not None and hydrogenated is not None
    assert original.num_atoms() == 3
    assert hydrogenated.num_atoms() == 9

    fingerprints = batch.fingerprint_morgan_list(generator=ck.MorganParams(radius=2), n_jobs=1)
    assert len(fingerprints) == len(batch)
    assert fingerprints[0] is not None
    assert fingerprints[0].n_bits() == 2048
    assert fingerprints[1] is None
    assert batch.to_smiles_list() == before
