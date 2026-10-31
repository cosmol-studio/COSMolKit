"""Small regressions for public batch help and its executable examples."""

import pytest
import cosmolkit as ck


def test_batch_help_examples_match_return_and_error_behavior():
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC"], errors="keep", n_jobs=1)
    with pytest.raises(ck.BatchValidationError, match="query not started"):
        batch.to_smiles()
    clean = batch.with_valid_records()
    before = clean.to_smiles()
    assert before == ["CCO"]
    assert clean.to_smiles(canonical=False, n_jobs=1) == before
    assert [error.index for error in batch.errors] == [1]

    with pytest.raises(ck.BatchValidationError):
        _ = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC"])
    with pytest.raises(ck.BatchValidationError):
        _ = batch.with_hydrogens(errors="raise")
    assert batch.with_hydrogens().valid_mask == [True, False]
    with_hs = batch.with_hydrogens(errors="keep", n_jobs=1)
    assert len(with_hs) == len(batch)
    assert with_hs[1] is None
    assert [{"index": error.index, "operation": error.operation, "message": error.message} for error in with_hs.errors] == [
        {"index": error.index, "operation": error.operation, "message": error.message} for error in batch.errors
    ]
    original, hydrogenated = batch[0], with_hs[0]
    assert original is not None and hydrogenated is not None
    assert original.num_atoms == 3
    assert hydrogenated.num_atoms == 9

    with pytest.raises(ck.BatchValidationError, match="query not started"):
        batch.fingerprint_morgan(radius=2, n_jobs=1)
    fingerprints = clean.fingerprint_morgan(radius=2, n_jobs=1)
    assert len(fingerprints) == len(clean)
    assert fingerprints[0] is not None
    assert fingerprints[0].n_bits == 2048
    assert batch.with_valid_records().to_smiles() == before
