"""Small default/configuration regressions, not corpus parity tests."""

import pytest
import cosmolkit as ck


def test_morgan_defaults_do_not_change_when_execution_options_are_supplied():
    molecule = ck.Molecule.from_smiles("c1ccccc1")
    batch = ck.MoleculeBatch.from_smiles_list(["c1ccccc1"])
    expected = [389, 1088, 1232, 1873]
    assert molecule.fingerprint_morgan().on_bits == expected
    results = [
        batch.fingerprint_morgan()[0],
        batch.fingerprint_morgan(n_jobs=1)[0],
        batch.fingerprint_morgan(n_jobs=2)[0],
        batch.fingerprint_morgan(execution=ck.BatchQueryParams())[0],
        batch.fingerprint_morgan(ck.MorganFingerprintParams(), execution=ck.BatchQueryParams())[0],
    ]
    for result in results:
        assert result is not None
        assert result.on_bits == expected
    for radius in (2, 3):
        options = ck.MorganFingerprintParams(radius=radius)
        result = batch.fingerprint_morgan(radius=radius)[0]
        assert result is not None
        assert result.on_bits == molecule.fingerprint_morgan_with_params(options, None).on_bits


@pytest.mark.parametrize("family", ["morgan", "atom_pair"])
def test_provenance_defaults_allow_omitted_collect_flag_and_worker_count(family: str):
    batch = ck.MoleculeBatch.from_smiles_list(["c1ccccc1"])
    if family == "morgan":
        options = ck.MorganFingerprintParams()
        reference = batch.fingerprint_morgan_with_output_with_params(
            options, True, ck.BatchQueryParams(n_jobs=1)
        )[0]
        results = [
            batch.fingerprint_morgan_with_output()[0],
            batch.fingerprint_morgan_with_output(n_jobs=1)[0],
            batch.fingerprint_morgan_with_output(n_jobs=2)[0],
            batch.fingerprint_morgan_with_output(execution=ck.BatchQueryParams())[0],
            batch.fingerprint_morgan_with_output(params=options)[0],
        ]
        without_output = batch.fingerprint_morgan_with_output(collect_additional_output=False)[0]
        with pytest.raises(TypeError, match="mutually exclusive"):
            # Deliberately invalid call: also verify runtime rejects it.
            batch.fingerprint_morgan_with_output(params=options, radius=options.radius)  # pyright: ignore[reportCallIssue]
    else:
        options = ck.AtomPairFingerprintParams()
        reference = batch.fingerprint_atom_pair_with_output_with_params(
            options, True, ck.BatchQueryParams(n_jobs=1)
        )[0]
        results = [
            batch.fingerprint_atom_pair_with_output()[0],
            batch.fingerprint_atom_pair_with_output(n_jobs=1)[0],
            batch.fingerprint_atom_pair_with_output(n_jobs=2)[0],
            batch.fingerprint_atom_pair_with_output(execution=ck.BatchQueryParams())[0],
            batch.fingerprint_atom_pair_with_output(params=options)[0],
        ]
        without_output = batch.fingerprint_atom_pair_with_output(collect_additional_output=False)[0]
        with pytest.raises(TypeError, match="mutually exclusive"):
            # Deliberately invalid call: also verify runtime rejects it.
            batch.fingerprint_atom_pair_with_output(params=options, min_distance=options.min_distance)  # pyright: ignore[reportCallIssue]
    assert reference is not None
    expected = reference.fingerprint.on_bits
    expected_counts = reference.additional_output.atom_counts
    for result in results:
        assert result is not None
        assert result.fingerprint.on_bits == expected
        assert result.additional_output.atom_counts == expected_counts
    assert without_output is not None
    assert without_output.fingerprint.on_bits == expected
    with pytest.raises(ck.BatchFingerprintOutputError):
        _ = without_output.additional_output
