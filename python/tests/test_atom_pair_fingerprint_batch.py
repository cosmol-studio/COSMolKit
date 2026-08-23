from concurrent.futures import ThreadPoolExecutor

import cosmolkit
import pytest


SMILES = ["C", "CCO", "c1ccccc1O", "C[C@H](O)F"]


def test_atom_pair_batch_all_result_forms_match_ordered_scalar_calls():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(SMILES)
    kwargs = {
        "n_bits": 256,
        "include_chirality": True,
        "count_simulation": True,
        "count_bounds": [1, 3, 5],
        "num_bits_per_feature": 2,
    }

    explicit = batch.fingerprint_atom_pair_list(**kwargs, n_jobs=4, progress_bar=False)
    assert [value.on_bits() for value in explicit] == [
        cosmolkit.Molecule.from_smiles(smiles)
        .fingerprint_atom_pair(**kwargs)
        .on_bits()
        for smiles in SMILES
    ]

    sparse_count = batch.fingerprint_atom_pair_sparse_count_list(
        **kwargs, n_jobs=3, progress_bar=False
    )
    assert [value.nonzero_elements() for value in sparse_count] == [
        cosmolkit.Molecule.from_smiles(smiles)
        .fingerprint_atom_pair_sparse_count(**kwargs)
        .nonzero_elements()
        for smiles in SMILES
    ]

    count = batch.fingerprint_atom_pair_count_list(
        **kwargs, n_jobs=2, progress_bar=False
    )
    assert [value.nonzero_elements() for value in count] == [
        cosmolkit.Molecule.from_smiles(smiles)
        .fingerprint_atom_pair_count(**kwargs)
        .nonzero_elements()
        for smiles in SMILES
    ]

    sparse_bits = batch.fingerprint_atom_pair_sparse_bits_list(
        **kwargs, n_jobs=4, progress_bar=False
    )
    assert [value.on_bits() for value in sparse_bits] == [
        cosmolkit.Molecule.from_smiles(smiles)
        .fingerprint_atom_pair_sparse_bits(**kwargs)
        .on_bits()
        for smiles in SMILES
    ]

    outputs = batch.fingerprint_atom_pair_with_output_list(
        **kwargs, n_jobs=4, progress_bar=False
    )
    assert [value.fingerprint().on_bits() for value in outputs] == [
        value.on_bits() for value in explicit
    ]
    assert [value.additional_output().atom_counts() for value in outputs] == [
        cosmolkit.Molecule.from_smiles(smiles)
        .fingerprint_atom_pair_with_output(**kwargs)
        .additional_output()
        .atom_counts()
        for smiles in SMILES
    ]


def test_atom_pair_batch_keeps_invalid_positions_and_reports_operation_indices():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CC", "not-a-smiles", "CCC"], errors="keep"
    )
    values = batch.fingerprint_atom_pair_list(n_jobs=2, progress_bar=False)
    assert [value is not None for value in values] == [True, False, True]
    assert [(error.index(), error.operation()) for error in batch.errors()] == [
        (1, "batch.from_smiles_list")
    ]

    with pytest.raises(Exception, match="batch.atom_pair_fingerprint"):
        batch.fingerprint_atom_pair_list(use_2d=False, n_jobs=2, progress_bar=False)


def test_atom_pair_batch_thread_counts_progress_defaults_and_repeats_are_stable():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(SMILES)
    configured = batch.with_parallel_jobs(4).with_progress_bar(False)
    assert configured.parallel_jobs() == 4
    assert configured.progress_bar() is False
    baseline = configured.fingerprint_atom_pair_list(n_jobs=1)
    for n_jobs in (1, 2, 4):
        for _ in range(3):
            assert [value.on_bits() for value in configured.fingerprint_atom_pair_list(
                n_jobs=n_jobs
            )] == [value.on_bits() for value in baseline]

    with pytest.raises(ValueError, match="n_jobs"):
        configured.fingerprint_atom_pair_list(n_jobs=0)


def test_atom_pair_batch_is_safe_during_concurrent_mixed_family_calls():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(SMILES)
    expected = [
        value.on_bits()
        for value in batch.fingerprint_atom_pair_list(n_jobs=2, progress_bar=False)
    ]

    def atom_pair_call():
        return [
            value.on_bits()
            for value in batch.fingerprint_atom_pair_list(n_jobs=4, progress_bar=False)
        ]

    def morgan_call():
        return [
            value.on_bits()
            for value in batch.fingerprint_morgan_list(n_jobs=3, progress_bar=False)
        ]

    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(atom_pair_call) for _ in range(3)]
        morgan = executor.submit(morgan_call)
        assert [future.result() for future in futures] == [expected, expected, expected]
        assert len(morgan.result()) == len(SMILES)

    assert [
        value.on_bits()
        for value in batch.fingerprint_atom_pair_list(n_jobs=2, progress_bar=False)
    ] == expected
