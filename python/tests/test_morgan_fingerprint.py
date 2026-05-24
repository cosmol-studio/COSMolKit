import cosmolkit
import pytest


def assert_morgan_is_deterministic(smiles, **kwargs):
    first = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan(**kwargs)
    second = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan(**kwargs)
    assert first.n_bits() == kwargs.get("n_bits", 2048)
    assert first.on_bits() == second.on_bits()
    return first


def test_morgan_fingerprint_default_and_tanimoto_are_self_consistent():
    ethanol = assert_morgan_is_deterministic("CCO", radius=2, n_bits=2048)

    propanol = cosmolkit.Molecule.from_smiles("CCCO").fingerprint_morgan()

    assert ethanol.on_bits()
    assert ethanol.tanimoto(ethanol) == pytest.approx(1.0)
    assert 0.0 <= ethanol.tanimoto(propanol) <= 1.0


def test_morgan_fingerprint_advanced_generators_and_counts_are_supported():
    feature_fp = assert_morgan_is_deterministic(
        "N[C@@H](C)C(=O)O",
        radius=2,
        n_bits=512,
        include_chirality=True,
        atom_invariants_generator="feature",
        num_bits_per_feature=2,
    )
    count_fp = assert_morgan_is_deterministic(
        "c1ccccc1O",
        radius=3,
        n_bits=256,
        atom_invariants_generator="connectivity",
        atom_invariants_include_ring_membership=False,
        bond_invariants_generator="morgan",
        bond_invariants_use_bond_types=False,
        count_simulation=True,
        count_bounds=[1, 2, 4, 8],
    )
    assert feature_fp.on_bits()
    assert count_fp.on_bits()


def test_morgan_fingerprint_custom_invariants_and_root_atoms_are_supported():
    smiles = "CC(C)O"
    fp = assert_morgan_is_deterministic(
        smiles,
        radius=2,
        n_bits=256,
        from_atoms=[0],
        custom_atom_invariants=[1, 2, 3, 4],
        custom_bond_invariants=[7, 8, 9],
        only_nonzero_invariants=True,
        include_redundant_environments=True,
    )
    assert fp.on_bits()


def test_morgan_additional_output_matches_returned_fingerprint():
    smiles = "CC(C)O"
    result = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan_with_output(
        radius=2, n_bits=256
    )
    output = result.additional_output()

    assert result.fingerprint().on_bits()
    assert output.atom_counts()
    assert len(output.atom_to_bits()) == cosmolkit.Molecule.from_smiles(smiles).num_atoms()
    assert set(output.bit_info_map()).issubset(set(result.fingerprint().on_bits()))
    assert set(output.atoms_per_bit()).issubset(set(result.fingerprint().on_bits()))


def test_morgan_batch_api_preserves_order_and_none_for_invalid_records():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "not-a-smiles", "CCCO"], errors="keep"
    )
    values = batch.fingerprint_morgan_list(n_bits=256, n_jobs=2)
    assert [value is not None for value in values] == [True, False, True]
    assert values[0] is not None
    assert values[0].on_bits() == cosmolkit.Molecule.from_smiles("CCO").fingerprint_morgan(
        n_bits=256
    ).on_bits()

    with_output = batch.fingerprint_morgan_with_output_list(n_bits=256, n_jobs=2)
    assert with_output[0] is not None
    assert with_output[0].additional_output().atom_counts()
    assert with_output[1] is None


def test_morgan_binding_rejects_unknown_generators():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    with pytest.raises(ValueError, match="atom_invariants_generator"):
        mol.fingerprint_morgan(atom_invariants_generator="unknown")
    with pytest.raises(ValueError, match="bond_invariants_generator"):
        mol.fingerprint_morgan(bond_invariants_generator="unknown")
