import cosmolkit
import pytest


def assert_morgan_is_deterministic(smiles, settings, call=None, *, atom_invariants=None, bond_invariants=None):
    generator = cosmolkit.MorganFingerprintGenerator(
        params=settings, atom_invariants=atom_invariants, bond_invariants=bond_invariants
    )
    first = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan_with_generator(generator, params=call)
    second = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan_with_generator(generator, params=call)
    assert first.n_bits() == settings.fp_size
    assert first.on_bits() == second.on_bits()
    return first


def test_morgan_fingerprint_default_and_tanimoto_are_self_consistent():
    ethanol = assert_morgan_is_deterministic("CCO", cosmolkit.MorganParams(radius=2, fp_size=2048))

    propanol = cosmolkit.Molecule.from_smiles("CCCO").fingerprint_morgan()

    assert ethanol.on_bits()
    assert ethanol.tanimoto(ethanol) == pytest.approx(1.0)
    assert 0.0 <= ethanol.tanimoto(propanol) <= 1.0


def test_morgan_fingerprint_advanced_generators_and_counts_are_supported():
    feature_fp = assert_morgan_is_deterministic(
        "N[C@@H](C)C(=O)O",
        cosmolkit.MorganParams(radius=2, fp_size=512, include_chirality=True, bits_per_feature=2),
        atom_invariants=cosmolkit.MorganAtomInvariantsGenerator.features(),
    )
    count_fp = assert_morgan_is_deterministic(
        "c1ccccc1O",
        cosmolkit.MorganParams(radius=3, fp_size=256, count_simulation=True, count_bounds=[1, 2, 4, 8]),
        atom_invariants=cosmolkit.MorganAtomInvariantsGenerator.connectivity(False),
        bond_invariants=cosmolkit.MorganBondInvariantsGenerator(use_bond_types=False),
    )
    assert feature_fp.on_bits()
    assert count_fp.on_bits()


def test_morgan_fingerprint_custom_invariants_and_root_atoms_are_supported():
    smiles = "CC(C)O"
    fp = assert_morgan_is_deterministic(
        smiles,
        cosmolkit.MorganParams(radius=2, fp_size=256, only_nonzero_invariants=True, include_redundant_environments=True),
        cosmolkit.MorganCallParams(from_atoms=[0], custom_atom_invariants=[1, 2, 3, 4], custom_bond_invariants=[7, 8, 9]),
    )
    assert fp.on_bits()


def test_morgan_additional_output_matches_returned_fingerprint():
    smiles = "CC(C)O"
    generator = cosmolkit.MorganFingerprintGenerator(params=cosmolkit.MorganParams(radius=2, fp_size=256))
    output = cosmolkit.FingerprintAdditionalOutput()
    output.allocate_atom_counts()
    output.allocate_atom_to_bits()
    output.allocate_bit_info_map()
    output.allocate_atoms_per_bit()
    fingerprint = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan_with_generator(
        generator, output=output
    )

    assert fingerprint.on_bits()
    assert output.atom_counts()
    assert len(output.atom_to_bits()) == cosmolkit.Molecule.from_smiles(smiles).num_atoms()
    assert set(output.bit_info_map()).issubset(set(fingerprint.on_bits()))
    assert set(output.atoms_per_bit()).issubset(set(fingerprint.on_bits()))


def test_morgan_batch_api_preserves_order_and_none_for_invalid_records():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "not-a-smiles", "CCCO"], errors="keep"
    )
    settings = cosmolkit.MorganParams(fp_size=256)
    execution = cosmolkit.BatchQueryParams(n_jobs=2)
    values = batch.fingerprint_morgan_list_with_generator_params(
        settings, None, None, cosmolkit.MorganCallParams(), execution
    )
    assert [value is not None for value in values] == [True, False, True]
    assert values[0] is not None
    assert values[0].on_bits() == cosmolkit.Molecule.from_smiles("CCO").fingerprint_morgan_with_generator(
        cosmolkit.MorganFingerprintGenerator(params=settings)
    ).on_bits()

    with_output = batch.fingerprint_morgan_with_output_list_with_generator_params(
        settings, None, None, cosmolkit.MorganCallParams(), True, execution
    )
    assert with_output[0] is not None
    assert with_output[0].additional_output().atom_counts()
    assert with_output[1] is None


def test_morgan_binding_rejects_unknown_generators():
    with pytest.raises(TypeError, match="MorganAtomInvariantsGenerator"):
        cosmolkit.MorganFingerprintGenerator(atom_invariants="unknown")
    with pytest.raises(TypeError, match="MorganBondInvariantsGenerator"):
        cosmolkit.MorganFingerprintGenerator(bond_invariants="unknown")
