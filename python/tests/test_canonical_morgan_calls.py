"""Canonical call/metadata transport; expected chemistry stays in source tests."""

import pytest

import cosmolkit as ck


@pytest.mark.parametrize(
    "method",
    [
        "fingerprint_morgan",
        "fingerprint_morgan_sparse",
        "fingerprint_morgan_count",
        "fingerprint_morgan_sparse_count",
    ],
)
def test_all_forms_fill_the_original_additional_output_and_preserve_input(method):
    molecule = ck.Molecule.from_smiles("CC(C)O")
    before = molecule.to_smiles()
    output = ck.FingerprintAdditionalOutput()
    output.allocate_atom_counts()
    output.allocate_atom_to_bits()
    output.allocate_bit_info_map()
    params = ck.MorganFingerprintParams(generator=ck.MorganParams(radius=2, fp_size=256))
    result = getattr(molecule, method + "_with_params")(params, output)
    ids = set(result.on_bits() if hasattr(result, "on_bits") else result.nonzero_elements())
    assert ids
    assert len(output.atom_counts()) == molecule.num_atoms()
    assert len(output.atom_to_bits()) == molecule.num_atoms()
    assert set(output.bit_info_map()) <= {value & 0xFFFFFFFF for value in ids}
    assert all(output.atom_counts())
    assert output.bit_paths() is None
    assert output.atoms_per_bit() is None
    assert molecule.to_smiles() == before
    # Detached read results cannot mutate the uniquely owned canonical output.
    snapshot = output.atom_counts()
    snapshot[0] = 0
    assert output.atom_counts()[0] > 0


def test_optional_empty_roots_are_distinct_and_call_parameters_are_immutable():
    molecule = ck.Molecule.from_smiles("CCO")
    params = ck.MorganFingerprintParams(from_atoms=[])
    assert params.from_atoms == []
    assert params.ignore_atoms is None
    assert molecule.fingerprint_morgan_with_params(params, None).on_bits() == []
    assert molecule.fingerprint_morgan_with_params(ck.MorganFingerprintParams(), None).on_bits()
    params.conformer_id = 7
    assert params.conformer_id == 7
    roots = params.from_atoms
    roots.append(0)
    assert params.from_atoms == []


def test_generator_errors_keep_the_source_cause_and_input():
    molecule = ck.Molecule.from_smiles("CCO")
    before = molecule.to_smiles()
    params = ck.MorganFingerprintParams(
        generator=ck.MorganParams(count_simulation=True, count_bounds=[])
    )
    with pytest.raises(ck.MorganReadError) as caught:
        molecule.fingerprint_morgan_with_params(params, None)
    assert caught.value.domain == "fingerprints"
    assert caught.value.kind == "Generator"
    assert caught.value.__cause__ is not None
    assert "bad count bounds provided" in str(caught.value)
    assert molecule.to_smiles() == before


def test_feature_and_custom_invariant_options_reach_the_owner():
    molecule = ck.Molecule.from_smiles("CCO")
    custom = ck.MorganFingerprintParams(
        from_atoms=[0],
        custom_atom_invariants=[1, 2, 3],
        custom_bond_invariants=[7, 8],
        generator=ck.MorganParams(radius=0),
    )
    assert molecule.fingerprint_morgan_sparse_count_with_params(custom, None).nonzero_elements() == {1: 1}
    features = ck.MorganFingerprintParams(invariants=ck.MorganInvariants.features())
    connectivity = ck.MorganFingerprintParams(invariants=ck.MorganInvariants.connectivity())
    assert molecule.fingerprint_morgan_sparse_count_with_params(features, None).nonzero_elements() != molecule.fingerprint_morgan_sparse_count_with_params(connectivity, None).nonzero_elements()
