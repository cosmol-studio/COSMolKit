import struct
from concurrent.futures import ThreadPoolExecutor

import cosmolkit
import pytest


def f64_bits(value: float) -> str:
    return struct.pack(">d", value).hex()


def ordered_f64_sum(values: list[float]) -> float:
    result = 0.0
    for value in values:
        result += value
    return result


def test_descriptor_bindings_preserve_exact_core_values_and_input():
    molecule = cosmolkit.Molecule.from_smiles("C=C")
    before = molecule.to_smiles()

    assert f64_bits(molecule.molecular_weight()) == "403c0dd2f1a9fbe6"
    assert f64_bits(molecule.exact_molecular_weight()) == "403c080349021ee0"
    assert molecule.molecular_formula() == "C2H4"
    assert molecule.num_hbd() == 0
    assert molecule.num_hba() == 0
    assert f64_bits(molecule.fraction_csp3()) == "0000000000000000"

    totals = molecule.crippen_descriptors()
    logp, molar_refractivity = totals.logp, totals.molar_refractivity
    assert f64_bits(logp) == "3fe9ab9f559b3d08"
    assert f64_bits(molar_refractivity) == "4026820c49ba5e36"
    assert f64_bits(molecule.tpsa()) == "0000000000000000"
    assert molecule.num_aromatic_rings() == 0

    for mode in (
        cosmolkit.RotatableBondsOptions.Default,
        cosmolkit.RotatableBondsOptions.NonStrict,
        cosmolkit.RotatableBondsOptions.Strict,
        cosmolkit.RotatableBondsOptions.StrictLinkages,
    ):
        assert molecule.num_rotatable_bonds_with_params(mode) == 0

    assert f64_bits(molecule.qed()) == "3fd60c3c2ca10e0d"
    assert molecule.to_smiles() == before


def test_descriptor_binding_options_and_errors_are_explicit():
    ethene_smiles = "C=C"
    ethene = cosmolkit.Molecule.from_smiles(ethene_smiles)
    assert f64_bits(ethene.molecular_weight_with_params(True)) == (
        "403805a1cac08312"
    )
    assert f64_bits(ethene.exact_molecular_weight_with_params(True)) == (
        "4038000000000000"
    )

    expected_crippen = {
        False: ("3fd3da5119ce075f", "401c1a9fbe76c8b4"),
        True: ("3fe9ab9f559b3d08", "4026820c49ba5e36"),
    }
    for include_hs in (False, True):
        for force in (False, True):
            # The option matrix measures each branch from the same fresh input.
            # Cache-order behavior is covered separately by the core parity tests.
            branch_molecule = cosmolkit.Molecule.from_smiles(ethene_smiles)
            totals = branch_molecule.crippen_descriptors_with_params(include_hs, force)
            logp, molar_refractivity = totals.logp, totals.molar_refractivity
            assert (f64_bits(logp), f64_bits(molar_refractivity)) == (
                expected_crippen[include_hs]
            )

    sulfur = cosmolkit.Molecule.from_smiles("F[C@@H]1O[C@H](Cl)S1")
    for force in (False, True):
        for include_sandp, expected_bits in (
            (False, "402275c28f5c28f6"),
            (True, "404143d70a3d70a4"),
        ):
            assert f64_bits(
                sulfur.tpsa_with_params(include_sandp, force)
            ) == expected_bits

    deuterated_water = cosmolkit.Molecule.from_smiles("[2H]O")
    assert deuterated_water.molecular_formula_with_params(True, True) == "HDO"
    assert deuterated_water.molecular_formula_with_params(True, False) == "H[2H]O"

    # Canonical enum spellings use the same selector as native enum inputs.
    assert deuterated_water.num_rotatable_bonds_with_params("strict") == (
        deuterated_water.num_rotatable_bonds_with_params(cosmolkit.RotatableBondsOptions.Strict)
    )
    with pytest.raises(ValueError, match="RotatableBondsOptions"):
        _ = deuterated_water.num_rotatable_bonds_with_params(
            "unknown",
        )


def _high_feasibility_descriptor_snapshot(
    molecule: cosmolkit.Molecule,
) -> tuple[object, ...]:
    return (
        f64_bits(molecule.chi_0()),
        f64_bits(molecule.chi_1()),
        f64_bits(molecule.chi_3_v()),
        f64_bits(molecule.kappa_2()),
        molecule.lipinski_hba(),
        molecule.num_heteroatoms(),
        tuple(molecule.mqns()),
        f64_bits(molecule.labute_asa()),
        tuple(f64_bits(value) for value in molecule.slogp_vsa()),
        tuple(f64_bits(value) for value in molecule.smr_vsa()),
    )


def test_high_feasibility_descriptor_bindings_compose_without_mutating_input():
    molecule = cosmolkit.Molecule.from_smiles("CC(O)c1ccncc1")
    before = molecule.to_smiles()
    expected = _high_feasibility_descriptor_snapshot(molecule)

    slogp_vsa = molecule.slogp_vsa()
    smr_vsa = molecule.smr_vsa()
    assert len(slogp_vsa) == 12
    assert len(smr_vsa) == 10
    assert f64_bits(molecule.slogp_vsa_1()) == f64_bits(slogp_vsa[0])
    assert f64_bits(molecule.slogp_vsa_12()) == f64_bits(
        slogp_vsa[11]
    )
    assert f64_bits(molecule.smr_vsa_1()) == f64_bits(smr_vsa[0])
    assert f64_bits(molecule.smr_vsa_10()) == f64_bits(smr_vsa[9])

    bins = [-0.2, 0.0, 0.25, 0.25, 0.8]
    assert len(molecule.slogp_vsa_with_params(bins, True)) == 6
    assert len(molecule.smr_vsa_with_params(bins, True)) == 6

    alpha, alpha_contributions = (
        molecule.hall_kier_alpha_with_contributions()
    )
    assert len(alpha_contributions) == molecule.num_atoms()
    assert f64_bits(ordered_f64_sum(alpha_contributions)) == f64_bits(alpha)

    contributions = molecule.labute_asa_contributions_with_params(True, True)
    asa, atom_contributions, hydrogen_contribution = (
        contributions.asa, contributions.atom_contributions, contributions.hydrogen_contribution
    )
    assert len(atom_contributions) == molecule.num_atoms()
    assert f64_bits(
        ordered_f64_sum(atom_contributions) + hydrogen_contribution
    ) == f64_bits(asa)

    clone = molecule.sanitize()
    _ = clone.chi_n_v_with_params(5, True)
    _ = clone.smr_vsa_with_params(None, True)
    assert _high_feasibility_descriptor_snapshot(clone) == expected
    assert _high_feasibility_descriptor_snapshot(molecule) == expected
    assert molecule.to_smiles() == before


def test_high_feasibility_descriptor_bindings_are_parallel_read_deterministic():
    molecule = cosmolkit.Molecule.from_smiles("CC(O)c1ccncc1")
    expected = _high_feasibility_descriptor_snapshot(molecule)
    with ThreadPoolExecutor(max_workers=8) as executor:
        actual = list(
            executor.map(
                _high_feasibility_descriptor_snapshot,
                [molecule] * 32,
            )
        )
    assert actual == [expected] * 32
