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

    assert f64_bits(cosmolkit.calc_mol_wt(molecule)) == "403c0dd2f1a9fbe6"
    assert f64_bits(cosmolkit.calc_exact_mol_wt(molecule)) == "403c080349021ee0"
    assert cosmolkit.calc_mol_formula(molecule) == "C2H4"
    assert cosmolkit.calc_num_hbd(molecule) == 0
    assert cosmolkit.calc_num_hba(molecule) == 0
    assert f64_bits(cosmolkit.calc_fraction_csp3(molecule)) == "0000000000000000"

    logp, molar_refractivity = cosmolkit.calc_crippen_descriptors(molecule)
    assert f64_bits(logp) == "3fe9ab9f559b3d08"
    assert f64_bits(molar_refractivity) == "4026820c49ba5e36"
    assert f64_bits(cosmolkit.calc_tpsa(molecule)) == "0000000000000000"
    assert cosmolkit.calc_num_aromatic_rings(molecule) == 0

    for mode in ("default", "non_strict", "strict", "strict_linkages"):
        assert cosmolkit.calc_num_rotatable_bonds(molecule, mode=mode) == 0

    assert f64_bits(cosmolkit.calc_qed(molecule)) == "3fd60c3c2ca10e0d"
    assert molecule.to_smiles() == before


def test_descriptor_binding_options_and_errors_are_explicit():
    ethene_smiles = "C=C"
    ethene = cosmolkit.Molecule.from_smiles(ethene_smiles)
    assert f64_bits(cosmolkit.calc_mol_wt(ethene, only_heavy=True)) == (
        "403805a1cac08312"
    )
    assert f64_bits(cosmolkit.calc_exact_mol_wt(ethene, only_heavy=True)) == (
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
            logp, molar_refractivity = cosmolkit.calc_crippen_descriptors(
                branch_molecule,
                include_hs=include_hs,
                force=force,
            )
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
                cosmolkit.calc_tpsa(
                    sulfur,
                    force=force,
                    include_sandp=include_sandp,
                )
            ) == expected_bits

    deuterated_water = cosmolkit.Molecule.from_smiles("[2H]O")
    assert cosmolkit.calc_mol_formula(
        deuterated_water,
        separate_isotopes=True,
        abbreviate_h_isotopes=True,
    ) == "HDO"
    assert cosmolkit.calc_mol_formula(
        deuterated_water,
        separate_isotopes=True,
        abbreviate_h_isotopes=False,
    ) == "H[2H]O"

    with pytest.raises(ValueError, match="unsupported rotatable-bond mode"):
        _ = cosmolkit.calc_num_rotatable_bonds(
            deuterated_water,
            mode="unknown",  # pyright: ignore[reportArgumentType]
        )


def _high_feasibility_descriptor_snapshot(
    molecule: cosmolkit.Molecule,
) -> tuple[object, ...]:
    return (
        f64_bits(cosmolkit.calc_chi_0(molecule)),
        f64_bits(cosmolkit.calc_chi_1(molecule)),
        f64_bits(cosmolkit.calc_chi_3v(molecule)),
        f64_bits(cosmolkit.calc_kappa_2(molecule)),
        cosmolkit.calc_lipinski_hba(molecule),
        cosmolkit.calc_num_heteroatoms(molecule),
        tuple(cosmolkit.calc_mqns(molecule)),
        f64_bits(cosmolkit.calc_labute_asa(molecule)),
        tuple(f64_bits(value) for value in cosmolkit.calc_slogp_vsa(molecule)),
        tuple(f64_bits(value) for value in cosmolkit.calc_smr_vsa(molecule)),
    )


def test_high_feasibility_descriptor_bindings_compose_without_mutating_input():
    molecule = cosmolkit.Molecule.from_smiles("CC(O)c1ccncc1")
    before = molecule.to_smiles()
    expected = _high_feasibility_descriptor_snapshot(molecule)

    slogp_vsa = cosmolkit.calc_slogp_vsa(molecule)
    smr_vsa = cosmolkit.calc_smr_vsa(molecule)
    assert len(slogp_vsa) == 12
    assert len(smr_vsa) == 10
    assert f64_bits(cosmolkit.calc_slogp_vsa_1(molecule)) == f64_bits(slogp_vsa[0])
    assert f64_bits(cosmolkit.calc_slogp_vsa_12(molecule)) == f64_bits(
        slogp_vsa[11]
    )
    assert f64_bits(cosmolkit.calc_smr_vsa_1(molecule)) == f64_bits(smr_vsa[0])
    assert f64_bits(cosmolkit.calc_smr_vsa_10(molecule)) == f64_bits(smr_vsa[9])

    bins = [-0.2, 0.0, 0.25, 0.25, 0.8]
    assert len(cosmolkit.calc_slogp_vsa(molecule, bins=bins, force=True)) == 6
    assert len(cosmolkit.calc_smr_vsa(molecule, bins=bins, force=True)) == 6

    alpha, alpha_contributions = (
        cosmolkit.calc_hall_kier_alpha_with_contributions(molecule)
    )
    assert len(alpha_contributions) == molecule.num_atoms()
    assert f64_bits(ordered_f64_sum(alpha_contributions)) == f64_bits(alpha)

    asa, atom_contributions, hydrogen_contribution = (
        cosmolkit.calc_labute_asa_contributions(molecule, force=True)
    )
    assert len(atom_contributions) == molecule.num_atoms()
    assert f64_bits(
        ordered_f64_sum(atom_contributions) + hydrogen_contribution
    ) == f64_bits(asa)

    clone = molecule.sanitize()
    _ = cosmolkit.calc_chi_nv(clone, 5, force=True)
    _ = cosmolkit.calc_smr_vsa(clone, force=True)
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
