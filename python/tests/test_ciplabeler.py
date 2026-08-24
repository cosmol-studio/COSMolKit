from __future__ import annotations

import cosmolkit
import pytest


def test_full_assignment_exposes_typed_atom_and_bond_state() -> None:
    molecule = cosmolkit.Molecule.from_smiles("C[C@H](F)Cl")
    labeled = molecule.with_cip_labels()

    assert labeled is not molecule
    assert not molecule.cip_computed()
    assert labeled.cip_computed()
    center = labeled.atoms()[1]
    assert center.cip_descriptor() in {"R", "S"}
    assert center.cip_neighbor_order() == [3, 2, 0]
    assert labeled.bonds()[0].cip_descriptor() is None


def test_selected_atom_and_bond_assignment_are_contextual_and_value_style() -> None:
    molecule = cosmolkit.Molecule.from_smiles("F/C=C/F")
    atom_selected = molecule.with_cip_labels(atoms=[1])
    bond_selected = molecule.with_cip_labels(bonds=[1])
    both_selected = molecule.with_cip_labels(atoms=[1], bonds=[1])

    assert atom_selected.atoms()[1].cip_descriptor() is None
    assert atom_selected.bonds()[1].cip_descriptor() is None
    assert bond_selected.bonds()[1].cip_descriptor() in {"E", "Z"}
    assert bond_selected.atoms()[1].cip_descriptor() is None
    assert both_selected.bonds()[1].cip_descriptor() in {"E", "Z"}
    assert molecule.atoms()[1].cip_descriptor() is None


def test_in_place_assignment_and_empty_dispatch_match_public_contract() -> None:
    molecule = cosmolkit.Molecule.from_smiles("C[C@H](F)Cl")
    assert molecule.assign_cip_labels_() is None
    assert molecule.cip_computed()
    assert molecule.atoms()[1].cip_descriptor() in {"R", "S"}

    empty = cosmolkit.Molecule.from_smiles("C[C@H](F)Cl")
    empty.assign_cip_labels_(atoms=[], bonds=[])
    assert empty.cip_computed()
    assert empty.atoms()[1].cip_descriptor() in {"R", "S"}


def test_recursion_limit_errors_are_structured_and_do_not_mutate_source() -> None:
    molecule = cosmolkit.Molecule.from_smiles(
        "Cc1ccc(S(=O)(=O)O)cc1.O=C(CNC1CC1)NC[C@H]1CC[C@]2(CC1)OO[C@]1(O2)C2CC3CC(C2)CC1C3"
    )
    before = molecule.mol_to_binary()
    with pytest.raises(ValueError, match="Max Iterations Exceeded"):
        molecule.assign_cip_labels_(max_recursive_iterations=1)
    assert molecule.mol_to_binary() == before
