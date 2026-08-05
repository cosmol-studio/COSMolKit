from __future__ import annotations

import pickle

import cosmolkit
import numpy as np
import pytest


TETRAHEDRAL_COORDINATES = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
)
OPPOSITE_TETRAHEDRAL_COORDINATES = TETRAHEDRAL_COORDINATES[[0, 1, 3, 2]]
T_SHAPED_COORDINATES = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
)


def _center_tag(molecule: cosmolkit.Molecule) -> cosmolkit.ChiralTag:
    return molecule.atoms()[0].chiral_tag()


def _tetrahedral_molecule() -> cosmolkit.Molecule:
    return cosmolkit.Molecule.from_smiles("C(F)(Cl)Br")


def test_default_and_specific_conformer_selection_are_distinct() -> None:
    molecule = (
        _tetrahedral_molecule()
        .with_added_3d_conformer(TETRAHEDRAL_COORDINATES)
        .with_added_3d_conformer(OPPOSITE_TETRAHEDRAL_COORDINATES)
    )

    default_result = molecule.with_chiral_tags_from_structure()
    first_result = molecule.with_chiral_tags_from_structure(conf_id=0)
    second_result = molecule.with_chiral_tags_from_structure(conf_id=1)

    assert _center_tag(default_result) == _center_tag(first_result)
    assert _center_tag(first_result) != cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert _center_tag(second_result) != cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert _center_tag(first_result) != _center_tag(second_result)
    for conformer_index in range(molecule.num_conformers()):
        assert np.array_equal(
            default_result.coordinates_3d(conformer_index),
            molecule.coordinates_3d(conformer_index),
        )
        assert np.array_equal(
            second_result.coordinates_3d(conformer_index),
            molecule.coordinates_3d(conformer_index),
        )


def test_replace_modes_preserve_or_assign_existing_tag() -> None:
    tagged = cosmolkit.Molecule.from_smiles("[C@@H](F)(Cl)Br").with_only_3d_conformer(
        TETRAHEDRAL_COORDINATES
    )
    original_tag = _center_tag(tagged)

    preserved = tagged.with_chiral_tags_from_structure(replace_existing_tags=False)
    replaced = tagged.with_chiral_tags_from_structure(replace_existing_tags=True)

    assert original_tag != cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert _center_tag(preserved) == original_tag
    assert _center_tag(replaced) != cosmolkit.ChiralTag.CHI_UNSPECIFIED


def test_value_and_in_place_forms_have_explicit_mutation_semantics() -> None:
    source = _tetrahedral_molecule().with_only_3d_conformer(TETRAHEDRAL_COORDINATES)
    source_coordinates = source.coordinates_3d().copy()

    value_result = source.with_chiral_tags_from_structure()

    assert value_result is not source
    assert _center_tag(source) == cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert _center_tag(value_result) != cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert np.array_equal(source.coordinates_3d(), source_coordinates)
    assert np.array_equal(value_result.coordinates_3d(), source_coordinates)

    assert source.assign_chiral_tags_from_structure_() is None
    assert _center_tag(source) == _center_tag(value_result)
    assert np.array_equal(source.coordinates_3d(), source_coordinates)


def test_non_3d_conformer_is_a_source_defined_noop() -> None:
    molecule = _tetrahedral_molecule().with_only_3d_conformer(
        TETRAHEDRAL_COORDINATES,
        is_3d=False,
    )
    coordinates = molecule.coordinates_3d().copy()

    result = molecule.with_chiral_tags_from_structure()

    assert _center_tag(molecule) == cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert _center_tag(result) == cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert np.array_equal(result.coordinates_3d(), coordinates)


def test_no_conformer_is_a_source_defined_noop() -> None:
    molecule = _tetrahedral_molecule()
    before = molecule.mol_to_binary()

    result = molecule.with_chiral_tags_from_structure()

    assert molecule.mol_to_binary() == before
    assert result.mol_to_binary() == before


def test_assignment_survives_pickle_and_binary_workflow_boundaries() -> None:
    assigned = (
        _tetrahedral_molecule()
        .with_only_3d_conformer(TETRAHEDRAL_COORDINATES)
        .with_chiral_tags_from_structure()
    )

    restored_values = [
        pickle.loads(pickle.dumps(assigned, protocol=pickle.HIGHEST_PROTOCOL)),
        cosmolkit.Molecule.mol_from_binary(assigned.mol_to_binary()),
    ]
    for restored in restored_values:
        assert restored.mol_to_binary() == assigned.mol_to_binary()
        assert restored.to_smiles() == assigned.to_smiles()
        assert np.array_equal(restored.coordinates_3d(), TETRAHEDRAL_COORDINATES)


def test_environment_exact_zero_disables_nontetrahedral_assignment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    molecule = cosmolkit.Molecule.from_smiles("[P](F)(Cl)Br").with_only_3d_conformer(
        T_SHAPED_COORDINATES
    )

    monkeypatch.delenv("RDK_ENABLE_NONTETRAHEDRAL_STEREO", raising=False)
    enabled = molecule.with_chiral_tags_from_structure()
    monkeypatch.setenv("RDK_ENABLE_NONTETRAHEDRAL_STEREO", "0")
    disabled = molecule.with_chiral_tags_from_structure()

    assert _center_tag(enabled) == cosmolkit.ChiralTag.CHI_SQUAREPLANAR
    assert _center_tag(disabled) == cosmolkit.ChiralTag.CHI_UNSPECIFIED


def test_missing_specific_conformer_is_a_structured_value_error() -> None:
    molecule = _tetrahedral_molecule().with_only_3d_conformer(TETRAHEDRAL_COORDINATES)
    before = molecule.mol_to_binary()

    with pytest.raises(
        ValueError,
        match=r"^with_chiral_tags_from_structure: stereo error: Can't find conformation with ID: 17$",
    ):
        _ = molecule.with_chiral_tags_from_structure(conf_id=17)

    assert molecule.mol_to_binary() == before
