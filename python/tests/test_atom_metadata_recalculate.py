"""Explicit recalculation versus validity-checked cached metadata queries."""

import pytest
import cosmolkit as ck


def metadata_values(rows: list[ck.AtomMetadata]) -> list[tuple[int, int, int, int, int]]:
    return [
        (
            row.degree(),
            row.explicit_valence(),
            row.implicit_hydrogens(),
            row.total_hydrogens(),
            row.total_valence(),
        )
        for row in rows
    ]


def raw_carbon() -> ck.Molecule:
    builder = ck.MoleculeBuilder.new()
    _ = builder.add_atom(ck.AtomSpec(ck.Element.C))
    return builder.build()


def test_atom_metadata_default_recalculates_without_installing_cache() -> None:
    molecule = raw_carbon()
    expected = [(0, 0, 4, 4, 4)]
    assert metadata_values(molecule.atom_metadata()) == expected
    assert metadata_values(molecule.atom_metadata(recalculate=True)) == expected
    assert metadata_values(molecule.atom_metadata(True)) == expected
    with pytest.raises(ck.ValenceError) as caught:
        _ = molecule.atom_metadata(recalculate=False)
    assert caught.value.kind == "ExplicitValenceCacheNotInitialized"
    assert molecule.num_atoms() == 1
    assert molecule.num_bonds() == 0
    assert metadata_values(molecule.atom_metadata()) == expected


def test_atom_metadata_cached_query_uses_prepared_assignment_and_keeps_source_cold() -> None:
    source = raw_carbon()
    prepared = source.with_assigned_valence()
    expected = [(0, 0, 4, 4, 4)]
    assert metadata_values(prepared.atom_metadata(False)) == expected
    assert metadata_values(prepared.atom_metadata(recalculate=False)) == expected
    assert metadata_values(prepared.atom_metadata()) == expected
    with pytest.raises(ck.ValenceError):
        _ = source.atom_metadata(False)
    assert source.num_atoms() == prepared.num_atoms() == 1
    assert source.num_bonds() == prepared.num_bonds() == 0


def test_atom_metadata_empty_molecule_both_modes() -> None:
    molecule = ck.Molecule.new()
    assert molecule.atom_metadata() == []
    assert molecule.atom_metadata(True) == []
    assert molecule.atom_metadata(False) == []
