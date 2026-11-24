import cosmolkit as ck
import pytest


def test_mapped_components_keep_their_own_molecule_and_indices() -> None:
    source = ck.Molecule.from_smiles("C1.CC.C1")
    first, second = source.fragments_with_atom_mapping()
    assert [first.source_atom_indices, second.source_atom_indices] == [(0, 3), (1, 2)]

    cut = ck.BondCut.at_bond(0)
    owned = first.molecule
    owned.cut_bonds_([cut])
    assert owned.num_atoms == 4
    assert source.num_atoms == 4
    assert first.molecule.to_smiles() == "CC"

    source.cut_bonds_([cut])
    del source
    indices, molecule = second.source_atom_indices, second.molecule
    del second
    assert indices == (1, 2)
    assert molecule.to_smiles() == "CC"


def test_mapped_extraction_configuration_controls_fragment_sanitization() -> None:
    source = ck.Molecule.from_smiles("C.C(F)(F)(F)(F)F", sanitize=False, remove_hs=False)
    with pytest.raises(ck.OperationError) as caught:
        source.fragments_with_atom_mapping()
    assert caught.value.fragments.component_index == 1

    params = ck.FragmentExtractionParams(sanitize=False)
    for parts in (
        source.fragments_with_atom_mapping(params),
        source.fragments_with_atom_mapping(sanitize=False),
        source.fragments_with_atom_mapping_with_params(params),
    ):
        assert [p.source_atom_indices for p in parts] == [(0,), (1, 2, 3, 4, 5, 6)]
        assert parts[1].molecule.num_bonds == 5

    with pytest.raises(TypeError):
        source.fragments_with_atom_mapping(params, sanitize=False)


def test_retained_hydrogen_maps_to_source_atom_65536() -> None:
    # A shallow tree keeps SMILES nesting small while exercising a 17-bit index.
    tree = "C"
    for _ in range(15):
        tree = f"C({tree}){tree}"
    source = ck.Molecule.from_smiles(f"{tree}C.[H]", sanitize=False, remove_hs=False)
    carbons, hydrogen = source.fragments_with_atom_mapping(sanitize=False)
    assert carbons.molecule.num_atoms == 65536
    assert hydrogen.source_atom_indices == (65536,)
    assert hydrogen.molecule.to_smiles() == "[H]"
