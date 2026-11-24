import cosmolkit as ck
import pytest


def test_canonical_order_then_extract_a_path_with_coordinates() -> None:
    source = ck.Molecule.from_smiles("OCC")
    source.add_3d_conformer_([[1, 0, 0], [2, 0, 0], [3, 0, 0]])
    ranks = source.canonical_atom_ranks()
    order = sorted(range(source.num_atoms), key=ranks.__getitem__)
    assert order == [2, 0, 1]
    ordered = source.with_atoms_renumbered(order)
    assert list(ordered.shortest_path(0, 2)) == [0, 2]
    bond = ordered.bond_between_atoms(0, 2)
    assert bond == 1
    assert ordered.bond_between_atoms(0, 1) is None
    assert ordered.adjacency_matrix().dimension == 3
    assert ordered.adjacency_matrix().get(0, 2) == 1

    selection = ck.SubgraphSelection.bonds([bond])
    result = ordered.extract_subgraph(selection)
    retained = ordered.retain_subgraph_(selection)
    assert list(retained.source_atom_indices) == [0, 2]
    assert list(retained.source_bond_indices) == [1]
    extracted, other = result.molecule, result.molecule
    assert list(extracted.coordinates_3d()[:, 0]) == [3, 2]
    assert list(source.coordinates_3d()[:, 0]) == [1, 2, 3]
    extracted.remove_atom_(0)
    del result, ordered
    assert other.num_atoms == 2
    assert list(retained.source_atom_indices) == [0, 2]


def test_selection_errors_and_self_insertion_leave_owned_results_usable() -> None:
    m = ck.Molecule.from_smiles("OCC")
    assert list(m.canonical_atom_ranks_in_fragment([0, 2], break_ties=False)) == [2, -1, 1]
    assert m.atom_environment(1, 0).atom_distances == {0: 0, 1: 1}
    assert m.subgraphs_in_range(1, 2)[1] == [[0, 1]]
    assert m.paths_of_length(2, use_atom_indices=True) == [[0, 1], [1, 2]]
    with pytest.raises(TypeError):
        m.paths_of_length(1, ck.PathEnumerationParams(), use_atom_indices=True)

    with pytest.raises(ck.OperationError) as invalid:
        m.retain_subgraph_(ck.SubgraphSelection.atoms_and_bonds([0], [0]))
    assert isinstance(invalid.value.__cause__, ck.IndexedGraphError)
    assert m.num_atoms == 3
    with pytest.raises(ck.OperationError):
        m.renumber_atoms_([1, 1, 0])
    assert m.num_bonds == 2

    clone = m.with_molecule_appended(m)
    m.insert_molecule_(m)
    assert m.to_smiles() == clone.to_smiles()
    assert m.num_atoms == 6 and m.bond_between_atoms(3, 4) == 2
    assert m.shortest_path(0, 3) == []


def test_sanitization_and_stereo_validation_precede_in_place_changes() -> None:
    m = ck.Molecule.from_smiles("C.C(F)(F)(F)(F)F", sanitize=False, remove_hs=False)
    selection = ck.SubgraphSelection.atoms(list(range(1, 7)))
    assert m.extract_subgraph(selection).molecule.num_bonds == 5
    peer = m.copy()
    with pytest.raises(ck.OperationError):
        m.retain_subgraph_(selection, sanitize=True)
    assert peer.num_atoms == 7
    assert m.adjacency_matrix().dimension == m.num_atoms

    m = ck.Molecule.from_smiles("F/C=C/Cl")
    spec = ck.BondSpec.new(0, 1, ck.BondOrder.Double).with_stereo(ck.BondStereo.Trans)
    with pytest.raises(ck.OperationError):
        m.replace_bond_(1, spec.with_stereo_atoms(3, 0))
    m.replace_bond_(1, spec.with_stereo_atoms(0, 3).with_prop("label", ck.PropertyValue.Int(42)))
    assert m.to_smiles() == ck.Molecule.from_smiles("F/C=C/Cl").to_smiles()
